! WRF:MODEL_LAYER:PHYSICS
 
    MODULE module_stoch_les_inflow_perts

    USE module_model_constants    

#ifdef DM_PARALLEL
    USE module_dm
#endif

    
    CONTAINS

!!=======================================================================

    SUBROUTINE force_down_meso_pblh( m_pblh, pblh,                        &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte         )

!-----------------------------------------------------------------------
!
! Assign pbl height to m_pblh array to be forced down for nested LES
!
!
    IMPLICIT NONE

    INTEGER, INTENT( IN )  &
    :: ids, ide, jds, jde, kds, kde,  &
       ims, ime, jms, jme, kms, kme,  &
       its, ite, jts, jte, kts, kte

    REAL , DIMENSION( ims:ime, jms:jme ), INTENT( IN  ) :: pblh
    REAL , DIMENSION( ims:ime, jms:jme ), INTENT( OUT ) :: m_pblh
    INTEGER :: itimestep
    INTEGER :: i,j
         

    DO j = jts, jte
       DO i = its, ite

          m_pblh(i,j) = pblh(i,j)

       END DO
    END DO
 
!    print*,'in m_pblh, m_pblh(its,jts)',m_pblh(its,jts)   

    END SUBROUTINE force_down_meso_pblh

!=======================================================================
 
    SUBROUTINE calc_pert_t( les_pert_opt,                      &
                            m_pblh_opt,                        &
                            prttms, prtdt, prtnk,              &
                            prtz, prtseed, pert_t,             &
                            mpblh,                             &
                            t, u, v, rdz,                      &
                            dx, dt,                            &
                            ids, ide, jds, jde, kds, kde,      &
                            ims, ime, jms, jme, kms, kme,      &
                            its, ite, jts, jte, kts, kte       )

!-----------------------------------------------------------------------

    IMPLICIT NONE

#ifdef DM_PARALLEL
    INCLUDE 'mpif.h'
#endif

    INTEGER, INTENT( IN    )  &
    :: ids, ide, jds, jde, kds, kde,  &
       ims, ime, jms, jme, kms, kme,  &
       its, ite, jts, jte, kts, kte

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( INOUT ) :: t         ! potential temperature                 [k]
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( INOUT ) :: pert_t    ! potential temperature perturbation    [k]
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN    ) :: u         ! zonal wind component                  [m/s]
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN    ) :: v         ! meridional wind component             [m/s]
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN    ) :: rdz       ! inverse vertical grid spacing         [1/m]
    REAL, DIMENSION( ims:ime, jms:jme ),          INTENT( IN    ) :: mpblh    ! PBL height from messoscale PBL scheme [m]
    REAL, DIMENSION( kms:kme ),                   INTENT( INOUT ) :: prttms    ! time since last perturbation refresh  [s]
    REAL, DIMENSION( kms:kme ),                   INTENT( INOUT ) :: prtdt     ! time scale to refresh perturbations   [s]
    
    INTEGER, DIMENSION( kms:kme ),                INTENT( INOUT ) :: prtseed   ! 1D variable to store variable-length seed array

    INTEGER, INTENT( IN    ) :: les_pert_opt           ! which perturbation option to use
    INTEGER, INTENT( IN    ) :: m_pblh_opt             ! whether to use (1) mesoscale PBL height or not (0)
    REAL,    INTENT( INOUT ) :: prtz                   ! height over swhich to apply perturbations
    INTEGER, INTENT( INOUT ) :: prtnk                  ! number of vertical grid cells over which to apply perturbations
    REAL,    INTENT( IN    ) :: dt                     ! model time step
    REAL,    INTENT( IN    ) :: dx                     ! horizontal grid spacing

    REAL, DIMENSION( its:ite, kts:kte, jts:jte ) :: z  ! height above surface at midpoint levels   
    REAL, DIMENSION( its:ite, jts:jte ) :: pblh        ! height above surface at midpoint levels   
    
    INTEGER :: i, j, k, big_k, slab_k, m            ! loop variables
    INTEGER :: i_start, i_end, j_start, j_end       ! patch endpoints
    INTEGER :: i_seed                               ! size of seed
    INTEGER :: seedsum                              ! variable for summing the seed
    INTEGER :: ngc_h, ngc_v                         ! number of gridcells in i and j directions for each perturbation cell
    INTEGER :: ni, nj                               ! number of perturbation cells in i and j directions (for random number query)
    INTEGER :: ncells_h
    INTEGER :: n_slabs_k                            ! number of perturbations cells in horizontal and vertical directions
    INTEGER :: k_slab_start, k_slab_end             ! vertical indices over which to apply perturbations
    INTEGER :: north, south, east, west, nsew       !logic for identification of LBs to perturb
    INTEGER :: k_geo                                ! k index at which h_wg is computed
    REAL    :: h_geo                                ! height at which to estimate ug
    REAL    :: ugeosum, vgeosum
    REAL    :: ugeolbsum, vgeolbsum
    REAL    :: ugeolbavg, vgeolbavg, wsgeolbavg
    REAL    :: uslabsum, vslabsum 
    REAL    :: uslablbsum, vslablbsum
    REAL    :: uslablbavg, vslablbavg, wsslablbavg  !
    REAL    :: anglelbavg 
    REAL    :: num_pts_in_sum                       ! Variables for computing averages
    REAL    :: tpertmag                             ! Magnitude of thermal perturbation
    REAL    :: dz                                   ! Depth
    REAL    :: pblh_def 
    REAL    :: sf                                  ! scaling factor for perturbation profile
    REAL    :: sf2 
    REAL    :: ek_opt = 0.20                        ! optimal value of the perturbation Ekman #
    REAL    :: pio2 = piconst/2.0
    REAL    :: pio4 = piconst/4.0    
    REAL    :: max_pblh, max_pert_z
    REAL    :: min_slab_z
    
    INTEGER :: idum, jdum
    
    REAL    :: counter
    
    INTEGER, dimension( : ), allocatable :: seed              ! random number seed
    REAL, dimension( :, : ), allocatable :: pxs,pxe,pys,pye   ! random number arrays for each lateral boundary
  
#ifdef DM_PARALLEL
    INTEGER :: ierr
    INTEGER :: tag
    INTEGER :: master
    INTEGER :: status(MPI_STATUS_SIZE)
#endif

#ifdef DM_PARALLEL
    master = 0
    tag    = 0
#endif

! END DECLARATIONS


    ngc_h = 8      !perturbation cell width in horizontal directions in i,j-index values
!    ngc_v = 1      !perturbation cell thickness in vertical direction in k-index values
    ngc_v = 2      !perturbation cell thickness in vertical direction in k-index values
    !DJW modified to increase depth of slabs / cells / blocks
    ncells_h = 3   !number of perturbation cells in horizontal direction
    
    n_slabs_k = kde/ngc_v
    
    pblh_def  = 500.0 !1000.0  !default pbl height if not using m_pblh_opt
    
    ! LOOP VER ALL VERTICAL PERTURBATION SLABS TO SEE IF IT IS TIME TO START OR REFRESH
    
    DO slab_k = 1, n_slabs_k     !Outer k-loop over number of vertical perturbation slabs
       
       prttms(slab_k) = prttms(slab_k) + dt
       
       ! print*,slab_k, prttms(slab_k), prtdt(slab_k)

       
       IF ( prttms(slab_k) .GE. prtdt(slab_k) ) THEN !prdt should be zero at startup
          
          prttms(slab_k) = dt

          prtdt(slab_k) = 300.0 !CHECK EVERY 10 MINUTES IF PBLH HAS GROWN
          
          print*,'Computing new perturbations, slab_k = ',slab_k
          
          ! ONLY COMPUTE/APPLY PERTURBATIONS FOR SLABS WITH A GLOBAL Z-VALUE BELOW GLOBAL PBLH MAXIMUM
          
          i_start = its
          i_end   = MIN(ite,ide)
          j_start = jts
          j_end   = MIN(jte,jde)
          
          IF (m_pblh_opt .EQ. 0 ) THEN !SET PBLH TO A DESIRED VALUE
             
             DO j=j_start, j_end
                DO i=i_start, i_end
                   pblh(i,j) = pblh_def
                END DO
             END DO
             
          ENDIF
          
          IF (m_pblh_opt .EQ. 1 ) THEN !SET PBLH TO A DESIRED VALUE
             
             DO j=j_start, j_end
                DO i=i_start, i_end
                   pblh(i,j) = mpblh(i,j)
                END DO
             END DO
             
          ENDIF
          
          ! COMPUTE THE LOCAL HEIGHT ABOVE THE SURFACE AT EACH GRID POINT
          
          DO j=j_start, j_end
             DO i=i_start, i_end
                z(i,kts,j)= 1.0/rdz(i,kts,j)
                DO k=kts+1,kde-1  
                   z(i,k,j) = z(i,k-1,j) + 1.0/rdz(i,k,j)
                END DO
             END DO
          END DO
          
          k_slab_start = MIN( (slab_k - 1)*ngc_v + 1,kde-1 )         !LOWEST K-VALUE IN SLAB
          k_slab_end = MIN( k_slab_start + ngc_v - 1, kde-1 )  !HIGHEST K-VALUE IN SLAB    
          
          max_pblh = MAXVAL( pblh(i_start:i_end,j_start:j_end) )        
          min_slab_z = MINVAL(z(i_start:i_end,k_slab_start,j_start:j_end))
          
          !print*,'Here 1 '
          !print*,'k_slab_start,k_slab_end ',k_slab_start,k_slab_end
          !print*,'max_pblh,min_slab_z ',max_pblh,min_slab_z

          idum = 1
          jdum = 1
          
#ifdef DM_PARALLEL
          
          CALL wrf_dm_maxval_real(max_pblh,idum,jdum)
          
          CALL wrf_dm_minval_real(min_slab_z,idum,jdum)
          
#endif
          
          !print*,'Here 2 '
          !print*,'max_pblh,min_slab_z ',max_pblh,min_slab_z


          max_pert_z = max_pblh*0.666666
!          max_pert_z = 200.0
!          
           IF ( min_slab_z .LE. max_pert_z ) THEN  
!          IF ( min_slab_z .LE. max_pblh ) THEN  !IF ANY K VALUE IN THE SLAB IS BELOW THE MAXIMUM PBL HEIGHT
             !THEN THAT GRID POINT MAY REQUIRE A PERTURBATION
           
            
             ! NOW THAT WE HAVE PBLH, CALCULATE UG, VG DOMAIN AVERAGE
             
             h_geo = 1.10*max_pblh

             ugeosum = 0.0
             vgeosum = 0.0
             uslabsum = 0.0
             vslabsum = 0.0
             num_pts_in_sum = 0.0
             
             
             ! print*,'a,h_geo ',h_geo 
             
             DO j=j_start, j_end
                
                DO i=i_start, i_end
                   
                   DO k=kts,kde-2
                      
                      IF ( ( z(i,k,j) .LE. h_geo ) .AND. ( z(i,k+1,j) .GT. h_geo ) )  k_geo = k
                      
                   END DO
                   
                   ugeosum = ugeosum + u(i,k_geo,j)
                   vgeosum = vgeosum + v(i,k_geo,j)
                   uslabsum = uslabsum + u(i,k_slab_end,j)
                   vslabsum = vslabsum + v(i,k_slab_end,j)
                   
                   num_pts_in_sum = num_pts_in_sum + 1.0
                   
                END DO
                
             END DO
             
             !print*,' Here 3',k_geo,num_pts_in_sum 
             
             
             !print*,'ugeosum ',ugeosum
             !print*,'vgeosum ',vgeosum
             !print*,'uslabsum ',uslabsum
             !print*,'vslabsum ',vslabsum
             
             
!GOTO 200
             
#ifdef DM_PARALLEL
             
             ugeolbsum  = wrf_dm_sum_real(ugeosum)
             
             vgeolbsum  = wrf_dm_sum_real(vgeosum)
             
             uslablbsum  = wrf_dm_sum_real(uslabsum)
             
             vslablbsum  = wrf_dm_sum_real(vslabsum)
             
             counter = wrf_dm_sum_real(num_pts_in_sum)
             
#endif
             
             !print*,' Here 4 counter',counter
             
             ugeolbavg = ugeolbsum/counter
             
             vgeolbavg = vgeolbsum/counter
             
             uslablbavg = uslablbsum/counter
             
             vslablbavg = vslablbsum/counter
             
             
             wsgeolbavg = sqrt( ugeolbavg*ugeolbavg + vgeolbavg*vgeolbavg )
             
             wsslablbavg = sqrt( uslablbavg*uslablbavg + vslablbavg*vslablbavg )
             
             
             anglelbavg = atan( abs( uslablbavg )/abs( vslablbavg ) )
             
             IF (anglelbavg .GT. pio4 ) anglelbavg = pio2 - anglelbavg 
             
             !print*,'ugeolbavg ',ugeolbavg
             !print*,'vgeolbavg ',vgeolbavg
             !print*,'uslablbsum ',uslablbavg
             !print*,'vslablbavg ',vslablbavg
             
             !print*,'anglelbavg ',anglelbavg
             
             north = 0
             east = 0
             south = 0
             west = 0
             
             IF ( ugeolbsum .GT. 0.0 ) west = 1
             IF ( ugeolbsum .LT. 0.0 ) east = 1
             IF ( vgeolbsum .GT. 0.0 ) south = 1
             IF ( vgeolbsum .LT. 0.0 ) north = 1
             
             ! print*,'uwsgpsumdsum,vwsgpsumdsum',uwsgpsumdsum,vwsgpsumdsum
             ! print*,'n, s, e, w',north,south,east,west     
             
             !
             ! Uncomment to select your own inflow/outflow boundaries
             !
!              north = 0
!              east = 0
!              south = 0
!              west = 1
             
             ! print*,'Applying perturbations to the following boundary(ies):'
             ! IF (north .EQ. 1) print*,'north'
             ! IF (south .EQ. 1) print*,'south'
             ! IF (east .EQ. 1) print*,'east'
             ! IF (west .EQ. 1) print*,'west'
             
             nsew = 0
             IF  ( (west + east) .EQ. 2) nsew = 1
             IF  ( (north + south) .EQ. 2) nsew = 1
             IF  ( (north + south + east + west) .EQ. 1) nsew = 1
             IF  ( (north + south + east + west) .EQ. 4) nsew = 1
             IF  ( (north + south + east + west) .EQ. 3) nsew = 0
             
             IF (wsslablbavg .EQ. 0.0) THEN !Prevent division by zero
                
                print*,'something wrong in stoch_les_inflow_perts'
                STOP
                
             ENDIF
             
             tpertmag = (wsgeolbavg*wsgeolbavg)/(ek_opt*cp)
!             tpertmag = 0.25

             
            ! print*,'tpertmag',tpertmag
             
             ! prtdt = 0.5*(1.0/cos(anglespavgdavg))*ngc_h*ncells_h*dx/wsspavgdavg 
             
             prtdt(slab_k) = (1.0/cos(anglelbavg))*ngc_h*ncells_h*dx/wsslablbavg

!JDM             prtdt(slab_k) = 100.0
             
             ! print*,'anglespavgdavg',anglespavgdavg*180.0/piconst
             
             !print*,'prtdt',prtdt(slab_k)

             
             ! CALCULATE (NEW) PERTURBATIONS ===============================================================================
             
             ni = (ide-1)/ngc_h 
             nj = (jde-1)/ngc_h  
             
             ! print*,'ni,nj',ni,nj
             
             ! Allocate perturbation array for each lateral boundary 
             
             ALLOCATE( pxs(1:nj,1:ncells_h) )
             ALLOCATE( pxe(1:nj,1:ncells_h) )
             ALLOCATE( pys(1:ni,1:ncells_h) )
             ALLOCATE( pye(1:ni,1:ncells_h) )
             
             CALL RANDOM_SEED(size=i_seed)  ! 1. Get the size of the seed, each processor.
             
             ALLOCATE( seed(1:i_seed) )     ! 2. Allocate an array to hold the seed, each processor.
             
#ifdef DM_PARALLEL
             
             IF ( mytask .EQ. master ) THEN
                
                seedsum = 0                     ! 3a. Sum the stored seed (up to i_seed) to see if it has been assigned.
                DO k = 1,i_seed             !     If so, now master has the current prtseed value in seed.
                   seed(k) = prtseed(k)
                   seedsum = seedsum + seed(k)
                END DO
                
                IF (seedsum .EQ. 0 ) THEN       !3b. If the seed is empty, get the first seed, which is now only on master.
                   
                   print*,'calling random seed for first time'   
                   CALL RANDOM_SEED(get=seed)
                   
                ENDIF
                
             ENDIF
             
             CALL MPI_BCAST(seed,i_seed,MPI_REAL,master,MPI_COMM_WORLD,ierr) ! Broadcast seed from master to all
             
             ! print*,'seed',seed
             
             CALL RANDOM_SEED(put=seed) !Use same seed value so all processors receive idential random numbers
             
             CALL RANDOM_NUMBER(pxs)
             CALL RANDOM_NUMBER(pxe)
             CALL RANDOM_NUMBER(pys)
             CALL RANDOM_NUMBER(pye)
             
             ! CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
             
             IF ( mytask .EQ. master ) THEN! Must get a new seed and store into prtseed so the sequence at the next perturbation update begins with new seed. 
                
                CALL RANDOM_SEED(get=seed)
                
                DO k = 1,i_seed
                   prtseed(k) = seed(k)
                ENDDO!
                
             ENDIF
             
#endif
             
             
             ! DO k = kts,kte  ! Zero out. Probably not necessary but good religion. 
             !    DO j=j_start, j_end
             !       DO i=i_start, i_end
             !          pert_t(i,k,j) = 0.0
             !       END DO
             !    END DO
             ! END DO
             
             sf = 1.0
             sf2 = 4.0*max_pblh/27.0 !maximum value of sf(z) = z*(1-z/h)^2 

             IF ( west .EQ. 1 ) THEN
                
                IF (its .LE. ids + ngc_h*ncells_h) THEN 
                   
                   DO j = MAX(jts,jds + ngc_h*ncells_h*(west - south - nsew)), MIN(jte,jde - 1 - ngc_h*ncells_h*(west - north - nsew))
                      
                      ! print*,'j,((j-1)/ngc_h+1)',j,((j-1)/ngc_h+1)
                      
                      DO i = its, MIN(ite, ide-1)
                         
                         ! print*,'i,((i-1)/ngc_h+1)',i,((i-1)/ngc_h+1)
                         
                         m = (i-1)/ngc_h+1
                         
                         IF ( m .LE. ncells_h ) THEN
                            
                            DO k = k_slab_start, k_slab_end

!                               sf = MAX( (z(i,k,j)*(1.0-z(i,k,j)/max_pblh)**2)/sf2, 0.0 )
!                               sf = MAX( cos(piconst*z(i,k,j)/(2*max_pblh))**2, 0.0 )
!                               sf = MAX( (1.0-z(i,k,j)/max_pblh)**2, 0.0 )

                               pert_t(i,k,j) =  (pxs(((j-1)/ngc_h+1),m)-0.5)*sf*2.0*tpertmag !Output variable only, for debugging.
                               
                               t(i,k,j) = t(i,k,j) + (pxs(((j-1)/ngc_h+1),m)-0.5)*sf*2.0*tpertmag
                               
                            END DO ! k
                            
                         ENDIF
                         
                      END DO ! i
                      
                   END DO! j
                   
                ENDIF
                
             ENDIF
             
             IF ( east .EQ. 1 ) THEN
                
                IF (ite .GE. ide - 1 - ngc_h*ncells_h) THEN 
                   
                   DO j = MAX(jts,jds + ngc_h*ncells_h*(east - south - nsew)), MIN(jte,jde - 1 - ngc_h*ncells_h*(east - north - nsew))
                      
                      DO i = MIN(ite, ide-1), its, -1
                         
                         m = (ide-i-1)/ngc_h+1
                         
                         IF ( m .LE. ncells_h ) THEN
                            
                            DO k = k_slab_start, k_slab_end

!                               sf = MAX( (z(i,k,j)*(1.0-z(i,k,j)/max_pblh)**2)/sf2, 0.0 )
!                               sf = MAX( cos(piconst*z(i,k,j)/(2*max_pblh))**2, 0.0 )
!                               sf = MAX( (1.0-z(i,k,j)/max_pblh)**2, 0.0 )
                               
                               pert_t(i,k,j) = (pxe(((j-1)/ngc_h+1),m)-0.5)*sf*2.0*tpertmag
                               
                               t(i,k,j) = t(i,k,j) + (pxe(((j-1)/ngc_h+1),m)-0.5)*sf*2.0*tpertmag
                               
                            END DO ! k
                            
                         ENDIF
                         
                      END DO ! i
                      
                   END DO ! j
                   
                ENDIF
                
             ENDIF
             
             IF ( south .EQ. 1 ) THEN
                
                IF (jts .LE. jds + ngc_h*ncells_h) THEN 
                   
                   DO i = MAX(its,ids + ngc_h*ncells_h*(south - west - nsew)), MIN(ite,ide - 1 - ngc_h*ncells_h*(south - east - nsew))
                      
                      DO j = jts, MIN(jte, jde-1)
                         
                         m = (j-1)/ngc_h+1
                         
                         IF ( m .LE. ncells_h ) THEN
                            
                            DO k = k_slab_start, k_slab_end
                               
                               ! print*,'i,(i-1)/ngc_h+1,j,m,pys',i,(i-1)/ngc_h+1,j,m
                               ! print*,'k,n,pys',k,n,pys(((i-1)/ngc_h+1),n,m)

!                               sf = MAX( (z(i,k,j)*(1.0-z(i,k,j)/max_pblh)**2)/sf2, 0.0 )
!                                sf = MAX( cos(piconst*z(i,k,j)/(2*max_pblh))**2, 0.0 )
!                                sf = MAX( (1.0-z(i,k,j)/max_pblh)**2, 0.0 )
                                
                               pert_t(i,k,j) = (pys(((i-1)/ngc_h+1),m)-0.5)*sf*2.0*tpertmag
                               
                               t(i,k,j) = t(i,k,j) + (pys(((i-1)/ngc_h+1),m)-0.5)*sf*2.0*tpertmag
                               
                            END DO ! k
                            
                         ENDIF
                         
                      END DO ! i
                      
                   END DO ! j
                   
                ENDIF
                
             ENDIF
             
             IF ( north .EQ. 1 ) THEN
                
                IF (jte .GE. jde - 1 - ngc_h*ncells_h) THEN 
                   
                   DO i = MAX(its,ids + ngc_h*ncells_h*(north - west - nsew)), MIN(ite,ide - 1 - ngc_h*ncells_h*(north - east - nsew))
                      
                      DO j = MIN(jte, jde-1),jts,-1
                         
                         m = (jde-j-1)/ngc_h+1
                         
                         IF ( m .LE. ncells_h ) THEN
                            
                            DO k = k_slab_start, k_slab_end

!                               sf = MAX( (z(i,k,j)*(1.0-z(i,k,j)/max_pblh)**2)/sf2, 0.0 )
!                               sf = MAX( cos(piconst*z(i,k,j)/(2*max_pblh))**2, 0.0 )
!                               sf = MAX( (1.0-z(i,k,j)/max_pblh)**2, 0.0 )
                               
                               pert_t(i,k,j) = (pye(((i-1)/ngc_h+1),m)-0.5)*sf*2.0*tpertmag
                               
                               t(i,k,j) = t(i,k,j) + (pye(((i-1)/ngc_h+1),m)-0.5)*sf*2.0*tpertmag
                               
                            END DO ! k
                            
                         ENDIF
                         
                      END DO ! i
                      
                   END DO ! j
                   
                ENDIF
                
             ENDIF
             
             DEALLOCATE(seed)
             
             DEALLOCATE( pxs )
             DEALLOCATE( pxe )
             DEALLOCATE( pys )
             DEALLOCATE( pye )
             
200          CONTINUE
             

!             print*,'min_slab_z .GE. max_pblh ',slab_k
             
          END IF ! IF ( min_slab_z .GE. max_pblh ) 

!200          CONTINUE
          
!          print*,'prttms(slab_k) .GE. prtdt(slab_k) ',slab_k
          
       END IF !(  prttms(slab_k) .GE. prtdt(slab_k) ) THEN
              
!       print*,'Big outer k-loop  ',slab_k
   
    END DO !Big outer k-loop over number of vertical perturbation slabs
    
    !       DO k=kts,kte - 1  
    !!          sf(k) = (z(itilemid,k,jtilemid)*(1.0-z(itilemid,k,jtilemid)/prtz)**2)/sf2
    !          sf(k) = 1.0 !jdm TEST 09142017
    !          IF ( z(itilemid,k,jtilemid) .LE. prtz ) prtnk = k
    !          IF ( z(itilemid,k,jtilemid) .LE. h_wg ) k_wg = k
    !          IF ( z(itilemid,k,jtilemid) .GT. prtz ) sf(k) = 0.0
    !!          print*,'k,z(itilemid,k,jtilemid),sf(k)',k,z(itilemid,k,jtilemid),sf(k)
    !       END DO
    
   
!=================================================================================

END SUBROUTINE calc_pert_t

!=======================================================================



!=======================================================================
!=======================================================================
!=======================================================================

    END MODULE module_stoch_les_inflow_perts

!=======================================================================
!=======================================================================
