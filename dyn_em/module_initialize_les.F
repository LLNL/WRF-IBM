!IDEAL:MODEL_LAYER:INITIALIZATION
!

!  This MODULE holds the routines which are used to perform various initializations
!  for the individual domains.  

!  This MODULE CONTAINS the following routines:

!  initialize_field_test - 1. Set different fields to different constant
!                             values.  This is only a test.  If the correct
!                             domain is not found (based upon the "id")
!                             then a fatal error is issued.               

!-----------------------------------------------------------------------

MODULE module_initialize_ideal

   USE module_domain
   USE module_io_domain
   USE module_state_description
   USE module_model_constants
   USE module_bc
   USE module_timing
   USE module_configure
   USE module_init_utilities
! KAL used for ideal soil initialization
   USE module_date_time
! KAL this adds random perturbations
   USE module_ranarray
! KAL for the immersed boundary method 
   USE module_ibm
!DJW 140408 added USE module_soil_pre along with subroutines 
!           monthly_min_max and monthly_interp_to_date. module_soil_pre
!           is necessary since several variables are defined in it that
!           are needed for monthly_interp_to_date.
   USE module_soil_pre
#ifdef DM_PARALLEL
   USE module_dm
#endif


CONTAINS


!-------------------------------------------------------------------
! this is a wrapper for the solver-specific init_domain routines.
! Also dereferences the grid variables and passes them down as arguments.
! This is crucial, since the lower level routines may do message passing
! and this will get fouled up on machines that insist on passing down
! copies of assumed-shape arrays (by passing down as arguments, the 
! data are treated as assumed-size -- ie. f77 -- arrays and the copying
! business is avoided).  Fie on the F90 designers.  Fie and a pox.

   SUBROUTINE init_domain ( grid )

   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid 
   !  Local data.
   INTEGER :: idum1, idum2

   CALL set_scalar_indices_from_config ( head_grid%id , idum1, idum2 )

     CALL init_domain_rk( grid &
!
#include "actual_new_args.inc"
!
                        )

   END SUBROUTINE init_domain

!-------------------------------------------------------------------

   SUBROUTINE init_domain_rk ( grid &
!
# include "dummy_new_args.inc"
!
)
   IMPLICIT NONE

   !KAL added for reading NetCDF files of terrain
   INCLUDE "netcdf.inc"

   !  Input data.
   TYPE (domain), POINTER :: grid

# include "dummy_new_decl.inc"

   TYPE (grid_config_rec_type)              :: config_flags

   !  Local data
   INTEGER                             ::                       &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte, &
                                  i, j, k

   ! Local data

!   INTEGER, PARAMETER :: nl_max = 1000
!KAL change to read high resolution soundings
   INTEGER, PARAMETER :: nl_max = 27000
   REAL, DIMENSION(nl_max) :: zk, p_in, theta, rho, u, v, qv, pd_in
   INTEGER :: nl_in


   INTEGER :: icm,jcm, ii, im1, jj, jm1, loop, error, fid, nxc, nyc
   REAL    :: u_mean,v_mean, f0, p_surf, p_level, qvf, z_at_v, z_at_u
   REAL    :: z_scale, xrad, yrad, zrad, rad, delt, cof1, cof2
!   REAL, EXTERNAL :: interp_0
   REAL    :: hm, xa, xa1
   REAL    :: pi
!KAL added offset for terrain
   REAL    :: x_dist, y_dist, offset, &
              slope !RSA

!  stuff from original initialization that has been dropped from the Registry 
   REAL    :: vnu, xnu, xnus, dinit0, cbh, p0_temp, t0_temp, zd, zt
   REAL    :: qvf1, qvf2, pd_surf, theta_surf
   INTEGER :: it
   real :: thtmp, ptmp, temp(3)

   LOGICAL :: moisture_init
   LOGICAL :: stretch_grid, dry_sounding

   !KAL adding specified height iteration
   REAL                                 :: mub_ave     !mub for specified height vector
   REAL, DIMENSION(grid%sm32:grid%em32) :: z_specified, p_specified !z and p for specified heigts

   !KAL the next three lines are added for random seed generation 
   REAL, DIMENSION(grid%sd31:grid%ed31,grid%sd33:grid%ed33) :: useed
   INTEGER :: iseed
   REAL :: amp

  INTEGER :: xs , xe , ys , ye
  REAL :: mtn_ht
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
!  For LES, add randx
   real :: randx
   LOGICAL :: temp_pert, vel_pert !RSA flags for temperature/velocity perturbations


  !DJW
  INTEGER :: ns   ! number of soil levels for idealized terrain/soil initialization

!DJW added for specifying different eta levels for each domain
  INTEGER :: ks, ke, id
  LOGICAL :: vnest !DJW T if using vertical nesting, otherwise F
  LOGICAL :: all_fine
  CHARACTER*256 :: message

!KAL variables for reading NetCDF files of terrain  
  INTEGER                                         :: ncid, & ! identifier for netcdf file
                                                     vid,  & !                netcdf variable
                                                     status
  INTEGER                                         :: nx, ny,   &
                                                     itstart,  &
                                                     itend,    &
                                                     jtstart,  &
                                                     jtend
  REAL,DIMENSION(:,:),ALLOCATABLE                 :: terrain_ht
  !RSA Read Granite Mountain terrain and surface data from netcdf file
  INTEGER :: idstart, idend, jdstart, jdend
  REAL,DIMENSION(:,:),ALLOCATABLE :: landuse_read
  REAL,DIMENSION(:,:),ALLOCATABLE :: soiltype_read
  REAL,DIMENSION(:,:),ALLOCATABLE :: lat_read
  REAL,DIMENSION(:,:),ALLOCATABLE :: long_read

!KAL adding for reading the map variables from another netcdf file
  LOGICAL :: read_geog   
  REAL,DIMENSION(grid%sp31:grid%ep31,grid%sp33:grid%ep33) :: dummy
! KAL  variables for reading text files of terrain
!  LOGICAL					   :: end_of_file
!  REAL, DIMENSION(nl_max)                         :: terrain_ht
  INTEGER :: width, ista, iend !DJW added for ideal_terrain
  REAL :: height !DJW added for ideal_terrain

#ifdef DM_PARALLEL
#    include "data_calls.inc"
#endif


   SELECT CASE ( model_data_order )
         CASE ( DATA_ORDER_ZXY )
   kds = grid%sd31 ; kde = grid%ed31 ;
   ids = grid%sd32 ; ide = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   kms = grid%sm31 ; kme = grid%em31 ;
   ims = grid%sm32 ; ime = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   kts = grid%sp31 ; kte = grid%ep31 ;   ! note that tile is entire patch
   its = grid%sp32 ; ite = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XYZ )
   ids = grid%sd31 ; ide = grid%ed31 ;
   jds = grid%sd32 ; jde = grid%ed32 ;
   kds = grid%sd33 ; kde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   jms = grid%sm32 ; jme = grid%em32 ;
   kms = grid%sm33 ; kme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   jts = grid%sp32 ; jte = grid%ep32 ;   ! note that tile is entire patch
   kts = grid%sp33 ; kte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XZY )
   ids = grid%sd31 ; ide = grid%ed31 ;
   kds = grid%sd32 ; kde = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   kms = grid%sm32 ; kme = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   kts = grid%sp32 ; kte = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch

   END SELECT


!  stretch_grid = .true.
!  FOR LES, set stretch to false
   stretch_grid = .false.
   temp_pert = .false. !RSA turn on/off temperature perturbations
   ! vel_pert = .false. !RSA turn on/off velocity perturbations !RSA changed to namelist argument
   delt = 3.
!   z_scale = .50
   z_scale = .40
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = (jde-jds)/2

   CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )

! here we check to see if the boundary conditions are set properly

   CALL boundary_condition_check( config_flags, bdyzone, error, grid%id )

   moisture_init = .true.

    grid%itimestep=0

#ifdef DM_PARALLEL
   CALL wrf_dm_bcast_bytes( icm , IWORDSIZE )
   CALL wrf_dm_bcast_bytes( jcm , IWORDSIZE )
#endif

    CALL nl_set_mminlu(1, '    ')
    CALL nl_set_iswater(1,0)
!KAL change the central lat and lon
! original code here    
!    CALL nl_set_cen_lat(1,40.)
!    CALL nl_set_cen_lon(1,-105.)
!    CALL nl_set_truelat1(1,0.)
!    CALL nl_set_truelat2(1,0.)
!    CALL nl_set_moad_cen_lat (1,0.)
!    CALL nl_set_stand_lon (1,0.)
!    CALL nl_set_pole_lon (1,0.)
!    CALL nl_set_pole_lat (1,90.)
!    CALL nl_set_map_proj(1,0)

!KAL changed code here
! this is to match the neutral boundary layer simulation
! that tina ran for JAS publication
!    CALL nl_set_cen_lat(1,43.29)
!    CALL nl_set_cen_lon(1,-100.)
!    CALL nl_set_truelat1(1,0.)
!    CALL nl_set_truelat2(1,0.)
!    CALL nl_set_moad_cen_lat (1,0.)
!    CALL nl_set_stand_lon (1,0.)
!    CALL nl_set_pole_lon (1,0.)
!    CALL nl_set_pole_lat (1,90.)
!    CALL nl_set_map_proj(1,0)

!KAL for Juergs intercomparison case
    CALL nl_set_cen_lat(1,36.)
    CALL nl_set_cen_lon(1,0.)
    CALL nl_set_truelat1(1,36.)
    CALL nl_set_truelat2(1,36.)
    CALL nl_set_moad_cen_lat (1,36.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_map_proj(1,0)
    !KAL add the following
    CALL nl_set_gmt(1,6.)
    CALL nl_set_julyr(1,2007)
    CALL nl_set_julday(1,80) 

!KAL for OKC simulation
!   CALL nl_set_cen_lat(1,35.5)
!   CALL nl_set_cen_lon(1,97.5)
!   CALL nl_set_truelat1(1,35.5)
!   CALL nl_set_truelat2(1,35.5)
!   CALL nl_set_moad_cen_lat (1,35.5)
!   CALL nl_set_stand_lon (1,97.5)
!   CALL nl_set_pole_lon (1,0.)
!   CALL nl_set_pole_lat (1,90.)
!   CALL nl_set_map_proj(1,0)


!  here we initialize data we currently is not initialized 
!  in the input data

    DO j = jts, jte
      DO i = its, ite
         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%msfvy(i,j)    = 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.
         grid%e(i,j)        = 0.
!  for LES, include Coriolis force
!         grid%f(i,j)        = 1.e-4 
!KAL
! do not add a coriolis force unless there is geostrophic flow	 
!RSA added if pert_coriolis and check for 2d case ideal_terrain = 21)
         IF ( ((config_flags%u_geo .NE. 0.0 ) .OR. (config_flags%v_geo .NE. 0.0) &
                                  .OR. (config_flags%pert_coriolis .EQ. .TRUE.)) &
              .AND. (config_flags%ideal_terrain .NE. 21)) THEN
            !KAL for OKC simulation
            !grid%f(i,j)        = 0.0000844598
            !KAL for other simulations
            ! grid%f(i,j)        = 0.0001000
            !RSA for flat plate case
            grid%f(i,j)        = 1.e-4
            if ((i .eq. its) .and. (j .eq. jts)) write(*,*) "RSA coriolis is on, f=", grid%f(i,j)
         ELSE
            grid%f(i,j)        = 0.0
            if ((i .eq. its) .and. (j .eq. jts)) write(*,*) "RSA coriolis is off, f=", grid%f(i,j)
         ENDIF
!KAL
      END DO
   END DO

    DO j = jts, jte
    DO k = kts, kte
      DO i = its, ite
         grid%ww(i,k,j)     = 0.
      END DO
   END DO
   END DO

   grid%step_number = 0

! set up the grid

   !DJW Added code for specifying multiple domains' eta_levels.
   !First check to make sure that we've not specified more
   !eta_levels than the dimensionality of eta_levels can handle! This
   !issue will most likely cause a break sometime before we real this
   !check, however it doesn't hurt to include it. To increase max_eta,
   !go to frame/module_driver_constants.F.
   vnest = .FALSE.
   DO id=1,model_config_rec%max_dom
     IF (model_config_rec%vert_refine_method(id) .NE. 0) THEN
       vnest = .TRUE.
     ENDIF
   ENDDO

   IF (model_config_rec%eta_levels(1) .EQ. -1) THEN !we do not have eta_levels from namelist      
     !DJW start of original code to set eta levels
     IF (stretch_grid) THEN ! exponential stretch for eta (nearly constant dz)
       CALL wrf_debug(0, "module_initialize_les: eta_levels is not specified in the namelist, setting levels with stretched spacing in eta.")
       DO k=1, kde
         grid%znw(k) = (exp(-(k-1)/float(kde-1)/z_scale) - exp(-1./z_scale))/ &
                       (1.-exp(-1./z_scale))
       ENDDO
     ELSE
       CALL wrf_debug(0,"module_initialize_les: eta_levels is not specified in the namelist, setting levels with constant spacing in eta.")
       DO k=1, kde
         grid%znw(k) = 1. - float(k-1)/float(kde-1)
       ENDDO
     ENDIF
   ELSE !we have specified eta levels from the namelist
     CALL wrf_debug(0,"module_initialize_les: vertical nesting is enabled, using eta_levels specified in namelist.input")
     i = 1
     ke = model_config_rec%e_vert(1)
     DO WHILE (grid%id .GT. id)
       id = id+1
       ke = ks+model_config_rec%e_vert(id)
     ENDDO
     IF (ke .GT. max_eta) THEN
       CALL wrf_error_fatal("too many vertical levels, increase max_eta in frame/module_driver_constants.F and recompile")
     ENDIF
     !Now set the eta_levels to what we specified in the namelist. We've
     !packed all the domains' eta_levels into a 'vector' and now we need
     !to pull only the section of the vector associated with our domain
     !of interest, which is between indicies ks and ks+e_vert(id)-1.
     ks = 1
     IF (grid%id .NE. 1) THEN
       id = 1
       DO WHILE (grid%id .GT. id)
         id = id+1
         ks = ks+model_config_rec%e_vert(id-1)
       ENDDO
     ENDIF
     DO k=1,kde
       grid%znw(k) = model_config_rec%eta_levels(ks+k-1)
       write(message,'(A,I3,A,F6.3)') "eta_levels(",k,") = ",grid%znw(k)
       CALL wrf_debug(0, TRIM(message))
     ENDDO
     !Check the value of the first and last eta level for our domain,
     !then check that the vector of eta levels is only decreasing
     all_fine = .TRUE.
     IF (grid%znw(1) .NE. 1.0) THEN
       write(message,'(A)') "error with specified eta_levels, first level is not 1.0"
       all_fine = .FALSE.
     ENDIF
     IF (grid%znw(kde) .NE. 0.0) THEN
       write(message,'(A)') "error with specified eta_levels, last level is not 0.0"
       all_fine = .FALSE.
     ENDIF
     DO k=2,kde
       IF (grid%znw(k) .GT. grid%znw(k-1)) THEN
         write(message,'(A)') "eta_levels are not uniformly decreasing from 1.0 to 0.0"
         all_fine = .FALSE.
       ENDIF
     ENDDO
     IF (.NOT. all_fine) THEN
       ks = 1
       ke = ks+model_config_rec%e_vert(1)
       write(*,'(A,I3,A,I3,A)') "eta_levels(",ks,":",ke,") = 1.0, ..., 0.0,"
       id = 1
       DO WHILE (model_config_rec%max_dom .GT. id)
         id = id+1
         ks = ks+model_config_rec%e_vert(id-1)
         ke = ks+model_config_rec%e_vert(id)-1
         write(*,'(A,I3,A,I3,A)') "eta_levels(",ks,":",ke,") = 1.0, ..., 0.0,"
       ENDDO
       CALL wrf_error_fatal(TRIM(message))
     ENDIF
   ENDIF

! **** IMPORTANT - if you are using this your input sounding needs to be at a finer resolution than your specified grid ****
! **** there are unexplained errors if it is not, and the grid can be unbalanced and move *****

!DJW previously hardcoded sets of vertical levels have been removed since they can now be
!    specified in the namelist. Check out an old commit if you need to figure out what
!    levels were hardcoded in previous simulations.

   IF ( (model_config_rec%z_specified(1) .NE. -1.0) .AND. &
        (model_config_rec%z_specified(2) .NE. -1.0) ) THEN
      !DJW find the number of z_specified entries that will precede the
      !    section associated with our domain of interest.
      ks = 1
      i = 1
      DO WHILE (i .LE. model_config_rec%max_dom)
         IF (model_config_rec%grid_id(i) .EQ. grid%grid_id) THEN
            i = model_config_rec%max_dom+1
         ELSE
            ks = ks+model_config_rec%e_vert(i)
            i = i+1
         ENDIF
      ENDDO

      DO k=ks,ks+kde-1
         z_specified(k-ks+1) = model_config_rec%z_specified(k)
      ENDDO

      DO k=1,kde
         write(wrf_err_message,'(A,I1,A,I3,A,F8.2)') "[module_initialize_les]: grid_id=",grid%grid_id," z_specified(",k,")=",z_specified(k)
         CALL wrf_message(wrf_err_message)
      ENDDO

      !Check that the levels specified in the namelist are uniformly increasing
      DO k=1,kde
         IF (k .NE. 1) THEN
            IF (z_specified(k) .LE. z_specified(k-1)) THEN
               write(wrf_err_message,'(2(A,I3),A,I1)') "--- ERROR: z_specified(",k,") <= z_specified(",k-1,") for grid_id=",grid%grid_id
               CALL wrf_error_fatal( wrf_err_message )
            ENDIF
         ENDIF
      ENDDO

      !Check that the top-most specified height matches the requested ztop
      IF (z_specified(kme) .NE. config_flags%ztop) THEN
         write(wrf_err_message,'(2(A,F8.2),A)') "---ERROR: ztop (",config_flags%ztop,") is not consistent with z_specified (",z_specified(kme),")"
         CALL wrf_error_fatal ( wrf_err_message )
      ENDIF
      !KAL get the sounding so that density can be calulated at the specified heights
      ! mu is also needed get sounding is called again later (in the original code).
      !This is probably not needed, but I am leaving the original code there.
      dry_sounding = .TRUE.
      CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, nl_max, nl_in, theta_surf )

      !KAL find the pressure at the top for calculating mu 		      
      !  find ptop for the desired ztop (ztop is input from the namelist),
      !  and find surface pressure
      grid%p_top = interp_0( p_in, zk, config_flags%ztop, nl_in ) 
      p_surf = interp_0( p_in, zk, z_specified(1), nl_in )
      mub_ave = p_surf-grid%p_top

      DO k=1,kde 
         p_specified(k) = interp_0( p_in, zk, z_specified(k), nl_in )
      ENDDO

      grid%znw(1) = 1.0
      DO k=2,kde-1  
         grid%znw(k) = (p_specified(k)-grid%p_top)/mub_ave
      ENDDO  
      grid%znw(kde) = 0.0
      DO k=1,kde
         write(*,'(A,I1,A,I3,A,F5.3)') "DJW[module_initialize_les]: grid%id=",grid%id,"   znw(",k,") = ",grid%znw(k)
      ENDDO
   ENDIF !z_specified

!  dnw is delta eta on full w levels
!  rdnw is inverse delta eta on full w levels
!  znu is eta values on half levels
   DO k=1, kde-1
    grid%dnw(k) = grid%znw(k+1) - grid%znw(k)
    grid%rdnw(k) = 1./grid%dnw(k)
    grid%znu(k) = 0.5*(grid%znw(k+1)+grid%znw(k))
   ENDDO
!  dn is delta eta on half levels
!  fnp is the lower weight for vertical stretching
!  fnm is the upper weight for verical stretching   
   DO k=2, kde-1
    grid%dn(k) = 0.5*(grid%dnw(k)+grid%dnw(k-1))
    grid%rdn(k) = 1./grid%dn(k)
    grid%fnp(k) = .5* grid%dnw(k  )/grid%dn(k)
    grid%fnm(k) = .5* grid%dnw(k-1)/grid%dn(k)
   ENDDO

   cof1 = (2.*grid%dn(2)+grid%dn(3))/(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(2) 
   cof2 =     grid%dn(2)        /(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(3) 
   grid%cf1  = grid%fnp(2) + cof1
   grid%cf2  = grid%fnm(2) - cof1 - cof2
   grid%cf3  = cof2       

   grid%cfn  = (.5*grid%dnw(kde-1)+grid%dn(kde-1))/grid%dn(kde-1)
   grid%cfn1 = -.5*grid%dnw(kde-1)/grid%dn(kde-1)
   grid%rdx = 1./config_flags%dx
   grid%rdy = 1./config_flags%dy

!  get the sounding from the ascii sounding file, first get dry sounding and 
!  calculate base state

  dry_sounding = .true.
  IF ( wrf_dm_on_monitor() ) THEN
  write(6,*) ' getting dry sounding for base state '

  CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, nl_max, nl_in, theta_surf )
  ENDIF
  CALL wrf_dm_bcast_real( zk , nl_max )
  CALL wrf_dm_bcast_real( p_in , nl_max )
  CALL wrf_dm_bcast_real( pd_in , nl_max )
  CALL wrf_dm_bcast_real( theta , nl_max )
  CALL wrf_dm_bcast_real( rho , nl_max )
  CALL wrf_dm_bcast_real( u , nl_max )
  CALL wrf_dm_bcast_real( v , nl_max )
  CALL wrf_dm_bcast_real( qv , nl_max )
  CALL wrf_dm_bcast_integer ( nl_in , 1 ) 

  write(6,*) ' returned from reading sounding, nl_in is ',nl_in

!  find ptop for the desired ztop (ztop is input from the namelist),
!  and find surface pressure

  grid%p_top = interp_0( p_in, zk, config_flags%ztop, nl_in )

  IF ((config_flags%ibm_opt .GT. 0) .and. (config_flags%vrm_wrf .EQ. 0)) THEN
      write(*,'(A,I1,A,I3)') "Calling ibm_terrain for grid%id=",grid%id," with ideal_terrain=",config_flags%ideal_terrain
      !KAL flat surface used for ibm
      DO j=jts,jte
      DO i=its,ite
         IF ((i .NE. ide) .AND. (j .NE. jde)) grid%ht(i,j) = config_flags%ideal_terrain_ht
      ENDDO
      ENDDO
      CALL ibm_terrain ( config_flags,                 &
                         grid%ibm_ht_u, grid%ibm_ht_v, &
                         grid%ibm_ht_w, grid%ibm_ht_c, &
                         grid%ibm_z0,                  &
                         ids, ide, jds, jde, kds, kde, &
                         ims, ime, jms, jme, kms, kme, &
                         its, ite, jts, jte, kts, kte )
  ELSE !define ht   
      !KAL adding ideal_terrain options
      ! ideal_terrain 0 = none
      ! 1 = X m flat plate
      ! 2 = 2d witch of agnesi hill in x
      ! 21 = 2d gaussian hill in x
      ! 3 = 3d witch of agnesi hill
      ! 31 = 3d gaussian hill
      ! 4 = 2d ridge or isolated cube (4,5 are not included below because they will not work for terrain following coordinates)
      ! 5 = 3d cube
      ! 6 = 2d idealized t-rex valley
      ! 7 = 3d idealized t-rex Valley
      ! 8 = v Shaped valley and ridges
      ! 81 = Sinusoidal valley and ridges RSA
      ! 9 = Schaer idealized advection test
      ! 10 = read data from file
      ! 11 = 2d witch of agnesi hill in y
      ! 12 = Alternate 5 hill terrain Schaer/Klemp/Shaw
      ! 15 = RSA/JB Askervein read from file
     
      IF (config_flags%ideal_terrain .EQ. 1) THEN
          DO j=jts,jte
          DO i=its,ite
              IF ((i .NE. ide) .AND. (j .NE. jde)) grid%ht(i,j) = config_flags%ideal_terrain_ht !RSA
          ENDDO
          ENDDO

      ELSEIF (config_flags%ideal_terrain .EQ. 2) THEN
          !hm = 100.0
          !xa = 250.0
          !offset = 0.0
          !RSA for Lundquist et al. (2012) case
          hm = 350.0
          xa = 800.0
          offset = 0.0
          DO j=jts,jte
          DO i=its,ite
             IF ((i .NE. ide) .AND. (j .NE. jde)) THEN
               ! hill in the x direction
               ! grid%ht(i,j) = offset + hm/(1.+(float(i-icm)/xa)**2) 
               ! hill in the y direction
               ! grid%ht(i,j) = offset + hm/(1.+(float(j-jcm)/xa)**2)
               !hill in the x direction
               x_dist = (-ide/2.)*config_flags%dx + i*config_flags%dx
               grid%ht(i,j) = offset + hm/(1.0+(x_dist/xa)**2)
               !hill in the y direction
               ! y_dist = (-jde/2.)*config_flags%dy + j*config_flags%dy
               ! grid%ht(i,j) = offset + hm/(1.0+(y_dist/xa)**2)
             ENDIF
          ENDDO
          ENDDO

      ELSEIF (config_flags%ideal_terrain .EQ. 21) THEN
          !RSA 2d gaussian hill to roughly match Askervein transect A
          hm = 100.0
          xa = 150.0
          offset = 0.0
          DO j=jts,jte
          DO i=its,ite
             IF ((i .NE. ide) .AND. (j .NE. jde)) THEN
               x_dist = (-ide/2)*config_flags%dx + i*config_flags%dx
               grid%ht(i,j) = offset + hm*exp(-0.5*(x_dist/xa)**2)
             ENDIF
          ENDDO
          ENDDO

      ELSEIF (config_flags%ideal_terrain .EQ. 11) THEN
          hm = 1000.0
          xa = 40.0
          offset = 13.0
          DO j=jts,jte
          DO i=its,ite
             IF ((i .NE. ide) .AND. (j .NE. jde)) THEN
               ! hill in the y direction
               ! grid%ht(i,j) = offset + hm/(1.+(float(j-jcm)/xa)**2)
               x_dist = (-jde/2)*config_flags%dy + j*config_flags%dy
               grid%ht(i,j) = offset + hm/(1.0+(x_dist/xa)**2)
             ENDIF
          ENDDO
          ENDDO

      ELSEIF (config_flags%ideal_terrain .EQ. 3) THEN
         !RSA for Lundquist et al. (2012) case
         hm = 350.0
         xa = 800.0
         offset = 0.0
         DO j=jts,jte
         DO i=its,ite
            IF ((i .NE. ide) .AND. (j .NE. jde)) THEN
                ! use this for gaussian bump
                ! grid%ht(i,j) = hm/(1.+(float(i-icm)/xa)**2+(float(j-jcm)/xa)**2)
                 x_dist = (-ide/2.)*config_flags%dx + i*config_flags%dx
                 y_dist = (-jde/2.)*config_flags%dy + j*config_flags%dy
                 grid%ht(i,j) = offset + hm/(1.+(x_dist/xa)**2+(y_dist/xa)**2)
            ENDIF
         ENDDO
         ENDDO

      ELSEIF (config_flags%ideal_terrain .EQ. 31) THEN
         !RSA 3d gaussian hill roughly the scale of Askervein
         hm = 100.0
         xa = 150.0
         offset = 0.0
         DO j=jts,jte
         DO i=its,ite
            IF ((i .NE. ide) .AND. (j .NE. jde)) THEN
                ! use this for gaussian bump
                ! grid%ht(i,j) = hm/(1.+(float(i-icm)/xa)**2+(float(j-jcm)/xa)**2)
                 x_dist = (-ide/2.)*config_flags%dx + i*config_flags%dx
                 y_dist = (-jde/2.)*config_flags%dy + j*config_flags%dy
                 grid%ht(i,j) = offset + hm*exp(-0.5*((x_dist/xa)**2+(y_dist/xa)**2))
            ENDIF
         ENDDO
         ENDDO

      ELSEIF (config_flags%ideal_terrain .EQ. 4) THEN
         grid%ht(i,j)=0.0

      ELSEIF (config_flags%ideal_terrain .EQ. 5) THEN
         grid%ht(i,j)=0.0

      ELSEIF (config_flags%ideal_terrain .EQ. 6) THEN
         hm = 1500.0
         DO j=jts,jte-1
         DO i=its,ite-1  
            !calculate x with x=0 as the middle of the domain
            x_dist = (-ide/2)*config_flags%dx + i*config_flags%dx
            IF (ABS(x_dist)<=500.0) then
               grid%ht(i,j) = 0.0
            ELSEIF ((ABS(x_dist) .GT. 500.0) .AND. (ABS(x_dist) .LT. 9500.0)) THEN
               grid%ht(i,j) = hm*(0.5 - 0.5*COS((pi*(ABS(x_dist)-500.0))/9000.0))
            ELSEIF ((ABS(x_dist) .GE. 9500.0) .AND. (ABS(x_dist) .LE. 10500.0)) THEN
               grid%ht(i,j) = hm
            ELSEIF ((ABS(x_dist) .GT. 10500.0) .AND. (ABS(x_dist) .LT. 19500.0)) THEN
               grid%ht(i,j) = hm*(0.5 + 0.5*COS((pi*(ABS(x_dist)-10500.0))/9000.0))
            ELSEIF (ABS(x_dist) .GE. 19500.0) THEN
               grid%ht(i,j) = 0.0
            ENDIF
         ENDDO
         ENDDO

      ELSEIF (config_flags%ideal_terrain .EQ. 7) THEN
         hm = 1500.0
         DO j=jts,jte-1
         DO i=its,ite-1
            !calculate x with x=0 as the middle of the domain
            x_dist = (-ide/2)*config_flags%dx + i*config_flags%dx
            IF (ABS(x_dist) .LE. 500.0) THEN
               grid%ht(i,j) = 0.0
            ELSEIF ((ABS(x_dist) .GT. 500.0) .AND. (ABS(x_dist) .LT. 9500.0)) THEN
               grid%ht(i,j) = hm*(0.5 - 0.5*COS((pi*(ABS(x_dist)-500.0))/9000.0))
            ELSEIF ((ABS(x_dist) .GE. 9500.0) .AND. (ABS(x_dist) .LE. 10500.0)) THEN
               grid%ht(i,j) = hm
            ELSEIF ((ABS(x_dist) .GT. 10500.0) .AND. (ABS(x_dist) .LT. 19500.0)) THEN
               grid%ht(i,j) = hm*(0.5 + 0.5*COS((pi*(ABS(x_dist)-10500.0))/9000.0))
            ELSEIF (ABS(x_dist) .GE. 19500.0) THEN
               grid%ht(i,j) = 0.0
            ENDIF
            !calculate the y distance with y=0 as the middle of the domain
            y_dist = (-jde/2)*config_flags%dy + j*config_flags%dy
            IF (ABS(y_dist) .LE. 100000.0) THEN
               grid%ht(i,j) = grid%ht(i,j)
            ELSEIF ((ABS(y_dist) .GT. 100000.0) .AND. (ABS(y_dist) .LT. 109000.0)) THEN
               grid%ht(i,j) = grid%ht(i,j)*(0.5 + 0.5*COS((pi*(ABS(y_dist)-100000.0))/9000.0))
            ELSEIF (ABS(y_dist) .GE. 109000.0) THEN
               grid%ht(i,j) = 0.0
            ENDIF
         ENDDO
         ENDDO

      ELSEIF (config_flags%ideal_terrain .EQ. 8) THEN
         !RSA v-shaped valley case
         offset = 0.0
         slope = 0.2126 !RSA 12 degrees
         xa = 1000.0 !side ridges
         write(*,*) "ide=",ide,"ide/2=",ide/2
         write(*,*) "dx=",config_flags%dx
         DO j=jts,jte-1
         DO i=its,ite-1
            !RSA modified for intermediate resolution IBM cases
            !calculate x with x=0 as the middle of the domain
            x_dist = (-ide/2)*config_flags%dx - config_flags%dx/2. + i*config_flags%dx
            IF (ABS(x_dist) .LT. (ide/2)*config_flags%dx - xa) THEN
               grid%ht(i,j) = offset + slope*ABS(x_dist)
            ELSE
               grid%ht(i,j) = offset + slope*((ide/2)*config_flags%dx - xa)
            ENDIF
            if (j .eq. jts) write(*,*) "i=",i,"xdist=",x_dist,"ht=",grid%ht(i,j)
            ! x_dist = i*config_flags%dx - (config_flags%dx/2.)
            !grid%ht(i,j) = 3.*x_dist/100.
         ENDDO
         ENDDO

      ELSEIF (config_flags%ideal_terrain .EQ. 81) THEN
         !RSA sinusoidal valley case
         offset = 550.0
         xa = 1000.0 !side ridges
         write(*,*) "ide=",ide,"ide/2=",ide/2
         write(*,*) "dx=",config_flags%dx
         DO j=jts,jte-1
         DO i=its,ite-1
            !calculate x with x=0 as the middle of the domain
            x_dist = (-ide/2)*config_flags%dx - config_flags%dx/2. + i*config_flags%dx
            IF (ABS(x_dist) .LT. (ide/2)*config_flags%dx - xa) THEN
               grid%ht(i,j) = offset/2 - offset/2*COS(2*pi*ABS(x_dist)/((ide-1)*config_flags%dx-2*xa))
            ELSE
               grid%ht(i,j) = offset
            ENDIF
            if (j .eq. jts) write(*,*) "i=",i,"xdist=",x_dist,"ht=",grid%ht(i,j)
         ENDDO
         ENDDO

      ELSEIF (config_flags%ideal_terrain .EQ. 9) THEN
         ! constants from Schaer 2002
         !hm = 2000.
         hm = 3000.0
         xa = 25000.0
         xa1  = 8000.0
         DO j=jts,jte-1
         DO i=its,ite-1
            !calculate x with x=0 as the middle of the domain
            x_dist = (-ide/2)*config_flags%dx + i*config_flags%dx
            IF (ABS(x_dist) .LE. xa) THEN
               grid%ht(i,j) = hm*COS((pi*x_dist)/(2*xa))*COS((pi*x_dist)/(2*xa))*COS((pi*x_dist)/(xa1))*COS((pi*x_dist)/(xa1))
            ELSEIF (ABS(x_dist) .GT. xa) THEN
               grid%ht(i,j) = 0.0
            ENDIF
         ENDDO
         ENDDO

      ELSEIF (config_flags%ideal_terrain .EQ. 10) THEN
         !code for reading NetCDF files - terrain from a wrfinput file (right now Jan 2000 boulder windstorm)
         itend = MIN(ite,ide-1)
         jtend = MIN(jte,jde-1)
         status = NF_OPEN("wrfinput_HRRR_from_real", nf_nowrite, ncid)
         status = NF_INQ_VARID(ncid, "HGT", vid) !use for wrfinput data made from real	 
         ALLOCATE(terrain_ht(its:itend,jts:jtend)) 
         !the A in VARA stands for array
         IF (KIND(terrain_ht) .EQ. 4) THEN
            status = NF_GET_VARA_REAL(ncid, vid, [its, jts, 1], [itend-its+1,jtend-jts+1, 1], terrain_ht)
         ELSEIF (KIND(terrain_ht) .EQ. 8) THEN
            status = NF_GET_VARA_DOUBLE(ncid, vid, [its, jts, 1], [itend-its+1,jtend-jts+1, 1], terrain_ht)
         ENDIF
         status = NF_CLOSE(ncid)
         DO j=jts,jtend
         DO i=its,itend
            !grid%ht(i,j) = terrain_ht(i,j)
            grid%ht(i,j) = 0.0
         ENDDO
         ENDDO
         DEALLOCATE(terrain_ht)

      ELSEIF (config_flags%ideal_terrain .EQ. 12) THEN
         !constants from Schaer 2002, Klemp 2011 
         hm = 1000.0 ! 250 for Schaer and 1000 for Klemp
         xa = 5000.0
         xa1 = 4000.0
         DO j=jts,jte-1
         DO i=its,ite-1
            ! calculate x with x=0 as the middle of the domain
            x_dist = (-ide/2.0)*config_flags%dx + i*config_flags%dx
            grid%ht(i,j) = hm*EXP(-1.*(x_dist/xa)**2.)*COS((pi*x_dist)/(xa1))*COS((pi*x_dist)/(xa1))
         ENDDO
         ENDDO

      ELSEIF (config_flags%ideal_terrain .EQ. 15) THEN
         !RSA/JB for Askervein 
         status = NF_OPEN("Ask_hill.nc", nf_nowrite, ncid)
         status = NF_INQ_DIMID(ncid, 'nx', vid)
         status = NF_INQ_DIMLEN(ncid, vid, nx)
         status = NF_INQ_DIMID(ncid, 'ny', vid)
         status = NF_INQ_DIMLEN(ncid, vid, ny)
         status = NF_INQ_VARID(ncid, "HT", vid)
         ALLOCATE(terrain_ht(nx,ny))
         terrain_ht = 10.
         !the A in VARA stands for array
         if (KIND(terrain_ht) == 4) then
             status = NF_GET_VARA_REAL(ncid, vid, [1, 1], [nx,ny],  terrain_ht)
         elseif (KIND(terrain_ht) == 8) then
             status = NF_GET_VARA_DOUBLE(ncid, vid, [1, 1], [nx,ny],terrain_ht)
         endif
         status = NF_CLOSE(ncid)
         idstart = (ide-(nx/2))/2 + 1
         idend = idstart+(nx/2)
         jdstart = (jde-(ny/2))/2 + 1
         jdend = jdstart+(ny/2)
         itstart = max(idstart,its)
         itend = min(idend,ite)
         jtstart = max(jdstart,jts)
         jtend = min(jdend,jte)
         !RSA/JB this assumes dx,dy of the terrain is 1/2 that of the model
         DO j=jtstart,jtend
         DO i=itstart,itend
            if ((i /= idend) .AND. (j /= jdend)) grid%ht(i,j) = terrain_ht((2*i)-(2*idstart)+2,(2*j)-(2*jdstart)+2)
         ENDDO
         ENDDO
         DEALLOCATE(terrain_ht)

      ENDIF ! end of ideal_terrain options
  ENDIF ! end of use_ibm 

  DO j=jts,jte
  DO i=its,ite
    grid%phb(i,1,j) = g * grid%ht(i,j)
    grid%ph0(i,1,j) = g * grid%ht(i,j)
  ENDDO
  ENDDO

  DO J = jts, jte
  DO I = its, ite

    p_surf = interp_0( p_in, zk, grid%phb(i,1,j)/g, nl_in )
    grid%mub(i,j) = p_surf-grid%p_top

!  this is dry hydrostatic sounding (base state), so given grid%p (coordinate),
!  interp theta (from interp) and compute 1/rho from eqn. of state

    DO K = 1, kte-1
      p_level = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
      grid%pb(i,k,j) = p_level
      grid%t_init(i,k,j) = interp_0( theta, p_in, p_level, nl_in ) - t0

      !KAL hard code
      grid%t_init(i,k,j) = 0.0
      IF ((i .EQ. its) .AND. (j .EQ. jts) .AND. (k .EQ. 1)) THEN
         write(*,'(A)') "DJW: t_init is hardcoded to 0.0 in module_initialize_les"
      ENDIF
       
      grid%alb(i,k,j) = (r_d/p1000mb)*(grid%t_init(i,k,j)+t0)*(grid%pb(i,k,j)/p1000mb)**cvpm
    ENDDO

!  calc hydrostatic balance (alternatively we could interp the geopotential from the
!  sounding, but this assures that the base state is in exact hydrostatic balance with
!  respect to the model eqns.

    DO k  = 2,kte
      grid%phb(i,k,j) = grid%phb(i,k-1,j) - grid%dnw(k-1)*grid%mub(i,j)*grid%alb(i,k-1,j)
    ENDDO

  ENDDO
  ENDDO

  IF ( wrf_dm_on_monitor() ) THEN
    write(6,*) ' ptop is ',grid%p_top
    write(6,*) ' base state grid%mub(1,1), p_surf is ',grid%mub(1,1),grid%mub(1,1)+grid%p_top
  ENDIF

!  calculate full state for each column - this includes moisture.

  write(6,*) ' getting moist sounding for full state '
  dry_sounding = .false.
  CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, nl_max, nl_in, theta_surf )

  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)

!  At this point grid%p_top is already set. find the DRY mass in the column 
!  by interpolating the DRY pressure.  

   pd_surf = interp_0( pd_in, zk, grid%phb(i,1,j)/g, nl_in )

!  compute the perturbation mass and the full mass

    grid%mu_1(i,j) = pd_surf-grid%p_top - grid%mub(i,j)
    grid%mu_2(i,j) = grid%mu_1(i,j)
    grid%mu0(i,j) = grid%mu_1(i,j) + grid%mub(i,j)

! given the dry pressure and coordinate system, interp the potential
! temperature and qv

    do k=1,kde-1

      p_level = grid%znu(k)*(pd_surf - grid%p_top) + grid%p_top

      moist(i,k,j,P_QV) = interp_0( qv, pd_in, p_level, nl_in )
      grid%t_1(i,k,j)          = interp_0( theta, pd_in, p_level, nl_in ) - t0
 
      !KAL hard code
      grid%t_1(i,k,j) = 0.0
      IF ((i .EQ. its) .AND. (j .EQ. jts) .AND. (k .EQ. 1)) THEN
         write(*,'(A)') "DJW: t_1 and t_2 are hardcoded to 0.0 in module_initialize_les"
      ENDIF
      
      grid%t_2(i,k,j)          = grid%t_1(i,k,j)
      

    enddo

!  integrate the hydrostatic equation (from the RHS of the bigstep
!  vertical momentum equation) down from the top to get grid%p.
!  first from the top of the model to the top pressure

    k = kte-1  ! top level

    qvf1 = 0.5*(moist(i,k,j,P_QV)+moist(i,k,j,P_QV))
    qvf2 = 1./(1.+qvf1)
    qvf1 = qvf1*qvf2

!    grid%p(i,k,j) = - 0.5*grid%mu_1(i,j)/grid%rdnw(k)
    grid%p(i,k,j) = - 0.5*(grid%mu_1(i,j)+qvf1*grid%mub(i,j))/grid%rdnw(k)/qvf2
    qvf = 1. + rvovrd*moist(i,k,j,P_QV)
    grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
    grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)

!  down the column

    do k=kte-2,1,-1
      qvf1 = 0.5*(moist(i,k,j,P_QV)+moist(i,k+1,j,P_QV))
      qvf2 = 1./(1.+qvf1)
      qvf1 = qvf1*qvf2
      grid%p(i,k,j) = grid%p(i,k+1,j) - (grid%mu_1(i,j) + qvf1*grid%mub(i,j))/qvf2/grid%rdn(k+1)
      qvf = 1. + rvovrd*moist(i,k,j,P_QV)
      grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                  (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
      grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
    enddo

!  this is the hydrostatic equation used in the model after the
!  small timesteps.  In the model, grid%al (inverse density)
!  is computed from the geopotential.


    grid%ph_1(i,1,j) = 0.
    DO k  = 2,kte
      grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - &
                    grid%dnw(k-1)*((grid%mub(i,j)+grid%mu_1(i,j))*grid%al(i,k-1,j)+ &
                    grid%mu_1(i,j)*grid%alb(i,k-1,j)  )
                                                   
      grid%ph_2(i,k,j) = grid%ph_1(i,k,j) 
      grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
    ENDDO

    IF ( wrf_dm_on_monitor() ) THEN
    if((i==2) .and. (j==2)) then
     write(6,*) ' grid%ph_1 calc ',grid%ph_1(2,1,2),grid%ph_1(2,2,2),&
                              grid%mu_1(2,2)+grid%mub(2,2),grid%mu_1(2,2), &
                              grid%alb(2,1,2),grid%al(1,2,1),grid%rdnw(1)
    endif
    ENDIF

  ENDDO
  ENDDO

!#if 0

  IF (temp_pert) THEN
   write(*,*) "RSA using temperature perturbations in module_initialize_les"

!  thermal perturbation to kick off convection

  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',delt

! For LES, change the initial random perturbations
! For 2D test, call randx outside I-loop
! For 3D runs, call randx inside both I-J loops

  DO J = jts, min(jde-1,jte)
!   yrad = config_flags%dy*float(j-nyc)/10000.
    yrad = 0.
    DO I = its, min(ide-1,ite)
!     xrad = config_flags%dx*float(i-nxc)/10000.
      xrad = 0.
      call random_number (randx)
      randx = randx - 0.5
!     DO K = 1, kte-1
      DO K = 1, 4 

!  No bubbles for LES!
!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

!       zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
!                  +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
!       zrad = (zrad-1500.)/1500.
        zrad = 0.
        RAD=SQRT(xrad*xrad+yrad*yrad+zrad*zrad)
        IF(RAD <= 1.) THEN
!          grid%t_1(i,k,j)=grid%t_1(i,k,j)+delt*COS(.5*PI*RAD)**2
           grid%t_1(i,k,j)=grid%t_1(i,k,j)+ 0.1 *randx
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
        ENDIF
      ENDDO

!  rebalance hydrostatically

      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (grid%dnw(k-1))*(       &
                     (grid%mub(i,j)+grid%mu_1(i,j))*grid%al(i,k-1,j)+ &
                      grid%mu_1(i,j)*grid%alb(i,k-1,j)  )
                                                   
        grid%ph_2(i,k,j) = grid%ph_1(i,k,j) 
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

    ENDDO
  ENDDO

  ENDIF !temp_pert

!#endif

   IF ( wrf_dm_on_monitor() ) THEN
   write(6,*) ' grid%mu_1 from comp ', grid%mu_1(1,1)
   write(6,*) ' full state sounding from comp, ph, grid%p, grid%al, grid%t_1, qv '
   do k=1,kde-1
     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,1)+grid%phb(1,k,1), &
                                      grid%p(1,k,1)+grid%pb(1,k,1), grid%alt(1,k,1), &
                                      grid%t_1(1,k,1)+t0, moist(1,k,1,P_QV)
   enddo

   write(6,*) ' pert state sounding from comp, grid%ph_1, pp, alp, grid%t_1, qv '
   do k=1,kde-1
     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,1), &
                                      grid%p(1,k,1), grid%al(1,k,1), &
                                      grid%t_1(1,k,1), moist(1,k,1,P_QV)
   enddo
   ENDIF

! interp v

  DO J = jts, jte
  DO I = its, min(ide-1,ite)

    IF (j == jds) THEN
      z_at_v = grid%phb(i,1,j)/g
    ELSE IF (j == jde) THEN
      z_at_v = grid%phb(i,1,j-1)/g
    ELSE
      z_at_v = 0.5*(grid%phb(i,1,j)+grid%phb(i,1,j-1))/g
    END IF
    p_surf = interp_0( p_in, zk, z_at_v, nl_in )

    DO K = 1, kte-1
      p_level = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
      grid%v_1(i,k,j) = interp_0( v, p_in, p_level, nl_in )
      grid%v_2(i,k,j) = grid%v_1(i,k,j)
    ENDDO

  ENDDO
  ENDDO

! interp u

  DO J = jts, min(jde-1,jte)
  DO I = its, ite

    IF (i == ids) THEN
      z_at_u = grid%phb(i,1,j)/g
    ELSE IF (i == ide) THEN
      z_at_u = grid%phb(i-1,1,j)/g
    ELSE
      z_at_u = 0.5*(grid%phb(i,1,j)+grid%phb(i-1,1,j))/g
    END IF

    p_surf = interp_0( p_in, zk, z_at_u, nl_in )

    DO K = 1, kte-1
      p_level = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
      grid%u_1(i,k,j) = interp_0( u, p_in, p_level, nl_in )
      grid%u_2(i,k,j) = grid%u_1(i,k,j)
    ENDDO

  ENDDO
  ENDDO

!  set w

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte
  DO I = its, min(ide-1,ite)
    grid%w_1(i,k,j) = 0.
    grid%w_2(i,k,j) = 0.
  ENDDO
  ENDDO
  ENDDO

!  set a few more things

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte-1
  DO I = its, min(ide-1,ite)
    grid%h_diabatic(i,k,j) = 0.
  ENDDO
  ENDDO
  ENDDO

  IF ( wrf_dm_on_monitor() ) THEN
  DO k=1,kte-1
    grid%t_base(k) = grid%t_1(1,k,1)
    grid%qv_base(k) = moist(1,k,1,P_QV)
    grid%u_base(k) = grid%u_1(1,k,1)
    grid%v_base(k) = grid%v_1(1,k,1)
    grid%z_base(k) = 0.5*(grid%phb(1,k,1)+grid%phb(1,k+1,1)+grid%ph_1(1,k,1)+grid%ph_1(1,k+1,1))/g
  ENDDO
  ENDIF
  CALL wrf_dm_bcast_real( grid%t_base , kte )
  CALL wrf_dm_bcast_real( grid%qv_base , kte )
  CALL wrf_dm_bcast_real( grid%u_base , kte )
  CALL wrf_dm_bcast_real( grid%v_base , kte )
  CALL wrf_dm_bcast_real( grid%z_base , kte )

  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)
     !DJW original code from v3.8.1 is below
!     grid%tsk(I,J)=theta_surf * (p_surf/p1000mb)**rcp
!     grid%tmn(I,J)=grid%tsk(I,J)-0.5
     !DJW the method used in IBMv3.6.1 is below
     thtmp = interp_0( theta, zk, grid%ht(i,j), nl_in )
     ptmp = interp_0( p_in, zk, grid%ht(i,j), nl_in )    
     grid%tsk(I,J)=thtmp * (ptmp/p1000mb)**rcp
     grid%tmn(I,J)=grid%tsk(I,J)-0.5
  ENDDO
  ENDDO

  !KAL add perturbations in u
  IF (config_flags%init_vel_pert) THEN
     write(*,*) "RSA using velocity perturbations"
     iseed = -2
     amp = 0.01
     DO k=1,(kde/3)
         CALL ranary(ide,jde-1,ids,ide,jds,jde-1,iseed,amp,useed)
         DO i=its,ite
         DO j=jts,min(jde-1,jte)
            grid%u_1(i,k,j) = grid%u_1(i,k,j)+useed(i,j)
            grid%u_2(i,k,j) = grid%u_2(i,k,j)+useed(i,j)
         ENDDO
         ENDDO
     ENDDO
  ENDIF !init_vel_pert

  !KAL-----changes for land surface model---------------------------
  IF (config_flags%ideal_soil) THEN

     !ideal_terrain = 14 is for Granite Mountain
     IF (config_flags%ideal_terrain .NE. 14) THEN
         CALL nl_set_mminlu(1,'USGS')
  
!         DO k =1,kte-1
!             grid%t_base(k) = 0.
!             grid%qv_base(k) = 0.
!             grid%u_base(k) = 0.
!             grid%v_base(k) = 0.
!             grid%z_base(k) = 0.
!         ENDDO
 
  
         DO j=jts,min(jde-1,jte)
         DO i=its,min(ide-1,ite)
             grid%ivgtyp(I,J) = 10 !savanna
             grid%lu_index(I,J) = grid%ivgtyp(I,J)
             grid%isltyp(I,J) = 3 !sandy loam

             !KAL 2 meter values
             grid%q2(I,J) = interp_0(qv, zk, grid%ht(i,j), nl_in)
             grid%th2(I,J) = interp_0(theta, zk, grid%ht(i,j), nl_in)
             grid%t2(I,J)=grid%th2(I,J)*(((grid%p(I,1,J)+grid%pb(I,1,J))/p1000mb)**rcp)
             !KAL this is calculated at the begining of solve anyway
             grid%psfc(I,J) = grid%p(I,1,J)+grid%pb(I,1,J) 
     
             !KAL 10 meter values
             grid%u10(I,J) = interp_0(u, zk, grid%ht(I,J) + 10.0, nl_in)
             grid%v10(I,J) = interp_0(v, zk, grid%ht(I,J) + 10.0, nl_in)

             !land mask = 1 for land, 2 for water
             grid%landmask(I,J) = 1.0
             grid%xland(I,J) = 1.0
      
             grid%greenfrac = 0.10000001
             grid%vegfra(I,J) = 0.10000001
             grid%shdmin(I,J) = 0.10000001
             grid%shdmax(I,J) = 0.10000001
       
         ENDDO
         ENDDO   
         !for the trex simulations, the veg fraction is specified above, so this is not needed 
!         CALL monthly_interp_to_date ( grid%greenfrac, current_date, grid%vegfra, &
!                                       ids, ide, jds, jde, kds, kde, &
!                                       ims, ime, jms, jme, kms, kme, &
!                                       its, ite, jts, jte, kts, kte ) 

!         CALL monthly_min_max ( grid%greenfrac , grid%shdmin, grid%shdmax, &
!                                ids, ide, jds, jde, kds, kde, &
!                                ims, ime, jms, jme, kms, kme, &
!                                its, ite, jts, jte, kts, kte )

         DO j=jts, MIN(jte,jde-1)
         DO i=its, MIN(ite,ide-1)
             grid%vegfra(i,j) = grid%vegfra(i,j)*100.0 !vegetation fraction
             grid%shdmax(i,j) = grid%shdmax(i,j)*100.0 !annual max veg fraction
             grid%shdmin(i,j) = grid%shdmin(i,j)*100.0 !annual min veg fraction
             grid%snoalb(i,j) = 60.0 !average from owens valley data
             DO ns = 1,config_flags%num_soil_layers
!                 grid%tslb(i,ns,j) = grid%t2(i,j) !initialize soil temp as equal to atm
                 grid%tslb(i,ns,j)=grid%tsk(i,j)
                 grid%smois(i,ns,j) = 0.0868 !twenty percent saturation rate * porosity for sandy loam
             ENDDO
      
             grid%sst(i,j) = grid%tsk(i,j) !sst is sea surface temperature, tsk is surface skin temperature
             grid%albbck(i,j) = 0.27 !background albedo
!             grid%xlat(i,j) = 36.0
!             grid%xlong(i,j) = 0.0
             grid%xlat(i,j) = 35.44 !OKC
             grid%xlong(i,j) = 97.55 !OKC
         ENDDO
         ENDDO
  
         DO j=jts, MIN(jte,jde-1)
         DO i=its, ite 
!             grid%xlat_u(i,j) = 36.0
!             grid%xlong_u(i,j) = 0.0
             grid%xlat_u(i,j) = 35.44 !OKC
             grid%xlong_u(i,j) = 97.55 !OKC
         ENDDO
         ENDDO
  
         DO j=jts, jte
         DO i=its, MIN(ite,ide-1)    
!             grid%xlat_v(i,j) = 36.0
!             grid%xlong_v(i,j) = 0.0
             grid%xlat_v(i,j) = 35.44 !OKC
             grid%xlong_v(i,j) = 97.55 !OKC
         ENDDO
         ENDDO    

         grid%zs =  [0.05, 0.25, 0.70, 1.50] !from real initialization, depths of centers of soil layers
         grid%dzs = [0.10, 0.30, 0.60, 1.00] !from real initialization, thicknesses of soil layers

     ELSEIF (config_flags%ideal_terrain .EQ. 14) THEN !RSA Granite Mountain
         print *, 'RSA ideal soil, ideal_terrain = 14, GRANITE MOUNTAIN'
         CALL nl_set_mminlu(1,'USGS')

         !Read in available land surface data
         status = NF_OPEN("granite_mtn.nc", nf_nowrite, ncid) 
         status = NF_INQ_DIMID(ncid, 'nx', vid)
         status = NF_INQ_DIMLEN(ncid, vid, nx)
         status = NF_INQ_DIMID(ncid, 'ny', vid)
         status = NF_INQ_DIMLEN(ncid, vid, ny) 

         !LANDUSE
         status = NF_INQ_VARID(ncid, "LANDUSE", vid)
         ALLOCATE(landuse_read(nx,ny))
         landuse_read = 0. 
         IF (KIND(landuse_read) .EQ. 4) THEN
             status = NF_GET_VARA_REAL(ncid, vid, [1, 1], [nx,ny], landuse_read) 
             IF (status .NE. NF_NOERR) print *, NF_STRERROR(status) 
         ELSEIF (KIND(landuse_read) .EQ. 8) THEN    
             status = NF_GET_VARA_DOUBLE(ncid, vid, [1, 1], [nx,ny], landuse_read)
             IF (status .ne. NF_NOERR) print *, NF_STRERROR(status) 
         ENDIF

         !SOILTYPE
         status = NF_INQ_VARID(ncid, "SOILTYPE", vid)
         ALLOCATE(soiltype_read(nx,ny))
         soiltype_read = 0. 
         IF (KIND(soiltype_read) .EQ. 4) THEN
             status = NF_GET_VARA_REAL(ncid, vid, [1, 1], [nx,ny], soiltype_read) 
             IF (status .NE. NF_NOERR) print *, NF_STRERROR(status) 
         ELSEIF (KIND(soiltype_read) .EQ. 8) THEN    
             status = NF_GET_VARA_DOUBLE(ncid, vid, [1, 1], [nx,ny], soiltype_read)
             IF (status .NE. NF_NOERR) print *, NF_STRERROR(status) 
         ENDIF

         !LATITUDE
         status = NF_INQ_VARID(ncid, "LAT", vid)
         ALLOCATE(lat_read(nx,ny))
         lat_read = 0. 
         IF (KIND(lat_read) .EQ. 4) THEN
             status = NF_GET_VARA_REAL(ncid, vid, [1, 1], [nx,ny], lat_read) 
             IF (status .NE. NF_NOERR) print *, NF_STRERROR(status) 
         ELSEIF (KIND(lat_read) .EQ. 8) THEN    
             status = NF_GET_VARA_DOUBLE(ncid, vid, [1, 1], [nx,ny], lat_read)
             IF (status .NE. NF_NOERR) print *, NF_STRERROR(status) 
         ENDIF

         !LONGITUDE
         status = NF_INQ_VARID(ncid, "LONG", vid)
         ALLOCATE(long_read(nx,ny))
         long_read = 0. 
         IF (KIND(long_read) .EQ. 4) THEN
             status = NF_GET_VARA_REAL(ncid, vid, [1, 1], [nx,ny], long_read) 
             IF (status .NE. NF_NOERR) print *, NF_STRERROR(status) 
         ELSEIF (KIND(long_read) .EQ. 8) THEN    
             status = NF_GET_VARA_DOUBLE(ncid, vid, [1, 1], [nx,ny], long_read)
             IF (status .NE. NF_NOERR) print *, NF_STRERROR(status) 
         ENDIF

         status = NF_CLOSE(ncid)

         idstart = (ide-(nx/2))/2 + 1
         idend = idstart+(nx/2)
         jdstart = (jde-(ny/2))/2 + 1
         jdend = jdstart+(ny/2)
         itstart = max(idstart,its)
         itend = min(idend,ite)
         jtstart = max(jdstart,jts)
         jtend = min(jdend,jte)

         DO j=jtstart,jtend
         DO i=itstart,itend
             IF ((i .NE. idend) .AND. (j .NE. jdend)) grid%xlat(i,j) = lat_read((2*i)-(2*idstart)+2,(2*j)-(2*jdstart)+2)
             IF ((i .NE. idend) .AND. (j .NE. jdend)) grid%xlong(i,j) = long_read((2*i)-(2*idstart)+2,(2*j)-(2*jdstart)+2)
             IF ((i .NE. idend) .AND. (j .NE. jdend)) grid%ivgtyp(i,j) = landuse_read((2*i)-(2*idstart)+2,(2*j)-(2*jdstart)+2)
             IF ((i .NE. idend) .AND. (j .NE. jdend)) grid%lu_index(i,j) = landuse_read((2*i)-(2*idstart)+2,(2*j)-(2*jdstart)+2)
             IF ((i .NE. idend) .AND. (j .NE. jdend)) grid%isltyp(i,j) = soiltype_read((2*i)-(2*idstart)+2,(2*j)-(2*jdstart)+2)
         ENDDO
         ENDDO

         DEALLOCATE(landuse_read)
         DEALLOCATE(soiltype_read)
         DEALLOCATE(lat_read)
         DEALLOCATE(long_read)

         DO j=jts,MIN(jde-1,jte)
         DO i=its,MIN(ide-1,ite)
             !KAL 2 meter values
             grid%q2(I,J) = interp_0(qv, zk, grid%ht(i,j), nl_in)
             grid%th2(I,J) = interp_0(theta, zk, grid%ht(i,j), nl_in)
             grid%t2(I,J)=grid%th2(I,J)*(((grid%p(I,1,J)+grid%pb(I,1,J))/p1000mb)**rcp)
             !KAL this is calculated at the begining of solve anyway
             grid%psfc(I,J) = grid%p(I,1,J)+grid%pb(I,1,J)
             !KAL 10 meter values
             grid%u10(I,J) = interp_0(u, zk, grid%ht(I,J) + 10.0, nl_in)
             grid%v10(I,J) = interp_0(v, zk, grid%ht(I,J) + 10.0, nl_in)
             !land mask = 1 for land, 2 for water
             grid%landmask(I,J) = 1.0
             grid%xland(I,J) = 1.0
             grid%greenfrac = 0.10000001
             grid%vegfra(I,J) = 0.10000001
             grid%shdmin(I,J) = 0.10000001
             grid%shdmax(I,J) = 0.10000001
         ENDDO
         ENDDO   

         DO j=jts,MIN(jte,jde-1)
         DO i=its,MIN(ite,ide-1)
             grid%vegfra(i,j) = grid%vegfra(i,j)*100.0 !vegetation fraction
             grid%shdmax(i,j) = grid%shdmax(i,j)*100.0 !annual max veg fraction
             grid%shdmin(i,j) = grid%shdmin(i,j)*100.0 !annual min veg fraction
!             grid%snoalb(i,j) = 60.0 !average from owens valley data
             DO ns=1,config_flags%num_soil_layers
!                 grid%tslb(i,ns,j) = grid%t2(i,j)     ! initialize soil temp as equal to atm
                 grid%tslb(i,ns,j)=grid%tsk(i,j)
                 grid%smois(i,ns,j) = 0.15
             ENDDO
         ENDDO
         ENDDO

         grid%zs =  [0.05, 0.25, 0.70, 1.50] !from real initialization, depths of centers of soil layers
         grid%dzs = [0.10, 0.30, 0.60, 1.00] !from real initialization, thicknesses of soil layers

     ENDIF !ideal_terrain
 ENDIF !ideal_soil

 END SUBROUTINE init_domain_rk
 
!---------------------------------------------------------------------
   !KAL / DJW: added subroutine from module_initialize_real.F
   SUBROUTINE monthly_min_max ( field_in , field_min , field_max , &
                      ids , ide , jds , jde , kds , kde , &
                      ims , ime , jms , jme , kms , kme , &
                      its , ite , jts , jte , kts , kte )

   !  Plow through each month, find the max, min values for each i,j.
   
      IMPLICIT NONE

      INTEGER , INTENT(IN)        :: ids , ide , jds , jde , kds , kde , &
                                     ims , ime , jms , jme , kms , kme , &
                                     its , ite , jts , jte , kts , kte

      REAL , DIMENSION(ims:ime,12,jms:jme) , INTENT(IN)  :: field_in
      REAL , DIMENSION(ims:ime,   jms:jme) , INTENT(OUT) :: field_min , field_max

      !  Local vars

      INTEGER :: i , j , l
      REAL :: minner , maxxer

      DO j = jts , MIN(jde-1,jte)
         DO i = its , MIN(ide-1,ite)
            minner = field_in(i,1,j)
            maxxer = field_in(i,1,j)
            DO l = 2 , 12
               IF ( field_in(i,l,j) .LT. minner ) THEN
                  minner = field_in(i,l,j)
               END IF
               IF ( field_in(i,l,j) .GT. maxxer ) THEN
                  maxxer = field_in(i,l,j)
               END IF
            END DO
            field_min(i,j) = minner
            field_max(i,j) = maxxer
         END DO
      END DO
   
   END SUBROUTINE monthly_min_max

!---------------------------------------------------------------------
   !KAL / DJW: added subroutine from module_initialize_real.F
   SUBROUTINE monthly_interp_to_date ( field_in , date_str , field_out , &
                      ids , ide , jds , jde , kds , kde , &
                      ims , ime , jms , jme , kms , kme , &
                      its , ite , jts , jte , kts , kte )

   !  Linrarly in time interpolate data to a current valid time.  The data is
   !  assumed to come in "monthly", valid at the 15th of every month.
   
      IMPLICIT NONE

      INTEGER , INTENT(IN)        :: ids , ide , jds , jde , kds , kde , &
                                     ims , ime , jms , jme , kms , kme , &
                                     its , ite , jts , jte , kts , kte

      CHARACTER (LEN=24) , INTENT(IN) :: date_str
      REAL , DIMENSION(ims:ime,12,jms:jme) , INTENT(IN)  :: field_in
      REAL , DIMENSION(ims:ime,   jms:jme) , INTENT(OUT) :: field_out

      !  Local vars

      INTEGER :: i , j , l
      INTEGER , DIMENSION(0:13) :: middle
      INTEGER :: target_julyr , target_julday , target_date
      INTEGER :: julyr , julday , int_month , month1 , month2
      REAL :: gmt
      CHARACTER (LEN=4) :: yr
      CHARACTER (LEN=2) :: mon , day15


      WRITE(day15,FMT='(I2.2)') 15
      DO l = 1 , 12
         WRITE(mon,FMT='(I2.2)') l
         CALL get_julgmt ( date_str(1:4)//'-'//mon//'-'//day15//'_'//'00:00:00.0000' , julyr , julday , gmt )
         middle(l) = julyr*1000 + julday
      END DO

      l = 0
      middle(l) = middle( 1) - 31

      l = 13
      middle(l) = middle(12) + 31

      CALL get_julgmt ( date_str , target_julyr , target_julday , gmt )
      target_date = target_julyr * 1000 + target_julday
      find_month : DO l = 0 , 12
         IF ( ( middle(l) .LT. target_date ) .AND. ( middle(l+1) .GE. target_date ) ) THEN
            DO j = jts , MIN ( jde-1 , jte )
               DO i = its , MIN (ide-1 , ite )
                  IF ( skip_middle_points_t ( ids , ide , jds , jde , i , j , em_width , hold_ups ) ) CYCLE
                  int_month = l
                  IF ( ( int_month .EQ. 0 ) .OR. ( int_month .EQ. 12 ) ) THEN
                     month1 = 12
                     month2 =  1
                  ELSE
                     month1 = int_month
                     month2 = month1 + 1
                  END IF
                  field_out(i,j) =  ( field_in(i,month2,j) * ( target_date - middle(l)   ) + &
                                      field_in(i,month1,j) * ( middle(l+1) - target_date ) ) / &
                                    ( middle(l+1) - middle(l) )
               END DO
            END DO
            EXIT find_month
         END IF
      END DO find_month

   END SUBROUTINE monthly_interp_to_date
   
!---------------------------------------------------------------------------------------------------   

   SUBROUTINE init_module_initialize
   END SUBROUTINE init_module_initialize

!---------------------------------------------------------------------

!  test driver for get_sounding
!
!      implicit none
!      integer n
!      parameter(n = 1000)
!      real zk(n),p(n),theta(n),rho(n),u(n),v(n),qv(n),pd(n)
!      logical dry
!      integer nl,k
!
!      dry = .false.
!      dry = .true.
!      call get_sounding( zk, p, pd, theta, rho, u, v, qv, dry, n, nl )
!      write(6,*) ' input levels ',nl
!      write(6,*) ' sounding '
!      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
!      do k=1,nl
!        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), pd(k), theta(k), rho(k), u(k), v(k), qv(k)
!      enddo
!      end
!
!---------------------------------------------------------------------------

      subroutine get_sounding( zk, p, p_dry, theta, rho, &
                               u, v, qv, dry, nl_max, nl_in, th_surf )
      implicit none

      integer nl_max, nl_in
      real zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
           u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max)
      logical dry

      integer n
      !parameter(n=1000)
      parameter(n=27000) !KAL
      logical debug
      parameter( debug = .true.)

! input sounding data

      real p_surf, th_surf, qv_surf
      real pi_surf, pi(n)
      real h_input(n), th_input(n), qv_input(n), u_input(n), v_input(n)

! diagnostics

      real rho_surf, p_input(n), rho_input(n)
      real pm_input(n)  !  this are for full moist sounding

! local data

      real r
      parameter (r = r_d)
      integer k, it, nl
      real qvf, qvf1, dz

!  first, read the sounding

      write(*,'(A)') "DJW[module_initialize_les/get_sounding]: calling read_sounding"
      call read_sounding( p_surf, th_surf, qv_surf, &
                          h_input, th_input, qv_input, u_input, v_input,n, nl, debug )

      if(dry) then
       do k=1,nl
         qv_input(k) = 0.
       enddo
      endif

      if(debug) write(6,*) ' number of input levels = ',nl

        nl_in = nl
        if(nl_in .gt. nl_max ) then
          write(6,*) ' too many levels for input arrays ',nl_in,nl_max
          call wrf_error_fatal ( ' too many levels for input arrays ' )
        end if

!  compute diagnostics,
!  first, convert qv(g/kg) to qv(g/g)

      do k=1,nl
        qv_input(k) = 0.001*qv_input(k)
      enddo

      p_surf = 100.*p_surf  ! convert to pascals
      qvf = 1. + rvovrd*qv_input(1) 
      rho_surf = 1./((r/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm))
      pi_surf = (p_surf/p1000mb)**(r/cp)

      if(debug) then
        write(6,*) ' surface density is ',rho_surf
        write(6,*) ' surface pi is      ',pi_surf
      end if


!  integrate moist sounding hydrostatically, starting from the
!  specified surface pressure
!  -> first, integrate from surface to lowest level

          qvf = 1. + rvovrd*qv_input(1) 
          qvf1 = 1. + qv_input(1)
          rho_input(1) = rho_surf
          dz = h_input(1)
          do it=1,10
            pm_input(1) = p_surf &
                    - 0.5*dz*(rho_surf+rho_input(1))*g*qvf1
            rho_input(1) = 1./((r/p1000mb)*th_input(1)*qvf*((pm_input(1)/p1000mb)**cvpm))
          enddo

! integrate up the column

          do k=2,nl
            rho_input(k) = rho_input(k-1)
            dz = h_input(k)-h_input(k-1)
            qvf1 = 0.5*(2.+(qv_input(k-1)+qv_input(k)))
            qvf = 1. + rvovrd*qv_input(k)   ! qv is in g/kg here
 
            !KAL changed it loop from 10 to 20
            !KAL the atm was not fully balanced with 10 iterations
            do it=1,20
              pm_input(k) = pm_input(k-1) &
                      - 0.5*dz*(rho_input(k)+rho_input(k-1))*g*qvf1
              rho_input(k) = 1./((r/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm))
            enddo
          enddo

!  we have the moist sounding

!  next, compute the dry sounding using p at the highest level from the
!  moist sounding and integrating down.

        p_input(nl) = pm_input(nl)

          do k=nl-1,1,-1
            dz = h_input(k+1)-h_input(k)
            p_input(k) = p_input(k+1) + 0.5*dz*(rho_input(k)+rho_input(k+1))*g
          enddo


        do k=1,nl

          zk(k) = h_input(k)
          p(k) = pm_input(k)
          p_dry(k) = p_input(k)
          theta(k) = th_input(k)
          rho(k) = rho_input(k)
          u(k) = u_input(k)
          v(k) = v_input(k)
          qv(k) = qv_input(k)

        enddo

     if(debug) then
      write(6,*) ' sounding '
      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
      do k=1,nl
        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), p_dry(k), theta(k), rho(k), u(k), v(k), qv(k)
      enddo

     end if

      end subroutine get_sounding

!-------------------------------------------------------

      subroutine read_sounding( ps,ts,qvs,h,th,qv,u,v,n,nl,debug )
      implicit none
      integer n,nl
      real ps,ts,qvs,h(n),th(n),qv(n),u(n),v(n)
      logical end_of_file
      logical debug

      integer k

      open(unit=10,file='input_sounding',form='formatted',status='old')
      rewind(10)
      read(10,*) ps, ts, qvs
      if(debug) then
        write(6,*) ' input sounding surface parameters '
        write(6,*) ' surface pressure (mb) ',ps
        write(6,*) ' surface pot. temp (K) ',ts
        write(6,*) ' surface mixing ratio (g/kg) ',qvs
      end if

      end_of_file = .false.
      k = 0

      do while (.not. end_of_file)

        read(10,*,end=100) h(k+1), th(k+1), qv(k+1), u(k+1), v(k+1)
        k = k+1
        if(debug) write(6,'(1x,i3,5(1x,e10.3))') k, h(k), th(k), qv(k), u(k), v(k)
        go to 110
 100    end_of_file = .true.
 110    continue
      enddo

      nl = k

      close(unit=10,status = 'keep')

      end subroutine read_sounding

END MODULE module_initialize_ideal
