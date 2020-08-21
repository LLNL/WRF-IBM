SUBROUTINE init_domain_constants_em ( parent , nest )
   USE module_domain, ONLY : domain
   IMPLICIT NONE
   TYPE(domain)  :: parent , nest

   INTEGER iswater, islake, isice, isurban, isoilwater, map_proj, julyr, julday
   REAL    truelat1 , truelat2 , gmt , moad_cen_lat , stand_lon, pole_lat, pole_lon
   CHARACTER (LEN=256) :: char_junk

! single-value constants

   nest%p_top   = parent%p_top
   nest%save_topo_from_real   = parent%save_topo_from_real
   nest%cfn     = parent%cfn
   nest%cfn1    = parent%cfn1
   nest%rdx     = 1./nest%dx
   nest%rdy     = 1./nest%dy
!  nest%dts     = nest%dt/float(nest%time_step_sound)
   nest%dtseps  = parent%dtseps  ! used in height model only?
   nest%resm    = parent%resm    ! used in height model only?
   nest%zetatop = parent%zetatop ! used in height model only?
   nest%cf1     = parent%cf1
   nest%cf2     = parent%cf2
   nest%cf3     = parent%cf3
   nest%gmt     = parent%gmt
   nest%julyr   = parent%julyr
   nest%julday  = parent%julday
   nest%iswater = parent%iswater
   nest%isice   = parent%isice
   nest%isurban = parent%isurban
   nest%islake  = parent%islake
   nest%isoilwater = parent%isoilwater
   nest%mminlu  = trim(parent%mminlu)
   nest%tiso    = parent%tiso
   nest%tlp     = parent%tlp
   nest%p00     = parent%p00
   nest%t00     = parent%t00
   nest%tlp_strat= parent%tlp_strat
   nest%p_strat = parent%p_strat
!cyl: variables for trajectory /float
   nest%traj_k = parent%traj_k
   nest%traj_long = parent%traj_long
   nest%traj_lat = parent%traj_lat
   nest%this_is_an_ideal_run = parent%this_is_an_ideal_run
   nest%lake_depth_flag = parent%lake_depth_flag
   !DJW: added for vertically nested initialization of nests
   nest%use_baseparam_fr_nml = parent%use_baseparam_fr_nml

   CALL nl_get_mminlu ( 1, char_junk )
   CALL nl_get_iswater( 1, iswater )
   CALL nl_get_islake ( 1, islake )
   CALL nl_get_isice  ( 1, isice )
   CALL nl_get_isurban( 1, isurban )
   CALL nl_get_isoilwater(1, isoilwater )
   CALL nl_get_truelat1 ( 1 , truelat1 )
   CALL nl_get_truelat2 ( 1 , truelat2 )
   CALL nl_get_moad_cen_lat ( 1 , moad_cen_lat )
   CALL nl_get_stand_lon ( 1 , stand_lon )
   CALL nl_get_pole_lat ( 1 , pole_lat )
   CALL nl_get_pole_lon ( 1 , pole_lon )
   CALL nl_get_map_proj ( 1 , map_proj )
   CALL nl_get_gmt ( 1 , gmt)
   CALL nl_get_julyr ( 1 , julyr)
   CALL nl_get_julday ( 1 , julday)
   IF ( nest%id .NE. 1 ) THEN
     CALL nl_set_gmt (nest%id, gmt)
     CALL nl_set_julyr (nest%id, julyr)
     CALL nl_set_julday (nest%id, julday)
     CALL nl_set_iswater ( nest%id, iswater )
     CALL nl_set_islake  ( nest%id, islake )
     CALL nl_set_isice   ( nest%id, isice )
     CALL nl_set_isurban ( nest%id, isurban )
     CALL nl_set_isoilwater ( nest%id, isoilwater )
     CALL nl_set_mminlu ( nest%id, char_junk )
     CALL nl_set_truelat1 ( nest%id , truelat1 )
     CALL nl_set_truelat2 ( nest%id , truelat2 )
     CALL nl_set_moad_cen_lat ( nest%id , moad_cen_lat )
     CALL nl_set_stand_lon ( nest%id , stand_lon )
     CALL nl_set_pole_lat ( nest%id , pole_lat )
     CALL nl_set_pole_lon ( nest%id , pole_lon )
     CALL nl_set_map_proj ( nest%id , map_proj )
   END IF
   nest%gmt     = gmt
   nest%julday  = julday
   nest%julyr   = julyr
   nest%iswater = iswater
   nest%islake  = islake
   nest%isice   = isice
   nest%isoilwater = isoilwater
   nest%mminlu  = trim(char_junk)
   nest%truelat1= truelat1
   nest%truelat2= truelat2
   nest%moad_cen_lat= moad_cen_lat
   nest%stand_lon= stand_lon
   nest%pole_lat= pole_lat
   nest%pole_lon= pole_lon
   nest%map_proj= map_proj

   nest%step_number  = parent%step_number

! 1D constants (Z)

   IF (nest%vert_refine_method .EQ. 0) THEN
     nest%fnm     = parent%fnm
     nest%fnp     = parent%fnp
     nest%rdnw    = parent%rdnw
     nest%rdn     = parent%rdn
     nest%dnw     = parent%dnw
     nest%dn      = parent%dn
     nest%znu     = parent%znu
     nest%t_base  = parent%t_base
     nest%u_base  = parent%u_base
     nest%v_base  = parent%v_base
     nest%qv_base = parent%qv_base
     nest%z_base  = parent%z_base
   ENDIF
   nest%znw     = parent%znw
   nest%dzs       = parent%dzs
   nest%zs        = parent%zs

END SUBROUTINE init_domain_constants_em


!---------------------------------------------------------------------------------------------------

SUBROUTINE init_domain_vert_nesting ( parent, nest, use_baseparam_fr_nml )

!KAL this is a driver to initialize the vertical coordinates for the nest when vertical nesting is used
   USE module_domain
   IMPLICIT NONE
   TYPE(domain), POINTER                      :: parent, nest
   LOGICAL                                    :: use_baseparam_fr_nml
   !local
   REAL, DIMENSION(parent%e_vert)             :: znw_c

   INTERFACE

      SUBROUTINE vert_cor_vertical_nesting_integer(nest,znw_c,k_dim_c)
         USE module_domain
         TYPE(domain), POINTER :: nest
         integer , intent(in) :: k_dim_c
         real , dimension(k_dim_c), INTENT(IN) :: znw_c
      END SUBROUTINE vert_cor_vertical_nesting_integer

      SUBROUTINE vert_cor_vertical_nesting_arbitrary(nest,znw_c,kde_c,use_baseparam_fr_nml)
          USE module_domain
          TYPE(domain), POINTER :: nest
          INTEGER, INTENT(IN   ) :: kde_c
          REAL, DIMENSION(kde_c), INTENT(IN   ) :: znw_c
          LOGICAL, INTENT(IN   ) :: use_baseparam_fr_nml
      END SUBROUTINE vert_cor_vertical_nesting_arbitrary
      
      SUBROUTINE vert_cor_vertical_nesting_z_specified(nest)
          USE module_domain
          USE module_init_utilities, ONLY : interp_0
          TYPE(domain), POINTER :: nest
      END SUBROUTINE vert_cor_vertical_nesting_z_specified

      SUBROUTINE vert_cor_vertical_nesting_z_specified_real(nest,use_baseparam_fr_nml)
          USE module_domain
          TYPE(domain), POINTER :: nest
          LOGICAL, INTENT(IN   ) :: use_baseparam_fr_nml
      END SUBROUTINE vert_cor_vertical_nesting_z_specified_real

   END INTERFACE

   
   ! save the coarse grid values here
   
   znw_c = nest%znw(1:parent%e_vert)
   
   write(*,'(A,I0.2)') "DJW[nest_init_utils/init_domain_vert_nesting]: starting vnest initialization for d",nest%id
   ! calculate the nest (fine) grid values here
   ! one of these calls goes to integer refinement in the vertical direction, and one goes to arbitrary refinement.  Eventually the call to integer refinement will be obsolete. 
   IF (nest%vert_refine_method .EQ. 1) THEN  !if you are in this subroutine there is vertical nesting- (i.e. nest%e_vert /= parent%e_vert to enter this subroutine)
      CALL vert_cor_vertical_nesting_integer(nest,znw_c,parent%e_vert)
   ELSEIF (nest%vert_refine_method .EQ. 2) THEN
      CALL vert_cor_vertical_nesting_arbitrary(nest,znw_c,parent%e_vert,use_baseparam_fr_nml)
   ELSEIF (nest%vert_refine_method .EQ. 3) THEN
      CALL vert_cor_vertical_nesting_z_specified(nest)
   ELSEIF (nest%vert_refine_method .EQ. 4) THEN
      CALL vert_cor_vertical_nesting_z_specified_real(nest,use_baseparam_fr_nml)
   ENDIF
   write(*,'(A,I0.2)') "DJW[nest_init_utils/init_domain_vert_nesting]: finished vnest initialization for d",nest%id

END SUBROUTINE init_domain_vert_nesting

!-----------------------------------------------------------------------------------------

SUBROUTINE vert_cor_vertical_nesting_z_specified( nest )
   USE module_domain
   USE module_init_utilities, ONLY : interp_0
   IMPLICIT NONE
   !inputs
   TYPE(domain), POINTER                   :: nest
   !local
   INTEGER, PARAMETER                      :: nl_max=27000
   REAL, DIMENSION(nl_max)                 :: zk,          &
                                              p_in,        &
                                              theta,       &
                                              rho,         &
                                              u,           &
                                              v,           &
                                              qv,          &
                                              pd_in
   INTEGER                                 :: nl_in,       &
                                              k,           &
                                              ks,          &
                                              i
   LOGICAL                                 :: keepgoing
   REAL                                    :: p_surf,      &
                                              theta_surf,  &
                                              mub_ave,     &
                                              cof1,        &
                                              cof2
   REAL, DIMENSION(1:nest%e_vert)          :: p_specified, &
                                              z_specified

   CALL wrf_debug(100, "dyn_em/nest_init_utils/subroutine vert_cor_vertical_nesting_z_specified")

   !DJW this statement below would be nice to have but unfortunately the model
   !    is coded in a way that makes it near-impossible to set
   !    'this_is_an_ideal_run' correctly for restart runs. In other words, if
   !    you restart an ideal run, that variable is always going to be false. :'(
!   IF (.NOT.(parent%this_is_an_ideal_run)) THEN
!      CALL wrf_error_fatal( "--- ERROR: vert_refine_method=3 can only be used with ideal.exe" )
!   ENDIF

   !DJW find the number of z_specified entries that will precede the
   !    section associated with our domain of interest.
   ks = 1
   i = 1
   DO WHILE (i .LE. model_config_rec%max_dom)
      IF (model_config_rec%grid_id(i) .EQ. nest%grid_id) THEN
         i = model_config_rec%max_dom+1
      ELSE
         ks = ks+model_config_rec%e_vert(i)
         i = i+1
      ENDIF
   ENDDO

   DO k=ks,ks+nest%e_vert-1
      z_specified(k-ks+1) = model_config_rec%z_specified(k)
   ENDDO

   DO k=1,nest%e_vert
      write(wrf_err_message,'(A,I1,A,I3,A,F8.2)') "DJW[nest_init_utils/vert_cor_vertical_nesting_z_specified]: grid_id=",nest%id," z_specified(",k,")=",z_specified(k)
      CALL wrf_message( wrf_err_message )
   ENDDO

   !Check that the levels specified in the namelist are uniformly increasing
   DO k=1,nest%e_vert
      IF (k .NE. 1) THEN
         IF (z_specified(k) .LE. z_specified(k-1)) THEN
            write(wrf_err_message,'(2(A,I3),A,I1)') "--- ERROR: z_specified(",k,") <= z_specified(",k-1,") for grid_id=",nest%grid_id
            CALL wrf_error_fatal( wrf_err_message )
         ENDIF
      ENDIF
   ENDDO

   !Check that the top-most specified height matches the requested ztop
   IF (z_specified(nest%e_vert) .NE. nest%ztop) THEN
      write(wrf_err_message,'(A,F6.1,A,I3,A,F6.1)') "--- ERROR: nest%ztop ",nest%ztop," is not consistent with z_specified(",nest%e_vert,") ",z_specified(nest%e_vert)
      CALL wrf_error_fatal( TRIM(wrf_err_message) )
   ENDIF

   CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, .TRUE., nl_max, nl_in, theta_surf )
   p_surf = interp_0( p_in, zk, z_specified(1), nl_in )
   mub_ave = p_surf-nest%p_top
   DO k=1,nest%e_vert
      p_specified(k) = interp_0( p_in, zk, z_specified(k), nl_in )
   ENDDO
   nest%znw(1) = 1.0
   DO k=2,nest%e_vert-1
      nest%znw(k) = (p_specified(k)-nest%p_top)/mub_ave
   ENDDO
   nest%znw(nest%e_vert) = 0.0

   DO k=1,nest%e_vert-1
      nest%dnw(k) = nest%znw(k+1)-nest%znw(k)
      nest%rdnw(k) = 1.0/nest%dnw(k)
      nest%znu(k) = 0.5*(nest%znw(k+1)+nest%znw(k))
   ENDDO

   DO k=1,nest%e_vert-1
      write(*,'(A,I2.2,2(A,I3,A,F6.3))') "d",nest%id,"   znw(",k,")=",nest%znw(k),"   znu(",k,")=",nest%znu(k)
   ENDDO
   write(*,'(A,I2.2,A,I3,A,F6.3)') "d",nest%id,"   znw(",nest%e_vert,")=",nest%znw(nest%e_vert)

   DO k=2,nest%e_vert-1
      nest%dn(k) = 0.5*(nest%dnw(k)+nest%dnw(k-1))
      nest%rdn(k) = 1.0/nest%dn(k)
      nest%fnp(k) = 0.5*nest%dnw(k)/nest%dn(k)
      nest%fnm(k) = 0.5*nest%dnw(k-1)/nest%dn(k)
   ENDDO

   cof1 = (2.0*nest%dn(2)+nest%dn(3))/(nest%dn(2)+nest%dn(3))*nest%dnw(1)/nest%dn(2)
   cof2 = nest%dn(2)/(nest%dn(2)+nest%dn(3))*nest%dnw(1)/nest%dn(3)

   nest%cf1  = nest%fnp(2)+cof1
   nest%cf2  = nest%fnm(2)-cof1-cof2
   nest%cf3  = cof2

   nest%cfn  = (0.5*nest%dnw(nest%e_vert-1)+nest%dn(nest%e_vert-1))/nest%dn(nest%e_vert-1)
   nest%cfn1 = -0.5*nest%dnw(nest%e_vert-1)/nest%dn(nest%e_vert-1)

END SUBROUTINE vert_cor_vertical_nesting_z_specified

!-----------------------------------------------------------------------------------------

!this is a direct copy of a subroutine that is in ndown, but I couldn't link to the subroutine in ndown because it is compiled after this file
!so a dependecy on ndown will not work. Additionally, ndown is not compiled for ideal cases. The variable is named parent in ndown, but it is actually operating on the nest. It has been renamed to nest here.
      SUBROUTINE vert_cor_vertical_nesting_integer(nest,znw_c,k_dim_c)
         USE module_domain
   IMPLICIT NONE
         TYPE(domain), POINTER :: nest
         integer , intent(in) :: k_dim_c
         real , dimension(k_dim_c), INTENT(IN) :: znw_c

       integer :: kde_c , kde_n ,n_refine,ii,kkk,k
       real :: dznw_m,cof1,cof2
!KAL this subroutine recalculates the vertical coordinates for the nest when vertical nesting is used.  This routine is copied from ndown and allows integer refinement only.

!KAL znw is eta values on full w levels
!KAL everything else is set from znw
!KAL dnw is delta eta on w levels
!KAL rdn is inverse delta eta on w levels
!KAL fnp 

        kde_c = k_dim_c
        kde_n = nest%e_vert
!        n_refine = nest%vert_refine_fact
        n_refine = (kde_n-1)/(kde_c-1)
	
         kkk = 0
         do k = 1 , kde_c-1
         dznw_m = znw_c(k+1) - znw_c(k)
         do ii = 1,n_refine
         kkk = kkk + 1
         nest%znw(kkk) = znw_c(k) + float(ii-1)/float(n_refine)*dznw_m
         enddo
         enddo
         nest%znw(kde_n) = znw_c(kde_c)
         nest%znw(1) = znw_c(1)

      DO k=1, kde_n-1
         nest%dnw(k) = nest%znw(k+1) - nest%znw(k)
         nest%rdnw(k) = 1./nest%dnw(k)
         nest%znu(k) = 0.5*(nest%znw(k+1)+nest%znw(k))
      END DO

   write(*,'(A,I1)') "INTEGER REFINEMENT RATIO OF ",n_refine
   DO k=1,nest%e_vert-1
       IF (k .GT. kde_c) THEN
           write(*,'(A,I2.2,2(A,I3,A,F6.3))') "d",nest%id,"   znw(",k,")=",nest%znw(k),"   znu(",k,")=",nest%znu(k)
       ELSE
           write(*,'(A,I2.2,3(A,I3,A,F6.3))') "d",nest%id,"   znw(",k,")=",nest%znw(k),"   znu(",k,")=",nest%znu(k),"   znw_c(",k,")=",znw_c(k)
       ENDIF
   ENDDO
   write(*,'(A,I2.2,A,I3,A,F6.3)') "d",nest%id,"   znw(",nest%e_vert,")=",nest%znw(nest%e_vert)

      DO k=2, kde_n-1
         nest%dn(k) = 0.5*(nest%dnw(k)+nest%dnw(k-1))
         nest%rdn(k) = 1./nest%dn(k)
         nest%fnp(k) = .5* nest%dnw(k  )/nest%dn(k)
         nest%fnm(k) = .5* nest%dnw(k-1)/nest%dn(k)
      END DO

  cof1 = (2.*nest%dn(2)+nest%dn(3))/(nest%dn(2)+nest%dn(3))*nest%dnw(1)/nest%dn(2)
  cof2 =     nest%dn(2)        /(nest%dn(2)+nest%dn(3))*nest%dnw(1)/nest%dn(3)

      nest%cf1  = nest%fnp(2) + cof1
      nest%cf2  = nest%fnm(2) - cof1 - cof2
      nest%cf3  = cof2

      nest%cfn  = (.5*nest%dnw(kde_n-1)+nest%dn(kde_n-1))/nest%dn(kde_n-1)
      nest%cfn1 = -.5*nest%dnw(kde_n-1)/nest%dn(kde_n-1)
      
      ! the variables dzs and zs are kept from the parent domain.  These are the depths and thickness of the soil layers, which are not included in vertical nesting.

      END SUBROUTINE vert_cor_vertical_nesting_integer

!-----------------------------------------------------------------------------------------

SUBROUTINE vert_cor_vertical_nesting_arbitrary(nest,znw_c,kde_c,use_baseparam_fr_nml)
    USE module_domain
    USE module_model_constants
    IMPLICIT NONE
    TYPE(domain), POINTER :: nest
    INTEGER, INTENT(IN   ) :: kde_c
    REAL, DIMENSION(kde_c), INTENT(IN   ) :: znw_c
    LOGICAL, INTENT(IN   ) :: use_baseparam_fr_nml
    INTEGER :: k, kde_n, ks, ke, id
    REAL :: cof1, cof2
    REAL :: max_dz = 1000
    REAL :: p00, t00, a, tiso, a_strat, p_strat
    INTEGER :: ids, ide, jds, jde, kds, kde, &
               ims, ime, jms, jme, kms, kme, &
               its, ite, jts, jte, kts, kte
    REAL, DIMENSION(max_eta) :: eta_levels
    LOGICAL :: all_fine

    IF ( use_baseparam_fr_nml ) then
      CALL nl_get_base_pres        ( 1 , p00        )
      CALL nl_get_base_temp        ( 1 , t00        )
      CALL nl_get_base_lapse       ( 1 , a          )
      CALL nl_get_iso_temp         ( 1 , tiso       )
      CALL nl_get_base_lapse_strat ( 1 , a_strat    )
      CALL nl_get_base_pres_strat  ( 1 , p_strat    )
      IF ((t00 .LT. 100.0) .OR. (p00 .LT. 10000.0)) THEN
        WRITE(wrf_err_message,*) '--- ERROR: bad base state for T00 or P00 in namelist.input file'
        CALL wrf_error_fatal(TRIM(wrf_err_message))
      END IF
    ELSE
      p00     = nest%p00
      t00     = nest%t00
      a       = nest%tlp
      tiso    = nest%tiso
      a_strat = nest%tlp_strat
      p_strat = nest%p_strat
      IF ((t00 .LT. 100.0) .OR. (p00 .LT. 10000.0)) THEN
        WRITE(wrf_err_message,*) '--- ERROR: did not find base state parameters in nest. Add use_baseparam_fr_nml = .true. in &dynamics and rerun'
        CALL wrf_error_fatal(TRIM(wrf_err_message))
      ENDIF
    ENDIF

    kde_n = nest%e_vert

    !DJW Added code for specifying multiple domains' eta_levels
    IF ((nest%this_is_an_ideal_run) .AND. (model_config_rec%eta_levels(1) .EQ. -1.0)) THEN
      !DJW If we're running an ideal case and do not have levels set in the
      !namelist then we set znw using constant spacing in eta
      DO k=1,kde_n
        CALL wrf_debug(0, "nest_init_utils: eta_levels are not specified in the namelist, setting levels with constant spacing in eta.")
        nest%znw(k) = 1.0-(k-1)/FLOAT((kde_n-1))
      ENDDO
    ELSEIF (.NOT.(nest%this_is_an_ideal_run) .AND. (model_config_rec%eta_levels(1) .EQ. -1.0)) THEN
      write(*,'(A,I2,A)') "--- WARNING: eta_levels are not specified in the namelist for grid_id=",nest%grid_id,", using WRF's default levels."
      CALL get_ijk_from_grid( nest, &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte )
      write(*,'(A,F10.3)') "--- USING: nest%p_top   = ",nest%p_top
      write(*,'(A,F10.3)') "--- USING: g            = ",g
      write(*,'(A,F10.3)') "--- USING: cvpm         = ",cvpm
      write(*,'(A,F10.3)') "--- USING: r_d          = ",r_d
      write(*,'(A,F10.3)') "--- USING: cp           = ",cp
      write(*,'(A,F10.3)') "--- USING: p1000mb      = ",p1000mb
      write(*,'(A,F10.3)') "--- USING: t0           = ",t0
      write(*,'(A,F10.3)') "--- USING: p00          = ",p00
      write(*,'(A,F10.3)') "--- USING: t00          = ",t00
      write(*,'(A,F10.3)') "--- USING: a            = ",a
      write(*,'(A,F10.3)') "--- USING: tiso         = ",tiso
      write(*,'(A,F10.3)') "--- USING: a_strat      = ",a_strat
      write(*,'(A,F10.3)') "--- USING: p_strat      = ",p_strat
      CALL compute_eta ( nest%znw, &
                         eta_levels, max_eta, max_dz, &
                         nest%p_top, g, p00, cvpm, a, r_d, cp, &
                         t00, p1000mb, t0, tiso, p_strat, a_strat, &
                         ids, ide, jds, jde, kds, kde, &
                         ims, ime, jms, jme, kms, kme, &
                         its, ite, jts, jte, kts, kte )
    ELSE
      !DJW If we're in here then we are suppose to read in eta_levels
      !from the namelist. We do so and then check to make sure they make sense
      ks = 1
      ke = ks+model_config_rec%e_vert(1)
      IF (nest%id .NE. 1) THEN
        id = 1
        DO WHILE (nest%id .GT. id)
          id = id+1
          ks = ks+model_config_rec%e_vert(id-1)
          ke = ks+model_config_rec%e_vert(id)-1
        ENDDO
      ENDIF
      DO k=1,kde_n
        nest%znw(k) = model_config_rec%eta_levels(ks+k-1)
        write(*,'(A,I2.2,A,I3,A,F6.3)') "d",nest%id," eta_levels(",k,") = ",nest%znw(k)
      ENDDO
      !Check the value of the first and last eta level for our domain,
      !then check that the vector of eta levels is only decreasing
      all_fine = .TRUE.
      IF (nest%znw(1) .NE. 1.0) THEN
        write(wrf_err_message,'(A,I2,A)') "--- ERROR: first eta_level for grid_id=",nest%grid_id," is not 1.0. Check eta_levels specified in namelist."
        all_fine = .FALSE.
      ENDIF
      write(*,'(A,F10.3)') "--- USING: g            = ",g
      write(*,'(A,F10.3)') "--- USING: cvpm         = ",cvpm
      write(*,'(A,F10.3)') "--- USING: r_d          = ",r_d
      write(*,'(A,F10.3)') "--- USING: cp           = ",cp
      write(*,'(A,F10.3)') "--- USING: p1000mb      = ",p1000mb
      write(*,'(A,F10.3)') "--- USING: t0           = ",t0
      IF (nest%znw(kde_n) .NE. 0.0) THEN
        write(wrf_err_message,'(A,I2,A)') "--- ERROR: last eta_level for grid_id=",nest%grid_id," is not 0.0. Check eta_levels specified in namelist."
        all_fine = .FALSE.
      ENDIF
      DO k=2,kde_n
        IF (nest%znw(k) .GT. nest%znw(k-1)) THEN
          write(wrf_err_message,'(A,I2,A)') "--- ERROR: eta_levels for grid_id=",nest%grid_id," are not monotonically decreasing. Check eta_levels specified in namelist."
          all_fine = .FALSE.
        ENDIF
      ENDDO
      IF (.NOT. all_fine) THEN
        write(*,'(A)') "--- ERROR: check your namelist and compare to below,"
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
        CALL wrf_error_fatal( wrf_err_message )
      ENDIF
    ENDIF
    !DJW End of added code for specifying eta_levels

    DO k=1,kde_n-1
      nest%dnw(k) = nest%znw(k+1)-nest%znw(k)
      nest%rdnw(k) = 1./nest%dnw(k)
      nest%znu(k) = 0.5*(nest%znw(k+1)+nest%znw(k))
    ENDDO
    nest%znu(kde_n) = 0.0

    DO k=2,kde_n-1
      nest%dn(k) = 0.5*(nest%dnw(k)+nest%dnw(k-1))
      nest%rdn(k) = 1./nest%dn(k)
      nest%fnp(k) = .5* nest%dnw(k  )/nest%dn(k)
      nest%fnm(k) = .5* nest%dnw(k-1)/nest%dn(k)
    ENDDO

    cof1 = (2.*nest%dn(2)+nest%dn(3))/(nest%dn(2)+nest%dn(3))*nest%dnw(1)/nest%dn(2)
    cof2 =     nest%dn(2)            /(nest%dn(2)+nest%dn(3))*nest%dnw(1)/nest%dn(3)

    nest%cf1  = nest%fnp(2) + cof1
    nest%cf2  = nest%fnm(2) - cof1 - cof2
    nest%cf3  = cof2
    nest%cfn  = (.5*nest%dnw(kde_n-1)+nest%dn(kde_n-1))/nest%dn(kde_n-1)
    nest%cfn1 = -.5*nest%dnw(kde_n-1)/nest%dn(kde_n-1)

END SUBROUTINE vert_cor_vertical_nesting_arbitrary

!-----------------------------------------------------------------------------------------

SUBROUTINE vert_cor_vertical_nesting_z_specified_real( nest, use_baseparam_fr_nml )
    USE module_domain
    USE module_model_constants
    IMPLICIT NONE
    TYPE(domain), POINTER          :: nest
    LOGICAL, INTENT(IN   )         :: use_baseparam_fr_nml
    !local
    REAL                           :: p00, t00, a, tiso, a_strat,   &
                                      p_strat, cof1, cof2, p_surf
    INTEGER                        :: ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte, &
                                      i, k, ks, nz_specified
    REAL, DIMENSION(1:nest%e_vert) :: z_specified

    CALL wrf_debug(100, "dyn_em/nest_init_utils/subroutine vert_cor_vertical_nesting_z_specified_real")

    write(*,'(A,I1)') "DJW[nest_init_utils/vert_cor_vertical_nesting_z_specified_real]: at top, nest%id=",nest%grid_id

    IF ( use_baseparam_fr_nml ) then
      CALL nl_get_base_pres        ( 1 , p00        )
      CALL nl_get_base_pres_strat  ( 1 , p_strat    )
      CALL nl_get_iso_temp         ( 1 , tiso       )
      CALL nl_get_base_temp        ( 1 , t00        )
      CALL nl_get_base_lapse       ( 1 , a          )
      CALL nl_get_base_lapse_strat ( 1 , a_strat    )
      IF ((t00 .LT. 100.0) .OR. (p00 .LT. 10000.0)) THEN
        WRITE(wrf_err_message,*) '--- ERROR: bad base state for T00 or P00 in namelist.input file'
        CALL wrf_error_fatal(TRIM(wrf_err_message))
      END IF
    ELSE
      p00     = nest%p00
      p_strat = nest%p_strat
      tiso    = nest%tiso
      t00     = nest%t00
      a       = nest%tlp
      a_strat = nest%tlp_strat
      IF ((t00 .LT. 100.0) .OR. (p00 .LT. 10000.0)) THEN
        WRITE(wrf_err_message,*) '--- ERROR: did not find base state parameters in nest. Add use_baseparam_fr_nml = .true. in &dynamics and rerun'
        CALL wrf_error_fatal(TRIM(wrf_err_message))
      ENDIF
    ENDIF

    p_surf = p00*EXP(-t00/a+((t00/a)**2-2.0*g*375.0/a/r_d)**0.5)

    CALL get_ijk_from_grid( nest, &
                            ids, ide, jds, jde, kds, kde, &
                            ims, ime, jms, jme, kms, kme, &
                            its, ite, jts, jte, kts, kte )

    !DJW find the number of z_specified entries that will precede the
    !    section associated with our domain of interest.
    ks = 1
    i = 1
    DO WHILE (i .LE. model_config_rec%max_dom)
       IF (model_config_rec%grid_id(i) .EQ. nest%grid_id) THEN
          i = model_config_rec%max_dom+1
       ELSE
          ks = ks+model_config_rec%e_vert(i)
          i = i+1
       ENDIF
    ENDDO

    IF (model_config_rec%z_specified(ks) .NE. 0.0) THEN
       write(wrf_err_message,'(A,I3,A,F9.3,A)') "---ERROR: z_specified(",ks,")=",model_config_rec%z_specified(ks)," but it must be zero."
       CALL wrf_error_fatal(TRIM(wrf_err_message))
    ENDIF
    IF (model_config_rec%z_specified(ks+nest%e_vert-1) .NE. -1.0) THEN
       write(wrf_err_message,'(A,I3,A,F9.3,A)') "---ERROR: z_specified(",ks+nest%e_vert,")=",model_config_rec%z_specified(ks+nest%e_vert)," but it shouldn't be specified."
       CALL wrf_error_fatal(TRIM(wrf_err_message))
    ENDIF

    k = ks
    DO WHILE (k .LT. ks+nest%e_vert-1)
       k = k+1
       IF ( (model_config_rec%z_specified(k)   .EQ. -1.0) .AND. &
            (model_config_rec%z_specified(k-1) .NE. -1.0) ) THEN
          nz_specified = k-ks
       ENDIF
    ENDDO

    DO k=ks+1,ks+nest%e_vert
       IF (k-ks+1 .LE. nz_specified) THEN
          IF (model_config_rec%z_specified(k) .LE. model_config_rec%z_specified(k-1)) THEN
             write(wrf_err_message,'(A)') "---ERROR: z_specified is not monotonically increasing."
             CALL wrf_error_fatal(TRIM(wrf_err_message))
          ENDIF
       ENDIF
    ENDDO

    DO k=ks,ks+nest%e_vert-1
       z_specified(k-ks+1) = model_config_rec%z_specified(k)
    ENDDO

    DO k=kds,kde
       write(*,'(A,I1,A,I3,A,F7.1)') "d0",nest%grid_id," z_specified(",k,") = ",z_specified(k)
    ENDDO

    CALL compute_eta_z_specified ( nest%znw,                         &
                                   z_specified(1:nz_specified),      &
                                   nz_specified,                     &
                                   p_surf, p00, nest%p_top, p_strat, &
                                   tiso, t00, t0,                    &
                                   a, a_strat,                       &
                                   r_d, cp, p1000mb, cvpm, g,        &
                                   nest%grid_id,                     &
                                   ids, ide, jds, jde, kds, kde,     &
                                   ims, ime, jms, jme, kms, kme,     &
                                   its, ite, jts, jte, kts, kte )

    DO k=kds,kde
       IF (k .LT. kde-1) THEN
          write(*,'(A,I1,A,I3,A,F6.4,A,I3,A,F8.5,A,I3,A,I3,A,F6.4)') "d0",nest%grid_id," znw(",k,") = ",nest%znw(k),"   dnw(",k,")=",nest%znw(k+1)-nest%znw(k),"   dnw(",k+1,")/dnw(",k,")=",(nest%znw(k+2)-nest%znw(k+1))/(nest%znw(k+1)-nest%znw(k))
       ELSEIF (k .LT. kde) THEN
          write(*,'(A,I1,A,I3,A,F6.4,A,I3,A,F8.5)') "d0",nest%grid_id," znw(",k,") = ",nest%znw(k),"   dnw(",k,")=",nest%znw(k+1)-nest%znw(k)
       ELSE
          write(*,'(A,I1,A,I3,A,F6.4)') "d0",nest%grid_id," znw(",k,") = ",nest%znw(k)
       ENDIF
    ENDDO

    DO k=1,kde-1
      nest%dnw(k) = nest%znw(k+1)-nest%znw(k)
      nest%rdnw(k) = 1.0/nest%dnw(k)
      nest%znu(k) = 0.5*(nest%znw(k+1)+nest%znw(k))
    ENDDO
    nest%znu(kde) = 0.0

    DO k=2,kde-1
      nest%dn(k) = 0.5*(nest%dnw(k)+nest%dnw(k-1))
      nest%rdn(k) = 1.0/nest%dn(k)
      nest%fnp(k) = 0.5* nest%dnw(k  )/nest%dn(k)
      nest%fnm(k) = 0.5* nest%dnw(k-1)/nest%dn(k)
    ENDDO

    cof1 = (2.0*nest%dn(2)+nest%dn(3))/(nest%dn(2)+nest%dn(3))*nest%dnw(1)/nest%dn(2)
    cof2 =      nest%dn(2)            /(nest%dn(2)+nest%dn(3))*nest%dnw(1)/nest%dn(3)

    nest%cf1  = nest%fnp(2) + cof1
    nest%cf2  = nest%fnm(2) - cof1 - cof2
    nest%cf3  = cof2
    nest%cfn  = (0.5*nest%dnw(kde-1)+nest%dn(kde-1))/nest%dn(kde-1)
    nest%cfn1 = -0.5*nest%dnw(kde-1)/nest%dn(kde-1)

    write(*,'(A,I1)') "DJW[nest_init_utils/vert_cor_vertical_nesting_z_specified_real]: at bot, nest%id=",nest%grid_id

END SUBROUTINE vert_cor_vertical_nesting_z_specified_real

!-----------------------------------------------------------------------------------------

SUBROUTINE compute_eta_z_specified ( znw,                          &
                                     z_specified,                  &
                                     nz_specified,                 &
                                     p_surf, p00, p_top, p_strat,  &
                                     tiso, t00, t0,                &
                                     A, A_strat,                   &
                                     r_d, cp, p1000mb, cvpm, g,    &
                                     grid_id,                      &
                                     ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme, &
                                     its, ite, jts, jte, kts, kte )
    IMPLICIT NONE
    REAL, DIMENSION(kds:kde), INTENT(  OUT)          :: znw
    REAL, DIMENSION(kds:nz_specified), INTENT(IN   ) :: z_specified
    REAL, INTENT(IN   )                              :: p_surf,    &
                                                        p00,       &
                                                        p_top,     &
                                                        p_strat,   &
                                                        tiso,      &
                                                        t00,       &
                                                        t0,        &
                                                        A,         &
                                                        A_strat,   &
                                                        r_d,       &
                                                        cp,        &
                                                        p1000mb,   &
                                                        cvpm,      &
                                                        g
    INTEGER, INTENT(IN)                              :: nz_specified, &
                                                        grid_id,   &
                                                        ids, ide, jds, jde, kds, kde, &
                                                        ims, ime, jms, jme, kms, kme, &
                                                        its, ite, jts, jte, kts, kte
    !local
    REAL, DIMENSION(kds:kde-1)                      :: dz
    REAL                                            :: pb,         &
                                                       temp,       &
                                                       t_init,     &
                                                       stretch,    &
                                                       alb_max,    &
                                                       t_init_max, &
                                                       pb_max,     &
                                                       ztop
    INTEGER                                         :: k,          &
                                                       kk,         &
                                                       loop,       &
                                                       loop1
    CHARACTER(LEN=256)                              :: wrf_err_message
    REAL, DIMENSION(kds:kde)                        :: phb, alb
    INTEGER, PARAMETER                              :: prac_levels = 59
    REAL(KIND=8), DIMENSION(prac_levels)            :: znw_prac
    LOGICAL                                         :: debug

    CALL wrf_debug(100, "dyn_em/nest_init_utils/subroutine compute_eta_z_specified")
    debug = .TRUE.
!    debug = .FALSE.
    IF (debug) THEN
       write(wrf_err_message,'(A,I1)') "At the top of compute_eta_z_specified for grid_id=",grid_id
       CALL wrf_message(TRIM(wrf_err_message))
    ENDIF

    IF (nz_specified .GE. kde) THEN
       write(wrf_err_message, '(2(A,I3),A,I1)') "--- ERROR: nz_specified (",nz_specified,") >= kde (",kde,") for grid_id=",grid_id
       CALL wrf_error_fatal(wrf_err_message)
    ENDIF

    IF (nz_specified .LE. 1) THEN
       write(wrf_err_message, '(A,I3,A,I1)') "--- ERROR: nz_specified (",nz_specified,") <= 1 for grid_id=",grid_id
       CALL wrf_error_fatal(wrf_err_message)
    ENDIF

    DO k=kds,nz_specified-1
       dz(k) = z_specified(k+1)-z_specified(k)
    ENDDO

    !Start with equally spaced eta-levels
    DO k=kds,nz_specified-1
       znw(k) = 1.0-(k-1.0)/(kde-1.0)
    ENDDO

    !Iteratively set eta-levels closest to the ground using heights
    !specified in the namelist variable 'z_specified'
    DO loop=1,10
       DO k=kds,nz_specified-1
          pb = 0.5*(znw(k)+znw(k+1))*(p_surf-p_top)+p_top
          IF (pb .LT. p_strat) THEN
             temp = tiso+A_strat*LOG(pb/p_strat)
          ELSE
             temp = MAX(tiso, t00+A*LOG(pb/p00))
          ENDIF
          t_init = temp*(p00/pb)**(r_d/cp)-t0
          alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
          znw(k+1) = znw(k)-dz(k)*g/((p_surf-p_top)*alb(k))
       ENDDO !k-loop
    ENDDO !loop

    IF (.FALSE.) THEN
       !Set the remaining eta-levels so that there is a constant increase in the
       !width of eta levels. In other words znw(k+1)-znw(k) = X * (znw(k)-znw(k-1))
       !where X is some stretching constant that we are solving for below. There is
       !probably a better way to get this constant, but I'm solving the problem
       !iteratively.
       stretch = 1.0 !initial guess
       loop = 0 !track iterations so we don't get stuck
       znw(kde) = 1.0 !this is simply for the DO WHILE conditonal below
       DO WHILE (ABS(znw(kde)) .GT. 0.001)
          loop = loop+1
          !recalculate eta-levels using the most recent value of 'stretch'
          DO k=nz_specified,kde-1
             znw(k+1) = znw(k)-(znw(k-1)-znw(k))*stretch
          ENDDO !k-loop
          !The following will gradually increase the stretching between eta
          !levels, pushing znw(kde) gradually closer to zero. Once znw(kde)
          !overshoots and turns negative, it is either going to be close enough to
          !zero that we are happy with the result and move on, or we will slowly
          !decrease the stretching between levels at a slower rate than we first
          !approached zero using.
          IF (znw(kde) .LT. 0.0) THEN
             stretch = stretch+0.00001*znw(kde)/ABS(znw(kde))
          ELSE
             stretch = stretch+0.0001*znw(kde)/ABS(znw(kde))
          ENDIF
          !Throw an error if we loop too many times.
          IF (loop .GE. 10000) THEN
             write(wrf_err_message, '(A,I7,A,I3,A,F5.3)') "--- ERROR: reached loop limit (",loop,") in compute_eta_z_specified with znw(",kde,")=",znw(kde)
             CALL wrf_error_fatal(wrf_err_message)
          ENDIF
       ENDDO
       write(*,'(A,I5,A,I3,A,F6.4,A,F7.4)') "DJW[nest_init_utils/compute_eta_z_specified]: it took ",loop," iterations for znw(",kde,") to reach ",znw(kde)," with stretch=",stretch
       !After iterating above, znw(kde) will be near zero but it needs to be
       !exactly zero.
       znw(kde) = 0.0
    ELSE
       !Set the remaining eta-levels so that there is a constant increase in the
       !physical width as we move upwards. In other words, dz(k+1) = A * dz(k)
       !where A is some constant that we are solving for. To do this, we first
       !start with an initial guess for a ztop and then iterate to find
       !eta-levels that give us uniform stretching. Then, since our ztop is just
       !an estimate, our top eta-level will not be zero, so we iterate the last
       !steps until we settle on a ztop that results in znw(kde)~=0.0.
       !First guess a ztop based upon the "practical levels" that WRF uses.
       znw_prac = (/ 1.0000_8, 0.9930_8,  0.9830_8, 0.9700_8, 0.9540_8, 0.9340_8, 0.9090_8, 0.8800_8, &
                     0.8500_8, 0.8000_8,  0.7500_8, 0.7000_8, 0.6500_8, 0.6000_8, 0.5500_8, 0.5000_8, &
                     0.4500_8, 0.4000_8,  0.3500_8, 0.3000_8, 0.2500_8, 0.2000_8, 0.1500_8, 0.1000_8, &
                     0.0800_8, 0.0600_8,  0.0400_8, 0.0200_8, 0.0150_8, 0.0100_8, 0.0090_8, 0.0080_8, &
                     0.0070_8, 0.0060_8,  0.0050_8, 0.0040_8, 0.0035_8, 0.0030_8, 0.0028_8, 0.0026_8, &
                     0.0024_8, 0.0022_8,  0.0020_8, 0.0018_8, 0.0016_8, 0.0014_8, 0.0012_8, 0.0010_8, &
                     0.0009_8, 0.0008_8,  0.0007_8, 0.0006_8, 0.0005_8, 0.0004_8, 0.0003_8, 0.0002_8, &
                     0.0001_8, 0.00005_8, 0.0000_8 /)
       DO k=1,prac_levels-1
          pb = (znw_prac(k)+znw_prac(k+1))*0.5*(p_surf-p_top)+p_top
          temp = MAX(tiso, t00+A*LOG(pb/p00))
          IF (pb .LT. p_strat) THEN
             temp = tiso+A_strat*LOG(pb/p_strat)
          ENDIF
          t_init = temp*(p00/pb)**(r_d/cp)-t0
          alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
       ENDDO
       phb(1) = 0.0
       DO k=2,prac_levels
          phb(k) = phb(k-1)-(znw_prac(k)-znw_prac(k-1))*(p00-p_top)*alb(k-1)
       ENDDO
       ztop = phb(prac_levels)/g
       !Now that we have an estimate for ztop, estimate the height of the
       !specified levels.
       DO k=1,prac_levels-1
          pb = (znw_prac(k)+znw_prac(k+1))*0.5*(p_surf-p_top)+p_top
          temp = MAX(tiso, t00+A*LOG(pb/p00))
          IF (pb .LT. p_strat) THEN
             temp = tiso+A_strat*LOG(pb/p_strat)
          ENDIF
          t_init = temp*(p00/pb)**(r_d/cp)-t0
          alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
       ENDDO
       phb(1) = 0.0
       DO k=2,nz_specified
          phb(k) = phb(k-1)-(znw(k)-znw(k-1))*(p00-p_top)*alb(k-1)
       ENDDO
       !Next we loop until we refine our ztop so that znw(kde) approaches zero.
       znw(kde) = 1.0 !this is for the do-while conditional
       loop1 = 0 !a counter to keep track of our loops
       stretch = 1.0 !the coefficient we are iterating to find
       DO WHILE (ABS(znw(kde)) .GT. 0.001)
          loop1 = loop1+1
          !Loop until we find the correct stretching coefficient. When we are
          !finished looping, all of the non-specified levels should be gradually
          !stretching in dz by an equal amount as we move up the column.
          phb(kde) = g*ztop*2.0 !this is simply for the do-whle conditional
          loop = 0 !a counter to keep track of our loop (just in case...)
          DO WHILE (ABS(phb(kde)/g-ztop)/ztop .GT. 0.0001)
             loop = loop+1
             IF (debug) write(*,'(A,I5,A,F8.6)',ADVANCE='NO') "   loop=",loop,"   stretch=",stretch
             DO k=nz_specified,kde-1
                find_prac : DO kk=1,prac_levels
                   IF (znw_prac(kk) .LT. znw(k)) THEN
                      EXIT find_prac
                   ENDIF
                ENDDO find_prac
                pb = 0.5*(znw(k)+znw_prac(kk))*(p00-p_top)+p_top
                IF (pb .LT. p_strat) THEN
                   temp = tiso+A_strat*LOG(pb/p_strat)
                ELSE
                   temp = MAX(tiso, t00+A*LOG(pb/p00))
                ENDIF
                t_init = temp*(p00/pb)**(r_d/cp)-t0
                alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
                znw(k+1) = znw(k)-dz(k-1)*stretch*g/((p00-p_top)*alb(k))
                pb = 0.5*(znw(k)+znw(k+1))*(p00-p_top)+p_top
                IF (pb .LT. p_strat) THEN
                   temp = tiso+A_strat*LOG(pb/p_strat)
                ELSE
                   temp = MAX(tiso, t00+A*LOG(pb/p00))
                ENDIF
                t_init = temp*(p00/pb)**(r_d/cp)-t0
                alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
                znw(k+1) = znw(k)-dz(k-1)*stretch*g/((p00-p_top)*alb(k))
                pb = 0.5*(znw(k)+znw(k+1))*(p00-p_top)+p_top
                phb(k+1) = phb(k)-(znw(k+1)-znw(k))*(p00-p_top)*alb(k)
                dz(k) = (phb(k+1)-phb(k))/g
             ENDDO
             IF (debug) write(*,'(A,F8.1,A,F8.1)') "   ztop=",ztop,"   phb(kde)/g=",phb(kde)/g
             IF (phb(kde)/g .LT. ztop) THEN
                stretch = stretch+0.00333*SQRT(ABS(phb(kde)/g-ztop)/ztop)
             ELSE
                stretch = stretch-0.000001
             ENDIF
             !Throw an error if we loop too many times.
             IF (loop .GE. 10000) THEN
                write(wrf_err_message, '(A,I7,A,I3,2(A,F7.1))') "--- ERROR: reached loop limit (",loop,") in compute_eta_z_specified with phb(",kde,")/g=",phb(kde)/g," and ztop=",ztop
                CALL wrf_error_fatal(wrf_err_message)
             ENDIF
          ENDDO !loop
          IF (debug) write(*,'(2(A,I3),A,F8.5,A,F7.1)') "loop1=",loop1,"   znw(",kde,")=",znw(kde),"   ztop=",ztop
          IF (znw(kde) .GT. 0.0) THEN
             ztop = ztop+1000.0*SQRT( znw(kde))
          ELSE
             ztop = ztop-1000.0*SQRT(-znw(kde))
          ENDIF
          !Throw an error if we loop too many times. If you're converging on a
          !solution but getting this error, then up the threshhold in the
          !conditional below.
          IF (loop1 .GE. 100) THEN
             write(wrf_err_message, '(A,I3,A,I3,A,F8.6)') "--- ERROR: reached loop limit (",loop1,") in compute_eta_z_specified with znw(",kde,")=",znw(kde)
             CALL wrf_error_fatal(TRIM(wrf_err_message))
          ENDIF
       ENDDO !loop1
       write(wrf_err_message,'(A,I3,A,F6.3,A,I3,A,F5.3,A,I5,A)') "Converged with znw(",kde,")=",znw(kde)," after ",loop1," outer-loops and stretching coefficient = ",stretch," after ",loop," inner-loops."
       CALL wrf_message(TRIM(wrf_err_message))
       znw(kde) = 0.0
    ENDIF

    IF (debug) THEN
       !Calculate and write some information that is helpful for debugging
       DO k=kds,kde-1
          pb = 0.5*(znw(k)+znw(k+1))*(p_surf-p_top)+p_top
          IF (pb .LT. p_strat) THEN
             temp = tiso+A_strat*LOG(pb/p_strat)
          ELSE
             temp = MAX(tiso, t00+A*LOG(pb/p00))
          ENDIF
          t_init = temp*(p00/pb)**(r_d/cp)-t0
          alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
       ENDDO !k-loop
       phb(kds) = 0.0
       DO k=kds+1,kde
          phb(k) = phb(k-1)-(znw(k)-znw(k-1))*(p_surf-p_top)*alb(k-1)
       ENDDO
       DO k=kds,kde
          write(*,'(A,I3,A,F6.4,A,I3,A,F8.1)',ADVANCE='NO') "DJW[nest_init_utils/compute_eta_z_specified]: znw(",k,")=",znw(k),"   phb(",k,")/g=",phb(k)/g
          IF (k .LE. nz_specified) THEN
             write(*,'(A,I3,A,F6.1)',ADVANCE='NO') "   z_specified(",k,")=",z_specified(k)
          ELSE
             write(*,'(A)',ADVANCE='NO') "                          "
          ENDIF
          IF (k .GT. kds+1) THEN
             write(*,'(A,F6.4,A,F6.1,A,F6.4)') "   dznw=",znw(k-1)-znw(k),"   dz=",(phb(k)-phb(k-1))/g,"   dz(k)/dz(k-1)=",(phb(k)-phb(k-1))/(phb(k-1)-phb(k-2))
          ELSEIF (k .GT. kds) THEN
             write(*,'(A,F6.4,A,F6.1)') "   dznw=",znw(k-1)-znw(k),"   dz=",(phb(k)-phb(k-1))/g
          ELSE
             write(*,'(A)') ""
          ENDIF
       ENDDO
    ENDIF

END SUBROUTINE compute_eta_z_specified

!-----------------------------------------------------------------------------------------

SUBROUTINE compute_eta ( znw , &
                         eta_levels , max_eta , max_dz ,        &
                         p_top_def , g_def , p00_def ,          &
                         cvpm_def , a_def , r_d_def , cp_def ,  &
                         t00_def , p1000mb_def , t0_def ,       &
                         tiso_def , p_strat_def , a_strat_def , &
                         ids , ide , jds , jde , kds , kde ,    &
                         ims , ime , jms , jme , kms , kme ,    &
                         its , ite , jts , jte , kts , kte )

      !  Compute eta levels, either using given values from the namelist (hardly
      !  a computation, yep, I know), or assuming a constant dz above the PBL,
      !  knowing p_top and the number of eta levels.

      IMPLICIT NONE

      INTEGER , INTENT(IN)        :: ids , ide , jds , jde , kds , kde , &
                                     ims , ime , jms , jme , kms , kme , &
                                     its , ite , jts , jte , kts , kte
      REAL , INTENT(IN)           :: max_dz
      REAL , INTENT(IN)           :: p_top_def , g_def , p00_def , cvpm_def , &
                                     a_def , r_d_def , cp_def , t00_def ,     &
                                     p1000mb_def , t0_def , tiso_def ,        &
                                     p_strat_def , a_strat_def
      INTEGER , INTENT(IN)        :: max_eta
      REAL , DIMENSION (max_eta)  :: eta_levels

      REAL , DIMENSION (kts:kte) , INTENT(OUT) :: znw

      !  Local vars

      INTEGER :: k , kk
      REAL(KIND=8) :: mub , t_init , p_surf , pb, ztop, ztop_pbl , dz , temp
      REAL(KIND=8) , DIMENSION(kts:kte) :: dnw
      REAL(KIND=8) :: p_top , g , p00 , cvpm , &
                      a , r_d , cp , t00 ,     &
                      p1000mb , t0 , tiso ,        &
                      p_strat , a_strat

      INTEGER , PARAMETER :: prac_levels = 59
      INTEGER :: loop , loop1
      REAL(KIND=8) , DIMENSION(prac_levels) :: znw_prac , znu_prac , dnw_prac
      REAL(KIND=8) , DIMENSION(MAX(prac_levels,kde)) :: alb , phb
      REAL(KIND=8) :: alb_max, t_init_max, pb_max, phb_max

      CHARACTER(LEN=256) :: message

      !  Compute top of the atmosphere with some silly levels.  We just want to
      !  integrate to get a reasonable value for ztop.  We use the planned
      !  PBL-esque
      !  levels, and then just coarse resolution above that.  We know p_top,
      !  and we
      !  have the base state vars.

      p_top   = p_top_def
      g       = g_def
      p00     = p00_def
      cvpm    = cvpm_def
      a       = a_def
      r_d     = r_d_def
      cp      = cp_def
      t00     = t00_def
      p1000mb = p1000mb_def
      t0      = t0_def
      tiso    = tiso_def
      p_strat = p_strat_def
      a_strat = a_strat_def

      p_surf = p00

      znw_prac = (/ 1.0000_8 , 0.9930_8 , 0.9830_8 , 0.9700_8 , 0.9540_8 , 0.9340_8 , 0.9090_8 , 0.8800_8 , &
                    0.8500_8 , 0.8000_8 , 0.7500_8 , 0.7000_8 , 0.6500_8 , 0.6000_8 , 0.5500_8 , 0.5000_8 , &
                    0.4500_8 , 0.4000_8 , 0.3500_8 , 0.3000_8 , 0.2500_8 , 0.2000_8 , 0.1500_8 , 0.1000_8 , &
                    0.0800_8 , 0.0600_8 , 0.0400_8 , 0.0200_8 , &
                    0.0150_8 , 0.0100_8 , 0.0090_8 , 0.0080_8 , 0.0070_8 , 0.0060_8 , 0.0050_8 , 0.0040_8 , &
                    0.0035_8 , 0.0030_8 , &
                    0.0028_8 , 0.0026_8 , 0.0024_8 , 0.0022_8 , 0.0020_8 , &
                    0.0018_8 , 0.0016_8 , 0.0014_8 , 0.0012_8 , 0.0010_8 , &
                    0.0009_8 , 0.0008_8 , 0.0007_8 , 0.0006_8 , 0.0005_8 , 0.0004_8 , 0.0003_8 , &
                    0.0002_8 , 0.0001_8 , 0.00005_8, 0.0000_8 /)

      DO k = 1 , prac_levels - 1
        znu_prac(k) = ( znw_prac(k) + znw_prac(k+1) ) * 0.5_8
        dnw_prac(k) = znw_prac(k+1) - znw_prac(k)
      END DO

      DO k = 1, prac_levels-1
        pb = znu_prac(k)*(p_surf - p_top) + p_top
        temp = MAX ( tiso, t00 + A*LOG(pb/p00) )
        IF ( pb .LT. p_strat ) THEN
          temp = tiso + A_strat*LOG(pb/p_strat)
        END IF
        t_init = temp*(p00/pb)**(r_d/cp) - t0
        alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
      END DO

      !  Base state mu is defined as base state surface pressure minus p_top

      mub = p_surf - p_top

      !  Integrate base geopotential, starting at terrain elevation.

      phb(1) = 0._8
      DO k  = 2,prac_levels
        phb(k) = phb(k-1) - dnw_prac(k-1)*mub*alb(k-1)
      END DO

      !  So, now we know the model top in meters.  Get the average depth above the PBL
      !  of each of the remaining levels.  We are going for a constant delta z thickness.

      ztop     = phb(prac_levels) / g
      ztop_pbl = phb(8          ) / g
      dz = ( ztop - ztop_pbl ) / REAL ( kde - 8 )

      IF ( dz .GE. max_dz ) THEN
        WRITE (message,FMT='("With a requested ",F7.1," Pa model top, the model lid will be about ",F7.1," m.")') p_top, ztop
        CALL wrf_message ( message )
        WRITE (message,FMT='("With ",I3," levels above the PBL, the level thickness will be about ",F6.1," m.")') kde-8, dz
        CALL wrf_message ( message )
        WRITE (message,FMT='("Thicknesses greater than ",F7.1," m are not recommended.")') max_dz
        CALL wrf_message ( message )
        CALL wrf_error_fatal ( 'Add more levels to namelist.input for e_vert' )
      END IF

      !  Standard levels near the surface so no one gets in trouble.

      DO k = 1 , 8
        eta_levels(k) = znw_prac(k)
      END DO

      !  Using d phb(k)/ d eta(k) = -mub * alb(k), eqn 2.9
      !  Skamarock et al, NCAR TN 468.  Use full levels, so
      !  use twice the thickness.

      DO k = 8, kte-1-2

        find_prac : DO kk = 1 , prac_levels
          IF (znw_prac(kk) .LT. eta_levels(k) ) THEN
            EXIT find_prac
          END IF
        end do find_prac

        pb = 0.5*(eta_levels(k)+znw_prac(kk)) * (p_surf - p_top) + p_top

        temp = MAX ( tiso, t00 + A*LOG(pb/p00) )
        IF ( pb .LT. p_strat ) THEN
          temp = tiso + A_strat * LOG ( pb/p_strat )
        END IF
!       temp =             t00 + A*LOG(pb/p00)
        t_init = temp*(p00/pb)**(r_d/cp) - t0
        alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
        eta_levels(k+1) = eta_levels(k) - dz*g / ( mub*alb(k) )
        pb = 0.5*(eta_levels(k)+eta_levels(k+1)) * (p_surf - p_top) + p_top

        temp = MAX ( tiso, t00 + A*LOG(pb/p00) )
        IF ( pb .LT. p_strat ) THEN
          temp = tiso + A_strat * LOG ( pb/p_strat )
        END IF
!       temp =             t00 + A*LOG(pb/p00)
        t_init = temp*(p00/pb)**(r_d/cp) - t0
        alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
        eta_levels(k+1) = eta_levels(k) - dz*g / ( mub*alb(k) )
        pb = 0.5*(eta_levels(k)+eta_levels(k+1)) * (p_surf - p_top) + p_top

        phb(k+1) = phb(k) - (eta_levels(k+1)-eta_levels(k)) * mub*alb(k)
      END DO

      alb_max = alb(kte-1-2)
      t_init_max = t_init
      pb_max = pb
      phb_max = phb(kte-1)

      DO k = 1 , kte-1-2
        znw(k) = eta_levels(k)
      END DO
      znw(kte-2) = 0.000

      !  There is some iteration.  We want the top level, ztop, to be
      !  consistent with the delta z, and we want the half level values
      !  to be consistent with the eta levels.  The inner loop to 10 gets
      !  the eta levels very accurately, but has a residual at the top, due
      !  to dz changing.  We reset dz five times, and then things seem OK.

      DO loop1 = 1 , 5
        DO loop = 1 , 10
          DO k = 8, kte-1-2-1
            pb = (znw(k)+znw(k+1))*0.5 * (p_surf - p_top) + p_top
            temp = MAX ( tiso, t00 + A*LOG(pb/p00) )
            IF ( pb .LT. p_strat ) THEN
              temp = tiso + A_strat * LOG ( pb/p_strat )
            END IF
            t_init = temp*(p00/pb)**(r_d/cp) - t0
            alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
            znw(k+1) = znw(k) - dz*g / ( mub*alb(k) )
          END DO
          pb = pb_max
          t_init = t_init_max
          alb(kte-1-2) = alb_max
          znw(kte-2) = znw(kte-1-2) - dz*g / ( mub*alb(kte-1-2) )
          IF ( ( loop1 .EQ. 5 ) .AND. ( loop .EQ. 10 ) ) THEN
            print *,'Converged znw(kte) should be about 0.0 = ',znw(kte-2)
          END IF
          znw(kte-2) = 0.000
        END DO

        !  Here is where we check the eta levels values we just computed.

        DO k = 1, kde-1-2
          pb = (znw(k)+znw(k+1))*0.5 * (p_surf - p_top) + p_top
          temp = MAX ( tiso, t00 + A*LOG(pb/p00) )
          IF ( pb .LT. p_strat ) THEN
            temp = tiso + A_strat * LOG ( pb/p_strat )
          END IF
          t_init = temp*(p00/pb)**(r_d/cp) - t0
          alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
        END DO

        phb(1) = 0.
        DO k  = 2,kde-2
          phb(k) = phb(k-1) - (znw(k)-znw(k-1)) * mub*alb(k-1)
        END DO

        !  Reset the model top and the dz, and iterate.

        ztop = phb(kde-2)/g
        ztop_pbl = phb(8)/g
        dz = ( ztop - ztop_pbl ) / REAL ( (kde-2) - 8 )
      END DO

      IF ( dz .GT. max_dz ) THEN
        print *,'z (m)            = ',phb(1)/g
        do k = 2 ,kte-2
          print *,'z (m) and dz (m) = ',phb(k)/g,(phb(k)-phb(k-1))/g
        end do
        print *,'dz (m) above fixed eta levels = ',dz
        print *,'namelist max_dz (m) = ',max_dz
        print *,'namelist p_top (Pa) = ',p_top
        CALL wrf_debug ( 0, 'You need one of three things:' )
        CALL wrf_debug ( 0, '1) More eta levels to reduce the dz: e_vert' )
        CALL wrf_debug ( 0, '2) A lower p_top so your total height is reduced: p_top_requested')
        CALL wrf_debug ( 0, '3) Increase the maximum allowable eta thickness: max_dz')
        CALL wrf_debug ( 0, 'All are namelist options')
        CALL wrf_error_fatal ( 'dz above fixed eta levels is too large')
      END IF

      !  Add those 2 levels back into the middle, just above the 8 levels
      !  that semi define a boundary layer.  After we open up the levels,
      !  then we just linearly interpolate in znw.  So now levels 1-8 are
      !  specified as the fixed boundary layer levels given in this routine.
      !  The top levels, 12 through kte are those computed.  The middle
      !  levels 9, 10, and 11 are equi-spaced in znw, and are each 1/2 the
      !  the znw thickness of levels 11 through 12.

      DO k = kte-2 , 9 , -1
        znw(k+2) = znw(k)
      END DO

      znw( 9) = 0.75 * znw( 8) + 0.25 * znw(12)
      znw(10) = 0.50 * znw( 8) + 0.50 * znw(12)
      znw(11) = 0.25 * znw( 8) + 0.75 * znw(12)

      DO k = 8, kte-1
        pb = (znw(k)+znw(k+1))*0.5 * (p_surf - p_top) + p_top
        temp = MAX ( tiso, t00 + A*LOG(pb/p00) )
        IF ( pb .LT. p_strat ) THEN
          temp = tiso + A_strat * LOG ( pb/p_strat )
        END IF
        t_init = temp*(p00/pb)**(r_d/cp) - t0
        alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
        phb(k) = phb(k-1) - (znw(k)-znw(k-1)) * mub*alb(k-1)
      END DO
      phb(kte) = phb(kte-1) - (znw(kte)-znw(kte-1)) * mub*alb(kte-1)

      k=1
      WRITE (*,FMT='("Full level index = ",I4,"     Height = ",F7.1," m")') k,phb(1)/g
      do k = 2 ,kte
        WRITE (*,FMT='("Full level index = ",I4,"     Height = ",F7.1," m      Thickness = ",F6.1," m")') k,phb(k)/g,(phb(k)-phb(k-1))/g
      end do

END SUBROUTINE compute_eta

!-----------------------------------------------------------------------------------------

SUBROUTINE blend_terrain ( ter_interpolated , ter_input , &
                           ids , ide , jds , jde , kds , kde , &
                           ims , ime , jms , jme , kms , kme , &
                           ips , ipe , jps , jpe , kps , kpe )

   USE module_configure
   IMPLICIT NONE

   INTEGER , INTENT(IN)                       :: ids , ide , jds , jde , kds , kde , &
                                                 ims , ime , jms , jme , kms , kme , &
                                                 ips , ipe , jps , jpe , kps , kpe
   REAL , DIMENSION(ims:ime,kms:kme,jms:jme) , INTENT(IN)    :: ter_interpolated
   REAL , DIMENSION(ims:ime,kms:kme,jms:jme) , INTENT(INOUT) :: ter_input

   REAL , DIMENSION(ims:ime,kms:kme,jms:jme) :: ter_temp
   INTEGER :: i , j , k , spec_bdy_width
   REAL    :: r_blend_zones
   INTEGER blend_cell, blend_width

   !  The fine grid elevation comes from the horizontally interpolated
   !  parent elevation for the first spec_bdy_width row/columns, so we need
   !  to get that value.  We blend the coarse and fine in the next blend_width
   !  rows and columns.  After that, in the interior, it is 100% fine grid.

   CALL nl_get_spec_bdy_width ( 1, spec_bdy_width)
   CALL nl_get_blend_width ( 1, blend_width)

   !  Initialize temp values to the nest ter elevation.  This fills in the values
   !  that will not be modified below.

   DO j = jps , MIN(jpe, jde-1)
      DO k = kps , kpe
         DO i = ips , MIN(ipe, ide-1)
            ter_temp(i,k,j) = ter_input(i,k,j)
         END DO
      END DO
   END DO

   !  To avoid some tricky indexing, we fill in the values inside out.  This allows
   !  us to overwrite incorrect assignments.  There are replicated assignments, and
   !  there is much unnecessary "IF test inside of a loop" stuff.  For a large
   !  domain, this is only a patch; for a small domain, this is not a biggy.

   r_blend_zones = 1./(blend_width+1)
   DO j = jps , MIN(jpe, jde-1)
      DO k = kps , kpe
         DO i = ips , MIN(ipe, ide-1)
            DO blend_cell = blend_width,1,-1
               IF   ( ( i .EQ.       spec_bdy_width + blend_cell ) .OR.  ( j .EQ.       spec_bdy_width + blend_cell ) .OR. &
                      ( i .EQ. ide - spec_bdy_width - blend_cell ) .OR.  ( j .EQ. jde - spec_bdy_width - blend_cell ) ) THEN
                  ter_temp(i,k,j) = ( (blend_cell)*ter_input(i,k,j) + (blend_width+1-blend_cell)*ter_interpolated(i,k,j) ) &
                                    * r_blend_zones
               END IF
            ENDDO
            IF      ( ( i .LE.       spec_bdy_width     ) .OR.  ( j .LE.       spec_bdy_width     ) .OR. &
                      ( i .GE. ide - spec_bdy_width     ) .OR.  ( j .GE. jde - spec_bdy_width     ) ) THEN
               ter_temp(i,k,j) =      ter_interpolated(i,k,j)
            END IF
         END DO
      END DO
   END DO

   !  Set nest elevation with temp values.  All values not overwritten in the above
   !  loops have been previously set in the initial assignment.

   DO j = jps , MIN(jpe, jde-1)
      DO k = kps , kpe
         DO i = ips , MIN(ipe, ide-1)
            ter_input(i,k,j) = ter_temp(i,k,j)
         END DO
      END DO
   END DO

END SUBROUTINE blend_terrain

SUBROUTINE copy_3d_field ( ter_interpolated , ter_input , &
                           ids , ide , jds , jde , kds , kde , &
                           ims , ime , jms , jme , kms , kme , &
                           ips , ipe , jps , jpe , kps , kpe )

   IMPLICIT NONE

   INTEGER , INTENT(IN)                       :: ids , ide , jds , jde , kds , kde , &
                                                 ims , ime , jms , jme , kms , kme , &
                                                 ips , ipe , jps , jpe , kps , kpe
   REAL , DIMENSION(ims:ime,kms:kme,jms:jme) , INTENT(OUT) :: ter_interpolated
   REAL , DIMENSION(ims:ime,kms:kme,jms:jme) , INTENT(IN)  :: ter_input

   INTEGER :: i , j , k

   DO j = jps , MIN(jpe, jde-1)
      DO k = kps , kpe
         DO i = ips , MIN(ipe, ide-1)
            ter_interpolated(i,k,j) = ter_input(i,k,j)
         END DO
      END DO
   END DO

END SUBROUTINE copy_3d_field

SUBROUTINE adjust_tempqv ( mub, save_mub, znw, p_top, &
                           th, pp, qv,  &
                           ids , ide , jds , jde , kds , kde , &
                           ims , ime , jms , jme , kms , kme , &
                           ips , ipe , jps , jpe , kps , kpe )

   !USE module_configure
   !USE module_domain
   USE module_model_constants

   !USE module_bc
   !USE module_io_domain
   !USE module_state_description
   !USE module_timing
   !USE module_soil_pre
   IMPLICIT NONE

   INTEGER , INTENT(IN)                       :: ids , ide , jds , jde , kds , kde , &
                                                 ims , ime , jms , jme , kms , kme , &
                                                 ips , ipe , jps , jpe , kps , kpe
   REAL , DIMENSION(ims:ime,jms:jme) , INTENT(IN)    :: mub, save_mub
   REAL , DIMENSION(kms:kme) , INTENT(IN)    :: znw
   REAL , DIMENSION(ims:ime,kms:kme,jms:jme) , INTENT(INOUT) :: th, pp, qv

   REAL , DIMENSION(ims:ime,kms:kme,jms:jme) :: p_old, p_new, rh
   REAL :: es,dth,tc,e,dth1
   INTEGER :: i , j , k

   real p_top


! p_old = full pressure before terrain blending; also compute initial RH
! which is going to be conserved during terrain blending
   DO j = jps , MIN(jpe, jde-1)
      DO k = kps , kpe-1
         DO i = ips , MIN(ipe, ide-1)
            p_old(i,k,j) = 0.5*(znw(k+1)+znw(k))*save_mub(i,j) + p_top + pp(i,k,j)
            tc = (th(i,k,j)+300.)*(p_old(i,k,j)/1.e5)**(2./7.) - 273.15
            es = 610.78*exp(17.0809*tc/(234.175+tc))
            e = qv(i,k,j)*p_old(i,k,j)/(0.622+qv(i,k,j))
            rh(i,k,j) = e/es
         END DO
      END DO
   END DO

! p_new = full pressure after terrain blending; also compute temperature correction and convert RH back to QV
   DO j = jps , MIN(jpe, jde-1)
      DO k = kps , kpe-1
         DO i = ips , MIN(ipe, ide-1)
            p_new(i,k,j) = 0.5*(znw(k+1)+znw(k))*mub(i,j) + p_top + pp(i,k,j)
! 2*(g/cp-6.5e-3)*R_dry/g = -191.86e-3
            dth1 = -191.86e-3*(th(i,k,j)+300.)/(p_new(i,k,j)+p_old(i,k,j))*(p_new(i,k,j)-p_old(i,k,j))
            dth = -191.86e-3*(th(i,k,j)+0.5*dth1+300.)/(p_new(i,k,j)+p_old(i,k,j))*(p_new(i,k,j)-p_old(i,k,j))
            th(i,k,j) = th(i,k,j)+dth
            tc = (th(i,k,j)+300.)*(p_new(i,k,j)/1.e5)**(2./7.) - 273.15
            es = 610.78*exp(17.0809*tc/(234.175+tc))
            e = rh(i,k,j)*es
            qv(i,k,j) = 0.622*e/(p_new(i,k,j)-e)
         END DO
      END DO
   END DO


END SUBROUTINE adjust_tempqv

SUBROUTINE input_terrain_rsmas ( grid ,                        &
                           ids , ide , jds , jde , kds , kde , &
                           ims , ime , jms , jme , kms , kme , &
                           ips , ipe , jps , jpe , kps , kpe )

   USE module_domain, ONLY : domain
   IMPLICIT NONE
   TYPE ( domain ) :: grid

   INTEGER , INTENT(IN)                       :: ids , ide , jds , jde , kds , kde , &
                                                 ims , ime , jms , jme , kms , kme , &
                                                 ips , ipe , jps , jpe , kps , kpe

   LOGICAL, EXTERNAL ::  wrf_dm_on_monitor

   INTEGER :: i , j , k , myproc
   INTEGER, DIMENSION(256) :: ipath  ! array for integer coded ascii for passing path down to get_terrain
   CHARACTER*256 :: message, message2
   CHARACTER*256 :: rsmas_data_path

#if DM_PARALLEL
! Local globally sized arrays
   REAL , DIMENSION(ids:ide,jds:jde) :: ht_g, xlat_g, xlon_g
#endif

   CALL wrf_get_myproc ( myproc )

#if 0
CALL domain_clock_get ( grid, current_timestr=message2 )
WRITE ( message , FMT = '(A," HT before ",I3)' ) TRIM(message2), grid%id
write(30+myproc,*)ipe-ips+1,jpe-jps+1,trim(message)
do j = jps,jpe
do i = ips,ipe
write(30+myproc,*)grid%ht(i,j)
enddo
enddo
#endif

   CALL nl_get_rsmas_data_path(1,rsmas_data_path)
   do i = 1, LEN(TRIM(rsmas_data_path))
      ipath(i) = ICHAR(rsmas_data_path(i:i))
   enddo

#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )

   CALL wrf_patch_to_global_real ( grid%xlat , xlat_g , grid%domdesc, ' ' , 'xy' ,       &
                                   ids, ide-1 , jds , jde-1 , 1 , 1 , &
                                   ims, ime   , jms , jme   , 1 , 1 , &
                                   ips, ipe   , jps , jpe   , 1 , 1   )
   CALL wrf_patch_to_global_real ( grid%xlong , xlon_g , grid%domdesc, ' ' , 'xy' ,       &
                                   ids, ide-1 , jds , jde-1 , 1 , 1 , &
                                   ims, ime   , jms , jme   , 1 , 1 , &
                                   ips, ipe   , jps , jpe   , 1 , 1   )

   IF ( wrf_dm_on_monitor() ) THEN
     CALL get_terrain ( grid%dx/1000., xlat_g(ids:ide,jds:jde), xlon_g(ids:ide,jds:jde), ht_g(ids:ide,jds:jde), &
                        ide-ids+1,jde-jds+1,ide-ids+1,jde-jds+1, ipath, LEN(TRIM(rsmas_data_path)) )
     WHERE ( ht_g(ids:ide,jds:jde) < -1000. ) ht_g(ids:ide,jds:jde) = 0.
   ENDIF

   CALL wrf_global_to_patch_real ( ht_g , grid%ht , grid%domdesc, ' ' , 'xy' ,         &
                                   ids, ide-1 , jds , jde-1 , 1 , 1 , &
                                   ims, ime   , jms , jme   , 1 , 1 , &
                                   ips, ipe   , jps , jpe   , 1 , 1   )
#else

   CALL get_terrain ( grid%dx/1000., grid%xlat(ids:ide,jds:jde), grid%xlong(ids:ide,jds:jde), grid%ht(ids:ide,jds:jde), &
                       ide-ids+1,jde-jds+1,ide-ids+1,jde-jds+1, ipath, LEN(TRIM(rsmas_data_path)) )
   WHERE ( grid%ht(ids:ide,jds:jde) < -1000. ) grid%ht(ids:ide,jds:jde) = 0.

#endif

#if 0
CALL domain_clock_get ( grid, current_timestr=message2 )
WRITE ( message , FMT = '(A," HT after ",I3)' ) TRIM(message2), grid%id
write(30+myproc,*)ipe-ips+1,jpe-jps+1,trim(message)
do j = jps,jpe
do i = ips,ipe
write(30+myproc,*)grid%ht(i,j)
enddo
enddo
#endif

END SUBROUTINE input_terrain_rsmas

SUBROUTINE update_after_feedback_em ( grid  &
!
#include "dummy_new_args.inc"
!
                 )
!
! perform core specific updates, exchanges after
! model feedback  (called from med_feedback_domain) -John
!

! Driver layer modules
   USE module_domain, ONLY : domain, get_ijk_from_grid
   USE module_configure
   USE module_driver_constants
   USE module_machine
   USE module_tiles
#ifdef DM_PARALLEL
   USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, mytask
   USE module_comm_dm, ONLY : HALO_EM_FEEDBACK_sub
#else
   USE module_dm
#endif
   USE module_bc
! Mediation layer modules
! Registry generated module
   USE module_state_description

   IMPLICIT NONE

   !  Subroutine interface block.

   TYPE(domain) , TARGET         :: grid

   !  Definitions of dummy arguments
#include "dummy_new_decl.inc"

   INTEGER                         :: ids , ide , jds , jde , kds , kde , &
                                      ims , ime , jms , jme , kms , kme , &
                                      ips , ipe , jps , jpe , kps , kpe

  CALL wrf_debug( 500, "entering update_after_feedback_em" )

!  Obtain dimension information stored in the grid data structure.
  CALL get_ijk_from_grid (  grid ,                   &
                            ids, ide, jds, jde, kds, kde,    &
                            ims, ime, jms, jme, kms, kme,    &
                            ips, ipe, jps, jpe, kps, kpe    )

  CALL wrf_debug( 500, "before HALO_EM_FEEDBACK.inc in update_after_feedback_em" )
#ifdef DM_PARALLEL
#include "HALO_EM_FEEDBACK.inc"
#endif
  CALL wrf_debug( 500, "leaving update_after_feedback_em" )

END SUBROUTINE update_after_feedback_em

!-------------------------------------------------------

!DJW this is a direct copy from module_initialize_les, which is required for
!    subroutine vert_cor_vert_nesting_z_specified. Copying the subroutine was
!    necessary because of the order in which the modules are compiled.
      subroutine get_sounding( zk, p, p_dry, theta, rho, &
                               u, v, qv, dry, nl_max, nl_in, th_surf )
      USE module_model_constants
      implicit none

      integer nl_max, nl_in
      real zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
           u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max)
      logical dry

      integer n
      !KAL
      !parameter(n=1000)
      parameter(n=27000)
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
! KAL the atm was not fully balanced with 10 iterations
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
      write(6,*) ' sounding from subroutine get sounding'
      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
      do k=1,nl
        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), p_dry(k), theta(k), rho(k), u(k), v(k), qv(k)
      enddo

     end if

      end subroutine get_sounding

!-------------------------------------------------------

      subroutine read_sounding( ps,ts,qvs,h,th,qv,u,v,n,nl,debug )
      USE module_model_constants
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

!---------------------------------------------------
