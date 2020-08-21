!WRF:MEDIATION_LAYER:SOLVER

#define BENCH_START(A)
#define BENCH_END(A)

MODULE module_first_rk_step_part2

CONTAINS

  SUBROUTINE first_rk_step_part2 (   grid , config_flags              &
                             , moist , moist_tend               &
                             , chem  , chem_tend                &
                             , tracer, tracer_tend              &
                             , scalar , scalar_tend             &
                             , fdda3d, fdda2d                   &
                             , ru_tendf, rv_tendf               &
                             , rw_tendf, t_tendf                &
                             , ph_tendf, mu_tendf               &
                             , tke_tend                         &
                             , adapt_step_flag , curr_secs      &
                             , psim , psih , wspd , gz1oz0 , chklowq &
                             , cu_act_flag , hol , th_phy        &
                             , pi_phy , p_phy , t_phy     &
                             , dz8w , p8w , t8w           &
                             , nba_mij, n_nba_mij         & !JDM
                             , nba_rij, n_nba_rij         & !JDM
                             , ph_save, n, start                & !RSA
                             , ids, ide, jds, jde, kds, kde     &
                             , ims, ime, jms, jme, kms, kme     &
                             , ips, ipe, jps, jpe, kps, kpe     &
                             , imsx,imex,jmsx,jmex,kmsx,kmex    &
                             , ipsx,ipex,jpsx,jpex,kpsx,kpex    &
                             , imsy,imey,jmsy,jmey,kmsy,kmey    &
                             , ipsy,ipey,jpsy,jpey,kpsy,kpey    &
                             , k_start , k_end                  &
                            )
    USE module_state_description
    USE module_model_constants
    USE module_domain, ONLY : domain
    USE module_configure, ONLY : grid_config_rec_type, model_config_rec
#ifdef DM_PARALLEL
    USE module_dm, ONLY : local_communicator, mytask, ntasks, ntasks_x, ntasks_y, local_communicator_periodic, &
                          wrf_dm_maxval, wrf_err_message, local_communicator_x, local_communicator_y
    USE module_comm_dm, ONLY : halo_em_tke_c_sub,halo_em_tke_d_sub,halo_em_tke_e_sub            &
            ,halo_em_phys_pbl_sub,halo_em_phys_shcu_sub &
            ,halo_em_fdda_sub,halo_em_phys_diffusion_sub,halo_em_tke_3_sub &
            ,halo_em_tke_5_sub,halo_obs_nudge_sub,period_bdy_em_a1_sub,period_bdy_em_phy_bc_sub &
            ,period_bdy_em_fdda_bc_sub,period_bdy_em_chem_sub,halo_em_phys_cu_sub,halo_em_helicity_sub &
            ,halo_em_ibm_eddy_init_sub,period_bdy_em_ibm_eddy_init_sub & !RSA
            ,halo_em_ibm_eddy_inside_sub,period_bdy_em_ibm_eddy_inside_sub & !RSA
            ,halo_em_ibm_ss_recon_init_sub, period_bdy_em_ibm_ss_recon_init_sub & !RSA
            ,halo_em_ibm_ss_recon_inside_sub, period_bdy_em_ibm_ss_recon_inside_sub & !RSA
            ,halo_em_ibm_ss_extrap_init_sub, period_bdy_em_ibm_ss_extrap_init_sub & !RSA
            ,halo_em_ibm_ss_extrap_sub, period_bdy_em_ibm_ss_extrap_sub & !RSA
            ,halo_em_ibm_zero_ph_save_sub, period_bdy_em_ibm_zero_ph_save_sub !RSA
#endif

    USE module_driver_constants
    USE module_diffusion_em, ONLY : phy_bc, cal_deform_and_div, compute_diff_metrics, &
                                    vertical_diffusion_2, horizontal_diffusion_2, calculate_km_kh, &
                                    tke_rhs, cal_helicity, &
                                    assign_heat_flux !KAL
    USE module_em, ONLY : calculate_phy_tend
    USE module_fddaobs_driver, ONLY : fddaobs_driver
    USE module_bc, ONLY : set_physical_bc3d, set_physical_bc2d, &
                          set_physical_bc3d_int, set_physical_bc2d_int !RSA
    USE module_big_step_utilities_em, ONLY : zero_tend !RSA
    USE module_physics_addtendc, ONLY : update_phy_ten

    USE module_sfs_driver !JDM
    USE module_stoch, ONLY : update_stoch_ten, perturb_physics_tend,RAND_PERT_UPDATE
    USE module_ibm, ONLY : shear_stress_ibm_1D, & !, shear_stress_ibm_3D
                           vel_recon_ibm_eddy_init, vel_recon_ibm_eddy_part1, vel_recon_ibm_eddy_part2, & !RSA
                           shear_stress_recon_ibm_init, shear_stress_recon_ibm_part1, shear_stress_recon_ibm_part2, & !RSA
                           shear_stress_extrap_ibm_init, shear_stress_extrap_ibm !RSA
    USE module_stoch_les_inflow_perts !JDM

    IMPLICIT NONE

    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags

    INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde,     &
                           ims, ime, jms, jme, kms, kme,     &
                           ips, ipe, jps, jpe, kps, kpe,     &
                           imsx,imex,jmsx,jmex,kmsx,kmex,    &
                           ipsx,ipex,jpsx,jpex,kpsx,kpex,    &
                           imsy,imey,jmsy,jmey,kmsy,kmey,    &
                           ipsy,ipey,jpsy,jpey,kpsy,kpey



    LOGICAL ,INTENT(IN)                        :: adapt_step_flag
    REAL, INTENT(IN)                           :: curr_secs

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),INTENT(INOUT)   :: moist
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),INTENT(INOUT)   :: moist_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),INTENT(INOUT)   :: chem
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),INTENT(INOUT)   :: chem_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_fdda3d),INTENT(INOUT)  :: fdda3d
    REAL    ,DIMENSION(ims:ime,1:1,jms:jme,num_fdda2d),INTENT(INOUT)      :: fdda2d
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: psim
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: psih
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: wspd
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: gz1oz0
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: chklowq
    LOGICAL ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: cu_act_flag
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: hol

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: th_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: pi_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: p_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: dz8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: p8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t8w

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: ru_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: rv_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: rw_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: ph_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: tke_tend

    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: mu_tendf

    INTEGER , INTENT(IN)                          ::  k_start, k_end

!JDM
  INTEGER, INTENT(  IN ) :: n_nba_mij, n_nba_rij

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_nba_mij), INTENT(INOUT) & !KAL added intent
  :: nba_mij

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_nba_rij) &
  :: nba_rij

!RSA
  INTEGER, INTENT(  IN ) :: n, & !number of nearest neighbors for ibm
                            start !diagnostic for writing out ghost and neighbor points like in start_em
  REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: ph_save

! Local

    REAL, DIMENSION( ims:ime, jms:jme ) :: ht_loc
    REAL :: scale_factor
    INTEGER, DIMENSION( ims:ime, jms:jme ) :: shadowmask 
    INTEGER                             :: ij
    INTEGER  num_roof_layers
    INTEGER  num_wall_layers
    INTEGER  num_road_layers
    INTEGER  iswater
    INTEGER  rk_step 
    
! KAL adding defor12, 13, and 23 at the cell center (defor12_cc, defor13_cc, defor23_cc) 
! Smagorinsky k value is incorrect at the first grid point when there is a deformation at the surface, because D12 and D13 are hard-coded to be zero at the surface in cal_deform_and_div
    REAL, DIMENSION (ims:ime,kms:kme,jms:jme) :: defor12_cc, defor13_cc, defor23_cc   
    
    INTEGER :: i,j,k !KAL for writing/debugging   

#if ( WRF_DFI_RADAR == 1 )
    INTEGER  i_start,i_end,j_start,j_end,i,j,k
#endif

 ! initialize all tendencies to zero in order to update physics
 ! tendencies first (separate from dry dynamics).

   rk_step = 1

      IF ((grid%skebs_on==1).and.(grid%id .EQ. 1 )) then
          ! update and backtransform T
          CALL RAND_PERT_UPDATE(grid,'T',                                     &
                          grid%SPTFORCS,grid%SPTFORCC,                        &
                          grid%SPT_AMP,grid%ALPH_T,                           &
                          ips, ipe, jps, jpe, kps, kpe,                       &
                          ids, ide, jds, jde, kds, kde,                       &
                          ims, ime, jms, jme, kms, kme,                       &
                          k_start, k_end,                                     &
                          imsx,imex,jmsx,jmex,kmsx,kmex,                      &
                          ipsx,ipex,jpsx,jpex,kpsx,kpex,                      &
                          imsy,imey,jmsy,jmey,kmsy,kmey,                      &
                          ipsy,ipey,jpsy,jpey,kpsy,kpey,                      &
                          grid%num_stoch_levels,grid%num_stoch_levels,        &
                          grid%num_stoch_levels,grid%num_stoch_levels,        &
                          config_flags%restart, grid%iseedarr_skebs,          &
                          grid%DX,grid%DY,grid%skebs_vertstruc,               &
                          grid%rt_tendf_stoch,                                &
                          grid%VERTSTRUCC,grid%VERTSTRUCS,grid%VERTAMPUV      )
          ! Update streamfunction, backtransform U
           CALL RAND_PERT_UPDATE(grid,'U',                                    &
                           grid%SPSTREAMFORCS,grid%SPSTREAMFORCC,             &
                           grid%SPSTREAM_AMP,grid%ALPH_PSI,                   &
                           ips, ipe, jps, jpe, kps, kpe,                      &
                           ids, ide, jds, jde, kds, kde,                      &
                           ims, ime, jms, jme, kms, kme,                      &
                           k_start, k_end,                                    &
                           imsx,imex,jmsx,jmex,kmsx,kmex,                     &
                           ipsx,ipex,jpsx,jpex,kpsx,kpex,                     &
                           imsy,imey,jmsy,jmey,kmsy,kmey,                     &
                           ipsy,ipey,jpsy,jpey,kpsy,kpey,                     &
                           grid% num_stoch_levels,grid% num_stoch_levels,     &
                           grid% num_stoch_levels,grid% num_stoch_levels,     &
                           config_flags%restart, grid%iseedarr_skebs,         &
                           grid%DX,grid%DY,grid%skebs_vertstruc,              &
                           grid%ru_tendf_stoch,                               &
                           grid%VERTSTRUCC,grid%VERTSTRUCS,grid%VERTAMPUV     )
          ! Don't update streamfunction, backtransform V
           CALL RAND_PERT_UPDATE(grid,'V',                                    &
                           grid%SPSTREAMFORCS,grid%SPSTREAMFORCC,             &
                           grid%SPSTREAM_AMP,grid%ALPH_PSI,                   &
                           ips, ipe, jps, jpe, kps, kpe,                      &
                           ids, ide, jds, jde, kds, kde,                      &
                           ims, ime, jms, jme, kms, kme,                      &
                           k_start, k_end,                                    &
                           imsx,imex,jmsx,jmex,kmsx,kmex,                     &
                           ipsx,ipex,jpsx,jpex,kpsx,kpex,                     &
                           imsy,imey,jmsy,jmey,kmsy,kmey,                     &
                           ipsy,ipey,jpsy,jpey,kpsy,kpey,                     &
                           grid% num_stoch_levels,grid% num_stoch_levels,     &
                           grid% num_stoch_levels,grid% num_stoch_levels,     &
                           config_flags%restart, grid%iseedarr_skebs,         &
                           grid%DX,grid%DY,grid%skebs_vertstruc,              &
                           grid%rv_tendf_stoch,                               &
                           grid%VERTSTRUCC,grid%VERTSTRUCS,grid%VERTAMPT      )
       ENDIF !skebs_on 

     if ((grid%sppt_on==1).and.(grid%id .EQ. 1 )) then
          CALL RAND_PERT_UPDATE(grid,'T',                                     &
                          grid%SPPTFORCS,grid%SPPTFORCC,                      &
                          grid%SPPT_AMP,grid%ALPH_SPPT,                       &
                          ips, ipe, jps, jpe, kps, kpe,                       &
                          ids, ide, jds, jde, kds, kde,                       &
                          ims, ime, jms, jme, kms, kme,                       &
                          k_start, k_end,                                     &
                          imsx,imex,jmsx,jmex,kmsx,kmex,                      &
                          ipsx,ipex,jpsx,jpex,kpsx,kpex,                      &
                          imsy,imey,jmsy,jmey,kmsy,kmey,                      &
                          ipsy,ipey,jpsy,jpey,kpsy,kpey,                      &
                          grid%num_stoch_levels,grid%num_stoch_levels,        &
                          grid%num_stoch_levels,grid%num_stoch_levels,        &
                          config_flags%restart, grid%iseedarr_sppt,           &
                          grid%DX,grid%DY,grid%sppt_vertstruc,                &
                          grid%rstoch,                                        &
                          grid%VERTSTRUCC,grid%VERTSTRUCS,grid%VERTAMPT       )
       ENDIF !sppt_on

      if ((grid%rand_perturb_on==1).and.(grid%id .EQ. 1 )) then
           CALL RAND_PERT_UPDATE(grid,'T',                                     &
                           grid%SPFORCS,grid%SPFORCC,                          &
                           grid%SP_AMP,grid%ALPH_RAND,                         &
                           ips, ipe, jps, jpe, kps, kpe,                       &
                           ids, ide, jds, jde, kds, kde,                       &
                           ims, ime, jms, jme, kms, kme,                       &
                           k_start, k_end,                                     &
                           imsx,imex,jmsx,jmex,kmsx,kmex,                      &
                           ipsx,ipex,jpsx,jpex,kpsx,kpex,                      &
                           imsy,imey,jmsy,jmey,kmsy,kmey,                      &
                           ipsy,ipey,jpsy,jpey,kpsy,kpey,                      &
                           grid%num_stoch_levels,grid%num_stoch_levels,        &
                           grid%num_stoch_levels,grid%num_stoch_levels,        &
                           config_flags%restart, grid%iseedarr_rand_pert,      &
                           grid%DX,grid%DY,grid%rand_pert_vertstruc,           &
                           grid%RAND_PERT,                                     &
                           grid%VERTSTRUCC,grid%VERTSTRUCS,grid%VERTAMPT       )
       ENDIF !rand_perturb_on

! calculate_phy_tend
!KAL the subroutine calculate_phy_tend couples the physics tendencies with mu

BENCH_START(cal_phy_tend)
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )

      DO ij = 1 , grid%num_tiles

        CALL wrf_debug ( 200 , ' call calculate_phy_tend' )
        CALL calculate_phy_tend (config_flags,grid%mut,grid%muu,grid%muv,pi_phy, &
                     grid%rthraten,                                    &
                     grid%rublten,grid%rvblten,grid%rthblten,          &
                     grid%rqvblten,grid%rqcblten,grid%rqiblten,        &
                     grid%rucuten,grid%rvcuten,grid%rthcuten,          &
                     grid%rqvcuten,grid%rqccuten,grid%rqrcuten,        &
                     grid%rqicuten,grid%rqscuten,                      &
                     grid%rushten,grid%rvshten,grid%rthshten,          &
                     grid%rqvshten,grid%rqcshten,grid%rqrshten,        &
                     grid%rqishten,grid%rqsshten,grid%rqgshten,        &
                     grid%RUNDGDTEN,grid%RVNDGDTEN,grid%RTHNDGDTEN,grid%RQVNDGDTEN, &
                     grid%RMUNDGDTEN,                                  &
                     scalar, scalar_tend, num_scalar,                  &
                     tracer, tracer_tend, num_tracer,                  &
                     ids,ide, jds,jde, kds,kde,                        &
                     ims,ime, jms,jme, kms,kme,                        &
                     grid%i_start(ij), min(grid%i_end(ij),ide-1),      &
                     grid%j_start(ij), min(grid%j_end(ij),jde-1),      &
                     k_start    , min(k_end,kde-1)                     )

      ENDDO
      !$OMP END PARALLEL DO
BENCH_END(cal_phy_tend)

! tke diffusion

       IF(config_flags%diff_opt .eq. 2 .OR. config_flags%diff_opt .eq. 1) THEN

BENCH_START(comp_diff_metrics_tim)
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles
           CALL wrf_debug ( 200 , ' call compute_diff_metrics ' )
           CALL compute_diff_metrics ( config_flags, grid%ph_2, grid%phb, grid%z, grid%rdz, grid%rdzw, &
                                       grid%zx, grid%zy, grid%rdx, grid%rdy,                      &
                                       ids, ide, jds, jde, kds, kde,          &
                                       ims, ime, jms, jme, kms, kme,          &
                                       grid%i_start(ij), grid%i_end(ij),      &
                                       grid%j_start(ij), grid%j_end(ij),      &
                                       k_start    , k_end                    )
         ENDDO
         !$OMP END PARALLEL DO
BENCH_END(comp_diff_metrics_tim)

!start of calls for Cell Perturbation Method
IF ( ( config_flags%m_pblh_opt .EQ. 1 ) .AND. ( config_flags%km_opt .EQ. 4 ) ) THEN   
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles
           CALL force_down_meso_pblh( grid%m_pblh, grid%pblh,            &   
                                      ids, ide, jds, jde, kds, kde,      &
                                      ims, ime, jms, jme, kms, kme,      &
                                      grid%i_start(ij), grid%i_end(ij),  &
                                      grid%j_start(ij), grid%j_end(ij),  &
                                      k_start    , k_end                  )
         ENDDO
         !$OMP END PARALLEL DO
ENDIF

IF ( config_flags%les_pert_opt .GT. 0 ) THEN   
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles

           CALL calc_pert_t ( config_flags%les_pert_opt,                  &
                              config_flags%m_pblh_opt,                    &
                              grid%prttms, grid%prtdt, grid%prtnk,        &
                              grid%prtz, grid%prtseed, grid%pert_t,       &  
                              grid%m_pblh,                                & 
                              grid%t_2, grid%u_2, grid%v_2, grid%rdz,     &
                              grid%dx, grid%dt,                           & 
                              ids, ide, jds, jde, kds, kde,               &
                              ims, ime, jms, jme, kms, kme,               &
                              grid%i_start(ij), grid%i_end(ij),           &
                              grid%j_start(ij), grid%j_end(ij),           &
                              k_start    , k_end                          )

         ENDDO
         !$OMP END PARALLEL DO
ENDIF
!end of calls for Cell Perturbation Method

#ifdef DM_PARALLEL
#  include "HALO_EM_TKE_C.inc"
#  include "PERIOD_BDY_EM_A1.inc"
#endif

BENCH_START(tke_diff_bc_tim)
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )

         DO ij = 1 , grid%num_tiles

           CALL wrf_debug ( 200 , ' call bc for diffusion_metrics ' )
           CALL set_physical_bc3d( grid%rdzw , 'w', config_flags,           &
                                   ids, ide, jds, jde, kds, kde,       &
                                   ims, ime, jms, jme, kms, kme,       &
                                   ips, ipe, jps, jpe, kps, kpe,       &
                                   grid%i_start(ij), grid%i_end(ij),   &
                                   grid%j_start(ij), grid%j_end(ij),   &
                                   k_start    , k_end                 )
           CALL set_physical_bc3d( grid%rdz , 'w', config_flags,            &
                                   ids, ide, jds, jde, kds, kde,       &
                                   ims, ime, jms, jme, kms, kme,       &
                                   ips, ipe, jps, jpe, kps, kpe,       &
                                   grid%i_start(ij), grid%i_end(ij),   &
                                   grid%j_start(ij), grid%j_end(ij),   &
                                   k_start    , k_end                 )
           CALL set_physical_bc3d( grid%z , 'w', config_flags,              &
                                   ids, ide, jds, jde, kds, kde,       &
                                   ims, ime, jms, jme, kms, kme,       &
                                   ips, ipe, jps, jpe, kps, kpe,       &
                                   grid%i_start(ij), grid%i_end(ij),   &
                                   grid%j_start(ij), grid%j_end(ij),   &
                                   k_start    , k_end                 )
           CALL set_physical_bc3d( grid%zx , 'e', config_flags,             &
                                   ids, ide, jds, jde, kds, kde,       &
                                   ims, ime, jms, jme, kms, kme,       &
                                   ips, ipe, jps, jpe, kps, kpe,       &
                                   grid%i_start(ij), grid%i_end(ij),   &
                                   grid%j_start(ij), grid%j_end(ij),   &
                                   k_start    , k_end                 )
           CALL set_physical_bc3d( grid%zy , 'f', config_flags,             &
                                   ids, ide, jds, jde, kds, kde,       &
                                   ims, ime, jms, jme, kms, kme,       &
                                   ips, ipe, jps, jpe, kps, kpe,       &
                                   grid%i_start(ij), grid%i_end(ij),   &
                                   grid%j_start(ij), grid%j_end(ij),   &
                                   k_start    , k_end                 )
           CALL set_physical_bc2d( grid%ustm, 't', config_flags,            &
                                   ids, ide, jds, jde,                 &
                                   ims, ime, jms, jme,                 &
                                   ips, ipe, jps, jpe,                 &
                                   grid%i_start(ij), grid%i_end(ij),   &
                                   grid%j_start(ij), grid%j_end(ij)   )
          CALL set_physical_bc2d( grid%ust, 't', config_flags,              &
                                   ids, ide, jds, jde,                 &
                                   ims, ime, jms, jme,                 &
                                   ips, ipe, jps, jpe,                 &
                                   grid%i_start(ij), grid%i_end(ij),   &
                                   grid%j_start(ij), grid%j_end(ij)   )
 
         ENDDO
         !$OMP END PARALLEL DO
BENCH_END(tke_diff_bc_tim)

!KAL the deforXX_surf variables were added for the no-slip boundary condition with native terrain following coordinates

BENCH_START(deform_div_tim)

         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )

         DO ij = 1 , grid%num_tiles

           CALL wrf_debug ( 200 , ' call cal_deform_and_div' )
           CALL cal_deform_and_div ( config_flags,grid%u_2,grid%v_2,grid%w_2,grid%div,  &
                                     grid%defor11,grid%defor22,grid%defor33,            &
                                     grid%defor12,grid%defor13,grid%defor23,            &
                                     grid%defor11_surf,                                 & !KAL
                                     grid%defor22_surf,                                 & !KAL
                                     grid%defor12_surfu,                                & !KAL
                                     grid%defor12_surfv,                                & !KAL
                                     defor12_cc,                                        & !KAL
                                     defor13_cc,                                        & !KAL
                                     defor23_cc,                                        & !KAL      
                                     nba_rij, n_nba_rij,                                & !JDM
                                     grid%u_base, grid%v_base,grid%msfux,grid%msfuy,    &
                                     grid%msfvx,grid%msfvy,grid%msftx,grid%msfty,       &
                                     grid%rdx, grid%rdy, grid%dn, grid%dnw, grid%rdz,   &
                                     grid%rdzw,grid%fnm,grid%fnp,grid%cf1,grid%cf2,     &
                                     grid%cf3,grid%zx,grid%zy,            &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     grid%i_start(ij), grid%i_end(ij),    &
                                     grid%j_start(ij), grid%j_end(ij),    &
                                     k_start    , k_end                  )
         ENDDO
         !$OMP END PARALLEL DO
BENCH_END(deform_div_tim)

! Updraft helicity between output times

#ifdef DM_PARALLEL
#     include "HALO_EM_HELICITY.inc"
#endif

       IF ( ( config_flags%nwp_diagnostics .eq. 1 ) .OR. &
            ( ( config_flags%afwa_diag_opt .eq. 1 ) .AND. ( config_flags%afwa_severe_opt .EQ. 1 ) ) ) THEN
BENCH_START(helicity_tim)

       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )

       DO ij = 1 , grid%num_tiles

          CALL wrf_debug ( 200 , ' call cal_helicity' )
          CALL cal_helicity ( config_flags,grid%u_2,grid%v_2,grid%w_2,  &
                              grid%uh,                             &
                              grid%up_heli_max,                    &
                              grid%ph_2,grid%phb,                  &
                              grid%msfux,grid%msfuy,               &
                              grid%msfvx,grid%msfvy,               &
                              grid%ht,                             &
                              grid%rdx, grid%rdy, grid%dn, grid%dnw, grid%rdz, grid%rdzw,   &
                              grid%fnm,grid%fnp,grid%cf1,grid%cf2,grid%cf3,grid%zx,grid%zy, &
                              ids, ide, jds, jde, kds, kde,        &
                              ims, ime, jms, jme, kms, kme,        &
                              grid%i_start(ij), grid%i_end(ij),    &
                              grid%j_start(ij), grid%j_end(ij),    &
                              k_start    , k_end                  )
       ENDDO
       !$OMP END PARALLEL DO
BENCH_END(helicity_tim)
       ENDIF

#ifdef DM_PARALLEL
#     include "HALO_EM_TKE_D.inc"
#endif

! calculate tke, kmh, and kmv
BENCH_START(calc_tke_tim)
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles

           CALL wrf_debug ( 200 , ' call calculate_km_kh' )
           CALL calculate_km_kh( grid%ibm_ht_w,                                         & !KAL pull this out later 
                                 grid%top_s, grid%u_2,grid%v_2,                         & !KAL pull this out later
                                 config_flags,grid%dt,grid%dampcoef,grid%zdamp,         &
                                 config_flags%damp_opt,                                 &
                                 grid%xkmh,grid%xkmv,grid%xkhh,grid%xkhv,grid%bn2,      &
                                 grid%khdif,grid%kvdif,grid%div,                        &
                                 grid%defor11,grid%defor22,grid%defor33,grid%defor12,   &
                                 grid%defor13,grid%defor23,                             &
                                 defor12_cc, defor13_cc,defor23_cc,                     &
                                 grid%tke_2,p8w,t8w,th_phy,                             &
                                 t_phy,p_phy,moist,grid%dn,grid%dnw,                    &
                                 grid%dx,grid%dy,grid%rdz,grid%rdzw,                    &
                                 grid%z, grid%ht,                                       & !KAL added these for Mason and Thompson near surface model
                                 config_flags%mix_isotropic,                            &
                                 num_moist,                                             &
                                 grid%cf1, grid%cf2, grid%cf3, grid%warm_rain,          &
                                 grid%mix_upper_bound,                                  &
                                 grid%msftx, grid%msfty,                                &
                                 grid%zx, grid%zy,                                      & 
                                 ids,ide, jds,jde, kds,kde,                             &
                                 ims,ime, jms,jme, kms,kme,                             &
                                 grid%i_start(ij), grid%i_end(ij),                      &
                                 grid%j_start(ij), grid%j_end(ij),                      &
                                 k_start    , k_end                          )
         ENDDO
       !$OMP END PARALLEL DO
BENCH_END(calc_tke_tim)

#ifdef DM_PARALLEL
#     include "HALO_EM_TKE_E.inc"
#endif

         !RSA If using either eddy viscosity reconstruction or shear stress reconstruction, 
         !need to zero out ph_save because it is not needed here. Note that this
         !does not affect the rest of the code because ph_save is zeroed out
         !anyway in solve_em right after the call to module_first_rk_step_part2.
         IF ( ((config_flags%ibm_opt .EQ. 3) .AND. (config_flags%ibm_eddy .EQ. 1)) .OR. &
               (config_flags%ibm_opt .EQ. 5) ) THEN

BENCH_START(ibm_zero_ph_save_tim)
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles

               CALL zero_tend ( ph_save,                                               &
                                ids, ide, jds, jde, kds, kde,                          &
                                ims, ime, jms, jme, kms, kme,                          &
                                grid%i_start(ij), grid%i_end(ij),                      &
                                grid%j_start(ij), grid%j_end(ij),                      &
                                k_start    , k_end                          )

         ENDDO
         !$OMP END PARALLEL DO
BENCH_END(ibm_zero_ph_save_tim)

!ph_save
#ifdef DM_PARALLEL
#  include "HALO_EM_IBM_ZERO_PH_SAVE.inc"
#  include "PERIOD_BDY_EM_IBM_ZERO_PH_SAVE.inc"
#endif

         ENDIF

         !RSA option to reconstruct the eddy viscosity and diffusivity above the immersed boundary
         IF ((config_flags%ibm_opt .EQ. 3) .AND. (config_flags%ibm_eddy .EQ. 1)) THEN

BENCH_START(ibm_eddy_init_tim)
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles

            CALL set_physical_bc3d( ph_save, 'w', config_flags,       &
                                    ids, ide, jds, jde, kds, kde,     &
                                    ims, ime, jms, jme, kms, kme,     &
                                    ips, ipe, jps, jpe, kps, kpe,     &
                                    grid%i_start(ij), grid%i_end(ij), &
                                    grid%j_start(ij), grid%j_end(ij), &
                                    k_start, k_end )
            CALL vel_recon_ibm_eddy_init ( config_flags,                       &
                                           grid%ibm_ht_w,                      &
                                           grid%phb, ph_save, grid%ph_2,       &
                                           grid%ht_s,                          &
                                           grid%top_s,                         &
                                           grid%inside_s,                      &
                                           grid%inside_all_s,                  &
                                           ids, ide, jds, jde, kds, kde,       &
                                           ims, ime, jms, jme, kms, kme,       &
                                           grid%i_start(ij),grid%i_end(ij),    &
                                           grid%j_start(ij),grid%j_end(ij),    &
                                           k_start,k_end )

         ENDDO
         !$OMP END PARALLEL DO
BENCH_END(ibm_eddy_init_tim)

!ht_s, top_s, inside_s
#ifdef DM_PARALLEL
#  include "HALO_EM_IBM_EDDY_INIT.inc"
#  include "PERIOD_BDY_EM_IBM_EDDY_INIT.inc"
#endif

BENCH_START(ibm_eddy_tim)
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles
            CALL set_physical_bc3d( grid%ht_s, 't', config_flags,     &
                                    ids, ide, jds, jde, kds, kde,     &
                                    ims, ime, jms, jme, kms, kme,     &
                                    ips, ipe, jps, jpe, kps, kpe,     &
                                    grid%i_start(ij), grid%i_end(ij), &
                                    grid%j_start(ij), grid%j_end(ij), &
                                    k_start, k_end )
            CALL set_physical_bc2d_int( grid%top_s, 't', config_flags,    &
                                        ids, ide, jds, jde,               &
                                        ims, ime, jms, jme,               &
                                        ips, ipe, jps, jpe,               &
                                        grid%i_start(ij), grid%i_end(ij), &
                                        grid%j_start(ij), grid%j_end(ij) )

            CALL vel_recon_ibm_eddy_part1 (config_flags,                                &
                                           grid%top_s, grid%inside_s,                   &
                                           ids, ide, jds, jde, kds, kde,                &
                                           ims, ime, jms, jme, kms, kme,                &
                                           grid%i_start(ij), grid%i_end(ij),            &
                                           grid%j_start(ij), grid%j_end(ij),            &
                                           k_start    , k_end                           )

!inside_s
#ifdef DM_PARALLEL
#  include "HALO_EM_IBM_EDDY_INSIDE.inc"
#  include "PERIOD_BDY_EM_IBM_EDDY_INSIDE.inc"
#endif

            CALL set_physical_bc3d_int( grid%inside_s, 't', config_flags, &
                                        ids, ide, jds, jde, kds, kde,     &
                                        ims, ime, jms, jme, kms, kme,     &
                                        ips, ipe, jps, jpe, kps, kpe,     &
                                        grid%i_start(ij), grid%i_end(ij), &
                                        grid%j_start(ij), grid%j_end(ij), &
                                        k_start, k_end               )

            CALL vel_recon_ibm_eddy_part2( config_flags, 't', start,             &
                                           grid%ht_s, grid%inside_s,             &
                                           grid%inside_all_s,                    &
                                           grid%top_s,                           &
                                           grid%xkmh, grid%xkmv,                 &
                                           grid%xkhh, grid%xkhv,                 &
                                           grid%ibm_ht_u, grid%ibm_ht_v,         &
                                           grid%ibm_ht_w, grid%ibm_ht_c,         &
                                           n,                                    &
                                           grid%msftx, grid%msfty, grid%rdzw,    & !DJW 190930
                                           grid%dx, grid%dy, grid%dt,            & !DJW 190930
                                           grid%mix_upper_bound,                 & !DJW 190930
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           grid%i_start(ij),grid%i_end(ij),      &
                                           grid%j_start(ij),grid%j_end(ij),      &
                                           k_start, k_end )

         ENDDO
         !$OMP END PARALLEL DO
BENCH_END(ibm_eddy_tim)

#ifdef DM_PARALLEL
#     include "HALO_EM_TKE_E.inc"
#endif

         ENDIF !RSA end reconstruction of eddy viscosity and diffusivity

       ENDIF

#ifdef DM_PARALLEL
#      include "PERIOD_BDY_EM_PHY_BC.inc"
       IF ( config_flags%grid_fdda .eq. 1) THEN
#      include "PERIOD_BDY_EM_FDDA_BC.inc"
       ENDIF
#      include "PERIOD_BDY_EM_CHEM.inc"
#endif

BENCH_START(phy_bc_tim)
       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )

       DO ij = 1 , grid%num_tiles

         CALL wrf_debug ( 200 , ' call phy_bc' )
         CALL phy_bc (config_flags,grid%div,grid%defor11,grid%defor22,grid%defor33,            &
                      grid%defor12,grid%defor13,grid%defor23,      &
                      grid%xkmh,grid%xkmv,grid%xkhh,grid%xkhv,     &
                      grid%tke_2,grid%rho,                         &
                      grid%rublten, grid%rvblten,                  &
                      grid%rucuten, grid%rvcuten,                  &
                      grid%rushten, grid%rvshten,                  &
                      ids, ide, jds, jde, kds, kde,                &
                      ims, ime, jms, jme, kms, kme,                &
                      ips, ipe, jps, jpe, kps, kpe,                &
                      grid%i_start(ij), grid%i_end(ij),            &
                      grid%j_start(ij), grid%j_end(ij),            &
                      k_start    , k_end                           )
       ENDDO
       !$OMP END PARALLEL DO
BENCH_END(phy_bc_tim)

!JDM
IF ( ( config_flags%sfs_opt .GT. 0 ) .AND. ( config_flags%diff_opt .eq. 2 ) ) THEN

 CALL sfs_driver( grid, config_flags,     &
                  nba_mij, n_nba_mij,     & 
                  nba_rij, n_nba_rij      ) 

ENDIF

#ifdef DM_PARALLEL
!-----------------------------------------------------------------------
!
! MPP for some physics tendency, km, kh, deformation, and divergence
!
!                                                         * * * * * * *
!                                            * * * * *    * * * * * * *
!               *                     *      * * * * *    * * * * * * *
!             * + *      * + *        +      * * + * *    * * * + * * *
!               *                     *      * * * * *    * * * * * * *
!                                            * * * * *    * * * * * * *
!                                                         * * * * * * *
!
! (for PBL)
! rublten                  x
! rvblten                             x
!
! (for Cumulus)
! rucuten                  x
! rvcuten                             x
!
! (for Shallow Cumulus)
! rushten                  x
! rvshten                             x
!
! (for FDDA)
! rundgdten     x
! rvndgdten     x
!
! (for TKE3)
! tke_2                                          x               
! (for TKE5)
! tke_2                                                         x
!
! (for diff_opt >= 1)
! defor11                  x
! defor22                             x
! defor12       x
! defor13                  x
! defor23                             x
! div           x
! xkmv          x
! xkmh          x
! xkhv          x
! xkhh          x
! tke           x
!
!-----------------------------------------------------------------------
       IF ( config_flags%bl_pbl_physics .ge. 1 ) THEN
#      include "HALO_EM_PHYS_PBL.inc"
       ENDIF
       IF ( config_flags%shcu_physics .gt. 1 ) THEN
#      include "HALO_EM_PHYS_SHCU.inc"
       ENDIF
       IF ( config_flags%cu_physics == SASSCHEME      .or.   &
            config_flags%cu_physics == TIEDTKESCHEME  .or.   &
            config_flags%cu_physics == NTIEDTKESCHEME .or.   &
            config_flags%cu_physics == CAMZMSCHEME    .or.   &
            config_flags%cu_physics == MESO_SAS       .or.   &
            config_flags%cu_physics == NSASSCHEME ) THEN
#      include "HALO_EM_PHYS_CU.inc"
       ENDIF
       IF ( config_flags%grid_fdda .ge. 1) THEN
#      include "HALO_EM_FDDA.inc"
       ENDIF
       IF ( config_flags%diff_opt .ge. 1 ) THEN
#      include "HALO_EM_PHYS_DIFFUSION.inc"
       ENDIF

       IF      ( config_flags%h_mom_adv_order <= 4 ) THEN
#       include "HALO_EM_TKE_3.inc"
       ELSE IF ( config_flags%h_mom_adv_order <= 6 ) THEN
#       include "HALO_EM_TKE_5.inc"
       ELSE
         WRITE(wrf_err_message,*)'solve_em: invalid h_mom_adv_order = ',config_flags%h_mom_adv_order
         CALL wrf_error_fatal(TRIM(wrf_err_message))
       ENDIF
#endif

!KAL the routine update_phy_ten adds the physics tendencies to the 'regular' tendencies
! meaning: ru_tendf = ru_tendf + grid%rubltendf

BENCH_START(update_phy_ten_tim)
       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )

       DO ij = 1 , grid%num_tiles

         CALL wrf_debug ( 200 , ' call update_phy_ten' )
#if ( WRF_DFI_RADAR == 1 )
         if (config_flags%cu_physics .gt. 0) then
           i_start = grid%i_start(ij)
           i_end   = min( grid%i_end(ij),ide-1 )
           j_start = grid%j_start(ij)
           j_end   = min( grid%j_end(ij),jde-1 )
           if (grid%dfi_stage == DFI_FWD ) &
                 CALL wrf_debug ( 200 , ' Zero out cu_physics' )
           DO j = j_start, j_end
           DO k = k_start, min( k_end,kde-1 ) - 1
           DO i = i_start, i_end
             if (grid%dfi_stage ==DFI_FWD  &
             .and. grid%dfi_tten_rad(i,k,j) >= 1.0e-7 .and.  &
                   grid%dfi_tten_rad(i,k,j) <= 10.) then
! zero out cu-param temp tendency
                grid%rthcuten(i,k,j) = 0.0
             endif
           ENDDO
           ENDDO
           ENDDO
         ENDIF
#endif
         CALL update_phy_ten(ph_tendf,t_tendf, ru_tendf, rv_tendf,moist_tend ,&
                           scalar_tend, mu_tendf,                           &
                           grid%rthraten,grid%rthblten,grid%rthcuten,grid%rthshten, &
                           grid%rublten,grid%rucuten,grid%rushten,          &
                           grid%rvblten,grid%rvcuten,grid%rvshten,          &
                           grid%rqvblten,grid%rqcblten,grid%rqiblten,       &
                           grid%rqniblten,                                  & !CAMUWPBL scheme
                           grid%rqvcuten,grid%rqccuten,grid%rqrcuten,       &
                           grid%rqicuten,grid%rqscuten,                     &
                           grid%rqcncuten,grid%rqincuten,                   & !BSINGH - Added two CU tends
                           grid%rqvshten,grid%rqcshten,grid%rqrshten,       &
                           grid%rqishten,grid%rqsshten,grid%rqgshten,       &
                           grid%rqcnshten,grid%rqinshten,                   &!BSINGH - Added two SHCU tends
                           grid%RUNDGDTEN,                                  &
                           grid%RVNDGDTEN,grid%RTHNDGDTEN,grid%RPHNDGDTEN,  &
                           grid%RQVNDGDTEN,grid%RMUNDGDTEN,                 &
                           grid%rthfrten,grid%rqvfrten,                     &  ! fire
                           num_moist,num_scalar,config_flags,rk_step,       &
                           grid%adv_moist_cond,                             &
                           ids, ide, jds, jde, kds, kde,                    &
                           ims, ime, jms, jme, kms, kme,                    &
                           grid%i_start(ij), grid%i_end(ij),                &
                           grid%j_start(ij), grid%j_end(ij),                &
                           k_start, k_end                               )

       END DO
       !$OMP END PARALLEL DO
BENCH_END(update_phy_ten_tim)

      IF (grid%skebs_on==1) then
          !$OMP PARALLEL DO   &
          !$OMP PRIVATE ( ij )
          DO ij = 1 , grid%num_tiles
               CALL wrf_debug ( 200 , ' call update_stoch_ten' )
                CALL update_stoch_ten(ru_tendf, rv_tendf, t_tendf,&
                               grid%ru_tendf_stoch,                          &
                               grid%rv_tendf_stoch,                          &
                               grid%rt_tendf_stoch,                          &
                               grid%mu_2 , grid%mub,                         &
                               ids, ide, jds, jde, kds, kde,                 &
                               ims, ime, jms, jme, kms, kme,                 &
                               grid%i_start(ij), grid%i_end(ij),             &
                               grid%j_start(ij), grid%j_end(ij),             &
                               k_start, k_end,                               &
                               grid%num_stoch_levels,grid%num_stoch_levels   )

           ENDDO
           !$OMP END PARALLEL DO
      ENDIF !skebs_on

      IF (grid%sppt_on==1) then
          !$OMP PARALLEL DO   &
          !$OMP PRIVATE ( ij )
          ! JB comment:  P_QV is in moist_tend(ims,kms,jms,2)
          DO ij = 1 , grid%num_tiles
                 call perturb_physics_tend(grid%gridpt_stddev_sppt,          &
                        grid%stddev_cutoff_sppt,grid%rstoch,                 &
                        ru_tendf,rv_tendf,t_tendf,moist_tend(ims,kms,jms,2), &
                        ids, ide, jds, jde, kds, kde,                        &
                        ims, ime, jms, jme, kms, kme,                        &
                        grid%i_start(ij), grid%i_end(ij),                    &
                        grid%j_start(ij), grid%j_end(ij),                    &
                        k_start, k_end,                                      &
                        grid%num_stoch_levels,grid%num_stoch_levels          )
           ENDDO
          !$OMP END PARALLEL DO
  ENDIF

#ifdef PLANET
       ! do rayleigh (and zonal-average newtonian) damping during
       ! first iteration of RK loop only

       IF ( (config_flags%damp_opt == 101) .OR. &
            (config_flags%damp_opt == 103)      ) THEN
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles
           CALL damptop( grid%u_2, grid%v_2, grid%t_2, &
                         grid%mut, grid%muu, grid%muv, &
                         pi_phy,                                &
                         t_tendf, ru_tendf, rv_tendf, P2SI,     &
                         ids, ide, jds, jde, kds, kde,          &
                         ims, ime, jms, jme, kms, kme,          &
                         grid%i_start(ij), grid%i_end(ij),      &
                         grid%j_start(ij), grid%j_end(ij),      &
                         k_start, k_end                         )
         END DO
         !$OMP END PARALLEL DO
       END IF
#endif
!KAL this may need more if then statements.  
! Need to avoid writing over hfx if it has been calculated by a physics routine	 
     IF (config_flags%sf_surface_physics == 0) THEN ! if land-surface model is not used, assign a surface heat flux
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles
            CALL assign_heat_flux( config_flags,                      &
                                   moist(ims,kms,jms,P_QV), grid%rho, &
                                   grid%hfx, grid%itimestep, grid%dt, &
                                   ids, ide, jds, jde, kds, kde,      &
                                   ims, ime, jms, jme, kms, kme,      &
                                   grid%i_start(ij), grid%i_end(ij),  &
                                   grid%j_start(ij), grid%j_end(ij),  &
                                   k_start, k_end )
         ENDDO
         !$OMP END PARALLEL DO
     ENDIF

       IF( config_flags%diff_opt .eq. 2 .and. config_flags%km_opt .eq. 2 ) THEN

BENCH_START(tke_rhs_tim)
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles

           CALL tke_rhs  ( tke_tend,grid%bn2,                           &
                         config_flags,grid%defor11,grid%defor22,      &
                         grid%defor33,                                &
                         grid%defor12,grid%defor13,grid%defor23,      &
                         grid%u_2,grid%v_2,grid%w_2,grid%div,         &
                         grid%tke_2,grid%mut,                         &
                         th_phy,p_phy,p8w,t8w,grid%z,grid%fnm,        & 
                         grid%fnp,grid%cf1,grid%cf2,grid%cf3,         &     
                         grid%msftx,grid%msfty,grid%xkmh,             &
                         grid%xkmv,grid%xkhv,grid%rdx,grid%rdy,       &
                         grid%dx,grid%dy,grid%dt,grid%zx,grid%zy,     &
                         grid%rdz,grid%rdzw,grid%dn,                  &
                         grid%dnw,config_flags%mix_isotropic,         &
                         grid%hfx, grid%qfx, moist(ims,kms,jms,P_QV), &
                         grid%ustm, grid%rho,                         &
                         ids, ide, jds, jde, kds, kde,                &
                         ims, ime, jms, jme, kms, kme,                &
                         grid%i_start(ij), grid%i_end(ij),            &
                         grid%j_start(ij), grid%j_end(ij),            &
                         k_start    , k_end                           )

         ENDDO
         !$OMP END PARALLEL DO
BENCH_END(tke_rhs_tim)

       ENDIF

! calculate vertical diffusion first and then horizontal
! (keep this order)

       IF(config_flags%diff_opt .eq. 2) THEN

         IF (config_flags%bl_pbl_physics .eq. 0) THEN
!            IF ((config_flags%ibm_opt .GT. 0) .AND. (config_flags%ibm_stress_opt .EQ. 3)) THEN
!               CALL shear_stress_ibm_3D
!            ENDIF

             !RSA modify "defor" variables for IBM shear stress reconstruction option
             !These will then be used in vertical_diffusion_2 and horizontal_diffusion_2 
             !to calculate the stress divergence terms in the u, v, and w tendencies.
             !Note that defor is set to the stress divided by xkmv here such that when defor 
             !is multiplied by xkmv/h in the diffusion subroutines, the correct stress is recovered
             !(mix_isotropic=1 such that xkmv = xkmh).
             IF (config_flags%ibm_opt .EQ. 5) THEN

                !First reconstruct stress at points immediately above the immersed boundary.
                !i.e., the same points that are reconstructed in the VRM routines.
BENCH_START(shear_stress_recon_ibm_init_tim)
                !$OMP PARALLEL DO   &
                !$OMP PRIVATE ( ij )
                DO ij = 1 , grid%num_tiles

                   CALL set_physical_bc3d( ph_save, 'w', config_flags,       &
                                           ids, ide, jds, jde, kds, kde,     &
                                           ims, ime, jms, jme, kms, kme,     &
                                           ips, ipe, jps, jpe, kps, kpe,     &
                                           grid%i_start(ij), grid%i_end(ij), &
                                           grid%j_start(ij), grid%j_end(ij), &
                                           k_start, k_end )
                   CALL shear_stress_recon_ibm_init ( config_flags,                                &
                                                      grid%ibm_ht_u, grid%ibm_ht_v,                &
                                                      grid%ibm_ht_w, grid%ibm_ht_c,                &
                                                      grid%phb, ph_save, grid%ph_2,                &
                                                      grid%ht_u, grid%ht_v, grid%ht_w,             &
                                                      grid%ht_tauii, grid%ht_tau12,                & 
                                                      grid%ht_tau13, grid%ht_tau23,                &
                                                      grid%top_u, grid%top_v, grid%top_w,          &
                                                      grid%top_tauii, grid%top_tau12,              &
                                                      grid%top_tau13, grid%top_tau23,              &
                                                      grid%inside_u, grid%inside_v, grid%inside_w, &
                                                      grid%inside_tauii, grid%inside_tau12,        &
                                                      grid%inside_tau13, grid%inside_tau23,        &
                                                      grid%inside_all_tauii, grid%inside_all_tau12,&
                                                      grid%inside_all_tau13, grid%inside_all_tau23,&
                                                      ids, ide, jds, jde, kds, kde,                &
                                                      ims, ime, jms, jme, kms, kme,                &
                                                      grid%i_start(ij), grid%i_end(ij),            &
                                                      grid%j_start(ij), grid%j_end(ij),            &
                                                      k_start    , k_end                           )

                ENDDO
                !$OMP END PARALLEL DO
BENCH_END(shear_stress_recon_ibm_init_tim)

!ht_x, top_x, inside_x
#ifdef DM_PARALLEL
#  include "HALO_EM_IBM_SS_RECON_INIT.inc"
#  include "PERIOD_BDY_EM_IBM_SS_RECON_INIT.inc"
#endif

BENCH_START(shear_stress_recon_ibm_tim)
                !$OMP PARALLEL DO   &
                !$OMP PRIVATE ( ij )
                DO ij = 1 , grid%num_tiles

                   CALL set_physical_bc3d( grid%ht_u, 'u', config_flags,     &
                                           ids, ide, jds, jde, kds, kde,     &
                                           ims, ime, jms, jme, kms, kme,     &
                                           ips, ipe, jps, jpe, kps, kpe,     &
                                           grid%i_start(ij), grid%i_end(ij), &
                                           grid%j_start(ij), grid%j_end(ij), &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%ht_v, 'v', config_flags,     &
                                           ids, ide, jds, jde, kds, kde,     &
                                           ims, ime, jms, jme, kms, kme,     &
                                           ips, ipe, jps, jpe, kps, kpe,     &
                                           grid%i_start(ij), grid%i_end(ij), &
                                           grid%j_start(ij), grid%j_end(ij), &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%ht_w, 'w', config_flags,     &
                                           ids, ide, jds, jde, kds, kde,     &
                                           ims, ime, jms, jme, kms, kme,     &
                                           ips, ipe, jps, jpe, kps, kpe,     &
                                           grid%i_start(ij), grid%i_end(ij), &
                                           grid%j_start(ij), grid%j_end(ij), &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%ht_tauii, 't', config_flags,     &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%ht_tau12, 'd', config_flags,     &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%ht_tau13, 'e', config_flags,     &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%ht_tau23, 'f', config_flags,     &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc2d_int( grid%top_u, 'u', config_flags,    &
                                               ids, ide, jds, jde,               &
                                               ims, ime, jms, jme,               &
                                               ips, ipe, jps, jpe,               &
                                               grid%i_start(ij), grid%i_end(ij), &
                                               grid%j_start(ij), grid%j_end(ij) )
                   CALL set_physical_bc2d_int( grid%top_v, 'v', config_flags,    &
                                               ids, ide, jds, jde,               &
                                               ims, ime, jms, jme,               &
                                               ips, ipe, jps, jpe,               &
                                               grid%i_start(ij), grid%i_end(ij), &
                                               grid%j_start(ij), grid%j_end(ij) )
                   CALL set_physical_bc2d_int( grid%top_w, 'w', config_flags,    &
                                               ids, ide, jds, jde,               &
                                               ims, ime, jms, jme,               &
                                               ips, ipe, jps, jpe,               &
                                               grid%i_start(ij), grid%i_end(ij), &
                                               grid%j_start(ij), grid%j_end(ij) )
                   CALL set_physical_bc2d_int( grid%top_tauii, 't', config_flags,    &
                                               ids, ide, jds, jde,                   &
                                               ims, ime, jms, jme,                   &
                                               ips, ipe, jps, jpe,                   &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij) )
                   CALL set_physical_bc2d_int( grid%top_tau12, 'd', config_flags,    &
                                               ids, ide, jds, jde,                   &
                                               ims, ime, jms, jme,                   &
                                               ips, ipe, jps, jpe,                   &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij) )
                   CALL set_physical_bc2d_int( grid%top_tau13, 'e', config_flags,    &
                                               ids, ide, jds, jde,                   &
                                               ims, ime, jms, jme,                   &
                                               ips, ipe, jps, jpe,                   &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij) )
                   CALL set_physical_bc2d_int( grid%top_tau23, 'f', config_flags,    &
                                               ids, ide, jds, jde,                   &
                                               ims, ime, jms, jme,                   &
                                               ips, ipe, jps, jpe,                   &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij) )
                   CALL set_physical_bc3d_int( grid%inside_u, 'u', config_flags, &
                                               ids, ide, jds, jde, kds, kde,     &
                                               ims, ime, jms, jme, kms, kme,     &
                                               ips, ipe, jps, jpe, kps, kpe,     &
                                               grid%i_start(ij), grid%i_end(ij), &
                                               grid%j_start(ij), grid%j_end(ij), &
                                               k_start, k_end )
                   CALL set_physical_bc3d_int( grid%inside_v, 'v', config_flags, &
                                               ids, ide, jds, jde, kds, kde,     &
                                               ims, ime, jms, jme, kms, kme,     &
                                               ips, ipe, jps, jpe, kps, kpe,     &
                                               grid%i_start(ij), grid%i_end(ij), &
                                               grid%j_start(ij), grid%j_end(ij), &
                                               k_start, k_end )
                   CALL set_physical_bc3d_int( grid%inside_w, 'w', config_flags, &
                                               ids, ide, jds, jde, kds, kde,     &
                                               ims, ime, jms, jme, kms, kme,     &
                                               ips, ipe, jps, jpe, kps, kpe,     &
                                               grid%i_start(ij), grid%i_end(ij), &
                                               grid%j_start(ij), grid%j_end(ij), &
                                               k_start, k_end )

                   CALL shear_stress_recon_ibm_part1 ( config_flags,                                &
                                                       grid%top_tauii, grid%top_tau12,              &
                                                       grid%top_tau13, grid%top_tau23,              &
                                                       grid%inside_tauii, grid%inside_tau12,        &
                                                       grid%inside_tau13, grid%inside_tau23,        &
                                                       ids, ide, jds, jde, kds, kde,                &
                                                       ims, ime, jms, jme, kms, kme,                &
                                                       grid%i_start(ij), grid%i_end(ij),            &
                                                       grid%j_start(ij), grid%j_end(ij),            &
                                                       k_start    , k_end                           )

!inside_taux
#ifdef DM_PARALLEL
#  include "HALO_EM_IBM_SS_RECON_INSIDE.inc"
#  include "PERIOD_BDY_EM_IBM_SS_RECON_INSIDE.inc"
#endif

                   CALL set_physical_bc3d_int( grid%inside_tauii, 't', config_flags, &
                                               ids, ide, jds, jde, kds, kde,         &
                                               ims, ime, jms, jme, kms, kme,         &
                                               ips, ipe, jps, jpe, kps, kpe,         &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij),     &
                                               k_start, k_end )
                   CALL set_physical_bc3d_int( grid%inside_tau12, 'd', config_flags, &
                                               ids, ide, jds, jde, kds, kde,         &
                                               ims, ime, jms, jme, kms, kme,         &
                                               ips, ipe, jps, jpe, kps, kpe,         &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij),     &
                                               k_start, k_end )
                   CALL set_physical_bc3d_int( grid%inside_tau13, 'e', config_flags, &
                                               ids, ide, jds, jde, kds, kde,         &
                                               ims, ime, jms, jme, kms, kme,         &
                                               ips, ipe, jps, jpe, kps, kpe,         &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij),     &
                                               k_start, k_end )
                   CALL set_physical_bc3d_int( grid%inside_tau23, 'f', config_flags, &
                                               ids, ide, jds, jde, kds, kde,         &
                                               ims, ime, jms, jme, kms, kme,         &
                                               ips, ipe, jps, jpe, kps, kpe,         &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij),     &
                                               k_start, k_end )

                   CALL shear_stress_recon_ibm_part2 ( config_flags, 't', start,                    &
                                                       grid%u_2, grid%v_2, grid%w_2,                &
                                                       grid%xkmv, grid%fnm, grid%fnp,               &
                                                       grid%defor11, grid%defor22, grid%defor33,    &
                                                       grid%defor12, grid%defor13, grid%defor23,    &
                                                       grid%ht_u, grid%ht_v, grid%ht_w,             &
                                                       grid%ht_tauii,                               &
                                                       grid%inside_u, grid%inside_v, grid%inside_w, &
                                                       grid%inside_all_u, grid%inside_all_v,        & 
                                                       grid%inside_all_w,                           &
                                                       grid%inside_tauii, grid%top_tauii,           &
                                                       grid%ibm_ht_u, grid%ibm_ht_v,                &
                                                       grid%ibm_ht_w, grid%ibm_ht_c,                &
                                                       n,                                           &
                                                       ids, ide, jds, jde, kds, kde,                &
                                                       ims, ime, jms, jme, kms, kme,                &
                                                       grid%i_start(ij), grid%i_end(ij),            &
                                                       grid%j_start(ij), grid%j_end(ij),            &
                                                       k_start, k_end )
                   CALL shear_stress_recon_ibm_part2 ( config_flags, 'd', start,                    &
                                                       grid%u_2, grid%v_2, grid%w_2,                &
                                                       grid%xkmv, grid%fnm, grid%fnp,               &
                                                       grid%defor11, grid%defor22, grid%defor33,    &
                                                       grid%defor12, grid%defor13, grid%defor23,    &
                                                       grid%ht_u, grid%ht_v, grid%ht_w,             &
                                                       grid%ht_tau12,                               &
                                                       grid%inside_u, grid%inside_v, grid%inside_w, &
                                                       grid%inside_all_u, grid%inside_all_v,        & 
                                                       grid%inside_all_w,                           &
                                                       grid%inside_tau12, grid%top_tau12,           &
                                                       grid%ibm_ht_u, grid%ibm_ht_v,                &
                                                       grid%ibm_ht_w, grid%ibm_ht_c,                &
                                                       n,                                           &
                                                       ids, ide, jds, jde, kds, kde,                &
                                                       ims, ime, jms, jme, kms, kme,                &
                                                       grid%i_start(ij), grid%i_end(ij),            &
                                                       grid%j_start(ij), grid%j_end(ij),            &
                                                       k_start, k_end )
                   CALL shear_stress_recon_ibm_part2 ( config_flags, 'e', start,                    &
                                                       grid%u_2, grid%v_2, grid%w_2,                &
                                                       grid%xkmv, grid%fnm, grid%fnp,               &
                                                       grid%defor11, grid%defor22, grid%defor33,    &
                                                       grid%defor12, grid%defor13, grid%defor23,    &
                                                       grid%ht_u, grid%ht_v, grid%ht_w,             &
                                                       grid%ht_tau13,                               &
                                                       grid%inside_u, grid%inside_v, grid%inside_w, &
                                                       grid%inside_all_u, grid%inside_all_v,        & 
                                                       grid%inside_all_w,                           &
                                                       grid%inside_tau13, grid%top_tau13,           &
                                                       grid%ibm_ht_u, grid%ibm_ht_v,                &
                                                       grid%ibm_ht_w, grid%ibm_ht_c,                &
                                                       n,                                           &
                                                       ids, ide, jds, jde, kds, kde,                &
                                                       ims, ime, jms, jme, kms, kme,                &
                                                       grid%i_start(ij), grid%i_end(ij),            &
                                                       grid%j_start(ij), grid%j_end(ij),            &
                                                       k_start, k_end )
                   CALL shear_stress_recon_ibm_part2 ( config_flags, 'f', start,                    &
                                                       grid%u_2, grid%v_2, grid%w_2,                &
                                                       grid%xkmv, grid%fnm, grid%fnp,               &
                                                       grid%defor11, grid%defor22, grid%defor33,    &
                                                       grid%defor12, grid%defor13, grid%defor23,    &
                                                       grid%ht_u, grid%ht_v, grid%ht_w,             &
                                                       grid%ht_tau23,                               &
                                                       grid%inside_u, grid%inside_v, grid%inside_w, &
                                                       grid%inside_all_u, grid%inside_all_v,        & 
                                                       grid%inside_all_w,                           &
                                                       grid%inside_tau23, grid%top_tau23,           &
                                                       grid%ibm_ht_u, grid%ibm_ht_v,                &
                                                       grid%ibm_ht_w, grid%ibm_ht_c,                &
                                                       n,                                           &
                                                       ids, ide, jds, jde, kds, kde,                &
                                                       ims, ime, jms, jme, kms, kme,                &
                                                       grid%i_start(ij), grid%i_end(ij),            &
                                                       grid%j_start(ij), grid%j_end(ij),            &
                                                       k_start, k_end )

                ENDDO
                !$OMP END PARALLEL DO
BENCH_END(shear_stress_recon_ibm_tim)

                !Now extrapolate stress to points immediately below the immersed boundary.
                !i.e., ghost points in the IBM no slip routines.
BENCH_START(shear_stress_extrap_ibm_init_tim)
                !$OMP PARALLEL DO   &
                !$OMP PRIVATE ( ij )
                DO ij = 1 , grid%num_tiles

                   CALL shear_stress_extrap_ibm_init ( config_flags,                                &
                                                       grid%ibm_ht_u, grid%ibm_ht_v,                &
                                                       grid%ibm_ht_w, grid%ibm_ht_c,                &
                                                       grid%ht_tauii, grid%ht_tau12,                & 
                                                       grid%ht_tau13, grid%ht_tau23,                &
                                                       grid%top_tauii, grid%top_tau12,              &
                                                       grid%top_tau13, grid%top_tau23,              &
                                                       grid%inside_tauii, grid%inside_tau12,        &
                                                       grid%inside_tau13, grid%inside_tau23,        &
                                                       ids, ide, jds, jde, kds, kde,                &
                                                       ims, ime, jms, jme, kms, kme,                &
                                                       grid%i_start(ij), grid%i_end(ij),            &
                                                       grid%j_start(ij), grid%j_end(ij),            &
                                                       k_start    , k_end                           )

                ENDDO
                !$OMP END PARALLEL DO
BENCH_END(shear_stress_extrap_ibm_init_tim)

!deforx, top_x, inside_x
#ifdef DM_PARALLEL
#  include "HALO_EM_IBM_SS_EXTRAP_INIT.inc"
#  include "PERIOD_BDY_EM_IBM_SS_EXTRAP_INIT.inc"
#endif

BENCH_START(shear_stress_extrap_ibm_tim)
                !$OMP PARALLEL DO   &
                !$OMP PRIVATE ( ij )
                DO ij = 1 , grid%num_tiles

                   CALL set_physical_bc3d( grid%defor11, 't', config_flags,      &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%defor22, 't', config_flags,      &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%defor33, 't', config_flags,      &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%defor12, 'd', config_flags,      &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%defor13, 'e', config_flags,      &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%defor23, 'f', config_flags,      &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc2d_int( grid%top_tauii, 't', config_flags,    &
                                               ids, ide, jds, jde,                   &
                                               ims, ime, jms, jme,                   &
                                               ips, ipe, jps, jpe,                   &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij) )
                   CALL set_physical_bc2d_int( grid%top_tau12, 'd', config_flags,    &
                                               ids, ide, jds, jde,                   &
                                               ims, ime, jms, jme,                   &
                                               ips, ipe, jps, jpe,                   &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij) )
                   CALL set_physical_bc2d_int( grid%top_tau13, 'e', config_flags,    &
                                               ids, ide, jds, jde,                   &
                                               ims, ime, jms, jme,                   &
                                               ips, ipe, jps, jpe,                   &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij) )
                   CALL set_physical_bc2d_int( grid%top_tau23, 'f', config_flags,    &
                                               ids, ide, jds, jde,                   &
                                               ims, ime, jms, jme,                   &
                                               ips, ipe, jps, jpe,                   &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij) )
                   CALL set_physical_bc3d_int( grid%inside_tauii, 't', config_flags, &
                                               ids, ide, jds, jde, kds, kde,         &
                                               ims, ime, jms, jme, kms, kme,         &
                                               ips, ipe, jps, jpe, kps, kpe,         &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij),     &
                                               k_start, k_end )
                   CALL set_physical_bc3d_int( grid%inside_tau12, 'd', config_flags, &
                                               ids, ide, jds, jde, kds, kde,         &
                                               ims, ime, jms, jme, kms, kme,         &
                                               ips, ipe, jps, jpe, kps, kpe,         &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij),     &
                                               k_start, k_end )
                   CALL set_physical_bc3d_int( grid%inside_tau13, 'e', config_flags, &
                                               ids, ide, jds, jde, kds, kde,         &
                                               ims, ime, jms, jme, kms, kme,         &
                                               ips, ipe, jps, jpe, kps, kpe,         &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij),     &
                                               k_start, k_end )
                   CALL set_physical_bc3d_int( grid%inside_tau23, 'f', config_flags, &
                                               ids, ide, jds, jde, kds, kde,         &
                                               ims, ime, jms, jme, kms, kme,         &
                                               ips, ipe, jps, jpe, kps, kpe,         &
                                               grid%i_start(ij), grid%i_end(ij),     &
                                               grid%j_start(ij), grid%j_end(ij),     &
                                               k_start, k_end )

                   CALL shear_stress_extrap_ibm ( config_flags, 't', start,                    &
                                                  grid%xkmv, grid%fnm, grid%fnp,               &
                                                  grid%defor11, grid%defor22, grid%defor33,    &
                                                  grid%defor12, grid%defor13, grid%defor23,    &
                                                  grid%ht_tauii,                               &
                                                  grid%inside_tauii,                           &
                                                  grid%inside_all_tauii,                       &
                                                  grid%top_tauii,                              &
                                                  grid%ibm_ht_u, grid%ibm_ht_v,                &
                                                  grid%ibm_ht_w, grid%ibm_ht_c,                &
                                                  n,                                           &
                                                  ids, ide, jds, jde, kds, kde,                &
                                                  ims, ime, jms, jme, kms, kme,                &
                                                  grid%i_start(ij), grid%i_end(ij),            &
                                                  grid%j_start(ij), grid%j_end(ij),            &
                                                  k_start, k_end )
                   CALL shear_stress_extrap_ibm ( config_flags, 'd', start,                    &
                                                  grid%xkmv, grid%fnm, grid%fnp,               &
                                                  grid%defor11, grid%defor22, grid%defor33,    &
                                                  grid%defor12, grid%defor13, grid%defor23,    &
                                                  grid%ht_tau12,                               &
                                                  grid%inside_tau12,                           &
                                                  grid%inside_all_tau12,                       &
                                                  grid%top_tau12,                              &
                                                  grid%ibm_ht_u, grid%ibm_ht_v,                &
                                                  grid%ibm_ht_w, grid%ibm_ht_c,                &
                                                  n,                                           &
                                                  ids, ide, jds, jde, kds, kde,                &
                                                  ims, ime, jms, jme, kms, kme,                &
                                                  grid%i_start(ij), grid%i_end(ij),            &
                                                  grid%j_start(ij), grid%j_end(ij),            &
                                                  k_start, k_end )
                   CALL shear_stress_extrap_ibm ( config_flags, 'e', start,                    &
                                                  grid%xkmv, grid%fnm, grid%fnp,               &
                                                  grid%defor11, grid%defor22, grid%defor33,    &
                                                  grid%defor12, grid%defor13, grid%defor23,    &
                                                  grid%ht_tau13,                               &
                                                  grid%inside_tau13,                           &
                                                  grid%inside_all_tau13,                       &
                                                  grid%top_tau13,                              &
                                                  grid%ibm_ht_u, grid%ibm_ht_v,                &
                                                  grid%ibm_ht_w, grid%ibm_ht_c,                &
                                                  n,                                           &
                                                  ids, ide, jds, jde, kds, kde,                &
                                                  ims, ime, jms, jme, kms, kme,                &
                                                  grid%i_start(ij), grid%i_end(ij),            &
                                                  grid%j_start(ij), grid%j_end(ij),            &
                                                  k_start, k_end )
                   CALL shear_stress_extrap_ibm ( config_flags, 'f', start,                    &
                                                  grid%xkmv, grid%fnm, grid%fnp,               &
                                                  grid%defor11, grid%defor22, grid%defor33,    &
                                                  grid%defor12, grid%defor13, grid%defor23,    &
                                                  grid%ht_tau23,                               &
                                                  grid%inside_tau23,                           &
                                                  grid%inside_all_tau23,                       &
                                                  grid%top_tau23,                              &
                                                  grid%ibm_ht_u, grid%ibm_ht_v,                &
                                                  grid%ibm_ht_w, grid%ibm_ht_c,                &
                                                  n,                                           &
                                                  ids, ide, jds, jde, kds, kde,                &
                                                  ims, ime, jms, jme, kms, kme,                &
                                                  grid%i_start(ij), grid%i_end(ij),            &
                                                  grid%j_start(ij), grid%j_end(ij),            &
                                                  k_start, k_end )

                ENDDO
                !$OMP END PARALLEL DO
BENCH_END(shear_stress_extrap_ibm_tim)

!deforx
#ifdef DM_PARALLEL
#  include "HALO_EM_IBM_SS_EXTRAP.inc"
#  include "PERIOD_BDY_EM_IBM_SS_EXTRAP.inc"
#endif

             ENDIF !RSA end shear stress reconstruction and extrapolation (ibm_opt 5)

BENCH_START(vert_diff_tim)
           !$OMP PARALLEL DO   &
           !$OMP PRIVATE ( ij )
           DO ij = 1 , grid%num_tiles

                IF (config_flags%ibm_opt .EQ. 5) THEN
                   !RSA add BC calls for deforx because they were just modified above
                   CALL set_physical_bc3d( grid%defor11, 't', config_flags,      &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%defor22, 't', config_flags,      &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%defor33, 't', config_flags,      &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%defor12, 'd', config_flags,      &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%defor13, 'e', config_flags,      &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                   CALL set_physical_bc3d( grid%defor23, 'f', config_flags,      &
                                           ids, ide, jds, jde, kds, kde,         &
                                           ims, ime, jms, jme, kms, kme,         &
                                           ips, ipe, jps, jpe, kps, kpe,         &
                                           grid%i_start(ij), grid%i_end(ij),     &
                                           grid%j_start(ij), grid%j_end(ij),     &
                                           k_start, k_end )
                ENDIF !RSA

             CALL wrf_debug ( 200 , ' call vertical_diffusion_2 ' )
             CALL vertical_diffusion_2( ru_tendf, rv_tendf, rw_tendf,            &
                                      t_tendf, tke_tend,                         &
                                      moist_tend, num_moist,                      &
                                      chem_tend, num_chem,                       &
                                      scalar_tend, num_scalar,                     &
                                      tracer_tend, num_tracer,                     &
                                      grid%u_2, grid%v_2,                                  &
                                      grid%t_2,grid%u_base,grid%v_base,grid%t_base,grid%qv_base,          &
                                      grid%mut,grid%tke_2,config_flags, &
                                      grid%defor13,grid%defor23,grid%defor33,                   &
                                      nba_mij, num_nba_mij,          & !JDM
                                      grid%div, moist, chem, scalar,tracer,         &
                                      grid%xkmv, grid%xkhv, grid%xkmh, config_flags%km_opt,       & ! xkmh added Oct2013
                                      grid%fnm, grid%fnp, grid%dn, grid%dnw, grid%rdz, grid%rdzw, &
                                      grid%hfx, grid%qfx, grid%ustm, grid%rho,   &
                                      ids, ide, jds, jde, kds, kde,              &
                                      ims, ime, jms, jme, kms, kme,              &
                                      grid%i_start(ij), grid%i_end(ij),          &
                                      grid%j_start(ij), grid%j_end(ij),          &
                                      k_start, k_end                             )

!KAL  added the calls below to the ibm routine for the shear stress reconstruction.
! Here the ghost points are determined, and saved throughout the small steps until 
! rk_order = 3.  If the ibm implementation is something other than stress reconstruction
! then the ghost points are calculated below- after the velocities are calculated.
! Here ph_2 is passed, but below ph_1 is passed- they are the same values. 

!        IF ((config_flags%ibm_opt == 4).OR.(config_flags%ibm_opt == 5)) THEN
!	         CALL ghost_pts_ibm ( config_flags, phb, ph_2,      &
!	                              ibm_ht_u, ibm_ht_v, ibm_ht_w, & 
!	                              ids, ide, jds, jde, kds, kde, &			 	
!                                 ims, ime, jms, jme, kms, kme, &
!			                      grid%i_start(ij),grid%i_end(ij),grid%j_start(ij),grid%j_end(ij),k_start,k_end)
!            CALL shear_stress_ibm ( config_flags, phb, ph_2,      &
!	                                 u_2, v_2, ru_tendf, rv_tendf, & 
!				                     mut, defor13, defor23,        &
!				                     xkmv, fnm, fnp,               &
!				                     ibm_ht_u, ibm_ht_v, ibm_ht_w, &  
!			                         ids, ide, jds, jde, kds, kde, &			 	
!                                    ims, ime, jms, jme, kms, kme, &
!			                         grid%i_start(ij),grid%i_end(ij),grid%j_start(ij),grid%j_end(ij),k_start,k_end)
!         ENDIF !end of setting shear stress for ibm 
!KAL      

           ENDDO
           !$OMP END PARALLEL DO
BENCH_END(vert_diff_tim)

         ENDIF !bl_pbl_physics=0
!
BENCH_START(hor_diff_tim)
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles

           CALL wrf_debug ( 200 , ' call horizontal_diffusion_2' )
           CALL horizontal_diffusion_2( t_tendf, ru_tendf, rv_tendf, rw_tendf, &
                                      tke_tend,                              &
                                      moist_tend, num_moist,                  &
                                      chem_tend, num_chem,                   &
                                      scalar_tend, num_scalar,                 &
                                      tracer_tend, num_tracer,                 &
                                      grid%t_2, th_phy,                           &
                                      grid%mut, grid%tke_2, config_flags,              &
                                      grid%defor11, grid%defor22, grid%defor12,             &
                                      grid%defor13, grid%defor23,   &
                                      grid%defor11_surf,                          & !KAL
                                      grid%defor22_surf,                          &
                                      grid%defor12_surfu,                         &
                                      grid%defor12_surfv,                         & 
                                      nba_mij, num_nba_mij,         & !JDM
                                      grid%div,                     &
                                      moist, chem, scalar,tracer,               &
                                      grid%msfux,grid%msfuy, grid%msfvx,grid%msfvy, grid%msftx,  &
                                      grid%msfty, grid%xkmh, grid%xkhh, config_flags%km_opt,     &
                                      grid%rdx, grid%rdy, grid%rdz, grid%rdzw,                   &
                                      grid%fnm, grid%fnp, grid%cf1, grid%cf2, grid%cf3,          &
                                      grid%zx, grid%zy, grid%dn, grid%dnw, grid%rho,             &
                                      ids, ide, jds, jde, kds, kde,          &
                                      ims, ime, jms, jme, kms, kme,          &
                                      grid%i_start(ij), grid%i_end(ij),      &
                                      grid%j_start(ij), grid%j_end(ij),      &
                                      k_start    , k_end                    )

           IF ((config_flags%ibm_opt .GT. 0) .AND. (config_flags%ibm_stress_opt .EQ. 1)) THEN
              CALL shear_stress_ibm_1D( config_flags,                       &
                                        grid%mut,                           &
                                        grid%defor13, grid%defor23,         &
                                        nba_mij, num_nba_mij,               &
                                        grid%xkmv, grid%fnm, grid%fnp,      &
                                        grid%rho,                           & !DJW needed for cal_titau_* in v3.8.1
                                        grid%u_2, grid%v_2,                 &
                                        grid%ht_u, grid%ht_v, grid%ht_w,    &
                                        grid%ibm_ht_u, grid%ibm_ht_v,       &
                                        grid%top_u, grid%top_v, grid%top_w, &
                                        ru_tendf, rv_tendf,                 &
                                        ids, ide, jds, jde, kds, kde,       &
                                        ims, ime, jms, jme, kms, kme,       &
                                        grid%i_start(ij), grid%i_end(ij),   &
                                        grid%j_start(ij), grid%j_end(ij),   &
                                        k_start, k_end )
           ENDIF
 
         ENDDO
         !$OMP END PARALLEL DO
BENCH_END(hor_diff_tim)
       ENDIF

       IF ( grid%obs_nudge_opt .EQ. 1 .AND. grid%xtime <= grid%fdda_end ) THEN
# ifdef DM_PARALLEL
#       include "HALO_OBS_NUDGE.inc"
#endif
!***********************************************************************
! This section for obs nudging
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )

         DO ij = 1 , grid%num_tiles

           CALL fddaobs_driver (grid%grid_id, model_config_rec%grid_id, &
                   model_config_rec%parent_id, config_flags%restart,    &
                   config_flags,                                        &
                   grid%obs_nudge_opt,                                  &
                   grid%obs_ipf_errob,                                  &
                   grid%obs_ipf_nudob,                                  &
                   grid%fdda_start,                                     &
                   grid%fdda_end,                                       &
                   grid%obs_nudge_wind,                                 &
                   grid%obs_nudge_temp,                                 &
                   grid%obs_nudge_mois,                                 &
                   grid%obs_nudge_pstr,                                 &
                   grid%obs_coef_wind,                                  &
                   grid%obs_coef_temp,                                  &
                   grid%obs_coef_mois,                                  &
                   grid%obs_coef_pstr,                                  &             
                   grid%obs_rinxy,                                      &
                   grid%obs_rinsig,                                     &
                   grid%obs_npfi,                                       &
                   grid%obs_ionf,                                       &
                   grid%obs_prt_max,                                    &
                   grid%obs_prt_freq,                                   &
                   grid%obs_idynin,                                     &
                   grid%obs_dtramp,                                     &
                   grid%parent_grid_ratio,                              &
                   grid%max_dom, grid%itimestep,                        &
                   grid%xtime,                                          &
                   grid%dt, grid%gmt, grid%julday, grid%fdob,           &
                   grid%max_obs,                                        &
                   model_config_rec%nobs_ndg_vars,                      &
                   model_config_rec%nobs_err_flds,                      &
                   grid%fdob%nstat, grid%fdob%varobs, grid%fdob%errf,   &
                   grid%dx, grid%KPBL,grid%HT,                          &
                   grid%mut, grid%muu, grid%muv,                        &
                   grid%msftx, grid%msfty, grid%msfux, grid%msfuy, grid%msfvx, grid%msfvy, &
                   p_phy, t_tendf, t0,                                  &
                   grid%u_2, grid%v_2, grid%t_2,                        &
                   moist(ims,kms,jms,P_QV),                             &
                   grid%pb, grid%p_top, grid%p, grid%phb, grid%ph_2,    &
                   grid%uratx, grid%vratx, grid%tratx,                  &
                   ru_tendf, rv_tendf,                                  &
                   moist_tend(ims,kms,jms,P_QV), grid%obs_savwt,        &
                   grid%regime, grid%pblh, grid%z_at_w, grid%z,         &
                   ids,ide, jds,jde, kds,kde,                           &
                   ims,ime, jms,jme, kms,kme,                           &
                   grid%i_start(ij), min(grid%i_end(ij),ide-1),         &
                   grid%j_start(ij), min(grid%j_end(ij),jde-1),         &
                   k_start    , min(k_end,kde-1)                     )
 
         ENDDO
         !$OMP END PARALLEL DO
       ENDIF  ! obs_nudge_opt .eq. 1
! 
!***********************************************************************

  END SUBROUTINE first_rk_step_part2

END MODULE module_first_rk_step_part2

