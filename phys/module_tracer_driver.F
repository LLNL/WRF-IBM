!WRF
!
MODULE module_tracer_driver

   USE module_configure, ONLY: grid_config_rec_type
   USE module_model_constants
!   USE module_state_description

CONTAINS

!=======================================================================
!             tracer_init.F
!
!  This subroutine initializes a passive tracer puff for
!
!  Branko Kosovic, March 12, 2010
!
SUBROUTINE tracer_init( scalar,                                 &
                        src_mask, q_src,                        &
                        ids, ide, jds, jde, kds, kde,           &
                        ims, ime, jms, jme, kms, kme,           &
                        its, ite, jts, jte, kts, kte             )
!-----------------------------------------------------------------------
! Begin declarations.

   IMPLICIT NONE

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT)   :: scalar
   INTEGER, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT( IN ) :: src_mask
   REAL, INTENT(IN   )                                       :: q_src
   INTEGER ,        INTENT(IN   ) ::        ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte

! local variables
   INTEGER                                               :: i,j,k
!-----------------------------------------------------------------------   
! executable starts here
  CALL wrf_debug(100,'phys/module_tracer_driver.F/subroutine tracer_init')   
!-----------------------------------------------------------------------


!    CALL wrf_error_fatal('wrong something')

   DO j=jts,MIN(jte,jde-1)
   DO k=kts,MIN(kte,kde-1)
   DO i=its,MIN(ite,ide-1)
      scalar(i,k,j) = q_src*float(src_mask(i,k,j))
   ENDDO
   ENDDO
   ENDDO

END SUBROUTINE tracer_init
!-----------------------------------------------------------------------

!=======================================================================
!             tracer_source_init.F
!
!  This subroutine initializes source of a passive tracer
!
!  Branko Kosovic, March 12, 2010
!
SUBROUTINE tracer_source_init( config_flags,                      &
                               src_mask,                          &
                               ids, ide, jds, jde, kds, kde,      &
                               ims, ime, jms, jme, kms, kme,      &
                               its, ite, jts, jte, kts, kte        )
!-----------------------------------------------------------------------
! Begin declarations.

   IMPLICIT NONE

   TYPE( grid_config_rec_type ), INTENT( IN )  :: config_flags
   INTEGER, DIMENSION(ims:ime,kms:kme,jms:jme) :: src_mask
   INTEGER ,        INTENT(IN   ) ::        ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte

! local variables
   INTEGER                                     :: i,j,k
   INTEGER                                     :: i_start, i_end, &
                                                  j_start, j_end, &
                                                  k_start, k_end
   INTEGER                                     :: is_src, ie_src, &
                                                  js_src, je_src, &
                                                  ks_src, ke_src
!-----------------------------------------------------------------------
! executable starts here
    CALL wrf_debug(100,'phys/module_tracer_driver.F/subroutine tracer_source_init')   
!-----------------------------------------------------------------------

    write(*,'(2(A,I3))') "DJW[module_tracer_driver/tracer_source_init]: config_flags%is_src=",config_flags%is_src,"   config_flags%ie_src=",config_flags%ie_src
    write(*,'(2(A,I3))') "DJW[module_tracer_driver/tracer_source_init]: config_flags%js_src=",config_flags%js_src,"   config_flags%je_src=",config_flags%je_src
    write(*,'(2(A,I3))') "DJW[module_tracer_driver/tracer_source_init]: config_flags%ks_src=",config_flags%ks_src,"   config_flags%ke_src=",config_flags%ke_src

    i_start = its
    i_end = MIN(ite,ide-1)
    j_start = jts
    j_end = MIN(jte,jde-1)
    k_start = kts
    k_end = MIN(kte,kde-1)

    write(*,'(6(A,I3))') "DJW[module_tracer_driver/tracer_source_init]: ids    =",ids,"   ide   =",ide,"   kds    =",kds,"   kde  =",kde,"   jds    =",jds,"   jde  =",jde
    write(*,'(6(A,I3))') "DJW[module_tracer_driver/tracer_source_init]: ims    =",ims,"   ime   =",ime,"   kms    =",kms,"   kme  =",kme,"   jms    =",jms,"   jme  =",jme
    write(*,'(6(A,I3))') "DJW[module_tracer_driver/tracer_source_init]: its    =",its,"   ite   =",ite,"   kts    =",kts,"   kte  =",kte,"   jts    =",jts,"   jte  =",jte
    write(*,'(6(A,I3))') "DJW[module_tracer_driver/tracer_source_init]: i_start=",i_start,"   i_ende=",i_end,"   k_start=",k_start,"   k_end=",k_end,"   j_start=",j_start,"   j_end=",j_end

    is_src = config_flags%is_src
    ie_src = config_flags%ie_src
    js_src = config_flags%js_src
    je_src = config_flags%je_src
    ks_src = config_flags%ks_src
    ke_src = config_flags%ke_src

    DO j=jts,jte
    DO k=kts,kte
    DO i=its,ite
       src_mask(i,k,j) = 0
    ENDDO
    ENDDO
    ENDDO

    IF ((.NOT. (is_src .GT. i_end)) .AND. (.NOT. (ie_src .LT. i_start)) .AND. &
        (.NOT. (js_src .GT. j_end)) .AND. (.NOT. (je_src .LT. j_start)) .AND. &
        (.NOT. (ks_src .GT. k_end)) .AND. (.NOT. (ke_src .LT. k_start)) ) THEN

       IF (is_src .LT. i_start) THEN
          is_src = i_start
       ENDIF
       IF (ie_src .GT. i_end) THEN
          ie_src = i_end
       ENDIF
       IF (js_src .LT. j_start) THEN
          js_src = j_start
       ENDIF
       IF (je_src .GT. j_end) THEN
          je_src = j_end
       ENDIF
       IF (ks_src .LT. k_start) THEN
          ks_src = k_start
       ENDIF
       IF (ke_src .GT. k_end) THEN
          ke_src = k_end
       ENDIF

       write(*,'(2(A,I3))') "DJW[module_tracer_driver/tracer_source_init]: is_src=",is_src,"   ie_src=",ie_src
       write(*,'(2(A,I3))') "DJW[module_tracer_driver/tracer_source_init]: js_src=",js_src,"   je_src=",je_src
       write(*,'(2(A,I3))') "DJW[module_tracer_driver/tracer_source_init]: ks_src=",ks_src,"   ke_src=",ke_src

       IF ((js_src .NE. je_src) .AND. &
           (ks_src .NE. ke_src) .AND. &
           (is_src .NE. ie_src))  THEN
          DO j=js_src, je_src
          DO k=ks_src, ke_src
          DO i=is_src, ie_src
             src_mask(i,k,j) = 1
          ENDDO
          ENDDO
          ENDDO
       ENDIF

       IF ((js_src .NE. je_src) .AND. &
           (ks_src .NE. ke_src) .AND. &
           (is_src .EQ. ie_src))  THEN
          DO j=js_src, je_src
          DO k=ks_src, ke_src
             src_mask(is_src,k,j) = 1
          ENDDO
          ENDDO
       ENDIF

       IF ((js_src .NE. je_src) .AND. &
           (ks_src .EQ. ke_src) .AND. &
           (is_src .NE. ie_src))  THEN
          DO j=js_src, je_src
          DO i=is_src, ie_src
             src_mask(i,ks_src,j) = 1
          ENDDO
          ENDDO
       ENDIF

       IF ((js_src .EQ. je_src) .AND. &
           (ks_src .NE. ke_src) .AND. &
           (is_src .NE. ie_src))  THEN
          DO k=ks_src, ke_src
          DO i=is_src, ie_src
             src_mask(i,k,js_src) = 1
          ENDDO
          ENDDO
       ENDIF

       IF ((js_src .EQ. je_src) .AND. &
           (ks_src .EQ. ke_src) .AND. &
           (is_src .NE. ie_src))  THEN
          DO i=is_src, ie_src
             src_mask(i,ks_src,js_src) = 1
          ENDDO
       ENDIF

       IF ((js_src .EQ. je_src) .AND. &
           (ks_src .NE. ke_src) .AND. &
           (is_src .EQ. ie_src))  THEN
          DO k=ks_src, ke_src
             src_mask(is_src,k,js_src) = 1
          ENDDO
       ENDIF

       IF ((js_src .NE. je_src) .AND. &
           (ks_src .EQ. ke_src) .AND. &
           (is_src .EQ. ie_src))  THEN
          DO j=js_src, je_src
             src_mask(is_src,ks_src,j) = 1
          ENDDO
       ENDIF

       IF ((js_src .EQ. je_src) .AND. &
           (ks_src .EQ. ke_src) .AND. &
           (is_src .EQ. ie_src))  THEN
          src_mask(is_src,ks_src,js_src) = 1
       ENDIF

    ELSE

!    Do nothing - source outside of the tile

    ENDIF


END SUBROUTINE tracer_source_init

!=======================================================================
!             tracer_source_driver.F
!
!  This subroutine adds a source term to a scalar tendency for 
!  a passive tracer
!
!  Branko Kosovic, March 12, 2010
!
SUBROUTINE tracer_source_driver( config_flags,                           &
                                 tracer, rho, dz8w,                      &
                                 src_mask,                               &
                                 dt,                                     &
                                 ids, ide, jds, jde, kds, kde,           &
                                 ims, ime, jms, jme, kms, kme,           &
                                 its, ite, jts, jte, kts, kte             )
!-----------------------------------------------------------------------
! Begin declarations.

   IMPLICIT NONE
   
   TYPE (grid_config_rec_type),  INTENT( IN )                :: config_flags
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme)                  :: tracer
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT( IN )    :: rho
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT( IN )    :: dz8w
   INTEGER, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT( IN ) :: src_mask
   REAL, INTENT(IN   )                                       :: dt
   INTEGER ,        INTENT(IN   ) ::        ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte
! local variables
   REAL                                                      :: q_src, &
                                                                dx, dy
   INTEGER                                                   :: i, j, k
!-----------------------------------------------------------------------
! executable starts here
   CALL wrf_debug(100,'phys/module_tracer_driver.F/subroutine tracer_source_driver')   
!-----------------------------------------------------------------------
   q_src=config_flags%q_src
   dx = config_flags%dx
   dy = config_flags%dy
   
   DO j=jts,MIN(jte,jde-1)
   DO k=kts,MIN(kte,kde-1)
   DO i=its,MIN(ite,ide-1)
      tracer(i,k,j) = tracer(i,k,j) + &
          q_src*float(src_mask(i,k,j))/rho(i,k,j)*dt/dz8w(i,k,j)/dx/dy
   ENDDO
   ENDDO
   ENDDO

END SUBROUTINE tracer_source_driver

!=======================================================================
!             tracer_nonperiodic.F
!
!  This subroutine imposes Neuman BCs to a passive tracer 
!  in X and Y directions, d(trc)/dx=0 and d(trc)/dy=0
!  
!  Branko Kosovic, March 29, 2010
!
SUBROUTINE tracer_nonperiodic( tracer,                                 &
                               ids, ide, jds, jde, kds, kde,           &
                               ims, ime, jms, jme, kms, kme,           &
                               its, ite, jts, jte, kts, kte             )
!-----------------------------------------------------------------------
! Begin declarations.

   IMPLICIT NONE

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme)                  :: tracer
   INTEGER ,        INTENT(IN   ) ::        ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte

! local variables
   INTEGER                                                   :: i,j,k
!-----------------------------------------------------------------------
! executable starts here
  CALL wrf_debug(100,'phys/module_tracer_driver.F/subroutine tracer_nonperiodic')   
!-----------------------------------------------------------------------

!    CALL wrf_error_fatal('wrong something')

    IF (its .eq. ids) THEN
       DO j=jts,MIN(jte,jde-1)
       DO k=kts,MIN(kte,kde-1)
          tracer(its+1,k,j)=tracer(its+2,k,j)
          tracer(its,k,j)=tracer(its+2,k,j)
       ENDDO
       ENDDO
    ENDIF

    IF (ite .eq. ide) THEN
       DO j=jts,MIN(jte,jde-1)
       DO k=kts,MIN(kte,kde-1)
          tracer(ite-1,k,j)=tracer(ite-2,k,j)
          tracer(ite,k,j)=tracer(ite-2,k,j)
       ENDDO
       ENDDO
    ENDIF

    IF (jts .eq. jds) THEN
       DO k=kts,MIN(kte,kde-1)
       DO i=its,MIN(ite,ide-1)
          tracer(i,k,jts+1)=tracer(i,k,jts+2)
          tracer(i,k,jts)=tracer(i,k,jts+2)
       ENDDO
       ENDDO
    ENDIF

    IF (jte .eq. jde) THEN
       DO k=kts,MIN(kte,kde-1)
       DO i=its,MIN(ite,ide-1)
          tracer(i,k,jte-1)=tracer(i,k,jte-2)
          tracer(i,k,jte)=tracer(i,k,jte-2)
       ENDDO
       ENDDO
    ENDIF

END SUBROUTINE tracer_nonperiodic

END MODULE module_tracer_driver
