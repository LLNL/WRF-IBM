MODULE module_ranarray
!KAL 10/13/2006 modified this ARPS routine for use in WRF 
! it is called from the ideal 2-d hill case

CONTAINS

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RANARY                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE ranary(nx,ny,nx1,nx2,ny1,ny2,iseed,amplit,rantp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate a 2-D array of machine-independent random numbers
!  between -amplit and +amplit with average value equal to zero.
!  The input parameter, iseed, must be a negative integer number
!  which needs to be defined only once in the calling subroutine
!  program.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: V.Wong and X.Song
!  7/20/1992.
!
!  Modification history:
!  7/25/1992. (MX)
!  Modified to allow assignment on a portion of the array.
!
!  12/11/1992 (MX)
!  Bug fix to the index bounds in loop 10.
!
!  9/1/94 (J. Levit &amp; Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx         Number of random numbers in the x-direction
!    ny         Number of random numbers in the y-direction
!    nx1,nx2    Begin and end 1st array indicies specifying a
!               subdomain in which the average is set to zero
!    ny1,ny2    Begin and end 2nd array indicies specifying a
!               subdomain in which the average is set to zero
!    iseed      an arbitrary negative integer as a seed for a
!               sequence of random numbers
!    amplit     The generated numbers stay within the range
!               [-amplit,amplit]
!
!  OUTPUT:
!
!    rantp      The array for storing random numbers
!
!  WORK ARRAY:
!
!    tem1      Temporary working array
!
!
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!
  IMPLICIT NONE        ! Force explicit declarations
!KAL taking out mpi stuff
!  INCLUDE 'mp.inc'
  INTEGER, PARAMETER :: nproc_x = 1
  INTEGER, PARAMETER :: nproc_y = 1
!KAL
  
  INTEGER :: i,j
  INTEGER :: nx,ny        ! Number of random numbers in each direction
  INTEGER :: nx1,nx2      ! Begin and end 1st array indecies specifying
                          ! a subdomain in which the average is set to
                          ! zero
  INTEGER :: ny1,ny2      ! Begin and end 2nd array indecies specifying
                          ! a subdomain in which the average is set to
                          ! zero
  INTEGER :: iseed        ! The seed for random number generation
  REAL :: amplit          ! The generated numbers stay within
                          ! [-amplit, amplit]
  REAL :: rantp (nx,ny)   ! Output array for storing the numbers
!KAL changing for wrf 
  REAL :: tem1 (nx,ny) 
!  REAL :: tem1 ((nx-3)*nproc_x+3,(ny-3)*nproc_y+3)   
                          ! For serial runs, this is just (nx,ny)
                          ! Temporary working array for storing the
                          ! numbers
!KAL 6-4-2007
! I am commenting out the type of ran3.  Allthough this has compiled with
! the intel compiler, I am getting an undefined reference error pointing
! towards ran3 with the pgi compiler
! according to the LC people, if the function is inside of the module
! that is calling it then it only needs to be typed (as in real) in the
! funtion.  If it is typed in the module that is calling the function,
! the the compiler is expecting the function to be external to the 
! module. 
! REAL :: ran3            ! The function to generate random number.
!KAL
!               
!-----------------------------------------------------------------------
!               
!  Miscellaneous local variables:
!  
!-----------------------------------------------------------------------   
  REAL :: wv,ranave,rantpmax  
!-----------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Initialize tem1 to zero
   tem1 = 0.0
!-----------------------------------------------------------------------
!
!  Generate random numbers between 0 and 1, for a given iseed.
!
!-----------------------------------------------------------------------
!KAL commenting out mpi/parallel stuff
!  IF (myproc == 0) THEN
    ranave=0.0
    DO j=ny1,ny2
      DO i=nx1,nx2        
        tem1(i,j)=ran3(iseed)
        ranave=ranave+tem1(i,j)
      END DO              
    END DO 
                 
    ranave=ranave/(FLOAT( (nx2-nx1+1)*(ny2-ny1+1) ))
! 
!-----------------------------------------------------------------------
!  Adjust the random numbers so that the new plane average value
!  equals zero and the random numbers stay within the range
!  [-amplit, amplit].
!-----------------------------------------------------------------------
!
    rantpmax=0.0
    DO j=ny1,ny2
      DO i=nx1,nx2
        rantpmax=AMAX1(rantpmax,ABS(ranave-tem1(i,j)))
      END DO
    END DO

    wv=amplit/rantpmax
    
    DO j=ny1,ny2
      DO i=nx1,nx2
        tem1(i,j)=(tem1(i,j)-ranave)*wv
      END DO
    END DO

!  ENDIF !end of processor 0 section

!KAL
!  CALL mpisplit2d(tem1,nx,ny,rantp)
!KAL adding below
    DO j=ny1,ny2
      DO i=nx1,nx2
        rantp(i,j)=tem1(i,j)
      END DO
    END DO
!KAL
  RETURN

END SUBROUTINE ranary

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  FUNCTION RAN3                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

  FUNCTION ran3(iseed)
!
!-----------------------------------------------------------------------
!  
!  PURPOSE:
!  
!  Generates a random number between 0 and 1 by feeding
!  a negative integer iseed.
!
!  Reference: "Seminumerical Algorithms" by Donald Knuth
!   
!-----------------------------------------------------------------------
!       
!  INPUT :
!   
!    iseed      an arbitrary negative integer as a seed for a
!               sequence of random numbers
!
!  OUTPUT:
!     
!    ran3       A random number between 0 and 1.
!     
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE        ! Force explicit declarations
  INTEGER :: iseed        ! The seed for random number generation
  REAL :: ran3            ! The function to generate random number.
!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables: 
!
!-----------------------------------------------------------------------
!
  INTEGER :: mbig,mseed,mz,k,ii,inext,inextp,i,iff,mk,mj
  REAL :: fac
  PARAMETER (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
  INTEGER :: ma(55)
  SAVE iff, inext, inextp, ma
  DATA iff /0/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Beginning of executable code...
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!----------------------------------------------------------------------
!  Initialize the sequence of random numbers between 0 and 1,
!  using iseed.
!----------------------------------------------------------------------
!
  IF (iseed < 0.OR.iff == 0) THEN
    iff=1
    mj=mseed-IABS(iseed)
    mj=MOD(mj,mbig)
    ma(55)=mj
    mk=1

    DO i=1,54
      ii=MOD(21*i,55)
      ma(ii)=mk
      mk=mj-mk
      IF (mk < mz) mk=mk+mbig
      mj=ma(ii)
    END DO

    DO k=1,4
      DO i=1,55
        ma(i)=ma(i)-ma(1+MOD(i+30,55))
        IF (ma(i) < mz) ma(i)=ma(i)+mbig
      END DO
    END DO

    inext=0
    inextp=31
    iseed=1
 
  END IF
!----------------------------------------------------------------------
! 
!  Start to generate a random number.
!
!----------------------------------------------------------------------

  inext=inext+1 
  IF (inext == 56) inext=1
  inextp=inextp+1
  IF (inextp == 56) inextp=1
  mj=ma(inext)-ma(inextp)
  IF (mj < mz) mj=mj+mbig
  ma(inext)=mj
  ran3=mj*fac

  RETURN 

  END FUNCTION ran3

END MODULE module_ranarray
