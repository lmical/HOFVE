!===============================================================================!
MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE InitializeParameters
  MODULE PROCEDURE InitializeParameters
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: InitializeParameters
!-------------------------------------------------------------------------------!
!
!
!
!===============================================================================!
CONTAINS
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE InitializeParameters()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: PI
USE MOD_FiniteVolume2D_vars,ONLY: CFL
USE MOD_FiniteVolume2D_vars,ONLY: TEnd
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
USE MOD_FiniteVolume2D_vars,ONLY: timescheme
USE MOD_FiniteVolume2D_vars,ONLY: WhichRiemannSolver
USE MOD_FiniteVolume2D_vars,ONLY: WhichOutput
USE MOD_FiniteVolume2D_vars,ONLY: nOutputFiles
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: VarNameVisu
USE MOD_FiniteVolume2D_vars,ONLY: GravitationalPotentialFlag
USE MOD_FiniteVolume2D_vars,ONLY: maxTimeSteps
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Gravity
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
CHARACTER(LEN=255) :: ErrorMessage
INTEGER            :: iarg, nargs  
REAL               :: iarg_real
CHARACTER(len=32)  :: arg
CHARACTER(LEN=80) :: NameTest
!-------------------------------------------------------------------------------!

PRINT*, "--------------------------"
PRINT*, "Initializing parameters   "
PRINT*, "--------------------------"

InitialCondition = 100  

nargs = command_argument_COUNT()
IF (nargs > 0) THEN
   CALL get_command_ARGUMENT(1, arg)
   READ(arg, *) iarg
   InitialCondition = iarg
END IF

SELECT CASE(InitialCondition)
!*------------------------------------------
!*NB: Due to the flag, a test for SW can have the same number as a test for Euler
!*But try to avoid it
!*------------------------------------------
#ifdef SW
  !*------------------------------------------
  !*[1] Unsteady smooth vortex for SW
  !*------------------------------------------
  CASE(1) !*UNSTEADY SMOOTH VORTEX
    NameTest="Unsteady smooth vortex for SW (NB: to be run with SW flag)"
    TEnd    = 0.1
    Gravity = 9.81
    Kappa   = 0.5*Gravity
    Gmm     = 2.0
    nElemsX = 512
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/3.0,3.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
    GravitationalPotentialFlag = 0       
#else
  !*------------------------------------------
  !*[2] Steady isentropic smooth vortex
  !*------------------------------------------
  CASE(2) 
    NameTest="Steady isentropic vortex"
    TEnd    = 0.1
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/-10.0,-10.0/)
    MESH_X1 = (/10.0,10.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
    GravitationalPotentialFlag = 0       
  !*------------------------------------------
  !*[3] Unsteady isentropic smooth vortex
  !*------------------------------------------
  CASE(3) 
    NameTest="Unsteady isentropic vortex"
    TEnd    = 0.1
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/-10.0,-10.0/)
    MESH_X1 = (/10.0,10.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
    GravitationalPotentialFlag = 0       
  !*------------------------------------------
  !*[4] Advection of smooth density
  !*------------------------------------------
  CASE(4) 
    NameTest="Advection of smooth density"
    TEnd    = 0.1
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
    GravitationalPotentialFlag = 0       

  !*------------------------------------------
  !*[892] Smooth periodic IC with the purpose of verifying conservation
  !*------------------------------------------
  CASE(892) 
    TEnd    = 0.1
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
    GravitationalPotentialFlag = 0       
#endif
  CASE DEFAULT
    ErrorMessage = "Initial condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

nargs = command_argument_COUNT()
IF (nargs == 2) THEN
   CALL get_command_ARGUMENT(2, arg)
   READ(arg, *) iarg
   nElemsX = iarg
   nElemsY = nElemsX
ELSE IF (nargs >2) THEN
   CALL get_command_ARGUMENT(2, arg)
   READ(arg, *) iarg
   nElemsX = iarg
   CALL get_command_ARGUMENT(3, arg)
   READ(arg, *) iarg
   nElemsY = iarg
END IF


CFL      = 0.5

IF (nargs >3) THEN
   CALL get_command_ARGUMENT(4, arg)
   READ(arg, *) iarg_real
   CFL = iarg_real
END IF


maxTimeSteps = 100000

!*---------------------------------------------
!*SPACE DISCRETIZATION (RECONSTRUCTION) LEGEND
!*---------------------------------------------
!* 1=First order FV
!* 2=MUSCL
!* 3=WENO3
!* 4=WENO5
!*---------------------------------------------
!* Different MINMOD limiters
!* 20=2 = MUSCL
!* 21   = k  !*OK
!* 22   = CO !*Ok only for k=1
!* 23   = VL !*Best
!* 24   = M  !*Best
!* 25   = VA !*Not 100% coded actually
!*---------------------------------------------

Reconstruction    = 4
ReconstructionFix = Reconstruction

!*---------------------------------------------
!*TIME SCHEME LEGEND
!*---------------------------------------------
!*  1 digit
!*  1 explicit euler, 2 SSPRK2, 3 SSPRK3, 4 SSPRK64,  5 RK65
!* 
!*  2 digits
!*  1* DeC   
!*  2* mPDeC 
!*---------------------------------------------

timescheme = 15

WhichRiemannSolver = 1 !* 1 Rusanov, 2 Exact

WhichOutput  = 0 ! 0 Nothing, 1 Octave, 2 Tecplot, 3 Both
nOutputFiles = 4

VarNameVisu(1) = "Density"
VarNameVisu(2) = "VelocityX"
VarNameVisu(3) = "VelocityY"
VarNameVisu(4) = "Pressure"
VarNameVisu(5) = "Gravitational_Potential"

PRINT*, "--------------------------"
PRINT*, "Test              = ", InitialCondition, TRIM(NameTest)
PRINT*, "Reconstruction    = ", Reconstruction
PRINT*, "ReconstructionFix = ", ReconstructionFix
PRINT*, "Time Scheme       = ", timescheme
SELECT CASE(WhichRiemannSolver)
  CASE(1)
    PRINT*, "Riemann solver: Rusanov"
  CASE(2)
    PRINT*, "Riemann solver: Exact Riemann solver"
  CASE DEFAULT
    PRINT*, "Wrong Riemann solver"
    PRINT*, WhichRiemannSolver
    STOP
END SELECT
PRINT*, "nElemsX = ", nElemsX, ", nElemsY = ", nElemsY 
PRINT*, "CFL = ", CFL
PRINT*, "--------------------------"

#ifdef SW
  PRINT*, "--------------------------"
  PRINT*, "SW setting"
  PRINT*, "p     = K*rho**gamma"
  PRINT*, "K     = ", Kappa
  PRINT*, "gamma = ", Gmm
  PRINT*, "--------------------------"
  !*Equivalence if K=g/2.0 and gamma=2.0
#endif

#ifdef PATANKAR
  PRINT*, "--------------------------"
  PRINT*, "Patankar active"
  PRINT*, "--------------------------"
#endif

#ifdef WELLBALANCED
  PRINT*, "--------------------------"
  PRINT*, "Well-balanced active"
  PRINT*, "--------------------------"
#endif




!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeParameters
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
