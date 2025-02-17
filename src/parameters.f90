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
USE MOD_FiniteVolume2D_vars,ONLY: source_flag
USE MOD_FiniteVolume2D_vars,ONLY: maxTimeSteps
USE MOD_FiniteVolume2D_vars,ONLY: LinearWeightsOnly
USE MOD_FiniteVolume2D_vars,ONLY: MStepsMax
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
    source_flag = 0       
#endif
#ifdef EqnEuler
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
    source_flag = 0       
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
    source_flag = 0       
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
    source_flag = 0       
  !*------------------------------------------
  !*[5] Advection of smooth density sin4
  !*------------------------------------------
  CASE(5) 
    NameTest="Advection of smooth density sin4"
    TEnd    = 0.1 !*1.0 !*
    Gmm     = 1.4
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
    source_flag = 0       

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
    source_flag = 0       
#endif

#ifdef EqnShallowWater
  CASE(1) !*UNSTEADY SMOOTH VORTEX
    TEnd    = 0.1
    Gravity = 9.81
    nElemsX = 120
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/3.0,3.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
    BathymetryFlag = 0

  CASE(2) !*LAKE AT REST
    TEnd    = 0.1
    Gravity = 9.81
    nElemsX = 10
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
    BathymetryFlag = 1

  CASE(20) !*LAKE AT REST PERTURBED
    TEnd    = 0.1
    Gravity = 9.81
    nElemsX = 10
    nElemsY = nElemsX
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/) !*PERIODIC BCs
    BathymetryFlag = 1

  CASE(21) ! *PERTURBATION ANALYSIS ON WET LAKE AT REST NON-SMOOTH
    TEnd    = 0.5
    Gravity = 9.8
    nElemsX = 400
    nElemsY = 120
    MESH_X0 = (/-5.0,-2.0/)
    MESH_X1 = (/5.0,2.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 2


  CASE(3) ! *WET-DRY LAKE AT REST
    TEnd    = 1.0
    Gravity = 9.8
    ! nElemsX = 400
    ! nElemsY = 120
    ! MESH_X0 = (/-5.0,-2.0/)
    ! MESH_X1 = (/5.0,2.0/)
    nElemsX = 25
    nElemsY = 25
    MESH_X0 = (/-5.0,-5.0/)
    MESH_X1 = (/5.0,5.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 2

  CASE(30) ! *PERTURBATION ANALYSIS ON WET-DRY LAKE AT REST
    TEnd    = 1.0
    Gravity = 9.8
    ! nElemsX = 400
    ! nElemsY = 120
    ! MESH_X0 = (/-5.0,-2.0/)
    ! MESH_X1 = (/5.0,2.0/)
    nElemsX = 50
    nElemsY = 50
    MESH_X0 = (/-5.0,-5.0/)
    MESH_X1 = (/5.0,5.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 2


  CASE(31) ! *PERTURBATION ANALYSIS ON WET-DRY LAKE AT REST 1D-LIKE
    TEnd    = 1.0
    Gravity = 9.8
    ! nElemsX = 400
    ! nElemsY = 120
    ! MESH_X0 = (/-5.0,-2.0/)
    ! MESH_X1 = (/5.0,2.0/)
    nElemsX = 50
    nElemsY = 5
    MESH_X0 = (/-5.0,-5.0/)
    MESH_X1 = (/5.0,5.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 21

  CASE(4) ! *CIRCULAR DAM BREAK 1
    TEnd    = 1.0
    Gravity = 9.8
    nElemsX = 100
    nElemsY = 100
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/40.0,40.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 0

  CASE(5) ! *CIRCULAR DAM BREAK 2
    TEnd    = 1.0
    Gravity = 9.8
    nElemsX = 200
    nElemsY = 200
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/50.0,50.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 0

  CASE(6) ! *WAVE OVER DRY ISLAND
    TEnd    = 1.0
    Gravity = 9.8
    nElemsX = 200!400
    nElemsY = 60 !120
    MESH_X0 = (/-5.0,-2.0/)
    MESH_X1 = (/5.0,2.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 2


  CASE(7) ! *WAVE TSUNAMI
    TEnd    = 3.0
    Gravity = 9.8
    nElemsX = 960 !240!480  !960
    nElemsY = 320 !80 !160  !320
    MESH_X0 = (/-5.0,-2.0/)
    MESH_X1 = (/7.0,2.0/)
    BoundaryConditionsType = (/4,2,4,3/)
    BathymetryFlag = 5

  CASE(8) ! *SOLITARY WAVE ON CONICAL ISLAND
    TEnd    = 40.0
    Gravity = 9.81
    nElemsX = 125
    nElemsY = 150 
    MESH_X0 = (/ 0.0 , 0.0/)
    MESH_X1 = (/ 25.0 , 30.0/)
    BoundaryConditionsType = (/2,2,2,3/)
    BathymetryFlag = 8



  CASE(10) ! *RUN UP BP04 : H_over_d=0.3
    TEnd    = 30.0
    Gravity = 1.0
    nElemsX = 400
    nElemsY = 5 
    MESH_X0 = (/ -10.0 , 0.0/)
    MESH_X1 = (/ 40.0 , (50.0/nElemsX)*nElemsY /)
    ! BoundaryConditionsType = (/1,2,1,4/)
    BoundaryConditionsType = (/4,2,4,4/)
    BathymetryFlag = 10


  CASE(11) ! *RUN UP BP04 : H_over_d=0.3
    TEnd    = 70.0
    Gravity = 1.0
    nElemsX = 400
    nElemsY = 5 
    MESH_X0 = (/ -10.0 , 0.0/)
    MESH_X1 = (/ 80.0 , (90.0/nElemsX)*nElemsY /)
    ! BoundaryConditionsType = (/1,2,1,4/)
    BoundaryConditionsType = (/4,2,4,4/)
    BathymetryFlag = 10

  !*----------------------------------------
  !*Lakes at rest 40-41-42-43
  !*with bump bathymetry 2
  !*domain [-5,5]x[-2,2]
  !*----------------------------------------
  !*->40 Wet unperturbed
  !*->41 Wet perturbed
  !*->42 Wet-Dry unperturbed
  !*->43 Wet-Dry perturbed
  !*----------------------------------------
  CASE(40)
    TEnd    = 1.0
    Gravity = 9.8
    nElemsX = 100
    nElemsY = 30
    MESH_X0 = (/-5.0,-2.0/)
    MESH_X1 = (/5.0,2.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 2
  CASE(41)
    TEnd    = 1.0
    Gravity = 9.8
    nElemsX = 100
    nElemsY = 30
    MESH_X0 = (/-5.0,-2.0/)
    MESH_X1 = (/5.0,2.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 2
  CASE(42)
    TEnd    = 1.0
    Gravity = 9.8
    nElemsX = 100
    nElemsY = 30
    MESH_X0 = (/-5.0,-2.0/)
    MESH_X1 = (/5.0,2.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 2
  CASE(43)
    TEnd    = 1.0
    Gravity = 9.8
    nElemsX = 100
    nElemsY = 40 
    MESH_X0 = (/-5.0,-2.0/)
    MESH_X1 = (/5.0,2.0/)
    BoundaryConditionsType = (/1,1,1,1/)
    BathymetryFlag = 2
  !*----------------------------------------




  
  CASE(150) ! supercritical flow (Kurganov)
    TEnd    = 50.0
    Gravity = 9.812
    nElemsY = 2
    nElemsX = 25*nElemsY
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/25.0,1.0/)
    BoundaryConditionsType = (/1,2,1,3/) ! inflow - trasmissive
    BathymetryFlag = 4

  CASE(350) ! supercritical flow (Kurganov)
    TEnd    = 50.0
    Gravity = 9.812
    nElemsY = 100
    nElemsX = 250
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/25.0,10.0/)
    BoundaryConditionsType = (/1,2,1,3/) ! inflow - trasmissive
    BathymetryFlag = 4  
  CASE(450) ! supercritical flow (Kurganov) perturbation
    TEnd    = 1.0
    Gravity = 9.812
    nElemsY = 100
    nElemsX = 250
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/25.0,10.0/)
    BoundaryConditionsType = (/1,2,1,3/) ! inflow - trasmissive
    BathymetryFlag = 4


  CASE(151) ! subcritical flow (Kurganov)
    TEnd    = 200.0
    Gravity = 9.812
    nElemsY = 10
    nElemsX = 25*nElemsY
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/25.0,1.0/)
    BoundaryConditionsType = (/1,3,1,3/) ! inflow - outflow
    BathymetryFlag = 4

  CASE(351) ! subcritical flow (Kurganov)
    TEnd    = 200.0
    Gravity = 9.812
    nElemsY = 100
    nElemsX = 250
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/25.0,10.0/)
    BoundaryConditionsType = (/1,3,1,3/) ! inflow - outflow
    BathymetryFlag = 4


  CASE(451) ! subcritical flow (Kurganov)
    TEnd    = 1.0
    Gravity = 9.812
    nElemsY = 100
    nElemsX = 250
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/25.0,10.0/)
    BoundaryConditionsType = (/1,3,1,3/) ! inflow - outflow
    BathymetryFlag = 4

  CASE(152) ! transcritical no shock
    TEnd    = 200.0
    Gravity = 9.812
    nElemsY = 10
    nElemsX = 25*nElemsY
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/25.0,1.0/)
    BoundaryConditionsType = (/1,3,1,3/) ! inflow - outflow
    BathymetryFlag = 4

  CASE(352) ! transcritical no shock
    TEnd    = 200.0
    Gravity = 9.812
    nElemsY = 100
    nElemsX = 250
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/25.0,10.0/)
    BoundaryConditionsType = (/1,3,1,3/) ! inflow - outflow
    BathymetryFlag = 4


  CASE(153) ! transcritical shock
    TEnd    = 200.0
    Gravity = 9.812
    nElemsY = 10
    nElemsX = 25*nElemsY
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/25.0,1.0/)
    BoundaryConditionsType = (/1,3,1,3/) ! inflow - outflow
    BathymetryFlag = 4

  CASE(250) ! supercritical flow (Kurganov) with perturbation
    TEnd    = 1.0
    Gravity = 9.812
    nElemsY = 2
    nElemsX = 25*nElemsY
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/25.0,1.0/)
    BoundaryConditionsType = (/1,2,1,3/) ! inflow - trasmissive
    BathymetryFlag = 4

  CASE(251) ! subcritical flow (Kurganov) with perturbation
    TEnd    = 1.0
    Gravity = 9.812
    nElemsY = 10
    nElemsX = 25*nElemsY
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/25.0,1.0/)
    BoundaryConditionsType = (/1,3,1,3/) ! inflow - outflow
    BathymetryFlag = 4

#endif

#ifdef EqnAcoustics

#endif
  CASE DEFAULT
    ErrorMessage = "Initial condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

nargs = command_argument_COUNT()
IF (nargs > 1) THEN
   CALL get_command_ARGUMENT(2, arg)
   READ(arg, *) iarg
   nElemsX = iarg
END IF

IF (nargs > 2) THEN
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
!* 4,5=WENO5
!* 7=WENO7
!* 9=WENO9
!* 11=WENO11
!* 13=WENO13
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
LinearWeightsOnly = .FALSE.

IF (nargs > 4) THEN
   CALL get_command_ARGUMENT(5, arg)
   READ(arg, *) iarg
   Reconstruction = iarg
   ReconstructionFix = Reconstruction
END IF


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

IF (nargs > 5) THEN
   CALL get_command_ARGUMENT(6, arg)
   READ(arg, *) iarg
   timescheme = iarg
END IF

WhichRiemannSolver = 1 !* 1 Rusanov, 2 Exact

IF (nargs > 6) THEN
   CALL get_command_ARGUMENT(7, arg)
   READ(arg, *) iarg
   WhichRiemannSolver = iarg
END IF


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
IF (Reconstruction .GE. 3) THEN
  PRINT*, "Linear Weights Only = ", LinearWeightsOnly
END IF
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



SELECT CASE (timescheme)
  CASE(12,21,22,-2)    !* 20-n DeCu; 21-n DeCdu 
    MstepsMax=2 !*Order 1,2
  CASE(13,-3,-4) 
    MstepsMax=3 !*Order 3,4
  CASE(15,-5) 
    MstepsMax=4 !*Order 5
  CASE(17,-7) 
    MstepsMax=5 !*Order 7
  CASE(19,-9) 
    MstepsMax=6 !*Order 9
  CASE(-11) 
    MstepsMax=7 !*Order 11
  CASE(-13) 
    MstepsMax=8 !*Order 13
  CASE DEFAULT
    MstepsMax=1
END SELECT 



!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeParameters
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
