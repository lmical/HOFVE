!===============================================================================!
MODULE MOD_Equation
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ExactFunction
  MODULE PROCEDURE ExactFunction
END INTERFACE

INTERFACE ExactFunctionWB
  MODULE PROCEDURE ExactFunctionWB
END INTERFACE

INTERFACE GlobalFluxTerms
  MODULE PROCEDURE GlobalFluxTerms
END INTERFACE

INTERFACE SourceTerms
  MODULE PROCEDURE SourceTerms
END INTERFACE

INTERFACE BoundaryConditions
  MODULE PROCEDURE BoundaryConditions
END INTERFACE

INTERFACE TimeStep
  MODULE PROCEDURE TimeStep
END INTERFACE

INTERFACE RiemannSolver
  MODULE PROCEDURE RiemannSolver
END INTERFACE

INTERFACE EvaluateFlux1D
  MODULE PROCEDURE EvaluateFlux1D
END INTERFACE

INTERFACE ConsToPrim
  MODULE PROCEDURE ConsToPrim
END INTERFACE

INTERFACE PrimToCons
  MODULE PROCEDURE PrimToCons
END INTERFACE

INTERFACE Gravitational_Potential
  MODULE PROCEDURE Gravitational_Potential
END INTERFACE

#ifdef GFWENO
INTERFACE RiemannSolverCorner
  MODULE PROCEDURE RiemannSolverCorner
END INTERFACE
#endif

!-------------------------------------------------------------------------------!
PUBLIC :: ExactFunction
PUBLIC :: ExactFunctionWB
PUBLIC :: SourceTerms
PUBLIC :: GlobalFluxTerms
PUBLIC :: BoundaryConditions
PUBLIC :: TimeStep
PUBLIC :: RiemannSolver
PUBLIC :: EvaluateFlux1D
PUBLIC :: ConsToPrim
PUBLIC :: PrimToCons
PUBLIC :: Gravitational_Potential
#ifdef GFWENO
PUBLIC :: RiemannSolverCorner 
#endif
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
SUBROUTINE ExactFunction(WhichInitialCondition,t,x,Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MESH_SX
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: MIN_POSITIVE_VAR
#ifdef EqnEuler
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#endif
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Gravity
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN) :: WhichInitialCondition
REAL,INTENT(IN)    :: t
REAL,INTENT(IN)    :: x(1:nDims)
REAL,INTENT(OUT)   :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: Prim(1:nVar)
REAL               :: xc(2), xm(2), r, r0, hl, hr, r2, r20
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!
!*OUR VARIABLES
REAL               :: Omega, Jamma, u_inf, v_inf, h_inf, DeltaH
INTEGER            :: power
REAL               :: ro_inf, p_inf, beta, delta_u, delta_v, delta_T
REAL               :: xmxc(1:2), x_wrt_BL(1:2), x_wrt_BL_bm(1:2), x_0(1:2), x_d(1:2)




Cons = 0.0
Prim = 0.0
SELECT CASE(WhichInitialCondition)
#ifdef SW
  !*------------------------------------------
  !*[1] Unsteady smooth vortex for SW
  !*------------------------------------------
CASE(1)
  u_inf = 2.
  v_inf = 3.
  H_inf=1.
  r0 = 1.

  xm(1) = MESH_X0(1)+0.5*MESH_SX(1)
  xm(2) = MESH_X0(2)+0.5*MESH_SX(2)
  xc(1) = MODULO( x(1)-u_inf*t-MESH_X0(1) , Mesh_SX(1) ) + MESH_X0(1)-xm(1)
  xc(2) = MODULO( x(2)-v_inf*t-MESH_X0(2) , Mesh_SX(2) ) + MESH_X0(2)-xm(2)
  r     = (xc(1)**2 + xc(2)**2)
  Omega = sqrt(2.*Gravity*hDerivSmoothAuxiliary(r))

  Prim(1) = H_inf
  Prim(2) = u_inf
  Prim(3) = v_inf

  IF (r .LT. 1) THEN
    Prim(1) = hSmoothAuxiliary(r)
    Prim(2) = Prim(2)+Omega*(+xc(2))
    Prim(3) = Prim(3)+Omega*(-xc(1))
  END IF

  Prim(4)= Kappa*Prim(1)**Gmm

  CALL PrimToCons(Prim,Cons)
#endif
#ifdef EqnEuler
  !*------------------------------------------
  !*[2] Steady isentropic vortex
  !*------------------------------------------
CASE(2)

  u_inf=0.0
  v_inf=0.0

  !*Center of the vortex
  xc=0.5*(MESH_X1+MESH_X0)

  !*Coordinates from the center of the vortex
  xmxc=x-xc

  !*Distance squared from the center of the vortex
  r2=xmxc(1)**2+xmxc(2)**2

  !*Vortex amplitude
  beta=5.0 !*5.0 0.1

  delta_u=beta/(2.0*PI)*EXP( 0.5*( 1.0-r2 ) )*( -xmxc(2) )
  delta_v=beta/(2.0*PI)*EXP( 0.5*( 1.0-r2 ) )*xmxc(1)
  delta_T=-(Gmm-1.0)*beta**2/(8.0*Gmm*Pi**2)*EXP( 1.0-r2 )

  Prim(1)=(1.0+delta_T)**( 1.0 / (Gmm-1.0) )
  Prim(2)=delta_u
  Prim(3)=delta_v
  Prim(4)=(1.0+delta_T)**( Gmm / (Gmm-1.0) )

  CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[3] Unsteady isentropic vortex
  !*------------------------------------------
CASE(3)

  u_inf=1.0
  v_inf=1.0

  !*Original center of the vortex before the moevement
  xc=0.5*(MESH_X1+MESH_X0)

  !*We want to get the initial position of x before the movement

  !*x with respect to bottom-left corner
  x_wrt_BL=x-MESH_X0

  !*x with respect to bottom-left corner before movement
  x_wrt_BL_bm(1)=x_wrt_BL(1)-u_inf*t
  x_wrt_BL_bm(2)=x_wrt_BL(2)-v_inf*t

  !*This is the position before movement modulo the length of the domain
  !*NB: MODULO RESULT IS ALWAYS POSITIVE
  x_wrt_BL_bm(1)=MODULO( x_wrt_BL_bm(1), MESH_SX(1) )
  x_wrt_BL_bm(2)=MODULO( x_wrt_BL_bm(2), MESH_SX(2) )

  !*This is the initial position
  x_0=MESH_X0+x_wrt_BL_bm

  !*Distance squared from the center of the vortex at the initial time
  x_d=x_0-xc

  r2=x_d(1)**2+x_d(2)**2

  !*Vortex amplitude
  beta=5.0 !*5.0 0.1

  delta_u=beta/(2.0*PI)*EXP( 0.5*( 1.0-r2 ) )*( -x_d(2) )
  delta_v=beta/(2.0*PI)*EXP( 0.5*( 1.0-r2 ) )*x_d(1)
  delta_T=-(Gmm-1.0)*beta**2/(8.0*Gmm*Pi**2)*EXP( 1.0-r2 )

  Prim(1)=(1.0+delta_T)**( 1.0 / (Gmm-1.0) )
  Prim(2)=u_inf+delta_u
  Prim(3)=v_inf+delta_v
  Prim(4)=(1.0+delta_T)**( Gmm / (Gmm-1.0) )

  CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[4] Advection of smooth density
  !*------------------------------------------
CASE(4)
  u_inf = 1.0
  v_inf = -0.5
  p_inf = 1.0
  Prim(1)=1.0+0.5*SIN( 4.0*Pi*( x(1)+x(2)-t*(u_inf+v_inf) ) )
  Prim(2)=u_inf
  Prim(3)=v_inf
  Prim(4)=p_inf

  CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[5] Advection of smooth density sin4
  !*------------------------------------------
CASE(5)
  u_inf = 1.0
  v_inf = -0.5
  p_inf = 1.0
  Prim(1) = 2.0+SIN(2.0*PI*( x(1)+x(2)-t*(u_inf+v_inf) ))**4
  Prim(2)=u_inf
  Prim(3)=v_inf
  Prim(4)=p_inf

  CALL PrimToCons(Prim,Cons)

  !*------------------------------------------
  !*[892] Smooth periodic IC with the purpose of verifying conservation
  !*------------------------------------------
CASE(892)
  Prim(1)=8.0+0.5*SIN(2.0*Pi*x(1))
  Prim(2)=0.5+COS(8.0*Pi*x(2))
  Prim(3)=0.5+SIN(4.0*Pi*x(2))
  Prim(4)=7.0+SIN(4.0*Pi*x(1))

  CALL PrimToCons(Prim,Cons)

#endif
#ifdef EqnAcoustics
#endif
CASE DEFAULT
  ErrorMessage = "Exact function not specified"
  WRITE(*,*) ErrorMessage
  STOP
END SELECT

!-------------------------------------------------------------------------------!
CONTAINS

#if defined (SW) || defined (EqnShallowWater)
REAL FUNCTION hSmoothAuxiliary(x)
  IMPLICIT NONE
  REAL, INTENT(IN) :: x

  hSmoothAuxiliary=1.-0.5*exp(-1./atan(1.-x)**3.)

END FUNCTION

REAL FUNCTION hDerivSmoothAuxiliary(x)
  IMPLICIT NONE
  REAL, INTENT(IN) :: x

  hDerivSmoothAuxiliary=3.*0.5*exp(1./atan(x - 1.)**3.)/(atan(x - 1.)**4*((x - 1.)**2 + 1.))

END FUNCTION
#endif

END SUBROUTINE ExactFunction
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ExactFunctionWB(WhichInitialCondition,x,Cons)
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI
USE MOD_FiniteVolume2D_vars,ONLY: MIN_POSITIVE_VAR
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN) :: WhichInitialCondition
REAL,INTENT(IN)    :: x(1:nDims)
REAL,INTENT(OUT)   :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: Prim(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage

SELECT CASE (WhichInitialCondition)

  CASE DEFAULT
    ErrorMessage = "Exact WB function not specified"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

END SUBROUTINE ExactFunctionWB
!===============================================================================!
!
!
#ifdef GFWENO

!===============================================================================!
SUBROUTINE GlobalFluxTerms(t)
!-------------------------------------------------------------------------------!
USE MOD_Reconstruction     ,ONLY: WENO1_SecondSweep
USE MOD_Reconstruction     ,ONLY: WENO3_SecondSweep
USE MOD_Reconstruction     ,ONLY: WENO5_SecondSweep
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Eta
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: Gravity
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP          
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP          
USE MOD_FiniteVolume2D_vars,ONLY: nGPs         
USE MOD_FiniteVolume2D_vars,ONLY: nVar         
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX      
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1      
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0      
USE MOD_FiniteVolume2D_vars,ONLY: RX
USE MOD_FiniteVolume2D_vars,ONLY: RY
USE MOD_FiniteVolume2D_vars,ONLY: FG
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: Ftemp(nVar,nGPs, nGPs)
REAL               :: Fe, error, 
REAL               :: FluxX(1:nVar,1:nGPs,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
REAL               :: FluxY(1:nVar,1:nGPs,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
REAL               :: FFX_interface(1:2,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
REAL               :: FFY_interface(1:2,1:nVar,-nGhosts:nElemsX+nGhosts+1,-nGhosts:nElemsY+nGhosts+1)
REAL               :: Utemp1X(0:nVar,nGPs,-2*nGhosts:nElemsX+2*nGhosts+1)
REAL               :: Utemp2X(0:nVar,nGPs,nGPs)
REAL               :: Utemp1Y(0:nVar,nGPs,-2*nGhosts:nElemsY+2*nGhosts+1)
REAL               :: Utemp2Y(0:nVar,nGPs,nGPs)
REAL               :: Vtemp(nVar,nGPs,nGPs), FluxX_int, FluxY_int
INTEGER            :: ii, iGP, iVar
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!
Ftemp = 0.
Utemp = 0.
Vtemp = 0.
FF    = 0.
FluxX = 0.
FluxY = 0.

SELECT CASE (Reconstruction)
  CASE(1)
    DO jj=-2*nGhosts,nElemsY+2*nGhosts+1
      DO ii=-nGhosts,nElemsX+nGhosts+1
        CALL WENO1_SecondSweep( Eta(ii,jj-nGhosts:jj+nGhosts) , Eta(ii,jj-nGhosts:jj+nGhosts), Utemp1(0,1:nGPs,ii) )
        CALL WENO1_SecondSweep( U(1,ii,jj-nGhosts:jj+nGhosts) , Eta(ii,jj-nGhosts:jj+nGhosts), Utemp1(1,1:nGPs,ii) )
        CALL WENO1_SecondSweep( U(2,ii,jj-nGhosts:jj+nGhosts) , U(2,ii,jj-nGhosts:jj+nGhosts), Utemp1(2,1:nGPs,ii) )
        CALL WENO1_SecondSweep( U(3,ii,jj-nGhosts:jj+nGhosts) , U(3,ii,jj-nGhosts:jj+nGhosts), Utemp1(3,1:nGPs,ii) )
      END DO

      DO ii=-nGhosts,nElemsX+nGhosts+1
        DO jGP=1,nGPs
          CALL WENO1_SecondSweep( Utemp1(1,jGP,ii-nGhosts:ii+nGhosts), Utemp1(0,jGP,ii-nGhosts:ii+nGhosts), Utemp2(1,1:nGPs,jGP) )
          CALL WENO1_SecondSweep( Utemp1(2,jGP,ii-nGhosts:ii+nGhosts), Utemp1(2,jGP,ii-nGhosts:ii+nGhosts), Utemp2(2,1:nGPs,jGP) )
          CALL WENO1_SecondSweep( Utemp1(3,jGP,ii-nGhosts:ii+nGhosts), Utemp1(3,jGP,ii-nGhosts:ii+nGhosts), Utemp2(3,1:nGPs,jGP) )
          DO iGP=1,nGPs
            CALL ConsToPrim( Utemp2(1:nVar,iGP, jGP) , Vtemp(1:nVar,iGP,jGP) )
            CALL EvaluateFlux2D_X( Vtemp(1:nVar,iGP, jGP) , Ftemp(1:nVar,iGP,jGP) )
            FluxX(1:nVar,jGP,ii,jj) = FluxX(1:nVar,jGP,ii,jj) + WeightsGP(iGP) * Ftemp(1:nVar,iGP,jGP) 
          END DO 
        END DO
      END DO  
    END DO

    DO ii=-2*nGhosts,nElemsX+2*nGhosts+1
      DO jj=-nGhosts,nElemsY+nGhosts+1
        CALL WENO1_SecondSweep( Eta(ii-nGhosts:ii+nGhosts,jj) , Eta(ii-nGhosts:ii+nGhosts,jj), Utemp1(0,1:nGPs,jj) )
        CALL WENO1_SecondSweep( U(1,ii-nGhosts:ii+nGhosts,jj) , Eta(ii-nGhosts:ii+nGhosts,jj), Utemp1(1,1:nGPs,jj) )
        CALL WENO1_SecondSweep( U(2,ii-nGhosts:ii+nGhosts,jj) , U(2,ii-nGhosts:ii+nGhosts,jj), Utemp1(2,1:nGPs,jj) )
        CALL WENO1_SecondSweep( U(3,ii-nGhosts:ii+nGhosts,jj) , U(3,ii-nGhosts:ii+nGhosts,jj), Utemp1(3,1:nGPs,jj) )
      END DO

      DO jj=-nGhosts,nElemsY+nGhosts+1
        DO iGP=1,nGPs
          CALL WENO1_SecondSweep( Utemp1(1,iGP,jj-nGhosts:jj+nGhosts), Utemp1(0,jGP,jj-nGhosts:jj+nGhosts), Utemp2(1,iGP, 1:nGPs) )
          CALL WENO1_SecondSweep( Utemp1(2,iGP,jj-nGhosts:jj+nGhosts), Utemp1(2,iGP,jj-nGhosts:jj+nGhosts), Utemp2(2,iGP, 1:nGPs) )
          CALL WENO1_SecondSweep( Utemp1(3,iGP,jj-nGhosts:jj+nGhosts), Utemp1(3,iGP,jj-nGhosts:jj+nGhosts), Utemp2(3,iGP, 1:nGPs) )
          DO jGP=1,nGPs
            CALL ConsToPrim( Utemp2(1:nVar,iGP, jGP) , Vtemp(1:nVar,iGP,jGP) )
            CALL EvaluateFlux2D_Y( Vtemp(1:nVar,iGP, jGP) , Ftemp(1:nVar,iGP,jGP) )
            FluxY(1:nVar,iGP,ii,jj) = FluxY(1:nVar,iGP,ii,jj) + WeightsGP(jGP) * Ftemp(1:nVar,iGP,jGP) 
          END DO 
        END DO
      END DO  
    END DO
  CASE(3)
    DO jj=-2*nGhosts,nElemsY+2*nGhosts+1
      DO ii=-nGhosts,nElemsX+nGhosts+1
        CALL WENO3_SecondSweep( Eta(ii,jj-nGhosts:jj+nGhosts) , Eta(ii,jj-nGhosts:jj+nGhosts), Utemp1(0,1:nGPs,ii) )
        CALL WENO3_SecondSweep( U(1,ii,jj-nGhosts:jj+nGhosts) , Eta(ii,jj-nGhosts:jj+nGhosts), Utemp1(1,1:nGPs,ii) )
        CALL WENO3_SecondSweep( U(2,ii,jj-nGhosts:jj+nGhosts) , U(2,ii,jj-nGhosts:jj+nGhosts), Utemp1(2,1:nGPs,ii) )
        CALL WENO3_SecondSweep( U(3,ii,jj-nGhosts:jj+nGhosts) , U(3,ii,jj-nGhosts:jj+nGhosts), Utemp1(3,1:nGPs,ii) )
      END DO

      DO ii=-nGhosts,nElemsX+nGhosts+1
        DO jGP=1,nGPs
          CALL WENO3_SecondSweep( Utemp1(1,jGP,ii-nGhosts:ii+nGhosts), Utemp1(0,jGP,ii-nGhosts:ii+nGhosts), Utemp2(1,1:nGPs,jGP) )
          CALL WENO3_SecondSweep( Utemp1(2,jGP,ii-nGhosts:ii+nGhosts), Utemp1(2,jGP,ii-nGhosts:ii+nGhosts), Utemp2(2,1:nGPs,jGP) )
          CALL WENO3_SecondSweep( Utemp1(3,jGP,ii-nGhosts:ii+nGhosts), Utemp1(3,jGP,ii-nGhosts:ii+nGhosts), Utemp2(3,1:nGPs,jGP) )
          DO iGP=1,nGPs
            CALL ConsToPrim( Utemp2(1:nVar,iGP, jGP) , Vtemp(1:nVar,iGP,jGP) )
            CALL EvaluateFlux2D_X( Vtemp(1:nVar,iGP, jGP) , Ftemp(1:nVar,iGP,jGP) )
            FluxX(1:nVar,jGP,ii,jj) = FluxX(1:nVar,jGP,ii,jj) + WeightsGP(iGP) * Ftemp(1:nVar,iGP,jGP) 
          END DO 
        END DO
      END DO  
    END DO

    DO ii=-2*nGhosts,nElemsX+2*nGhosts+1
      DO jj=-nGhosts,nElemsY+nGhosts+1
        CALL WENO3_SecondSweep( Eta(ii-nGhosts:ii+nGhosts,jj) , Eta(ii-nGhosts:ii+nGhosts,jj), Utemp1(0,1:nGPs,jj) )
        CALL WENO3_SecondSweep( U(1,ii-nGhosts:ii+nGhosts,jj) , Eta(ii-nGhosts:ii+nGhosts,jj), Utemp1(1,1:nGPs,jj) )
        CALL WENO3_SecondSweep( U(2,ii-nGhosts:ii+nGhosts,jj) , U(2,ii-nGhosts:ii+nGhosts,jj), Utemp1(2,1:nGPs,jj) )
        CALL WENO3_SecondSweep( U(3,ii-nGhosts:ii+nGhosts,jj) , U(3,ii-nGhosts:ii+nGhosts,jj), Utemp1(3,1:nGPs,jj) )
      END DO

      DO jj=-nGhosts,nElemsY+nGhosts+1
        DO iGP=1,nGPs
          CALL WENO3_SecondSweep( Utemp1(1,iGP,jj-nGhosts:jj+nGhosts), Utemp1(0,jGP,jj-nGhosts:jj+nGhosts), Utemp2(1,iGP, 1:nGPs) )
          CALL WENO3_SecondSweep( Utemp1(2,iGP,jj-nGhosts:jj+nGhosts), Utemp1(2,iGP,jj-nGhosts:jj+nGhosts), Utemp2(2,iGP, 1:nGPs) )
          CALL WENO3_SecondSweep( Utemp1(3,iGP,jj-nGhosts:jj+nGhosts), Utemp1(3,iGP,jj-nGhosts:jj+nGhosts), Utemp2(3,iGP, 1:nGPs) )
          DO jGP=1,nGPs
            CALL ConsToPrim( Utemp2(1:nVar,iGP, jGP) , Vtemp(1:nVar,iGP,jGP) )
            CALL EvaluateFlux2D_Y( Vtemp(1:nVar,iGP, jGP) , Ftemp(1:nVar,iGP,jGP) )
            FluxY(1:nVar,iGP,ii,jj) = FluxY(1:nVar,iGP,ii,jj) + WeightsGP(jGP) * Ftemp(1:nVar,iGP,jGP) 
          END DO 
        END DO
      END DO  
    END DO
  CASE(4)
    DO jj=-2*nGhosts,nElemsY+2*nGhosts+1
      DO ii=-nGhosts,nElemsX+nGhosts+1
        CALL WENO5_SecondSweep( Eta(ii,jj-nGhosts:jj+nGhosts) , Eta(ii,jj-nGhosts:jj+nGhosts), Utemp1(0,1:nGPs,ii) )
        CALL WENO5_SecondSweep( U(1,ii,jj-nGhosts:jj+nGhosts) , Eta(ii,jj-nGhosts:jj+nGhosts), Utemp1(1,1:nGPs,ii) )
        CALL WENO5_SecondSweep( U(2,ii,jj-nGhosts:jj+nGhosts) , U(2,ii,jj-nGhosts:jj+nGhosts), Utemp1(2,1:nGPs,ii) )
        CALL WENO5_SecondSweep( U(3,ii,jj-nGhosts:jj+nGhosts) , U(3,ii,jj-nGhosts:jj+nGhosts), Utemp1(3,1:nGPs,ii) )
      END DO

      DO ii=-nGhosts,nElemsX+nGhosts+1
        DO jGP=1,nGPs
          CALL WENO5_SecondSweep( Utemp1(1,jGP,ii-nGhosts:ii+nGhosts), Utemp1(0,jGP,ii-nGhosts:ii+nGhosts), Utemp2(1,1:nGPs,jGP) )
          CALL WENO5_SecondSweep( Utemp1(2,jGP,ii-nGhosts:ii+nGhosts), Utemp1(2,jGP,ii-nGhosts:ii+nGhosts), Utemp2(2,1:nGPs,jGP) )
          CALL WENO5_SecondSweep( Utemp1(3,jGP,ii-nGhosts:ii+nGhosts), Utemp1(3,jGP,ii-nGhosts:ii+nGhosts), Utemp2(3,1:nGPs,jGP) )
          DO iGP=1,nGPs
            CALL ConsToPrim( Utemp2(1:nVar,iGP, jGP) , Vtemp(1:nVar,iGP,jGP) )
            CALL EvaluateFlux2D_X( Vtemp(1:nVar,iGP, jGP) , Ftemp(1:nVar,iGP,jGP) )
            FluxX(1:nVar,jGP,ii,jj) = FluxX(1:nVar,jGP,ii,jj) + WeightsGP(iGP) * Ftemp(1:nVar,iGP,jGP) 
          END DO 
        END DO
      END DO  
    END DO

    DO ii=-2*nGhosts,nElemsX+2*nGhosts+1
      DO jj=-nGhosts,nElemsY+nGhosts+1
        CALL WENO5_SecondSweep( Eta(ii-nGhosts:ii+nGhosts,jj) , Eta(ii-nGhosts:ii+nGhosts,jj), Utemp1(0,1:nGPs,jj) )
        CALL WENO5_SecondSweep( U(1,ii-nGhosts:ii+nGhosts,jj) , Eta(ii-nGhosts:ii+nGhosts,jj), Utemp1(1,1:nGPs,jj) )
        CALL WENO5_SecondSweep( U(2,ii-nGhosts:ii+nGhosts,jj) , U(2,ii-nGhosts:ii+nGhosts,jj), Utemp1(2,1:nGPs,jj) )
        CALL WENO5_SecondSweep( U(3,ii-nGhosts:ii+nGhosts,jj) , U(3,ii-nGhosts:ii+nGhosts,jj), Utemp1(3,1:nGPs,jj) )
      END DO

      DO jj=-nGhosts,nElemsY+nGhosts+1
        DO iGP=1,nGPs
          CALL WENO5_SecondSweep( Utemp1(1,iGP,jj-nGhosts:jj+nGhosts), Utemp1(0,jGP,jj-nGhosts:jj+nGhosts), Utemp2(1,iGP, 1:nGPs) )
          CALL WENO5_SecondSweep( Utemp1(2,iGP,jj-nGhosts:jj+nGhosts), Utemp1(2,iGP,jj-nGhosts:jj+nGhosts), Utemp2(2,iGP, 1:nGPs) )
          CALL WENO5_SecondSweep( Utemp1(3,iGP,jj-nGhosts:jj+nGhosts), Utemp1(3,iGP,jj-nGhosts:jj+nGhosts), Utemp2(3,iGP, 1:nGPs) )
          DO jGP=1,nGPs
            CALL ConsToPrim( Utemp2(1:nVar,iGP, jGP) , Vtemp(1:nVar,iGP,jGP) )
            CALL EvaluateFlux2D_Y( Vtemp(1:nVar,iGP, jGP) , Ftemp(1:nVar,iGP,jGP) )
            FluxY(1:nVar,iGP,ii,jj) = FluxY(1:nVar,iGP,ii,jj) + WeightsGP(jGP) * Ftemp(1:nVar,iGP,jGP) 
          END DO 
        END DO
      END DO  
    END DO

  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

FluxX(2,:,:,:) = FluxX(2,:,:,:) + RX  ! (FluxX)_{jq,\bar i}= (Fx+Rx)_{jq,\bar i} 
FluxY(3,:,:,:) = FluxY(3,:,:,:) + RY
FFX = 0.
FFY = 0.
FFX_interface = 0.
FFY_interface = 0.

FG = 0.

DO ii=-nGhosts,nElemsX+nGhosts+1
  DO jj=-nGhosts,nElemsY+nGhosts+1
      DO iVar=1,nVar
         ! Build FFX
         CALL SourceInterpIntegralCoeff(FluxX(iVar,1:nGPs,ii,jj), FluxX_int)
         FFX(iVar,ii,jj) = FFX_interface(1,iVar,ii,jj) + FluxX_int * MESH_DX(2)  !\int^{y_j} FX + RX until cell average 

         FluxX_int = 0.
         DO jGP = 1,nGPs
            FluxX_int  = FluxX_int + WeightsGP(jGP) * FluxX(iVar,jGP,ii,jj) 
         END DO

        ! \int^{y_j+1/2} FX + RX until right interface
         FFX_interface(2,iVar,ii,jj) = FFX_interface(1,iVar,ii,jj) + MESH_DX(2)* FluxX_int
         FFX_interface(1,iVar,ii,jj+1) = FFX_interface(2,iVar,ii,jj)  ! no jump
         
         ! Build FFY
         CALL SourceInterpIntegralCoeff(FluxY(iVar,1:nGPs,ii,jj), FluxY_int)
         FFY(nVar,ii,jj) = FFY_interface(1,iVar,ii+1,jj) + FluxY_int * MESH_DX(1) !\int^{x_i} FY + RY until cell average

         FluxY_int = 0.
         DO iGP = 1,nGPs
            FluxY_int  = FluxY_int + WeightsGP(iGP) * FluxY(iVar,iGP,ii,jj) 
         END DO
         
         ! \int^{x_i+1/2} FY + RY until right interface
         FFY_interface(2,iVar,ii,jj) = FFY_interface(1,iVar,ii,jj) + MESH_DX(1)*FluxY_int
         FFY_interface(1,iVar,ii+1,jj) = FFY_interface(2,iVar,ii,jj)  ! no jump
      END DO
   END DO  
END DO

! 2D global flux computation
FG = FFX + FFY


!-------------------------------------------------------------------------------!
END SUBROUTINE GlobalFluxTerms
!===============================================================================!
!
!
!===============================================================================!
SUBROUTINE SourceTerms(t)
!-------------------------------------------------------------------------------!
  USE MOD_Reconstruction     ,ONLY: MUSCL
  USE MOD_Reconstruction     ,ONLY: WENO1_SecondSweep
  USE MOD_Reconstruction     ,ONLY: WENO3_SecondSweep
  USE MOD_Reconstruction     ,ONLY: WENO5_SecondSweep
  USE MOD_Reconstruction     ,ONLY: WENO1_FirstSweep
  USE MOD_Reconstruction     ,ONLY: WENO3_FirstSweep
  USE MOD_Reconstruction     ,ONLY: WENO5_FirstSweep
  USE MOD_FiniteVolume2D_vars,ONLY: S
  USE MOD_FiniteVolume2D_vars,ONLY: RX
  USE MOD_FiniteVolume2D_vars,ONLY: RY
  USE MOD_FiniteVolume2D_vars,ONLY: RX_interface
  USE MOD_FiniteVolume2D_vars,ONLY: RY_interface
  USE MOD_FiniteVolume2D_vars,ONLY: nGPs
  USE MOD_FiniteVolume2D_vars,ONLY: nVar
  USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
  USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
  USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
  USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
  USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
  USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
  USE MOD_FiniteVolume2D_vars,ONLY: U
  USE MOD_FiniteVolume2D_vars,ONLY: V
  USE MOD_FiniteVolume2D_vars,ONLY: V_reconstructed
  USE MOD_FiniteVolume2D_vars,ONLY: Bath
  USE MOD_FiniteVolume2D_vars,ONLY: Eta
  USE MOD_FiniteVolume2D_vars,ONLY: EtaL
  USE MOD_FiniteVolume2D_vars,ONLY: EtaR
  USE MOD_FiniteVolume2D_vars,ONLY: EtaB
  USE MOD_FiniteVolume2D_vars,ONLY: EtaT
  USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
  USE MOD_FiniteVolume2D_vars,ONLY: Gravity
  USE MOD_FiniteVolume2D_vars,ONLY: Bath_interfaceX
  USE MOD_FiniteVolume2D_vars,ONLY: Bath_interfaceY
  USE MOD_FiniteVolume2D_vars,ONLY: BathymetryFlag
  USE MOD_FiniteVolume2D_vars,ONLY: manning, coriolis
  !-------------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------------!
  ! >> FORMAL ARGUMENTS                                                           !
  !-------------------------------------------------------------------------------!
  REAL,INTENT(IN)  :: t
  !-------------------------------------------------------------------------------!
  ! >> LOCAL VARIABLES                                                            !
  !-------------------------------------------------------------------------------!
  REAL             :: Source_weights(1:nVar,1:nGPs,1:nGPs,-nGhosts:nElemsX+1+nGhosts,-nGhosts:nElemsY+1+nGhosts)
  REAL             :: Eta_quad(1:nGPs,1:nGPs,-nGhosts:nElemsX+1+nGhosts,-nGhosts:nElemsY+1+nGhosts)
  REAL             :: U_quad(1:nVar, 1:nGPs,1:nGPs,-nGhosts:nElemsX+1+nGhosts,-nGhosts:nElemsY+1+nGhosts)
  REAL             :: U_temp(0:nVar, 1:nGPs,-nGhosts:nElemsX+1+nGhosts,-nGhosts:nElemsY+1+nGhosts) ! bath in var=0, eta in var = 1, qx,qy in var=2,3
  REAL             :: S_int(1:nVar, 1:nGPs,-nGhosts:nElemsX+1+nGhosts,-nGhosts:nElemsY+1+nGhosts) ! bath in var=0, eta in var = 1, qx,qy in var=2,3
  REAL             :: Bath_quad(1:nGPs,1:nGPs,-nGhosts:nElemsX+1+nGhosts,-nGhosts:nElemsY+1+nGhosts)
  REAL             :: bath_interf_local
  REAL             :: BathDerivX_quad(1:nGPs,1:nGPs,-nGhosts:nElemsX+1+nGhosts,-nGhosts:nElemsY+1+nGhosts)
  REAL             :: BathDerivY_quad(1:nGPs,1:nGPs,-nGhosts:nElemsX+1+nGhosts,-nGhosts:nElemsY+1+nGhosts)
  REAL             :: B2X(1:nGPs,-nGhosts:nElemsX+1+nGhosts,-nGhosts:nElemsY+1+nGhosts)
  REAL             :: B2Y(1:nGPs,-nGhosts:nElemsX+1+nGhosts,-nGhosts:nElemsY+1+nGhosts)
  REAL             :: AveEta, jumpB2, jumpB
  INTEGER          :: ii, jj, iVar, iGP, jGP
  CHARACTER(LEN=255) :: ErrorMessage
  !-------------------------------------------------------------------------------!

  ! define cell averages of Eta
  Eta(:) = V(1,:) + Bath(:)

  ! compute source integral
  S = 0.0
  B2X = 0.0
  B2Y = 0.0
  BathDerivX_quad = 0.0
  BathDerivY_quad = 0.0
  Eta_quad = 0.0
  Bath_quad = 0.0
  Bath_interfaceX = 0.0
  Bath_interfaceY = 0.0
  Source_weights = 0.0

  !int_{iixjj} S(1:nVar,ii,jj) dxdy

  SELECT CASE (Reconstruction)

  CASE(1)


    DO jj=-nGhosts,nElemsY+nGhosts+1
      DO ii=-2*nGhosts,nElemsX+2*nGhosts+1
        CALL WENO1_SecondSweep( Bath(ii,jj-nGhosts:jj+nGhosts), Eta(ii,jj-nGhosts:jj+nGhosts), U_temp(0,1:nGPs,ii,jj) )
        CALL WENO1_SecondSweep( Eta(ii,jj-nGhosts:jj+nGhosts) , Eta(ii,jj-nGhosts:jj+nGhosts), U_temp(1,1:nGPs,ii,jj) )
        DO iVar=2:nVar
              CALL WENO1_SecondSweep( U(2,ii,jj-nGhosts:jj+nGhosts) , U(2,ii,jj-nGhosts:jj+nGhosts), U_temp(iVar,1:nGPs,ii,jj) )
        END DO
      END DO
    END DO

    DO jj=-nGhosts,nElemsY+nGhosts+1
        DO ii=-nGhosts,nElemsX+nGhosts+1
          DO jGP = 1:nGPs
            CALL WENO1_SecondSweep( U_temp(0,jGP, ii-nGhosts:ii+nGhosts,jj) , U_temp(1,jGP, ii-nGhosts:ii+nGhosts,jj), Bath_quad(1:nGPs,jGP,ii,jj) )
            CALL WENO1_SecondSweep( U_temp(1,jGP, ii-nGhosts:ii+nGhosts,jj) , U_temp(1,jGP, ii-nGhosts:ii+nGhosts,jj), Eta_quad(1:nGPs,jGP,ii,jj) )
            DO iVar=2:nVar
              CALL WENO1_SecondSweep( U_temp(iVar,jGP, ii-nGhosts:ii+nGhosts,jj) , U_temp(iVar,jGP, ii-nGhosts:ii+nGhosts,jj), U_quad(iVar,1:nGPs,jGP,ii,jj) )
            END DO
          END DO
        END DO
    END DO

  CASE(3)


    DO jj=-nGhosts,nElemsY+nGhosts+1
      DO ii=-2*nGhosts,nElemsX+2*nGhosts+1
        CALL WENO3_SecondSweep( Bath(ii,jj-nGhosts:jj+nGhosts), Eta(ii,jj-nGhosts:jj+nGhosts), U_temp(0,1:nGPs,ii,jj) )
        CALL WENO3_SecondSweep( Eta(ii,jj-nGhosts:jj+nGhosts) , Eta(ii,jj-nGhosts:jj+nGhosts), U_temp(1,1:nGPs,ii,jj) )
        DO iVar=2:nVar
              CALL WENO3_SecondSweep( U(2,ii,jj-nGhosts:jj+nGhosts) , U(2,ii,jj-nGhosts:jj+nGhosts), U_temp(iVar,1:nGPs,ii,jj) )
        END DO
      END DO
    END DO

    DO jj=-nGhosts,nElemsY+nGhosts+1
        DO ii=-nGhosts,nElemsX+nGhosts+1
          DO jGP = 1:nGPs
            CALL WENO3_SecondSweep( U_temp(0,jGP, ii-nGhosts:ii+nGhosts,jj) , U_temp(1,jGP, ii-nGhosts:ii+nGhosts,jj), Bath_quad(1:nGPs,jGP,ii,jj) )
            CALL WENO3_SecondSweep( U_temp(1,jGP, ii-nGhosts:ii+nGhosts,jj) , U_temp(1,jGP, ii-nGhosts:ii+nGhosts,jj), Eta_quad(1:nGPs,jGP,ii,jj) )
            DO iVar=2:nVar
              CALL WENO3_SecondSweep( U_temp(iVar,jGP, ii-nGhosts:ii+nGhosts,jj) , U_temp(iVar,jGP, ii-nGhosts:ii+nGhosts,jj), U_quad(iVar,1:nGPs,jGP,ii,jj) )
            END DO
          END DO
        END DO
    END DO


  CASE(4)


    DO jj=-nGhosts,nElemsY+nGhosts+1
      DO ii=-2*nGhosts,nElemsX+2*nGhosts+1
        CALL WENO5_SecondSweep( Bath(ii,jj-nGhosts:jj+nGhosts), Eta(ii,jj-nGhosts:jj+nGhosts), U_temp(0,1:nGPs,ii,jj) )
        CALL WENO5_SecondSweep( Eta(ii,jj-nGhosts:jj+nGhosts) , Eta(ii,jj-nGhosts:jj+nGhosts), U_temp(1,1:nGPs,ii,jj) )
        DO iVar=2:nVar
              CALL WENO5_SecondSweep( U(2,ii,jj-nGhosts:jj+nGhosts) , U(2,ii,jj-nGhosts:jj+nGhosts), U_temp(iVar,1:nGPs,ii,jj) )
        END DO
      END DO
    END DO

    DO jj=-nGhosts,nElemsY+nGhosts+1
        DO ii=-nGhosts,nElemsX+nGhosts+1
          DO jGP = 1:nGPs
            CALL WENO5_SecondSweep( U_temp(0,jGP, ii-nGhosts:ii+nGhosts,jj) , U_temp(1,jGP, ii-nGhosts:ii+nGhosts,jj), Bath_quad(1:nGPs,jGP,ii,jj) )
            CALL WENO5_SecondSweep( U_temp(1,jGP, ii-nGhosts:ii+nGhosts,jj) , U_temp(1,jGP, ii-nGhosts:ii+nGhosts,jj), Eta_quad(1:nGPs,jGP,ii,jj) )
            DO iVar=2:nVar
              CALL WENO5_SecondSweep( U_temp(iVar,jGP, ii-nGhosts:ii+nGhosts,jj) , U_temp(iVar,jGP, ii-nGhosts:ii+nGhosts,jj), U_quad(iVar,1:nGPs,jGP,ii,jj) )
            END DO
          END DO
        END DO
    END DO


  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE(*,*) ErrorMessage
    STOP
  END SELECT


  DO ii=-nGhosts,nElemsX+nGhosts+1
    DO jj=-nGhosts,nElemsY+nGhosts+1
      DO jGP=1,nGPs
          CALL Bath2Derivative(Bath_quad(1:nGPs,jGP, ii, jj),BathDerivX_quad(1:nGPs,jGP,ii,jj))
          BathDerivX_quad(1:nGPs,jGP,ii,jj)=BathDerivX_quad(1:nGPs,jGP,ii,jj)/MESH_DX(1)
      END DO

      DO iGP=1,nGPs
        CALL Bath2Derivative(Bath_quad(iGP,1:nGPs, ii, jj),BathDerivY_quad(iGP,1:nGPs,ii,jj))
        BathDerivY_quad(iGP,1:nGPs,ii,jj) = BathDerivY_quad(iGP,1:nGPs,ii,jj)/MESH_DX(2)
      END DO

      Source_weights(2,1:nGPs,1:nGPs,ii,jj) = Gravity*Eta_quad(1:nGPs,1:nGPs,ii,jj) * BathDerivX_quad(1:nGPs,1:nGPs,ii,jj) +&
        Gravity*manning**2*abs(U_quad(2,1:nGPs,1:nGPs,ii,jj))*U_quad(2,1:nGPs,1:nGPs,ii,jj)/(Eta_quad(1:nGPs,1:nGPs,ii,jj)-Bath_quad(1:nGPs,1:nGPs,ii,jj))**(7./3.) &
        + coriolis*U_quad(3,1:nGPs,1:nGPs,ii,jj)

      Source_weights(3,1:nGPs,1:nGPs,ii,jj) = Gravity*Eta_quad(1:nGPs,1:nGPs,ii,jj) * BathDerivY_quad(1:nGPs,1:nGPs,ii,jj) +&
        Gravity*manning**2*abs(U_quad(3,1:nGPs,1:nGPs,ii,jj))*U_quad(3,1:nGPs,1:nGPs,ii,jj)/(Eta_quad(1:nGPs,1:nGPs,ii,jj)-Bath_quad(1:nGPs,1:nGPs,ii,jj))**(7./3.) &
        - coriolis*U_quad(2,1:nGPs,1:nGPs,ii,jj)

    END DO
  END DO


  ! compute \sum_q { w_q * 1/2 * ( b(x_q)^2 - b(x_{i-1/2})^2 ) }
  B2X(1:nVar,:,:,:) = 0.0
  B2Y(1:nVar,:,:,:) = 0.0

  ! Bath_interfaceX(1,jGP,ii,jj) is the reconstruction at the left interface of the cell
  ! Bath_interfaceX(2,jGP,ii,jj) is the reconstruction at the right interface of the cell
  ! Bath_interfaceY(1,iGP,ii,jj) is the reconstruction at the bottom interface of the cell
  ! Bath_interfaceY(2,iGP,ii,jj) is the reconstruction at the top interface of the cell
  DO ii=-nGhosts,nElemsX+1+nGhosts
    DO jj=-nGhosts,nElemsX+1+nGhosts
      DO jGP=1,nGPs
        CALL Bath2Interfaces( Bath_quad(1:nGPs,jGP,ii,jj) , Bath_interfaceX(1:2,jGP, ii,jj) )
        CALL Bath2Interfaces( Eta_quad(1:nGPs,jGP,ii,jj) , Eta_interfaceX(1:2,jGP, ii,jj) )
      END DO
      DO iGP=1,nGPs
        CALL Bath2Interfaces( Bath_quad(iGP,1:nGPs,ii,jj) , Bath_interfaceY(1:2,iGP, ii,jj) )
        CALL Bath2Interfaces( Eta_quad(iGP,1:nGPs,ii,jj) ,  Eta_interfaceY(1:2,iGP, ii,jj) )
      END DO
      !
      DO jGP=1,nGPs
        DO iGP=1,nGPs
          B2X(jGP, ii,jj) = B2X(jGP, ii,jj) + WeightsGP(iGP)  * 0.5 * Bath_quad(iGP,jGP,ii,jj)**2 * Gravity
        END DO
        B2X(jGP, ii,jj) = B2X(jGP, ii,jj) - 0.5 * Bath_interfaceX(1,jGP,ii,jj)**2 * Gravity
      END DO

      DO iGP=1,nGPs
        DO jGP=1,nGPs
          B2Y(iGP, ii,jj) = B2Y(iGP, ii,jj) + WeightsGP(jGP)  * 0.5 * Bath_quad(iGP,jGP,ii,jj)**2 * Gravity
        END DO
        B2Y(iGP, ii,jj) = B2Y(iGP, ii,jj) - 0.5 * Bath_interfaceY(1,iGP,ii,jj)**2 * Gravity
      END DO

    END DO
  END DO

  DO ii=-nGhosts,nElemsX+1+nGhosts
    DO jj=-nGhosts,nElemsX+1+nGhosts
      DO iGP=1,nGPs
        DO jGP=1,nGPs
          S_int(2,jGP,ii,jj) = S(2,jGP,ii,jj) + WeightsGP(iGP) * Source_weights(2,iGP,jGP,ii,jj)
          S_int(3,iGP,ii,jj) = S(3,iGP,ii,jj) + WeightsGP(jGP) * Source_weights(3,iGP,jGP,ii,jj)
        ENDDO
      ENDDO
      DO jGP=1,nGPs
        S_int(2,jGP,ii,jj) = S(2,jGP,ii,jj) - 0.5 * Gravity * ( Bath_interfaceX(2,jGP,ii,jj)**2 - Bath_interfaceX(1,jGP,ii,jj)**2 ) / MESH_DX(1)
      ENDDO
      DO iGP=1,nGPs
        S_int(3,iGP,ii,jj) = S(3,iGP,ii,jj) - 0.5 * Gravity * ( Bath_interfaceY(2,iGP,ii,jj)**2 - Bath_interfaceY(1,iGP,ii,jj)**2 ) / MESH_DX(2)
      ENDDO
    ENDDO
  ENDDO

  RX = 0.
  RY = 0.
  RX_interface = 0.
  RY_interface = 0.



  DO ii=-nGhosts,nElemsX+1+nGhosts
    DO jj=-nGhosts,nElemsX+1+nGhosts
      DO jGP=1,nGPs
        CALL SourceInterpIntegralCoeff(Source_weights(2,1:nGPs,jGP,ii,jj), RX(jGP,ii,jj) )
        RX(jGP,ii,jj) = RX_interface(1,jGP,ii,jj) + RX(jGP,ii,jj) * MESH_DX(1) - B2X(jGP, ii,jj)
        RX_interface(2,jGP,ii,jj) = RX_interface(1,jGP,ii,jj) + MESH_DX(1)*S_int(2,jGP,ii,jj)
        jumpB  = Gravity * (Bath_interfaceX(1,jGP,ii+1,jj) - Bath_interfaceX(2,jGP,ii,jj))
        jumpB2 = 0.5* Gravity * (Bath_interfaceX(1,jGP,ii+1,jj)**2 - Bath_interfaceX(2,jGP,ii,jj)**2)
        AveEta = 0.5 * (Eta_interfaceX(1,jGP,ii+1,jj) +Eta_interfaceX(2,jGP,ii,jj))
        RX_interface(1,jGP,ii+1,jj) =RX_interface(2,jGP,ii,jj) +(AveEta*jumpB-jumpB2) 
      END DO
      DO iGP=1,nGPs
        CALL SourceInterpIntegralCoeff(Source_weights(3,jGP,1:nGPs,ii,jj), RY(iGP,ii,jj) )
        RY(iGP,ii,jj) = RY_interface(1,iGP,ii,jj) + RY(iGP,ii,jj) * MESH_DX(2) - B2Y(iGP, ii,jj)
        RY_interface(2,iGP,ii,jj) = RY_interface(1,iGP,ii,jj) + MESH_DX(2)*S_int(3,iGP,ii,jj)
        jumpB  = Gravity * (Bath_interfaceY(1,iGP,ii,jj+1) - Bath_interfaceY(2,iGP,ii,jj))
        jumpB2 = 0.5* Gravity * (Bath_interfaceY(1,iGP,ii,jj+1)**2 - Bath_interfaceY(2,iGP,ii,jj)**2)
        AveEta = 0.5 * (Eta_interfaceY(1,iGP,ii,jj+1) +Eta_interfaceY(2,iGP,ii,jj))
        RY_interface(1,iGP,ii,jj+1) =RY_interface(2,iGP,ii,jj) +(AveEta*jumpB-jumpB2) 
      END DO
    END DO
  END DO


!-------------------------------------------------------------------------------!
END SUBROUTINE SourceTerms
!===============================================================================!
!
!
!
!
!
!
!===============================================================================!
SUBROUTINE Bath2Interfaces(b_quad,BathInterface)
   !-------------------------------------------------------------------------------!
   USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
   USE MOD_FiniteVolume2D_vars,ONLY: nGPs
   !-------------------------------------------------------------------------------!
   IMPLICIT NONE
   !-------------------------------------------------------------------------------!
   ! >> FORMAL ARGUMENTS                                                           !
   !-------------------------------------------------------------------------------!
   REAL,INTENT(IN)  :: b_quad(1:nGPs)
   REAL,INTENT(OUT) :: BathInterface(2)
   CHARACTER(LEN=255) :: ErrorMessage
   !-------------------------------------------------------------------------------!

   !------------------------------------------------------------------!
   ! Coefficients interpolation polynomials in borders                !
   !------------------------------------------------------------------!

   SELECT CASE (Reconstruction)
    CASE(1,2)
      BathInterface(1) = b_quad(1)
      BathInterface(2) = b_quad(1)
    CASE(3)
      BathInterface(1) = +1.3660254037844386 * b_quad(1)-0.3660254037844387 * b_quad(2)
      BathInterface(2) = -0.3660254037844387 * b_quad(1)+1.3660254037844386 * b_quad(2)
    CASE(4)
      BathInterface(1) = +1.5267881254572668 * b_quad(1)-0.8136324494869273 * b_quad(2)+0.4007615203116504 * b_quad(3)-0.1139171962819899 * b_quad(4)
      BathInterface(2) = -0.1139171962819899 * b_quad(1)+0.4007615203116504 * b_quad(2)-0.8136324494869273 * b_quad(3)+1.5267881254572668 * b_quad(4)
    CASE DEFAULT
      ErrorMessage = "Reconstruction not implemented in Source Coefficients"
      WRITE(*,*) ErrorMessage
      STOP
   END SELECT



END SUBROUTINE Bath2Interfaces
!===============================================================================!
!
!
!
!
!===============================================================================!
SUBROUTINE SourceInterpIntegralCoeff(S_quad, R_ave)
!-------------------------------------------------------------------------------!
! This function provides the coeffienct to compute the cell average of the the integral
! of S_quad: meaning \sum_q w_q \int_{x_i-1/2}^xq L_\theta(x) S(x_theta) 
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: S_quad(1:nGPs)
REAL,INTENT(OUT) :: R_ave
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
CHARACTER(LEN=255) :: ErrorMessage


SELECT CASE (Reconstruction)
    CASE(0,1,2)
      R_ave = +0.5000000000000000 * S_quad(1)
    CASE(3)
      R_ave = +0.3943375672974064 * S_quad(1)+0.1056624327025936 * S_quad(2)
    CASE(4)
      R_ave = +0.1618513208623103 * S_quad(1)+0.2184655362953806 * S_quad(2)+0.1076070411358925 * S_quad(3)+0.0120761017064166 * S_quad(4)

    CASE DEFAULT
      ErrorMessage = "Reconstruction not implemented in Source Coefficients"
      WRITE(*,*) ErrorMessage
      STOP
  END SELECT


!-------------------------------------------------------------------------------!
END SUBROUTINE SourceInterpIntegralCoeff

 !===============================================================================!
 !
#else
!
!===============================================================================!
   SUBROUTINE SourceTerms(t)
!-------------------------------------------------------------------------------!
      USE MOD_Reconstruction     ,ONLY: MUSCL
      USE MOD_Reconstruction     ,ONLY: WENO3_SecondSweep
      USE MOD_Reconstruction     ,ONLY: WENO5_SecondSweep
      USE MOD_Reconstruction     ,ONLY: WENO7_SecondSweep
! USE MOD_Reconstruction     ,ONLY: WENO9_SecondSweep
! USE MOD_Reconstruction     ,ONLY: WENO11_SecondSweep
! USE MOD_Reconstruction     ,ONLY: WENO13_SecondSweep
      USE MOD_FiniteVolume2D_vars,ONLY: S
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: nDims
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
      USE MOD_FiniteVolume2D_vars,ONLY: MeshGP
      USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP
      USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
      USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
      USE MOD_FiniteVolume2D_vars,ONLY: U
      USE MOD_FiniteVolume2D_vars,ONLY: V
      USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
      USE MOD_FiniteVolume2D_vars,ONLY: source_flag
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: S_in_qp(1:nVar,nGPs,nGPs,nElemsX,nElemsY) !*Source in quadrature points
      REAL             :: Vtemp(1:nVar,-nGhosts:nElemsX+nGhosts+1,1:nElemsY,1:nGPs)
      REAL             :: Vtemp2(1:nVar,1:nGPs,1:nGPs,1:nElemsX,1:nElemsY)
      INTEGER          :: ii, jj, iVar, iGP, jGP
      CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

      S = 0.0 ! S(1:nVar,nElemsX,nElemsY)

!int_{iixjj} S(1:nVar,ii,jj) dxdy

      IF (source_flag .GT. 0) THEN
         SELECT CASE (Reconstruction)
          CASE(1,2,20,21,22,23,24,25)
            DO jj=1,nElemsY
               DO ii=1,nElemsX
                  S_in_qp(1:nVar,nGPs,nGPs,ii,jj) = SourceFunc( U(1:nVar,ii,jj) , MeshBary(:,ii,jj) )
               END DO
            END DO
          CASE(3)
            DO iVar=1,nVar
               DO jj=1,nElemsY
                  DO ii=-nGhosts,nElemsX+nGhosts+1
                     CALL WENO3_SecondSweep( U(iVar,ii,jj-nGhosts:jj+nGhosts) , Vtemp(iVar,ii,jj,1:nGPs) )
                  END DO
               END DO
            END DO
            DO jj=1,nElemsY
               DO ii=1,nElemsX
                  DO iGP=1,nGPs
                     DO iVar=1,nVar
                        CALL WENO3_SecondSweep( Vtemp(iVar,ii-nGhosts:ii+nGhosts,jj,iGP) , Vtemp2(iVar,1:nGPs,iGP,ii,jj) )
                     END DO
                     DO jGP=1,nGPs
                        S_in_qp(1:nVar,jGP,iGP,ii,jj) = SourceFunc( Vtemp2(1:nVar,jGP,iGP,ii,jj) , MeshGP(:,ii,jj,jGP,iGP)  )
                     END DO
                  END DO
               END DO
            END DO
          CASE(4,5)
            DO iVar=1,nVar
               DO jj=1,nElemsY
                  DO ii=-nGhosts,nElemsX+nGhosts+1
                     CALL WENO5_SecondSweep( U(iVar,ii,jj-nGhosts:jj+nGhosts) , Vtemp(iVar,ii,jj,1:nGPs) )
                  END DO
               END DO
            END DO
            DO jj=1,nElemsY
               DO ii=1,nElemsX
                  DO iGP=1,nGPs
                     DO iVar=1,nVar
                        CALL WENO5_SecondSweep( Vtemp(iVar,ii-nGhosts:ii+nGhosts,jj,iGP) , Vtemp2(iVar,1:nGPs,iGP,ii,jj) )
                     END DO
                     DO jGP=1,nGPs
                        S_in_qp(1:nVar,jGP,iGP,ii,jj) = SourceFunc( Vtemp2(1:nVar,jGP,iGP,ii,jj) , MeshGP(:,ii,jj,jGP,iGP)  )
                     END DO
                  END DO
               END DO
            END DO
          CASE(7)
            DO iVar=1,nVar
               DO jj=1,nElemsY
                  DO ii=-nGhosts,nElemsX+nGhosts+1
                     CALL WENO7_SecondSweep( U(iVar,ii,jj-nGhosts:jj+nGhosts) , Vtemp(iVar,ii,jj,1:nGPs) )
                  END DO
               END DO
            END DO
            DO jj=1,nElemsY
               DO ii=1,nElemsX
                  DO iGP=1,nGPs
                     DO iVar=1,nVar
                        CALL WENO7_SecondSweep( Vtemp(iVar,ii-nGhosts:ii+nGhosts,jj,iGP) , Vtemp2(iVar,1:nGPs,iGP,ii,jj) )
                     END DO
                     DO jGP=1,nGPs
                        S_in_qp(1:nVar,jGP,iGP,ii,jj) = SourceFunc( Vtemp2(1:nVar,jGP,iGP,ii,jj) , MeshGP(:,ii,jj,jGP,iGP)  )
                     END DO
                  END DO
               END DO
            END DO
          CASE DEFAULT
            ErrorMessage = "Reconstruction not implemented in Source"
            WRITE(*,*) ErrorMessage
            STOP
         END SELECT


         DO jj=1,nElemsY
            DO ii=1,nElemsX
               DO iGP=1,nGPs
                  DO jGP=1,nGPs
                     S(1:nVar,ii,jj) = S(1:nVar,ii,jj) + WeightsGP(iGP,jGP) * S_in_qp(1:nVar,iGP,jGP,ii,jj)
                  END DO
               END DO
            END DO
         END DO

      END IF
!-------------------------------------------------------------------------------!
   END SUBROUTINE SourceTerms
!===============================================================================!
!
#endif
!
!
!===============================================================================!
   FUNCTION SourceFunc(Q,X) RESULT(S) !*ALERT, it should be in conserved variables but I never tested it
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nDims
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      IMPLICIT NONE
      REAL, DIMENSION(1:nVar), INTENT(IN)  :: Q
      REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
      REAL, DIMENSION(1:nVar) :: S

#ifdef EqnEuler
      S(1) = 0.
      S(2) = -Q(1)*Gravitational_Potential_X(X)
      S(3) = -Q(1)*Gravitational_Potential_Y(X)
      S(4) = -(Q(2)*Gravitational_Potential_X(X)+Q(3)*Gravitational_Potential_Y(X))
#endif
#ifdef EqnShallowWater
      S(1) = 0.
      S(2) = 0. ! -Q(1)*Bathymetry_X(X)
      S(3) = 0. ! -Q(1)*Bathymetry_Y(X)
#endif
#ifdef EqnAcoustics
      S(1) = 0. ! mass source
      S(2) = 0. ! coriolis, friction
      S(3) = 0. ! coriolis, friction
#endif

!-------------------------------------------------------------------------------!
   END FUNCTION SourceFunc
!===============================================================================!
!
!
!
#ifdef EqnEuler
!===============================================================================!
   REAL FUNCTION Gravitational_Potential(X)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nDims
      USE MOD_FiniteVolume2D_vars,ONLY: PI
      USE MOD_FiniteVolume2D_vars,ONLY: source_flag
      IMPLICIT NONE
      REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
      REAL                            :: r2

      SELECT CASE (source_flag)

       CASE DEFAULT
         Gravitational_Potential = 0.
      END SELECT

!-------------------------------------------------------------------------------!
   END FUNCTION Gravitational_Potential
!===============================================================================!
!
!
!
!===============================================================================!
   REAL FUNCTION Gravitational_Potential_X(X)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nDims
      USE MOD_FiniteVolume2D_vars,ONLY: PI
      USE MOD_FiniteVolume2D_vars,ONLY: source_flag
      IMPLICIT NONE
      REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
      REAL                            :: r2

      SELECT CASE (source_flag)

       CASE DEFAULT
         Gravitational_Potential_X = 0.
      END SELECT

!-------------------------------------------------------------------------------!
   END FUNCTION Gravitational_Potential_X
!===============================================================================!
!
!
!
!===============================================================================!
   REAL FUNCTION Gravitational_Potential_Y(X)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nDims
      USE MOD_FiniteVolume2D_vars,ONLY: PI
      USE MOD_FiniteVolume2D_vars,ONLY: source_flag
      IMPLICIT NONE
      REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
      REAL                            :: r2

      SELECT CASE (source_flag)

       CASE DEFAULT
         Gravitational_Potential_Y = 0.
      END SELECT

!-------------------------------------------------------------------------------!
   END FUNCTION Gravitational_Potential_Y
!===============================================================================!
#endif
!
#ifdef EqnShallowWater
! Definition of bathymetry and derivatives
#endif
#ifdef EqnAcoustics
! Definition of coriolis, friction etc
#endif
!
!
!
!===============================================================================!
   SUBROUTINE BoundaryConditions(t)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
      USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
      USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
      USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
      USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
      USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
      USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
      USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
      USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
      USE MOD_FiniteVolume2D_vars,ONLY: U
      USE MOD_FiniteVolume2D_vars,ONLY: V
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)    :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      INTEGER            :: ii, jj
      INTEGER            :: idx_vx, idx_vy
      REAL               :: x0, xc, xt
      REAL               :: Prim_in(1:nVar), Prim_out(1:nVar)
      REAL               :: Cons_in(1:nVar), Cons_out(1:nVar)
      REAL               :: ConsRefState1(1:nVar), ConsRefState2(1:nVar)
      CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

#if defined(EqnEuler) || defined(EqnShallowWater) || defined(EqnAcoustics)
      idx_vx = 2
      idx_vy = 3
#endif

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
      SELECT CASE(BoundaryConditionsType(4))
       CASE(1) ! Periodic
         DO jj=1,nElemsY
            DO ii=0,nGhosts
               U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nElemsX-nGhosts+ii,jj)
            END DO
         END DO
       CASE(2) ! Transmissive
         DO jj=1,nElemsY
            DO ii=0,nGhosts
               U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
            END DO
         END DO
       CASE(3) ! Inflow
         Prim_in(1:nVar) = PrimRefState4(1:nVar)
         CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
         DO jj=1,nElemsY
            DO ii=0,nGhosts
               U(1:nVar,-nGhosts+ii,jj) = Cons_in(1:nVar)
            END DO
         END DO
       CASE(4) ! Outflow
         Prim_out(1:nVar) = PrimRefState4(1:nVar)
         CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
         DO jj=1,nElemsY
            DO ii=0,nGhosts
               U(1:nVar,-nGhosts+ii,jj) = Cons_out(1:nVar)
            END DO
         END DO
       CASE(5) ! Reflecting
         DO jj=1,nElemsY
            DO ii=0,nGhosts
               U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
               U(idx_vx,-nGhosts+ii,jj) =-U(idx_vx,nGhosts-ii+1,jj)
            END DO
         END DO
       CASE DEFAULT
         ErrorMessage = "Boundary condition not implemented"
         WRITE(*,*) ErrorMessage
         STOP
      END SELECT

!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
      SELECT CASE(BoundaryConditionsType(2))
       CASE(1) ! Periodic
         DO jj=1,nElemsY
            DO ii=1,nGhosts+1
               U(1:nVar,nElemsX+ii,jj) = U(1:nVar,ii,jj)
            END DO
         END DO
       CASE(2) ! Transmissive
         DO jj=1,nElemsY
            DO ii=1,nGhosts+1
               U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
            END DO
         END DO
       CASE(3) ! Inflow
         Prim_in(1:nVar) = PrimRefState2(1:nVar)
         CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
         DO jj=1,nElemsY
            DO ii=1,nGhosts+1
               U(1:nVar,nElemsX+ii,jj) = Cons_in(1:nVar)
            END DO
         END DO
       CASE(4) ! Outflow
         Prim_out(1:nVar) = PrimRefState2(1:nVar)
         CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
         DO jj=1,nElemsY
            DO ii=1,nGhosts+1
               U(1:nVar,nElemsX+ii,jj) = Cons_out(1:nVar)
            END DO
         END DO
       CASE(5) ! Reflecting
         DO jj=1,nElemsY
            DO ii=1,nGhosts+1
               U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
               U(idx_vx,nElemsX+ii,jj) =-U(idx_vx,nElemsX-ii+1,jj)
            END DO
         END DO
       CASE DEFAULT
         ErrorMessage = "Boundary condition not implemented"
         WRITE(*,*) ErrorMessage
         STOP
      END SELECT

!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
      SELECT CASE(BoundaryConditionsType(3))
       CASE(1) ! Periodic
         DO ii=1,nElemsX
            DO jj=1,nGhosts+1
               U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,jj)
            END DO
         END DO
       CASE(2) ! Transmissive
         DO ii=1,nElemsX
            DO jj=1,nGhosts+1
               U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
            END DO
         END DO
       CASE(3) ! Inflow
         Prim_in(1:nVar) = PrimRefState3(1:nVar)
         CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
         DO ii=1,nElemsX
            DO jj=1,nGhosts+1
               U(1:nVar,ii,nElemsY+jj) = Cons_in(1:nVar)
            END DO
         END DO
       CASE(4) ! Outflow
         Prim_out(1:nVar) = PrimRefState3(1:nVar)
         CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
         DO ii=1,nElemsX
            DO jj=1,nGhosts+1
               U(1:nVar,ii,nElemsY+jj) = Cons_out(1:nVar)
            END DO
         END DO
       CASE(5) ! Reflecting
         DO ii=1,nElemsX
            DO jj=1,nGhosts+1
               U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
               U(idx_vy,ii,nElemsY+jj) =-U(idx_vy,ii,nElemsY-jj+1)
            END DO
         END DO
       CASE DEFAULT
         ErrorMessage = "Boundary condition not implemented"
         WRITE(*,*) ErrorMessage
         STOP
      END SELECT

!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
      SELECT CASE(BoundaryConditionsType(1))
       CASE(1) ! Periodic
         DO ii=1,nElemsX
            DO jj=0,nGhosts
               U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nElemsY-nGhosts+jj)
            END DO
         END DO
       CASE(2) ! Transmissive
         DO ii=1,nElemsX
            DO jj=0,nGhosts
               U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
            END DO
         END DO
       CASE(3) ! Inflow
         Prim_in(1:nVar) = PrimRefState1(1:nVar)
         CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
         DO ii=1,nElemsX
            DO jj=0,nGhosts
               U(1:nVar,ii,-nGhosts+jj) = Cons_in(1:nVar)
            END DO
         END DO
       CASE(4) ! Outflow
         Prim_out(1:nVar) = PrimRefState1(1:nVar)
         CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
         DO ii=1,nElemsX
            DO jj=0,nGhosts
               U(1:nVar,ii,-nGhosts+jj) = Cons_out(1:nVar)
            END DO
         END DO
       CASE(5) ! Reflecting
         DO ii=1,nElemsX
            DO jj=0,nGhosts
               U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
               U(idx_vy,ii,-nGhosts+jj) =-U(idx_vy,ii,nGhosts-jj+1)
            END DO
         END DO
       CASE DEFAULT
         ErrorMessage = "Boundary condition not implemented"
         WRITE(*,*) ErrorMessage
         STOP
      END SELECT

!------------------------------!
! Upper Corners Boundary Conditions!
!------------------------------!
      SELECT CASE(BoundaryConditionsType(4))
       CASE(1) ! Periodic
         DO jj=-nGhosts,0
            DO ii=0,nGhosts
               U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nElemsX-nGhosts+ii,jj)
            END DO
         END DO

         DO jj=nElemsY+1,nElemsY+1+nGhosts
            DO ii=0,nGhosts
               U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nElemsX-nGhosts+ii,jj)
            END DO
         END DO

!  CASE DEFAULT
!    ErrorMessage = "Boundary condition not implemented"
!    WRITE(*,*) ErrorMessage
!    STOP
      END SELECT


!------------------------------!
! Right Corners Boundary Conditions    !
!------------------------------!
      SELECT CASE(BoundaryConditionsType(2))
       CASE(1) ! Periodic
         DO jj=-nGhosts,0
            DO ii=1,nGhosts+1
               U(1:nVar,nElemsX+ii,jj) = U(1:nVar,ii,jj)
            END DO
         END DO

         DO jj=nElemsY+1,nElemsY+1+nGhosts
            DO ii=1,nGhosts+1
               U(1:nVar,nElemsX+ii,jj) = U(1:nVar,ii,jj)
            END DO
         END DO
!  CASE DEFAULT
!    ErrorMessage = "Boundary condition not implemented"
!    WRITE(*,*) ErrorMessage
!    STOP
      END SELECT

      DO jj=-nGhosts,nElemsY+nGhosts+1
         DO ii=-nGhosts,nElemsX+nGhosts+1
            CALL ConsToPrim(U(1:nVar,ii,jj),V(1:nVar,ii,jj))
         END DO
      END DO

!-------------------------------------------------------------------------------!
   END SUBROUTINE BoundaryConditions
!===============================================================================!
!
!
!
!===============================================================================!
   FUNCTION TimeStep()
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: U
      USE MOD_FiniteVolume2D_vars,ONLY: CFL
      USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
      USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
      USE MOD_FiniteVolume2D_vars,ONLY: LambdaMaxX
      USE MOD_FiniteVolume2D_vars,ONLY: LambdaMaxY
      USE MOD_FiniteVolume2D_vars,ONLY: MIN_TIMESTEP
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL    :: TimeStep
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL    :: FastestWaveX, FastestWaveY
      REAL    :: Prim(1:nVar)
      INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

      LambdaMaxX = 0.0
      LambdaMaxY = 0.0
      TimeStep = HUGE(1.0)

      DO jj=1,nElemsY
         DO ii=1,nElemsX
            CALL ConsToPrim(U(1:nVar,ii,jj),Prim(1:nVar))
            CALL WaveSpeeds2D(Prim(1:nVar),FastestWaveX,FastestWaveY)
            LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveX))
            LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveY))
         END DO
      END DO
      TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)

      TimeStep = CFL*TimeStep

      IF (TimeStep .LT. MIN_TIMESTEP) THEN
         TimeStep = MIN_TIMESTEP
      END IF

!-------------------------------------------------------------------------------!
   END FUNCTION TimeStep
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE WaveSpeeds1D(Prim,slowest,fastest)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
#ifdef EqnShallowWater
      USE MOD_FiniteVolume2D_vars,ONLY: Gravity
#endif
#ifdef EqnEuler
      USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#endif
#ifdef SW
      USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)           :: Prim(1:nVar)
      REAL,INTENT(OUT),OPTIONAL :: slowest
      REAL,INTENT(OUT),OPTIONAL :: fastest
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

#ifdef EqnEuler
!-------------------------------------------------------------------------------!
      REAL                      :: ro, vdir, p
!-------------------------------------------------------------------------------!

      ro = Prim(1)
      vdir = Prim(2)
      p  = Prim(4)
#ifdef SW
      p  = Kappa*ro**Gmm
#endif

      IF(PRESENT(slowest)) THEN
         slowest = ABS(vdir) - SQRT(Gmm*p/ro)
      END IF

      IF(PRESENT(fastest)) THEN
         fastest = ABS(vdir) + SQRT(Gmm*p/ro)
      END IF
#endif

#ifdef EqnShallowWater
!-------------------------------------------------------------------------------!
      REAL                      :: h, vdir
!-------------------------------------------------------------------------------!

      h = Prim(1)
      vdir = Prim(2)

      IF(PRESENT(slowest)) THEN
         slowest = ABS(vdir) - SQRT(Gravity*h)
      END IF

      IF(PRESENT(fastest)) THEN
         fastest = ABS(vdir) + SQRT(Gravity*h)
      END IF
#endif

#ifdef EqnAcoustics

      IF(PRESENT(slowest)) THEN
         slowest=1.
      END IF
      IF(PRESENT(fastest)) THEN
         fastest = 1.
      END IF
#endif


!-------------------------------------------------------------------------------!
   END SUBROUTINE WaveSpeeds1D
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE WaveSpeeds2D(Prim,fastestx,fastesty)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
#ifdef EqnShallowWater
      USE MOD_FiniteVolume2D_vars,ONLY: Gravity
#endif
#ifdef EqnEuler
      USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#endif
#ifdef SW
      USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: Prim(1:nVar)
      REAL,INTENT(OUT) :: fastestx
      REAL,INTENT(OUT) :: fastesty
!-------------------------------------------------------------------------------!


#ifdef EqnEuler

!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL                      :: ro, vx, vy, p
!-------------------------------------------------------------------------------!
      ro = Prim(1)
      vx = Prim(2)
      vy = Prim(3)
      p  = Prim(4)
#ifdef SW
      p  = Kappa*ro**Gmm
#endif

      fastestx = ABS(vx) + SQRT(Gmm*p/ro)
      fastesty = ABS(vy) + SQRT(Gmm*p/ro)
#endif

#ifdef EqnShallowWater
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             ::h, vx, vy
!-------------------------------------------------------------------------------!
      h  = Prim(1)
      vx = Prim(2)
      vy = Prim(3)

      fastestx = ABS(vx) + SQRT(Gravity*h)
      fastesty = ABS(vy) + SQRT(Gravity*h)
#endif

#ifdef EqnAcoustics
      fastestx = 1.
      fastesty = 1.
#endif

!-------------------------------------------------------------------------------!
   END SUBROUTINE WaveSpeeds2D
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE ConsToPrim(Cons, Prim)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: MIN_POSITIVE_VAR, MIN_SPEED

#ifdef EqnEuler
      USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#endif
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: Cons(1:nVar)
      REAL,INTENT(OUT) :: Prim(1:nVar)

#ifdef EqnEuler
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: ro, rovx, rovy, Energy, rot
!-------------------------------------------------------------------------------!

      ro     = Cons(1)
      rovx   = Cons(2)
      rovy   = Cons(3)
      Energy = Cons(4)


#ifdef PATANKAR
      rot = ro + MIN_POSITIVE_VAR/ro
#else

      IF (ro .LT. MIN_POSITIVE_VAR) THEN
         rot = MIN_POSITIVE_VAR
         rovx = 0.
         rovy = 0.
      END IF


#endif


      Prim(1) = ro
      Prim(2) = rovx/rot
      Prim(3) = rovy/rot
      Prim(4) = (Gmm-1.0)*( Energy-0.5*ro*(Prim(2)**2+Prim(3)**2) )

#endif



#ifdef EqnShallowWater
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: h, qx, qy, ht
!-------------------------------------------------------------------------------!

      h     = Cons(1)
      qx    = Cons(2)
      qy    = Cons(3)

#ifdef PATANKAR
      ht = h + MIN_POSITIVE_VAR/h
#else

      IF (h .LT. MIN_POSITIVE_VAR) THEN
         ht = MIN_POSITIVE_VAR
         qx = 0.
         qy = 0.
      END IF

#endif


      Prim(1) = h
      Prim(2) = qx/ht
      Prim(3) = qy/ht

#endif
#ifdef EqnAcoustics
      Prim(1:nVar) = Cons(1:nVar)
#endif

!-------------------------------------------------------------------------------!
   END SUBROUTINE ConsToPrim
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE PrimToCons(Prim, Cons)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: MIN_POSITIVE_VAR, MIN_SPEED

#ifdef EqnEuler
      USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#endif
#ifdef SW
      USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: Prim(1:nVar)
      REAL,INTENT(OUT) :: Cons(1:nVar)

#ifdef EqnEuler
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: ro, vx, vy, p
!-------------------------------------------------------------------------------!

      ro  = Prim(1)
      vx  = Prim(2)
      vy  = Prim(3)
      p   = Prim(4)
#ifdef SW
      p   = Kappa*ro**Gmm
#endif

#ifdef PATANKAR

#else
      IF (ro .LT. MIN_POSITIVE_VAR) THEN
         ro = MIN_POSITIVE_VAR
         vx=0.
         vy=0.
      END IF
#endif

      Cons(1) = ro
      Cons(2) = ro*vx
      Cons(3) = ro*vy
      Cons(4) = p/(Gmm-1.0)+0.5*ro*(vx**2+vy**2)
#endif

#ifdef EqnShallowWater
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: h, vx, vy
!-------------------------------------------------------------------------------!

      h   = Prim(1)
      vx  = Prim(2)
      vy  = Prim(3)

#ifdef PATANKAR

#else
      IF (h .LT. MIN_POSITIVE_VAR) THEN
         h  = MIN_POSITIVE_VAR
         vx=0.
         vy=0.
      END IF
#endif

      Cons(1) = h
      Cons(2) = h*vx
      Cons(3) = h*vy
#endif

#ifdef EqnAcoustics
      Cons(1:nVar)=Prim(1:nVar)
#endif

!-------------------------------------------------------------------------------!
   END SUBROUTINE PrimToCons
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE EvaluateFlux1D(Prim,Flux)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: MIN_POSITIVE_VAR, MIN_SPEED
#ifdef EqnEuler
      USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#endif
#ifdef SW
      USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
#ifdef EqnShallowWater
      USE MOD_FiniteVolume2D_vars,ONLY: Gravity
#endif
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: Prim(1:nVar)
      REAL,INTENT(OUT) :: Flux(1:nVar)

#ifdef EqnEuler
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: ro, vx, vy, p, Energy
!-------------------------------------------------------------------------------!

      ro = Prim(1)
      vx = Prim(2)
      vy = Prim(3)
      p  = Prim(4)
#ifdef SW
      p  = Kappa*ro**Gmm
#endif

#ifdef PATANKAR

#else
      IF (ro .LT. MIN_POSITIVE_VAR) THEN
         ro = MIN_POSITIVE_VAR
         vx=0.
         vy=0.
      END IF
#endif

      Energy = p/(Gmm-1.0)+0.5*ro*(vx**2+vy**2)

      Flux(1) = ro*vx
      Flux(2) = ro*vx**2 + p
      Flux(3) = ro*vx*vy
      Flux(4) = vx*(Energy+p)
!*NB: Reference:
!*Eleuterio F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
!*3.2.4 The Split Three–Dimensional Riemann Problem
#endif

#ifdef EqnShallowWater
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: h, vx, vy
!-------------------------------------------------------------------------------!

      h  = Prim(1)
      vx = Prim(2)
      vy = Prim(3)

#ifdef PATANKAR

#else
      IF (h .LT. MIN_POSITIVE_VAR) THEN
         h  = MIN_POSITIVE_VAR
         vx=0.
         vy=0.
      END IF
#endif

      Flux(1) = h*vx
      Flux(2) = h*vx**2 + 0.5d0*Gravity*h**2
      Flux(3) = h*vx*vy
#endif

#ifdef EqnAcoustics
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: p, vx, vy
!-------------------------------------------------------------------------------!

      p  = Prim(1)
      vx = Prim(2)
      vy = Prim(3)

      Flux(1) = vx
      Flux(2) = p
      Flux(3) = 0.
#endif

!-------------------------------------------------------------------------------!
   END SUBROUTINE EvaluateFlux1D
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE EvaluateFlux2D_X(Prim,Flux)
!-------------------------------------------------------------------------------!
  USE MOD_FiniteVolume2D_vars,ONLY: nVar
  USE MOD_FiniteVolume2D_vars,ONLY: MIN_POSITIVE_VAR, MIN_SPEED
#ifdef EqnEuler
  USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#endif
#ifdef SW
  USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
#ifdef EqnShallowWater
  USE MOD_FiniteVolume2D_vars,ONLY: Gravity
#endif
!-------------------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
  REAL,INTENT(IN)  :: Prim(1:nVar)
  REAL,INTENT(OUT) :: Flux(1:nVar)

#ifdef EqnEuler
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
  REAL             :: ro, vx, vy, p, Energy
!-------------------------------------------------------------------------------!

  ro = Prim(1)
  vx = Prim(2)
  vy = Prim(3)
  p  = Prim(4)
#ifdef SW
  p  = Kappa*ro**Gmm
#endif

#ifdef PATANKAR

#else
  IF (ro .LT. MIN_POSITIVE_VAR) THEN
      ro = MIN_POSITIVE_VAR
      vx=0.
      vy=0.
  END IF
#endif

  Energy = p/(Gmm-1.0)+0.5*ro*(vx**2+vy**2)

  Flux(1) = ro*vx
  Flux(2) = ro*vx**2 + p
  Flux(3) = ro*vx*vy
  Flux(4) = vx*(Energy+p)
!*NB: Reference:
!*Eleuterio F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
!*3.2.4 The Split Three–Dimensional Riemann Problem
#endif

#ifdef EqnShallowWater
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
  REAL             :: h, vx, vy
!-------------------------------------------------------------------------------!

  h  = Prim(1)
  vx = Prim(2)
  vy = Prim(3)

#ifdef PATANKAR

#else
  IF (h .LT. MIN_POSITIVE_VAR) THEN
      h  = MIN_POSITIVE_VAR
      vx=0.
      vy=0.
  END IF
#endif

  Flux(1) = h*vx
  Flux(2) = h*vx**2 + 0.5d0*Gravity*h**2
  Flux(3) = h*vx*vy
#endif

#ifdef EqnAcoustics
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
  REAL             :: p, vx, vy
!-------------------------------------------------------------------------------!

  p  = Prim(1)
  vx = Prim(2)
  vy = Prim(3)

  Flux(1) = vx
  Flux(2) = p
  Flux(3) = 0.
#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE EvaluateFlux2D_X
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE EvaluateFlux2D_Y(Prim,Flux)
  !-------------------------------------------------------------------------------!
  USE MOD_FiniteVolume2D_vars,ONLY: nVar
  USE MOD_FiniteVolume2D_vars,ONLY: MIN_POSITIVE_VAR, MIN_SPEED
#ifdef EqnEuler
  USE MOD_FiniteVolume2D_vars,ONLY: Gmm
#endif
#ifdef SW
  USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
#ifdef EqnShallowWater
  USE MOD_FiniteVolume2D_vars,ONLY: Gravity
#endif
!-------------------------------------------------------------------------------!
  IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
  REAL,INTENT(IN)  :: Prim(1:nVar)
  REAL,INTENT(OUT) :: Flux(1:nVar)

#ifdef EqnEuler
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
  REAL             :: ro, vx, vy, p, Energy
!-------------------------------------------------------------------------------!

  ro = Prim(1)
  vx = Prim(2)
  vy = Prim(3)
  p  = Prim(4)
#ifdef SW
  p  = Kappa*ro**Gmm
#endif

#ifdef PATANKAR

#else
  IF (ro .LT. MIN_POSITIVE_VAR) THEN
      ro = MIN_POSITIVE_VAR
      vx=0.
      vy=0.
  END IF
#endif

  Energy = p/(Gmm-1.0)+0.5*ro*(vx**2+vy**2)

  Flux(1) = ro*vy
  Flux(2) = ro*vx*vy
  Flux(3) = ro*vy**2 + p
  Flux(4) = vy*(Energy+p)
!*NB: Reference:
!*Eleuterio F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
!*3.2.4 The Split Three–Dimensional Riemann Problem
#endif

#ifdef EqnShallowWater
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
  REAL             :: h, vx, vy
!-------------------------------------------------------------------------------!

  h  = Prim(1)
  vx = Prim(2)
  vy = Prim(3)

#ifdef PATANKAR

#else
  IF (h .LT. MIN_POSITIVE_VAR) THEN
      h  = MIN_POSITIVE_VAR
      vx=0.
      vy=0.
  END IF
#endif

  Flux(1) = h*vy
  Flux(2) = h*vx*vy
  Flux(3) = h*vy**2 + 0.5d0*Gravity*h**2
#endif

#ifdef EqnAcoustics
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
  REAL             :: p, vx, vy
!-------------------------------------------------------------------------------!

  p  = Prim(1)
  vx = Prim(2)
  vy = Prim(3)

  Flux(1) = vy
  Flux(2) = 0.
  Flux(3) = p
#endif

!-------------------------------------------------------------------------------!
END SUBROUTINE EvaluateFlux2D_Y
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE RiemannSolver(ConsL,ConsR,NormVect,TangVect,Flux)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: nDims
      USE MOD_FiniteVolume2D_vars,ONLY: Gmm
      USE MOD_FiniteVolume2D_vars,ONLY: WhichRiemannSolver
      USE exact_riemann_mod,      ONLY: exact_riemann
      USE exact_riemann_mod,      ONLY: sample
#ifdef SW
      USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: ConsL(1:nVar,1:nGPs)
      REAL,INTENT(IN)  :: ConsR(1:nVar,1:nGPs)
      REAL,INTENT(IN)  :: NormVect(1:nDims,1:nGPs)
      REAL,INTENT(IN)  :: TangVect(1:nDims,1:nGPs)
      REAL,INTENT(OUT) :: Flux(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: PrimLL(1:nVar,1:nGPs), PrimRR(1:nVar,1:nGPs)
      REAL             :: ConsLL(1:nVar,1:nGPs), ConsRR(1:nVar,1:nGPs)
      INTEGER          :: iGP
!-------------------------------------------------------------------------------!
! >> LOCAL FOR EXACT RIEMANN SOLVER                                             !
!-------------------------------------------------------------------------------!
      REAL, PARAMETER             :: s=0.0
      REAL                        :: al, ar
      REAL                        :: pl, pr
      REAL                        :: rho_star_l,rho_star_r
      REAL                        :: speedl, speedr
      REAL                        :: pm, um
      REAL                        :: u_norm_l, u_norm_r
      REAL                        :: u_tan_l,  u_tan_r
      REAL, DIMENSION(3+nDims)    :: vstar
      REAL, DIMENSION(2+nDims)    :: w
!-------------------------------------------------------------------------------!



      DO iGP=1,nGPs
         ! Rotating the vector quantities       !
         ConsLL(1,iGP) = ConsL(1,iGP)
         ConsLL(2,iGP) = NormVect(1,iGP)*ConsL(2,iGP) + NormVect(2,iGP)*ConsL(3,iGP)
         ConsLL(3,iGP) = TangVect(1,iGP)*ConsL(2,iGP) + TangVect(2,iGP)*ConsL(3,iGP)
         ConsLL(4,iGP) = ConsL(4,iGP)

         ConsRR(1,iGP) = ConsR(1,iGP)
         ConsRR(2,iGP) = NormVect(1,iGP)*ConsR(2,iGP) + NormVect(2,iGP)*ConsR(3,iGP)
         ConsRR(3,iGP) = TangVect(1,iGP)*ConsR(2,iGP) + TangVect(2,iGP)*ConsR(3,iGP)
         ConsRR(4,iGP) = ConsR(4,iGP)

         CALL ConsToPrim(ConsLL(1:nVar,iGP),PrimLL(1:nVar,iGP))
         CALL ConsToPrim(ConsRR(1:nVar,iGP),PrimRR(1:nVar,iGP))

         SELECT CASE(WhichRiemannSolver)
          CASE(1) !*Rusanov
            CALL RiemannSolverByRusanov(&
               ConsLL(1:nVar,iGP),ConsRR(1:nVar,iGP),&
               PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),Flux(1:nVar,iGP))
#ifdef EqnEuler
          CASE(2) !*Exact
            al=SQRT(Gmm*PrimLL(4,iGP)/PrimLL(1,iGP)) !*sound_ro_e_scal(ul(1),ul(2+ndim),eos)
            ar=SQRT(Gmm*PrimRR(4,iGP)/PrimRR(1,iGP)) !*sound_ro_e_scal(ur(1),ur(2+ndim),eos)
            pl=PrimLL(4,iGP) !*pres_ro_e_scal (ul(1),ul(2+ndim),eos)
            pr=PrimRR(4,iGP) !*pres_ro_e_scal (ur(1),ur(2+ndim),eos)
            u_norm_l=PrimLL(2,iGP) !*SUM(ul(2:1+ndim)*n_norm)
            u_norm_r=PrimRR(2,iGP) !*SUM(ur(2:1+ndim)*n_norm)
            u_tan_l=PrimLL(3,iGP)  !*-ul(2)*n_norm(2)+ul(3)*n_norm(1)
            u_tan_r=PrimRR(3,iGP)  !*-ur(2)*n_norm(2)+ur(3)*n_norm(1)


            !*CALL exact_riemann(Gmm,         ul(1),         ur(1), rho_star_l, rho_star_r, u_norm_l, u_norm_r, um,   pl, pr, pm,   al,ar, speedl, speedr)
            CALL exact_riemann(  Gmm, PrimLL(1,iGP), PrimRR(1,iGP), rho_star_l, rho_star_r, u_norm_l, u_norm_r, um,   pl, pr, pm,   al,ar, speedl, speedr)
            !*OUT       !*OUT                           !*OUT         !*OUT        !*OUT   !*OUT

            !*CALL sample(s,  vstar(3+ndim), vstar(2), vstar(1), ul(1),         ur(1),         u_norm_l, u_norm_r, um, pl, pr, pm, al, ar)
            CALL sample(  s, vstar(3+nDims), vstar(2), vstar(1), PrimLL(1,iGP), PrimRR(1,iGP), u_norm_l, u_norm_r, um, pl, pr, pm, al, ar)
            !*OUT          !*OUT     !*OUT

            !*I DO NOT NEED TO PASS TO THE INTERNAL ENERGY
            !*vstar(2+ndim)=e_ro_pres_scal(vstar(1),vstar(3+ndim),eos) ! energie interne



            w(1)=vstar(1)
            IF (um>0.0) THEN
               w(2)=vstar(2)
               w(3)=u_tan_l
            ELSE
               w(2)=vstar(2)
               w(3)=u_tan_r
            ENDIF
            w(4)=vstar(3+nDims)

            !*I DO NOT NEED TO ROTATE
            ! vstar(2)=w(2)*n_norm(1)-w(3)*n_norm(2)
            ! vstar(3)=w(2)*n_norm(2)+w(3)*n_norm(1)
            ! vstar: ici rho, u,v,eint,p

            CALL EvaluateFlux1D(w,Flux(1:nVar,iGP))
#endif
          CASE DEFAULT
            PRINT*, "Riemann Solver not defined"
            PRINT*, "Riemann Solver was", WhichRiemannSolver
            STOP
         END SELECT

         ! Rotating back the momentum components
         Flux(2:3,iGP) = NormVect(1:nDims,iGP)*Flux(2,iGP) &
            + TangVect(1:nDims,iGP)*Flux(3,iGP)
      END DO

!-------------------------------------------------------------------------------!
   END SUBROUTINE RiemannSolver
!===============================================================================!
!
!
#ifdef GFWENO
!
!
!===============================================================================!
SUBROUTINE RiemannSolverCorner(FluxBL,FluxBR,FluxTR,FluxTL,&
  ConsBL, ConsBR, ConsTR,ConsTL,&
  NumFluxBL, NumFluxBR, NumFluxTR,NumFluxTL)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: WhichRiemannSolver
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE exact_riemann_mod,      ONLY: exact_riemann
USE exact_riemann_mod,      ONLY: sample
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: ConsBL(1:nVar),ConsTL(1:nVar)
REAL,INTENT(IN)  :: ConsBR(1:nVar),ConsTR(1:nVar)
REAL,INTENT(IN)  :: FluxBL(1:nVar),FluxTL(1:nVar)
REAL,INTENT(IN)  :: FluxBR(1:nVar),FluxTR(1:nVar)
REAL,INTENT(OUT) :: NumFluxBL(1:nVar),NumFluxTL(1:nVar)
REAL,INTENT(OUT) :: NumFluxBR(1:nVar),NumFluxTR(1:nVar)

!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: Cons_aver(1:nVar)
REAL             :: Jump_Cons(1:nVar)
REAL             :: Jump_Flux(1:nVar)
REAL             :: Central_Flux(1:nVar)
INTEGER          :: iGP
!-------------------------------------------------------------------------------!
! >> LOCAL FOR EXACT RIEMANN SOLVER                                             !
!-------------------------------------------------------------------------------!
REAL, PARAMETER             :: C=0.50
REAL                        :: tau
REAL                        :: al, ar
REAL                        :: pl, pr
REAL                        :: rho_star_l,rho_star_r
REAL                        :: speedl, speedr
REAL                        :: s_max
REAL                        :: pm, um
REAL                        :: u_norm_l, u_norm_r
REAL                        :: u_tan_l,  u_tan_r
REAL, DIMENSION(3+nDims)    :: vstar
REAL, DIMENSION(2+nDims)    :: w
REAL, DIMENSION(nVar,nVar)  :: JX, JY
!-------------------------------------------------------------------------------!

Cons_aver(1:nVar) = 0.25*(ConsBL(1:nVar)+ConsBR(1:nVar)+ConsTR(1:nVar)+ConsTL(1:nVar))
Jump_Cons(1:nVar) = ConsBL(1:nVar)-ConsBR(1:nVar)-ConsTL(1:nVar)+ConsTR(1:nVar)

Central_Flux(1:nVar) = 0.25*(FluxBL(1:nVar)+FluxBR(1:nVar)+FluxTR(1:nVar)+FluxTL(1:nVar))
Jump_Flux(1:nVar) = FluxBL(1:nVar)-FluxBR(1:nVar)+FluxTR(1:nVar)-FluxTL(1:nVar)

SELECT CASE(WhichRiemannSolver)
  CASE(1) !*Rusanov
    s_max = problem.max_eigenvalue(Cons_aver)


    NumFluxBL(1:nVar) = +( Central_Flux(1:nVar)) + MESH_DX(1)* s_max* (ConsBL(1:nVar)-Cons_aver(1:nVar))   
    NumFluxBR(1:nVar) = -( Central_Flux(1:nVar)) + MESH_DX(1)* s_max* (ConsBR(1:nVar)-Cons_aver(1:nVar))  
    NumFluxTL(1:nVar) = -( Central_Flux(1:nVar)) + MESH_DX(1)* s_max* (ConsTL(1:nVar)-Cons_aver(1:nVar))   
    NumFluxTR(1:nVar) = +( Central_Flux(1:nVar)) + MESH_DX(1)* s_max* (ConsTR(1:nVar)-Cons_aver(1:nVar))  

#ifdef EqnShallowWater 
  CASE(3) !* SUPG
    
    JX = problem.JacobianX(q_aver)
    JY = problem.JacobianY(q_aver)
    
    s_max = problem.max_eigenvalue(Cons_aver)

    tau = C * 1./s_max   /4.

    NumFluxBL(1:nVar) = +Central_Flux(1:nVar)  + tau*(-MATMUL(JX,Jump_Flux)&
                                                      -MATMUL(JY,Jump_Flux))
    NumFluxBR(1:nVar) = -Central_Flux(1:nVar)  + tau*( MATMUL(JX,Jump_Flux)&
                                                      -MATMUL(JY,Jump_Flux))
    NumFluxTL(1:nVar) = -Central_Flux(1:nVar)  + tau*(-MATMUL(JX,Jump_Flux)&
                                                      +MATMUL(JY,Jump_Flux))
    NumFluxTR(1:nVar) = +Central_Flux(1:nVar)  + tau*( MATMUL(JX,Jump_Flux)&
                                                      +MATMUL(JY,Jump_Flux))
#endif
#ifdef EqnEuler
  CASE(2) !*Exact
    al=SQRT(Gmm*PrimLL(4,iGP)/PrimLL(1,iGP)) !*sound_ro_e_scal(ul(1),ul(2+ndim),eos)
    ar=SQRT(Gmm*PrimRR(4,iGP)/PrimRR(1,iGP)) !*sound_ro_e_scal(ur(1),ur(2+ndim),eos)
    pl=PrimLL(4,iGP) !*pres_ro_e_scal (ul(1),ul(2+ndim),eos)
    pr=PrimRR(4,iGP) !*pres_ro_e_scal (ur(1),ur(2+ndim),eos)
    u_norm_l=PrimLL(2,iGP) !*SUM(ul(2:1+ndim)*n_norm)
    u_norm_r=PrimRR(2,iGP) !*SUM(ur(2:1+ndim)*n_norm)
    u_tan_l=PrimLL(3,iGP)  !*-ul(2)*n_norm(2)+ul(3)*n_norm(1)
    u_tan_r=PrimRR(3,iGP)  !*-ur(2)*n_norm(2)+ur(3)*n_norm(1)


    !*CALL exact_riemann(Gmm,         ul(1),         ur(1), rho_star_l, rho_star_r, u_norm_l, u_norm_r, um,   pl, pr, pm,   al,ar, speedl, speedr)
    CALL exact_riemann(  Gmm, PrimLL(1,iGP), PrimRR(1,iGP), rho_star_l, rho_star_r, u_norm_l, u_norm_r, um,   pl, pr, pm,   al,ar, speedl, speedr)
    !*OUT       !*OUT                           !*OUT         !*OUT        !*OUT   !*OUT

    !*CALL sample(s,  vstar(3+ndim), vstar(2), vstar(1), ul(1),         ur(1),         u_norm_l, u_norm_r, um, pl, pr, pm, al, ar)
    CALL sample(  s, vstar(3+nDims), vstar(2), vstar(1), PrimLL(1,iGP), PrimRR(1,iGP), u_norm_l, u_norm_r, um, pl, pr, pm, al, ar)
    !*OUT          !*OUT     !*OUT

    !*I DO NOT NEED TO PASS TO THE INTERNAL ENERGY
    !*vstar(2+ndim)=e_ro_pres_scal(vstar(1),vstar(3+ndim),eos) ! energie interne



    w(1)=vstar(1)
    IF (um>0.0) THEN
        w(2)=vstar(2)
        w(3)=u_tan_l
    ELSE
        w(2)=vstar(2)
        w(3)=u_tan_r
    ENDIF
    w(4)=vstar(3+nDims)

    !*I DO NOT NEED TO ROTATE
    ! vstar(2)=w(2)*n_norm(1)-w(3)*n_norm(2)
    ! vstar(3)=w(2)*n_norm(2)+w(3)*n_norm(1)
    ! vstar: ici rho, u,v,eint,p

    CALL EvaluateFlux1D(w,Flux(1:nVar,iGP))
#endif
  CASE DEFAULT
    PRINT*, "Riemann Solver not defined"
    PRINT*, "Riemann Solver was", WhichRiemannSolver
    STOP
  END SELECT

  ! Rotating back the momentum components
  Flux(2:3,iGP) = NormVect(1:nDims,iGP)*Flux(2,iGP) &
    + TangVect(1:nDims,iGP)*Flux(3,iGP)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE RiemannSolverCorner
!===============================================================================!
!
!
#endif
!
!===============================================================================!
   SUBROUTINE RiemannSolverByRusanov(ConsL,ConsR,PrimL,PrimR,Flux)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: ConsL(1:nVar), ConsR(1:nVar)
      REAL,INTENT(IN)  :: PrimL(1:nVar), PrimR(1:nVar)
      REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: FluxL(1:nVar), FluxR(1:nVar)
      REAL             :: LambdaMax, fastestL, fastestR
!-------------------------------------------------------------------------------!

      CALL EvaluateFlux1D(PrimL,FluxL)
      CALL EvaluateFlux1D(PrimR,FluxR)
      CALL WaveSpeeds1D(PrimL,fastest=fastestL)
      CALL WaveSpeeds1D(PrimR,fastest=fastestR)

      LambdaMax = MAX(ABS(fastestL),ABS(fastestR))

      Flux = 0.5*((FluxL + FluxR) - LambdaMax*(ConsR - ConsL))

!-------------------------------------------------------------------------------!
   END SUBROUTINE RiemannSolverByRusanov
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Equation
!-------------------------------------------------------------------------------!
