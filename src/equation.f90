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

!-------------------------------------------------------------------------------!
PUBLIC :: ExactFunction
PUBLIC :: ExactFunctionWB
PUBLIC :: SourceTerms
PUBLIC :: BoundaryConditions
PUBLIC :: TimeStep
PUBLIC :: RiemannSolver
PUBLIC :: EvaluateFlux1D
PUBLIC :: ConsToPrim
PUBLIC :: PrimToCons
PUBLIC :: Gravitational_Potential
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
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
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
#else
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
  !*[892] Smooth periodic IC with the purpose of verifying conservation
  !*------------------------------------------
  CASE(892)
    Prim(1)=8.0+0.5*SIN(2.0*Pi*x(1))
    Prim(2)=0.5+COS(8.0*Pi*x(2))
    Prim(3)=0.5+SIN(4.0*Pi*x(2))
    Prim(4)=7.0+SIN(4.0*Pi*x(1))

    CALL PrimToCons(Prim,Cons)

#endif
  CASE DEFAULT
    ErrorMessage = "Exact function not specified"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
CONTAINS
   
#ifdef SW
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
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY
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
!
!===============================================================================!
SUBROUTINE SourceTerms(t)
!-------------------------------------------------------------------------------!
USE MOD_Reconstruction     ,ONLY: MUSCL
USE MOD_Reconstruction     ,ONLY: WENO3_SecondSweep 
USE MOD_Reconstruction     ,ONLY: WENO5_SecondSweep 
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
USE MOD_FiniteVolume2D_vars,ONLY: GravitationalPotentialFlag
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

IF (GravitationalPotentialFlag .GT. 0) THEN
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
    CASE(4)
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

    CASE DEFAULT
      ErrorMessage = "Reconstruction not implemented"
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

S(1) = 0.
S(2) = -Q(1)*Gravitational_Potential_X(X) 
S(3) = -Q(1)*Gravitational_Potential_Y(X) 
S(4) = -(Q(2)*Gravitational_Potential_X(X)+Q(3)*Gravitational_Potential_Y(X))

!-------------------------------------------------------------------------------!
END FUNCTION SourceFunc
!===============================================================================!
!
!
!
!===============================================================================!
REAL FUNCTION Gravitational_Potential(X)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI   
USE MOD_FiniteVolume2D_vars,ONLY: GravitationalPotentialFlag   
IMPLICIT NONE
REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
REAL                            :: r2

SELECT CASE (GravitationalPotentialFlag)

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
USE MOD_FiniteVolume2D_vars,ONLY: GravitationalPotentialFlag   
IMPLICIT NONE
REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
REAL                            :: r2

SELECT CASE (GravitationalPotentialFlag)

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
USE MOD_FiniteVolume2D_vars,ONLY: GravitationalPotentialFlag   
IMPLICIT NONE
REAL, DIMENSION(1:nDims) , INTENT(IN)  :: X
REAL                            :: r2

SELECT CASE (GravitationalPotentialFlag)

  CASE DEFAULT
   Gravitational_Potential_Y = 0.
END SELECT

!-------------------------------------------------------------------------------!
END FUNCTION Gravitational_Potential_Y
!===============================================================================!
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

idx_vx = 2
idx_vy = 3

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
    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)
  END DO
END DO

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
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
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
REAL                      :: ro, vx, vy, p
!-------------------------------------------------------------------------------!

ro = Prim(1)
vx = Prim(2)
vy = Prim(3)
p  = Prim(4)
#ifdef SW
p  = Kappa*ro**Gmm
#endif

IF(PRESENT(slowest)) THEN
  slowest = ABS(vx) - SQRT(Gmm*p/ro)
END IF

IF(PRESENT(fastest)) THEN
  fastest = ABS(vx) + SQRT(Gmm*p/ro)
END IF

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
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
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
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: ro, vx, vy, p
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
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY, MIN_SPEED
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Cons(1:nVar)
REAL,INTENT(OUT) :: Prim(1:nVar)
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
rot = ro + MIN_DENSITY/ro
#else

IF (ro .LT. MIN_DENSITY) THEN
  ro = MIN_DENSITY
  rovx = 0.
  rovy = 0.
END IF

rot = ro

#endif


Prim(1) = ro
Prim(2) = rovx/rot
Prim(3) = rovy/rot
Prim(4) = (Gmm-1.0)*( Energy-0.5*ro*(Prim(2)**2+Prim(3)**2) )

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
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY, MIN_SPEED
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
IF (ro .LT. MIN_DENSITY) THEN
 ro = MIN_DENSITY
 vx=0.
 vy=0.
END IF
#endif

Cons(1) = ro
Cons(2) = ro*vx
Cons(3) = ro*vy
Cons(4) = p/(Gmm-1.0)+0.5*ro*(vx**2+vy**2)

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
USE MOD_FiniteVolume2D_vars,ONLY: Gmm
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY, MIN_SPEED
#ifdef SW
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
#endif
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
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
IF (ro .LT. MIN_DENSITY) THEN
 ro = MIN_DENSITY
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

!-------------------------------------------------------------------------------!
END SUBROUTINE EvaluateFlux1D
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
