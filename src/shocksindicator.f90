!===============================================================================!
MODULE MOD_ShocksIndicator
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ShocksIndicatorX
  MODULE PROCEDURE ShocksIndicatorX
END INTERFACE

INTERFACE ShocksIndicatorY
  MODULE PROCEDURE ShocksIndicatorY
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: ShocksIndicatorX
PUBLIC :: ShocksIndicatorY
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
SUBROUTINE ShocksIndicatorX()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj
LOGICAL :: Shock
!-------------------------------------------------------------------------------!

Shock      = .FALSE.
Ind(1,:,:) = .FALSE.

DO jj=1,nElemsY
  DO ii=1,nElemsX
    Shock = JamesonIndicator(U(1:nVar,ii-1:ii+1,jj))
    Ind(1,ii-1:ii+1,jj) = Shock
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE ShocksIndicatorX
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ShocksIndicatorY()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: V
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj
LOGICAL :: Shock
!-------------------------------------------------------------------------------!

Shock      = .FALSE.
Ind(2,:,:) = .FALSE.

DO jj=1,nElemsY
  DO ii=1,nElemsX
    Shock = JamesonIndicator(U(1:nVar,ii,jj-1:jj+1))
    Ind(2,ii,jj-1:jj+1) = Shock
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE ShocksIndicatorY
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION JamesonIndicator(Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: Cons(1:nVar,-1:1)
LOGICAL         :: JamesonIndicator
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL            :: Indicator
REAL            :: IndicatorLimit
REAL            :: IndVar(-1:1)
!-------------------------------------------------------------------------------!

IndicatorLimit = 0.5E-03

JamesonIndicator = .FALSE.

IndVar    = Cons(1,:)
Indicator = ABS(IndVar(-1) + IndVar(+1) - 2.0*IndVar(0))
Indicator = Indicator/(ABS(IndVar(-1)) + ABS(IndVar(+1)) + ABS(2.0*IndVar(0)))
IF (Indicator .GT. IndicatorLimit) THEN
  JamesonIndicator = .TRUE.
END IF

!-------------------------------------------------------------------------------!
END FUNCTION JamesonIndicator
!===============================================================================!
!
!
!
!-------------------------------------------------------------------------------!
END MODULE MOD_ShocksIndicator
!===============================================================================!
