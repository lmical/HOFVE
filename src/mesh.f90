!===============================================================================!
MODULE MOD_Mesh
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE BuildMesh
  MODULE PROCEDURE BuildMesh
END INTERFACE
INTERFACE GlobalElem
  MODULE PROCEDURE GlobalElem
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: BuildMesh
PUBLIC :: GlobalElem
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
SUBROUTINE BuildMesh()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: MESH_SX
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: NormVectX
USE MOD_FiniteVolume2D_vars,ONLY: NormVectY
USE MOD_FiniteVolume2D_vars,ONLY: TangVectX
USE MOD_FiniteVolume2D_vars,ONLY: TangVectY
USE MOD_FiniteVolume2D_vars,ONLY: MeshNodes
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: MeshGP   
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGP   
USE MOD_FiniteVolume2D_vars,ONLY: WeightsGPBnd  
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! Local Variables
INTEGER :: ii, jj, iGP, jGP
REAL, DIMENSION(nGPs) :: quadWeights1D, quadNodes1D 
!-------------------------------------------------------------------------------!

MeshNodes = 0.0
MeshBary  = 0.0

Mesh_SX    = ABS(Mesh_X1-Mesh_X0)
Mesh_DX(1) = ABS(Mesh_SX(1))/(REAL(nElemsX))
Mesh_DX(2) = ABS(Mesh_SX(2))/(REAL(nElemsY))

DO jj=0,nElemsY
  DO ii=0,nElemsX
    MeshNodes(1:nDims,ii,jj) = Mesh_X0(1:2) + (/REAL(ii),REAL(jj)/)*Mesh_DX(1:2)
  END DO
END DO

DO jj=1,nElemsY
  DO ii=1,nElemsX
    MeshBary(1:nDims,ii,jj) = MeshNodes(1:nDims,ii-1,jj-1) + 0.5*Mesh_DX(1:2)
  END DO
END DO

DO iGP=1,nGPs
  !------------------------------!
  ! Normal vectors: x-direction  !
  !------------------------------!
  DO jj=1,nElemsY
    DO ii=0,nElemsX
      NormVectX(1:nDims,iGP,ii,jj) = (/1.0,0.0/)
    END DO
  END DO

  !------------------------------!
  ! Normal vectors: y-direction  !
  !------------------------------!
  DO jj=0,nElemsY
    DO ii=1,nElemsX
      NormVectY(1:nDims,iGP,ii,jj) = (/0.0,1.0/)
    END DO
  END DO

  !------------------------------!
  ! Tangent vectors: x-direction !
  !------------------------------!
  DO jj=1,nElemsY
    DO ii=0,nElemsX
      TangVectX(1:nDims,iGP,ii,jj) = (/0.0,1.0/)
    END DO
  END DO

  !------------------------------!
  ! Tangent vectors: y-direction !
  !------------------------------!
  DO jj=0,nElemsY
    DO ii=1,nElemsX
      TangVectY(1:nDims,iGP,ii,jj) = (/-1.0,0.0/)
    END DO
  END DO
END DO

  !------------------------------!
  !   MeshGP Quadrature          !
  !------------------------------!

SELECT CASE (nGPs)
  CASE(1)
    quadWeights1D(1) = 1.0
    quadNodes1D(1)  =  0.0 
  CASE(2)
    quadWeights1D = (/0.5,0.5 /)
    quadNodes1D  =  (/- 1./(2.*sqrt(3.)), 1./(2.*sqrt(3.)) /)
  CASE(3) !*NOT OK FOR WENO5 because nonlinear weights arise
    quadWeights1D = (/5./18.,4./9., 5./18. /)
    quadNodes1D  =  (/- 0.5*sqrt(3./5.), 0.,  0.5*sqrt(3./5.) /)
  CASE(4) !*4 points-> order 8
    quadWeights1D = (/  (18.-SQRT(30.))/72., (18.+SQRT(30.))/72., (18.+SQRT(30.))/72. , (18.-SQRT(30.))/72. /)
    quadNodes1D  =  (/  -0.5*SQRT(3./7.+2./7.*SQRT(6./5.)), -0.5*SQRT(3./7.-2./7.*SQRT(6./5.)),&
     0.5*SQRT(3./7.-2./7.*SQRT(6./5.)), 0.5*SQRT(3./7.+2./7.*SQRT(6./5.))  /)
  CASE(5) !*5 points-> order 10  !*NOT OK FOR WENO8 because nonlinear weights arise
    quadWeights1D = (/ 0.568888888888888888888888, 0.478628670499366468041291, 0.478628670499366468041291, 0.236926885056189087514264, 0.236926885056189087514264  /)
    quadWeights1D = 0.5*quadWeights1D
    quadNodes1D   = (/ 0.0, -0.5384693101056830910363, 0.53846931010568309103631, -0.9061798459386639927976, 0.90617984593866399279762  /)
    quadNodes1D   = 0.5*quadNodes1D
  CASE(6) !*6 points-> order 12
    quadWeights1D = (/ 0.360761573048138607569833, 0.360761573048138607569833, 0.467913934572691047389870, 0.467913934572691047389870, 0.171324492379170345040296, 0.171324492379170345040296  /)
    quadWeights1D = 0.5*quadWeights1D
    quadNodes1D   = (/ 0.6612093864662645136613, -0.661209386466264513661, -0.238619186083196908630, 0.2386191860831969086305, -0.932469514203152027812, 0.9324695142031520278123  /)
    quadNodes1D   = 0.5*quadNodes1D
  CASE(7) !*7 points-> order 14
    quadWeights1D = (/ 0.4179591836734693877551020, 0.3818300505051189449503697, 0.3818300505051189449503697, 0.2797053914892766679014677, 0.2797053914892766679014677, 0.1294849661688696932706114, 0.1294849661688696932706114/)
    quadWeights1D = 0.5*quadWeights1D
    quadNodes1D   = (/ 0.0, 0.40584515137739716690660641, -0.40584515137739716690660641, -0.741531185599394439863864, 0.74153118559939443986386477, -0.94910791234275852452618968, 0.94910791234275852452618968  /)
    quadNodes1D   = 0.5*quadNodes1D
  CASE(8) !*8 points-> order 16
    quadWeights1D = (/ 0.3626837833783619829651, 0.3626837833783619829651, 0.3137066458778872873379, 0.3137066458778872873379, 0.2223810344533744705443, 0.2223810344533744705443, 0.1012285362903762591525, 0.1012285362903762591525  /)
    quadWeights1D = 0.5*quadWeights1D
    quadNodes1D   = (/ -0.1834346424956498049394, 0.18343464249564980493947, -0.5255324099163289858177, 0.52553240991632898581773, -0.7966664774136267395915, 0.79666647741362673959155, -0.9602898564975362316835, 0.96028985649753623168356  /)
    quadNodes1D   = 0.5*quadNodes1D
  CASE DEFAULT
    PRINT*, "Quadrature not implemented"
    STOP
END SELECT


DO iGP = 1, nGPs
  WeightsGPBnd(iGP) = quadWeights1D(iGP)
  DO jGP = 1, nGPs
    WeightsGP(iGP,jGP) = quadWeights1D(iGP)* quadWeights1D(jGP)
    DO jj=1,nElemsY
      DO ii=1,nElemsX
        MeshGP(1:nDims,ii,jj,iGP,jGP) = (/ MeshBary(1,ii,jj) +quadNodes1D(iGP)*Mesh_DX(1) , MeshBary(2,ii,jj)  +quadNodes1D(jGP)*Mesh_DX(2) /)
      END DO
    END DO
  END DO
END DO



!-------------------------------------------------------------------------------!
END SUBROUTINE BuildMesh
!===============================================================================!
!
!
!
!===============================================================================!
INTEGER FUNCTION GlobalElem(ii,jj)
!-------------------------------------------------------------------------------!

USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! Local Variables
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!
GlobalElem = MODULO(jj-1,nElemsY)*NelemsX + MODULO(ii-1,NelemsX) + 1
!-------------------------------------------------------------------------------!
END FUNCTION GlobalElem
!===============================================================================!
!
!===============================================================================!
END MODULE MOD_Mesh
!-------------------------------------------------------------------------------!
