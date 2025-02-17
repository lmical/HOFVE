!===============================================================================!
MODULE MOD_Reconstruction
!===============================================================================!
   IMPLICIT NONE
!-------------------------------------------------------------------------------!
   PRIVATE
!-------------------------------------------------------------------------------!
   INTERFACE ReconstructionX
      MODULE PROCEDURE ReconstructionX
   END INTERFACE

   INTERFACE ReconstructionY
      MODULE PROCEDURE ReconstructionY
   END INTERFACE

   INTERFACE ReconstructionFixX
      MODULE PROCEDURE ReconstructionFixX
   END INTERFACE

   INTERFACE ReconstructionFixY
      MODULE PROCEDURE ReconstructionFixY
   END INTERFACE

#ifdef GFWENO

   INTERFACE ReconstructionXY_Global
      MODULE PROCEDURE ReconstructionXY_Global
   END INTERFACE
   INTERFACE WENO_Global_1D_Eta
      MODULE PROCEDURE WENO_Global_1D_Eta
   END INTERFACE
#endif

   INTERFACE MUSCL
      MODULE PROCEDURE MUSCL
   END INTERFACE

   INTERFACE kMUSCL
      MODULE PROCEDURE kMUSCL
   END INTERFACE

   INTERFACE COMUSCL
      MODULE PROCEDURE COMUSCL
   END INTERFACE

   INTERFACE VLMUSCL
      MODULE PROCEDURE VLMUSCL
   END INTERFACE

   INTERFACE M2MUSCL
      MODULE PROCEDURE M2MUSCL
   END INTERFACE

   INTERFACE VAMUSCL
      MODULE PROCEDURE VAMUSCL
   END INTERFACE

   INTERFACE WENO1_SecondSweep
      MODULE PROCEDURE WENO1_SecondSweep
   END INTERFACE

   INTERFACE WENO3_SecondSweep
      MODULE PROCEDURE WENO3_SecondSweep
   END INTERFACE

   INTERFACE WENO5_SecondSweep
      MODULE PROCEDURE WENO5_SecondSweep
   END INTERFACE

   INTERFACE WENO7_SecondSweep
      MODULE PROCEDURE WENO7_SecondSweep
   END INTERFACE

! INTERFACE WENO9_SecondSweep
   ! MODULE PROCEDURE WENO9_SecondSweep
! END INTERFACE
!
! INTERFACE WENO11_SecondSweep
   ! MODULE PROCEDURE WENO11_SecondSweep
! END INTERFACE
!
! INTERFACE WENO13_SecondSweep
   ! MODULE PROCEDURE WENO13_SecondSweep
! END INTERFACE
!-------------------------------------------------------------------------------!
   PUBLIC :: ReconstructionX
   PUBLIC :: ReconstructionY
   PUBLIC :: ReconstructionFixX
   PUBLIC :: ReconstructionFixY
   PUBLIC :: MUSCL
   PUBLIC :: WENO1_SecondSweep
   PUBLIC :: WENO3_SecondSweep
   PUBLIC :: WENO5_SecondSweep
   PUBLIC :: WENO7_SecondSweep
#ifdef GFWENO
   PUBLIC :: ReconstructionEtaGlobal
   PUBLIC :: ReconstructionXY_Global
   PUBLIC :: WENO_Global_1D_Eta
#endif
! PUBLIC :: WENO9_SecondSweep
! PUBLIC :: WENO11_SecondSweep
! PUBLIC :: WENO13_SecondSweep
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
   SUBROUTINE ReconstructionX()
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: U
      USE MOD_FiniteVolume2D_vars,ONLY: WM
      USE MOD_FiniteVolume2D_vars,ONLY: WP
      USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
      USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: Ind
      USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      INTEGER            :: ii, jj
      CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

      SELECT CASE (Reconstruction)
       CASE(1)
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               WM(1:nVar,nGPs,ii,jj) = U(1:nVar,ii,jj)
               WP(1:nVar,nGPs,ii,jj) = U(1:nVar,ii,jj)
            END DO
         END DO
       CASE(2,20) !*MUSCL
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               CALL MUSCL(&
                  U(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(1))
            END DO
         END DO
       CASE(21) !*k-MUSCL
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               CALL kMUSCL(&
                  U(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(1))
            END DO
         END DO
       CASE(22) !*MUSCL_CO
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               CALL COMUSCL(&
                  U(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(1))
            END DO
         END DO
       CASE(23) !*MUSCL_VL
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               CALL VLMUSCL(&
                  U(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(1))
            END DO
         END DO
       CASE(24) !*MUSCL_M
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               CALL M2MUSCL(&
                  U(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(1))
            END DO
         END DO
       CASE(25) !*MUSCL_VA
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               CALL VAMUSCL(&
                  U(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(1))
            END DO
         END DO
       CASE(3,4,5,7)
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               IF (.NOT. Ind(1,ii,jj)) THEN
                  CALL WENO_XDIR(&
                     U(1:nVar,-nGhosts+ii:ii+nGhosts,-nGhosts+jj:jj+nGhosts),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(1),&
                     Reconstruction)
               END IF
            END DO
         END DO
       CASE DEFAULT
         ErrorMessage = "Reconstruction not implemented"
         WRITE(*,*) ErrorMessage
         STOP
      END SELECT

!-------------------------------------------------------------------------------!
   END SUBROUTINE ReconstructionX
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE ReconstructionY()
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: U
      USE MOD_FiniteVolume2D_vars,ONLY: WM
      USE MOD_FiniteVolume2D_vars,ONLY: WP
      USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
      USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: Ind
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      INTEGER            :: ii, jj
      CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

      SELECT CASE (Reconstruction)
       CASE(1)
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               WM(1:nVar,nGPs,ii,jj) = U(1:nVar,ii,jj)
               WP(1:nVar,nGPs,ii,jj) = U(1:nVar,ii,jj)
            END DO
         END DO
       CASE(2,20) !*MUSCL
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               CALL MUSCL(&
                  U(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))
            END DO
         END DO
       CASE(21) !*k-MUSCL
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               CALL kMUSCL(&
                  U(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))
            END DO
         END DO
       CASE(22) !*MUSCL_CO
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               CALL COMUSCL(&
                  U(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))
            END DO
         END DO
       CASE(23) !*MUSCL_VL
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               CALL VLMUSCL(&
                  U(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))
            END DO
         END DO
       CASE(24) !*MUSCL_M
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               CALL M2MUSCL(&
                  U(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))
            END DO
         END DO
       CASE(25) !*MUSCL_VA
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               CALL VAMUSCL(&
                  U(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                  WM(1:nVar,1:nGPs,ii,jj),&
                  WP(1:nVar,1:nGPs,ii,jj),&
                  MESH_DX(2))
            END DO
         END DO
       CASE(3,4,5,7)
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               IF (.NOT. Ind(2,ii,jj)) THEN
                  CALL WENO_YDIR(&
                     U(1:nVar,-nGhosts+ii:ii+nGhosts,-nGhosts+jj:jj+nGhosts),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(2),&
                     Reconstruction)
               END IF
            END DO
         END DO
       CASE DEFAULT
         ErrorMessage = "Reconstruction not implemented"
         WRITE(*,*) ErrorMessage
         STOP
      END SELECT

!-------------------------------------------------------------------------------!
   END SUBROUTINE ReconstructionY
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE ReconstructionFixX()
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: U
      USE MOD_FiniteVolume2D_vars,ONLY: WM
      USE MOD_FiniteVolume2D_vars,ONLY: WP
      USE MOD_FiniteVolume2D_vars,ONLY: Ind
      USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
      USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
      USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      INTEGER            :: ii, jj, iGP
      CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

      IF (       (Reconstruction .EQ. 1)  .OR. (Reconstruction .EQ. 2)  &
      & .OR. (Reconstruction .EQ. 20) .OR. (Reconstruction .EQ. 21) &
      & .OR. (Reconstruction .EQ. 22) .OR. (Reconstruction .EQ. 23) &
      & .OR. (Reconstruction .EQ. 24) .OR. (Reconstruction .EQ. 25) ) THEN
         RETURN
      END IF


      SELECT CASE (ReconstructionFix)
       CASE(1)
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               IF (Ind(1,ii,jj)) THEN
                  DO iGP=1,nGPs
                     WM(1:nVar,iGP,ii,jj) = U(1:nVar,ii,jj)
                     WP(1:nVar,iGP,ii,jj) = U(1:nVar,ii,jj)
                  END DO
               END IF
            END DO
         END DO
       CASE(2,20)
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               IF (Ind(1,ii,jj)) THEN
                  CALL MUSCL(&
                     U(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(1))
               END IF
            END DO
         END DO
       CASE(21) !*k-MUSCL
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               IF (Ind(1,ii,jj)) THEN
                  CALL kMUSCL(&
                     U(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(1))
               END IF
            END DO
         END DO
       CASE(22) !*MUSCL_CO
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               IF (Ind(1,ii,jj)) THEN
                  CALL COMUSCL(&
                     U(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(1))
               END IF
            END DO
         END DO
       CASE(23) !*MUSCL_VL
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               IF (Ind(1,ii,jj)) THEN
                  CALL VLMUSCL(&
                     U(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(1))
               END IF
            END DO
         END DO
       CASE(24) !*MUSCL_M
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               IF (Ind(1,ii,jj)) THEN
                  CALL M2MUSCL(&
                     U(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(1))
               END IF
            END DO
         END DO
       CASE(25) !*MUSCL_VA
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               IF (Ind(1,ii,jj)) THEN
                  CALL VAMUSCL(&
                     U(1:nVar,-nGhosts+ii:ii+nGhosts,jj),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(1))
               END IF
            END DO
         END DO
       CASE(3,4,5,7)
         DO jj=1,nElemsY
            DO ii=0,nElemsX+1
               IF (Ind(1,ii,jj)) THEN
                  CALL WENO_XDIR(&
                     U(1:nVar,-nGhosts+ii:ii+nGhosts,-nGhosts+jj:jj+nGhosts),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(1),&
                     ReconstructionFix)
               END IF
            END DO
         END DO
       CASE DEFAULT
         ErrorMessage = "ReconstructionFix not implemented"
         WRITE(*,*) ErrorMessage
         STOP
      END SELECT

!-------------------------------------------------------------------------------!
   END SUBROUTINE ReconstructionFixX
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE ReconstructionFixY()
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: U
      USE MOD_FiniteVolume2D_vars,ONLY: WM
      USE MOD_FiniteVolume2D_vars,ONLY: WP
      USE MOD_FiniteVolume2D_vars,ONLY: Ind
      USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
      USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
      USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      INTEGER            :: ii, jj, iGP
      CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

      IF (       (Reconstruction .EQ. 1)  .OR. (Reconstruction .EQ. 2)  &
      & .OR. (Reconstruction .EQ. 20) .OR. (Reconstruction .EQ. 21) &
      & .OR. (Reconstruction .EQ. 22) .OR. (Reconstruction .EQ. 23) &
      & .OR. (Reconstruction .EQ. 24) .OR. (Reconstruction .EQ. 25) ) THEN
         RETURN
      END IF

      SELECT CASE (ReconstructionFix)
       CASE(1)
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               IF (Ind(2,ii,jj)) THEN
                  DO iGP=1,nGPs
                     WM(1:nVar,iGP,ii,jj) = U(1:nVar,ii,jj)
                     WP(1:nVar,iGP,ii,jj) = U(1:nVar,ii,jj)
                  END DO
               END IF
            END DO
         END DO
       CASE(2,20)
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               IF (Ind(2,ii,jj)) THEN
                  CALL MUSCL(&
                     U(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(2))
               END IF
            END DO
         END DO
       CASE(21) !*k-MUSCL
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               IF (Ind(2,ii,jj)) THEN
                  CALL kMUSCL(&
                     U(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(2))
               END IF
            END DO
         END DO
       CASE(22) !*MUSCL_CO
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               IF (Ind(2,ii,jj)) THEN
                  CALL COMUSCL(&
                     U(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(2))
               END IF
            END DO
         END DO
       CASE(23) !*MUSCL_VL
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               IF (Ind(2,ii,jj)) THEN
                  CALL VLMUSCL(&
                     U(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(2))
               END IF
            END DO
         END DO
       CASE(24) !*MUSCL_M
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               IF (Ind(2,ii,jj)) THEN
                  CALL M2MUSCL(&
                     U(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(2))
               END IF
            END DO
         END DO
       CASE(25) !*MUSCL_VA
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               IF (Ind(2,ii,jj)) THEN
                  CALL VAMUSCL(&
                     U(1:nVar,ii,-nGhosts+jj:jj+nGhosts),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(2))
               END IF
            END DO
         END DO
       CASE(3,4,5,7)
         DO jj=0,nElemsY+1
            DO ii=1,nElemsX
               IF (Ind(2,ii,jj)) THEN
                  CALL WENO_YDIR(&
                     U(1:nVar,-nGhosts+ii:ii+nGhosts,-nGhosts+jj:jj+nGhosts),&
                     WM(1:nVar,1:nGPs,ii,jj),&
                     WP(1:nVar,1:nGPs,ii,jj),&
                     MESH_DX(2),&
                     ReconstructionFix)
               END IF
            END DO
         END DO
       CASE DEFAULT
         ErrorMessage = "ReconstructionFix not implemented"
         WRITE(*,*) ErrorMessage
         STOP
      END SELECT

!-------------------------------------------------------------------------------!
   END SUBROUTINE ReconstructionFixY
!===============================================================================!
!
!
!
#ifdef GFWENO
!
!
!
!===============================================================================!
SUBROUTINE ReconstructionXY_Global()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: Eta
USE MOD_FiniteVolume2D_vars,ONLY: Ind
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: FG_reconstructed_corner
USE MOD_FiniteVolume2D_vars,ONLY: Cons_reconstructed_corner
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
REAL               :: tempL(nVar,0:nElemsX+1,-nGhosts:nElemsY+nGhosts+1)
REAL               :: tempR(nVar,0:nElemsX+1,-nGhosts:nElemsY+nGhosts+1)
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

! global flux reconstruction
SELECT CASE (Reconstruction)
  CASE(1)
    ! Reconstruct F+G global flux in corners
    DO ii=0,nElemsX+1 
      DO jj=-nGhosts,nElemsY+nGhosts+1
          DO iVar=1,nVar
            CALL WENO1_FirstSweep(&
                      FG(iVar,ii-nGhosts:ii-nGhosts,jj),tempL(iVar,ii,jj),tempR(iVar,ii,jj))
          END DO                  
      END DO
    END DO

    DO ii=0,nElemsX+1 
      DO jj=0,nElemsY+1
          DO iVar=1,nVar
            CALL WENO1_FirstSweep(&
                      tempL(iVar,ii,jj-nGhosts:jj+nGhosts),&
                      FG_reconstructed_corner(iVar,1,1,ii,jj),FG_reconstructed_corner(iVar,1,2,ii,jj))
            CALL WENO1_FirstSweep(&
                      tempR(iVar,ii,jj-nGhosts:jj+nGhosts),&
                      FG_reconstructed_corner(iVar,2,1,ii,jj),FG_reconstructed_corner(iVar,2,2,ii,jj))

          END DO                  
      END DO
    END DO

    ! Reconstruct conservative variables in corners
    ! I'm not sure if I should reconstruct Eta instead of h...
    DO ii=0,nElemsX+1 
      DO jj=-nGhosts,nElemsY+nGhosts+1
          DO iVar=1,nVar
            CALL WENO1_FirstSweep(&
                      U(iVar,ii-nGhosts:ii-nGhosts,jj),tempL(iVar,ii,jj),tempR(iVar,ii,jj))
          END DO                  
      END DO
    END DO

    DO ii=0,nElemsX+1 
      DO jj=0,nElemsY+1
          DO iVar=1,nVar
            CALL WENO1_FirstSweep(&
                      tempL(iVar,ii,jj-nGhosts:jj+nGhosts),&
                      Cons_reconstructed_corner(iVar,1,1,ii,jj),Cons_reconstructed_corner(iVar,1,2,ii,jj))
            CALL WENO1_FirstSweep(&
                      tempR(iVar,ii,jj-nGhosts:jj+nGhosts),&
                      Cons_reconstructed_corner(iVar,2,1,ii,jj),Cons_reconstructed_corner(iVar,2,2,ii,jj))

          END DO                  
      END DO
    END DO

 CASE(3)
    ! Reconstruct F+G global flux in corners
    DO ii=0,nElemsX+1 
      DO jj=-nGhosts,nElemsY+nGhosts+1
          DO iVar=1,nVar
            CALL WENO3_FirstSweep(&
                      FG(iVar,ii-nGhosts:ii-nGhosts,jj),tempL(iVar,ii,jj),tempR(iVar,ii,jj))
          END DO                  
      END DO
    END DO

    DO ii=0,nElemsX+1 
      DO jj=0,nElemsY+1
          DO iVar=1,nVar
            CALL WENO3_FirstSweep(&
                      tempL(iVar,ii,jj-nGhosts:jj+nGhosts),&
                      FG_reconstructed_corner(iVar,1,1,ii,jj),FG_reconstructed_corner(iVar,1,2,ii,jj))
            CALL WENO3_FirstSweep(&
                      tempR(iVar,ii,jj-nGhosts:jj+nGhosts),&
                      FG_reconstructed_corner(iVar,2,1,ii,jj),FG_reconstructed_corner(iVar,2,2,ii,jj))

          END DO                  
      END DO
    END DO

    ! Reconstruct conservative variables in corners
    ! I'm not sure if I should reconstruct Eta instead of h...
    DO ii=0,nElemsX+1 
      DO jj=-nGhosts,nElemsY+nGhosts+1
          DO iVar=1,nVar
            CALL WENO3_FirstSweep(&
                      U(iVar,ii-nGhosts:ii-nGhosts,jj),tempL(iVar,ii,jj),tempR(iVar,ii,jj))
          END DO                  
      END DO
    END DO

    DO ii=0,nElemsX+1 
      DO jj=0,nElemsY+1
          DO iVar=1,nVar
            CALL WENO3_FirstSweep(&
                      tempL(iVar,ii,jj-nGhosts:jj+nGhosts),&
                      Cons_reconstructed_corner(iVar,1,1,ii,jj),Cons_reconstructed_corner(iVar,1,2,ii,jj))
            CALL WENO3_FirstSweep(&
                      tempR(iVar,ii,jj-nGhosts:jj+nGhosts),&
                      Cons_reconstructed_corner(iVar,2,1,ii,jj),Cons_reconstructed_corner(iVar,2,2,ii,jj))

          END DO                  
      END DO
    END DO

 CASE(4)
    ! Reconstruct F+G global flux in corners
    DO ii=0,nElemsX+1 
      DO jj=-nGhosts,nElemsY+nGhosts+1
          DO iVar=1,nVar
            CALL WENO5_FirstSweep(&
                      FG(iVar,ii-nGhosts:ii-nGhosts,jj),tempL(iVar,ii,jj),tempR(iVar,ii,jj))
          END DO                  
      END DO
    END DO

    DO ii=0,nElemsX+1 
      DO jj=0,nElemsY+1
          DO iVar=1,nVar
            CALL WENO5_FirstSweep(&
                      tempL(iVar,ii,jj-nGhosts:jj+nGhosts),&
                      FG_reconstructed_corner(iVar,1,1,ii,jj),FG_reconstructed_corner(iVar,1,2,ii,jj))
            CALL WENO5_FirstSweep(&
                      tempR(iVar,ii,jj-nGhosts:jj+nGhosts),&
                      FG_reconstructed_corner(iVar,2,1,ii,jj),FG_reconstructed_corner(iVar,2,2,ii,jj))

          END DO                  
      END DO
    END DO

    ! Reconstruct conservative variables in corners
    ! I'm not sure if I should reconstruct Eta instead of h...
    DO ii=0,nElemsX+1 
      DO jj=-nGhosts,nElemsY+nGhosts+1
          DO iVar=1,nVar
            CALL WENO5_FirstSweep(&
                      U(iVar,ii-nGhosts:ii-nGhosts,jj),tempL(iVar,ii,jj),tempR(iVar,ii,jj))
          END DO                  
      END DO
    END DO

    DO ii=0,nElemsX+1 
      DO jj=0,nElemsY+1
          DO iVar=1,nVar
            CALL WENO5_FirstSweep(&
                      tempL(iVar,ii,jj-nGhosts:jj+nGhosts),&
                      Cons_reconstructed_corner(iVar,1,1,ii,jj),Cons_reconstructed_corner(iVar,1,2,ii,jj))
            CALL WENO5_FirstSweep(&
                      tempR(iVar,ii,jj-nGhosts:jj+nGhosts),&
                      Cons_reconstructed_corner(iVar,2,1,ii,jj),Cons_reconstructed_corner(iVar,2,2,ii,jj))

          END DO                  
      END DO
    END DO



  CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE ReconstructionXY_Global
!===============================================================================!
!
!
#endif
!===============================================================================!
   SUBROUTINE MUSCL(Q,WM,WP,dx)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: dx
      REAL,INTENT(IN)  :: Q(1:nVar,-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: WM(1:nVar,1:nGPs)
      REAL,INTENT(OUT) :: WP(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: sm, sp, slope
      INTEGER          :: iVar, iGP
!-------------------------------------------------------------------------------!

      DO iGP=1,nGPs
         DO iVar = 1, nVar
            sp    = (Q(iVar,+1) - Q(iVar,+0))/dx
            sm    = (Q(iVar,+0) - Q(iVar,-1))/dx
            slope = MINMOD(sm,sp)

            WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
            WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
         END DO
      END DO

!-------------------------------------------------------------------------------!
   END SUBROUTINE MUSCL
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE kMUSCL(Q,WM,WP,dx)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts

!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: dx
      REAL,INTENT(IN)  :: Q(1:nVar,-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: WM(1:nVar,1:nGPs)
      REAL,INTENT(OUT) :: WP(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: slope, a, b, s1, s2
      REAL, PARAMETER  :: k=1.5 !* 1<=k<=2
      INTEGER          :: iVar, iGP
!-------------------------------------------------------------------------------!

      DO iGP=1,nGPs
         DO iVar = 1, nVar
            a    = Q(iVar,+1) - Q(iVar,+0)
            b    = Q(iVar,+0) - Q(iVar,-1)

            s1=ABS( MINMOD( k*a , b   ) )
            s2=ABS( MINMOD( a   , k*b ) )

            !*SIGN(A,B) returns the value of A with the sign of B.
            !*NB: SIGN should work for reals
            slope = SIGN(1.,a)/dx*MAX( s1, s2 )

            WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
            WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
         END DO
      END DO

!-------------------------------------------------------------------------------!
   END SUBROUTINE kMUSCL
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE COMUSCL(Q,WM,WP,dx)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: dx
      REAL,INTENT(IN)  :: Q(1:nVar,-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: WM(1:nVar,1:nGPs)
      REAL,INTENT(OUT) :: WP(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: slope, a, b, s1
      REAL, PARAMETER  :: k=1. !* 1<=k<=2 !*NB: In my opinion it only works with k=1
      INTEGER          :: iVar, iGP
!-------------------------------------------------------------------------------!

      DO iGP=1,nGPs
         DO iVar = 1, nVar
            a    = Q(iVar,+1) - Q(iVar,+0)
            b    = Q(iVar,+0) - Q(iVar,-1)

            s1=MINMOD( k*a , b )
            slope = s1/dx

            WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
            WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
         END DO
      END DO

!-------------------------------------------------------------------------------!
   END SUBROUTINE COMUSCL
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE VLMUSCL(Q,WM,WP,dx)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: dx
      REAL,INTENT(IN)  :: Q(1:nVar,-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: WM(1:nVar,1:nGPs)
      REAL,INTENT(OUT) :: WP(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: slope, a, b
      INTEGER          :: iVar, iGP
!-------------------------------------------------------------------------------!

      DO iGP=1,nGPs
         DO iVar = 1, nVar
            a    = Q(iVar,+1) - Q(iVar,+0)
            b    = Q(iVar,+0) - Q(iVar,-1)

            slope = 0.

            IF ( (a*b) .GT. 0. ) THEN
               slope = 2.*a*b/(a+b)/dx
            ELSE
               slope = 0.
            END IF

            WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
            WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
         END DO
      END DO

!-------------------------------------------------------------------------------!
   END SUBROUTINE VLMUSCL
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE M2MUSCL(Q,WM,WP,dx)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: dx
      REAL,INTENT(IN)  :: Q(1:nVar,-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: WM(1:nVar,1:nGPs)
      REAL,INTENT(OUT) :: WP(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: slope, a, b, s1, s2
      INTEGER          :: iVar, iGP
!-------------------------------------------------------------------------------!

      DO iGP=1,nGPs
         DO iVar = 1, nVar
            a    = Q(iVar,+1) - Q(iVar,+0)
            b    = Q(iVar,+0) - Q(iVar,-1)

            s1=(a+b)/2.
            s2=2.*MINMOD(a,b)

            slope = MINMOD(s1,s2)/dx

            WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
            WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
         END DO
      END DO

!-------------------------------------------------------------------------------!
   END SUBROUTINE M2MUSCL
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE VAMUSCL(Q,WM,WP,dx)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: dx
      REAL,INTENT(IN)  :: Q(1:nVar,-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: WM(1:nVar,1:nGPs)
      REAL,INTENT(OUT) :: WP(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: slope, a, b, c
      INTEGER          :: iVar, iGP
!-------------------------------------------------------------------------------!

      DO iGP=1,nGPs
         DO iVar = 1, nVar
            a    = Q(iVar,+1) - Q(iVar,+0)
            b    = Q(iVar,+0) - Q(iVar,-1)

            PRINT*, "MUSCL_VA coded but I do not really know what the parameter c is"
            PRINT*, "Check before using it!"
            STOP

            slope = ( a*b + c**2 ) * (a+b) / (a**2 + b**2 + 2.*c**2 ) / dx

            WM(iVar,iGP) = Q(iVar,0) - 0.5*slope*dx
            WP(iVar,iGP) = Q(iVar,0) + 0.5*slope*dx
         END DO
      END DO

!-------------------------------------------------------------------------------!
   END SUBROUTINE VAMUSCL
!===============================================================================!
!
!
!
!===============================================================================!
   FUNCTION MINMOD(x,y)
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN) :: x, y
      REAL            :: MINMOD
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

      MINMOD = 0.5*(SIGN(1.0,x) + SIGN(1.0,y))*MIN(ABS(x),ABS(y))

!-------------------------------------------------------------------------------!
   END FUNCTION MINMOD
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE WENO_XDIR(V,WM,WP,dx,WhichReconstruction)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)    :: dx
      REAL,INTENT(IN)    :: V(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
      REAL,INTENT(OUT)   :: WM(1:nVar,1:nGPs)
      REAL,INTENT(OUT)   :: WP(1:nVar,1:nGPs)
      INTEGER,INTENT(IN) :: WhichReconstruction
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL               :: VtempM(1:nVar,-nGhosts:nGhosts)
      REAL               :: VtempP(1:nVar,-nGhosts:nGhosts)
      INTEGER            :: iVar, ii, jj, iGP
      CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

      SELECT CASE(WhichReconstruction)
       CASE(3)
         DO iVar=1,nVar
            DO jj=-nGhosts,nGhosts
               CALL WENO3_FirstSweep(&
                  V(iVar,-nGhosts:nGhosts,jj),VtempM(iVar,jj),VtempP(iVar,jj))
            END DO
            CALL WENO3_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
            CALL WENO3_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
         END DO
       CASE(4,5)
         DO iVar=1,nVar
            DO jj=-nGhosts,nGhosts
               CALL WENO5_FirstSweep(&
                  V(iVar,-nGhosts:nGhosts,jj),VtempM(iVar,jj),VtempP(iVar,jj))
            END DO
            CALL WENO5_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
            CALL WENO5_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
         END DO
       CASE(7)
         DO iVar=1,nVar
            DO jj=-nGhosts,nGhosts
               CALL WENO7_FirstSweep(&
                  V(iVar,-nGhosts:nGhosts,jj),VtempM(iVar,jj),VtempP(iVar,jj))
            END DO
            CALL WENO7_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
            CALL WENO7_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
         END DO
       CASE DEFAULT
         ErrorMessage = "Reconstruction not implemented"
         WRITE(*,*) ErrorMessage
         STOP
      END SELECT

!-------------------------------------------------------------------------------!
   END SUBROUTINE WENO_XDIR
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE WENO_YDIR(V,WM,WP,dy,WhichReconstruction)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nVar
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)    :: dy
      REAL,INTENT(IN)    :: V(1:nVar,-nGhosts:nGhosts,-nGhosts:nGhosts)
      REAL,INTENT(OUT)   :: WM(1:nVar,1:nGPs)
      REAL,INTENT(OUT)   :: WP(1:nVar,1:nGPs)
      INTEGER,INTENT(IN) :: WhichReconstruction
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL               :: VtempM(1:nVar,-nGhosts:nGhosts)
      REAL               :: VtempP(1:nVar,-nGhosts:nGhosts)
      INTEGER            :: iVar, ii, jj, iGP
      CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

      SELECT CASE(WhichReconstruction)
       CASE(3)
         DO iVar=1,nVar
            DO ii=-nGhosts,nGhosts
               CALL WENO3_FirstSweep(&
                  V(iVar,ii,-nGhosts:nGhosts),VtempM(iVar,ii),VtempP(iVar,ii))
            END DO
            CALL WENO3_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
            CALL WENO3_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
         END DO
       CASE(4,5)
         DO iVar=1,nVar
            DO ii=-nGhosts,nGhosts
               CALL WENO5_FirstSweep(&
                  V(iVar,ii,-nGhosts:nGhosts),VtempM(iVar,ii),VtempP(iVar,ii))
            END DO
            CALL WENO5_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
            CALL WENO5_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
         END DO
       CASE(7)
         DO iVar=1,nVar
            DO ii=-nGhosts,nGhosts
               CALL WENO7_FirstSweep(&
                  V(iVar,ii,-nGhosts:nGhosts),VtempM(iVar,ii),VtempP(iVar,ii))
            END DO
            CALL WENO7_SecondSweep(VtempM(iVar,-nGhosts:nGhosts),WM(iVar,1:nGPs))
            CALL WENO7_SecondSweep(VtempP(iVar,-nGhosts:nGhosts),WP(iVar,1:nGPs))
         END DO
       CASE DEFAULT
         ErrorMessage = "Reconstruction not implemented"
         WRITE(*,*) ErrorMessage
         STOP
      END SELECT

!-------------------------------------------------------------------------------!
   END SUBROUTINE WENO_YDIR
!===============================================================================!
!
!
   !
!
!===============================================================================!
SUBROUTINE WENO1_FirstSweep(Q,WM,WP)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: WM
REAL,INTENT(OUT) :: WP
WM = Q(0)
WP = Q(0)
!-------------------------------------------------------------------------------!
END SUBROUTINE WENO1_FirstSweep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WENO1_SecondSweep(Q,QCoeff,W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(IN)  :: QCoeff(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: W(1:nGPs)
W(1) = Q(0)
!-------------------------------------------------------------------------------!
END SUBROUTINE WENO1_SecondSweep
!===============================================================================!
!
!
!===============================================================================!
   SUBROUTINE WENO3_FirstSweep(Q,WM,WP)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
      USE MOD_FiniteVolume2D_vars,ONLY: LinearWeightsOnly
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: WM
      REAL,INTENT(OUT) :: WP
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: alpha1, alpha2
      REAL             :: beta1, beta2
      REAL             :: gamma1, gamma2
      REAL             :: omega1, omega2
      REAL             :: W1, W2
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
      beta1 = (Q(-1) - Q(+0))**2.0
      beta2 = (Q(+0) - Q(+1))**2.0


!------------------------------!
! WM: x_{i-1/2}                !
!------------------------------!

! Linear Weights
      gamma1 = 2.0/3.0
      gamma2 = 1.0/3.0

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2)
      omega2 = alpha2/(alpha1 + alpha2)

      W1 = 0.5*(    Q(-1) + Q(+0))
      W2 = 0.5*(3.0*Q(+0) - Q(+1))

      IF (LinearWeightsOnly) THEN
         WM = gamma1*W1 + gamma2*W2
      ELSE
         WM = omega1*W1 + omega2*W2
      END IF

!------------------------------!
! WP: x_{i+1/2}                !
!------------------------------!

! Linear Weights
      gamma1 = 1.0/3.0
      gamma2 = 2.0/3.0

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2)
      omega2 = alpha2/(alpha1 + alpha2)

! Reconstructed Polynomial
      W1 = 0.5*(-Q(-1) + 3.0*Q(+0))
      W2 = 0.5*( Q(+0) +     Q(+1))

      IF (LinearWeightsOnly) THEN
         WP = gamma1*W1 + gamma2*W2
      ELSE
         WP = omega1*W1 + omega2*W2
      END IF

!-------------------------------------------------------------------------------!
   END SUBROUTINE WENO3_FirstSweep
!===============================================================================!
!
!
!===============================================================================!
SUBROUTINE WENO3_SecondSweep(Q,QCoeff,W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(IN)  :: QCoeff(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: W(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2
REAL             :: beta1, beta2
REAL             :: gamma1, gamma2
REAL             :: omega1, omega2
REAL             :: W1, W2
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (QCoeff(-1) - QCoeff(+0))**2.0
beta2 = (QCoeff(+0) - QCoeff(+1))**2.0


!------------------------------!
! Point: x_{j-1/(2*sqrt(3))}   !
!------------------------------!

! Linear Weights
gamma1 = 0.5
gamma2 = 0.5

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

! Reconstructed Polynomial
W1   = (1.0/6.0)*(SQRT(3.0)*Q(-1) + 6.0*Q(+0) - SQRT(3.0)*Q(+0))
W2   = (1.0/6.0)*(SQRT(3.0)*Q(+0) + 6.0*Q(+0) - SQRT(3.0)*Q(+1))
W(1) = omega1*W1 + omega2*W2


!------------------------------!
! Point: x_{j+1/(2*sqrt(3))}   !
!------------------------------!

! Linear Weights
gamma1 = 0.5
gamma2 = 0.5

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

! Reconstructed Polynomial
W1   = (1.0/6.0)*(-SQRT(3.0)*Q(-1) + 6.0*Q(+0) + SQRT(3.0)*Q(+0))
W2   = (1.0/6.0)*(-SQRT(3.0)*Q(+0) + 6.0*Q(+0) + SQRT(3.0)*Q(+1))
W(2) = omega1*W1 + omega2*W2

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO3_SecondSweep
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE WENO5_FirstSweep(Q,WM,WP)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
      USE MOD_FiniteVolume2D_vars,ONLY: LinearWeightsOnly
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: WM
      REAL,INTENT(OUT) :: WP
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: alpha1, alpha2, alpha3
      REAL             :: beta1,  beta2,  beta3
      REAL             :: gamma1, gamma2, gamma3
      REAL             :: omega1, omega2, omega3
      REAL             :: W1, W2, W3
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
      beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
         + 11.0*Q(-2)*Q(+0) - 31.0*Q(-1)*Q(+0) + 10.0*Q(+0)*Q(+0))
      beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(+0) + 13.0*Q(+0)*Q(+0) &
         +  5.0*Q(-1)*Q(+1) - 13.0*Q(+0)*Q(+1) +  4.0*Q(+1)*Q(+1))
      beta3 = (1.0/3.0)*(10.0*Q(+0)*Q(+0) - 31.0*Q(+0)*Q(+1) + 25.0*Q(+1)*Q(+1) &
         + 11.0*Q(+0)*Q(+2) - 19.0*Q(+1)*Q(+2) +  4.0*Q(+2)*Q(+2))


!------------------------------!
! WM: x_{i-1/2}                !
!------------------------------!

! Linear Weights
      gamma1 = 3.0/10.0
      gamma2 = 3.0/5.0
      gamma3 = 1.0/10.0

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3)

      W1 = (1.0/6.0)*(    -Q(-2) + 5.0*Q(-1) + 2.0*Q(+0))
      W2 = (1.0/6.0)*( 2.0*Q(-1) + 5.0*Q(+0) -     Q(+1))
      W3 = (1.0/6.0)*(11.0*Q(+0) - 7.0*Q(+1) + 2.0*Q(+2))

      IF (LinearWeightsOnly) THEN
         WM = gamma1*W1 + gamma2*W2 + gamma3*W3
      ELSE
         WM = omega1*W1 + omega2*W2 + omega3*W3
      END IF

!------------------------------!
! WP: x_{i+1/2}                !
!------------------------------!

! Linear Weights
      gamma1 = 1.0/10.0
      gamma2 = 3.0/5.0
      gamma3 = 3.0/10.0

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3)

      W1 = (1.0/6.0)*(2.0*Q(-2) - 7.0*Q(-1) + 11.0*Q(+0))
      W2 = (1.0/6.0)*(   -Q(-1) + 5.0*Q(+0) +  2.0*Q(+1))
      W3 = (1.0/6.0)*(2.0*Q(+0) + 5.0*Q(+1) -      Q(+2))

      IF (LinearWeightsOnly) THEN
         WP = gamma1*W1 + gamma2*W2 + gamma3*W3
      ELSE
         WP = omega1*W1 + omega2*W2 + omega3*W3
      END IF
!-------------------------------------------------------------------------------!
   END SUBROUTINE WENO5_FirstSweep
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE WENO5_SecondSweep2nGPs(Q,W)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
      USE MOD_FiniteVolume2D_vars,ONLY: LinearWeightsOnly
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: W(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: alpha1, alpha2, alpha3
      REAL             :: beta1,  beta2,  beta3
      REAL             :: gamma1, gamma2, gamma3
      REAL             :: omega1, omega2, omega3
      REAL             :: W1, W2, W3
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
      beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
         + 11.0*Q(-2)*Q(+0) - 31.0*Q(-1)*Q(+0) + 10.0*Q(+0)*Q(+0))
      beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(+0) + 13.0*Q(+0)*Q(+0) &
         +  5.0*Q(-1)*Q(+1) - 13.0*Q(+0)*Q(+1) +  4.0*Q(+1)*Q(+1))
      beta3 = (1.0/3.0)*(10.0*Q(+0)*Q(+0) - 31.0*Q(+0)*Q(+1) + 25.0*Q(+1)*Q(+1) &
         + 11.0*Q(+0)*Q(+2) - 19.0*Q(+1)*Q(+2) +  4.0*Q(+2)*Q(+2))


!------------------------------!
! Point: x_{j-1/(2*sqrt(3))}   !
!------------------------------!

! Linear Weights
      gamma1 = (210.0 + SQRT(3.0))/1080.0
      gamma2 = 11.0/18.0
      gamma3 = (210.0 - SQRT(3.0))/1080.0

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
      W1     = -     SQRT(3.0)*Q(-2) &
         + 4.0*SQRT(3.0)*Q(-1) &
         +          12.0*Q(+0) &
         - 3.0*SQRT(3.0)*Q(+0)
      W2     = + 1.0*SQRT(3.0)*Q(-1) &
         +          12.0*Q(+0) &
         - 1.0*SQRT(3.0)*Q(+1)
      W3     = +          12.0*Q(+0) &
         + 3.0*SQRT(3.0)*Q(+0) &
         - 4.0*SQRT(3.0)*Q(+1) &
         +     SQRT(3.0)*Q(+2)

      IF (LinearWeightsOnly) THEN
         W(1)   = (gamma1*W1 + gamma2*W2 + gamma3*W3)/12.0
      ELSE
         W(1)   = (omega1*W1 + omega2*W2 + omega3*W3)/12.0
      END IF
!------------------------------!
! Point: x_{j+1/(2*sqrt(3))}   !
!------------------------------!

! Linear Weights
      gamma1 = (210.0 - SQRT(3.0))/1080.0
      gamma2 = 11.0/18.0
      gamma3 = (210.0 + SQRT(3.0))/1080.0

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
      W1     = + 1.0*SQRT(3.0)*Q(-2) &
         - 4.0*SQRT(3.0)*Q(-1) &
         +          12.0*Q(+0) &
         + 3.0*SQRT(3.0)*Q(+0)
      W2     = - 1.0*SQRT(3.0)*Q(-1) &
         +          12.0*Q(+0) &
         +     SQRT(3.0)*Q(+1)
      W3     = +          12.0*Q(+0) &
         - 3.0*SQRT(3.0)*Q(+0) &
         + 4.0*SQRT(3.0)*Q(+1) &
         -     SQRT(3.0)*Q(+2)

      IF (LinearWeightsOnly) THEN
         W(2)   = (gamma1*W1 + gamma2*W2 + gamma3*W3)/12.0
      ELSE
         W(2)   = (omega1*W1 + omega2*W2 + omega3*W3)/12.0
      END IF
!-------------------------------------------------------------------------------!
   END SUBROUTINE WENO5_SecondSweep2nGPs
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE WENO5_SecondSweep3nGPs(Q,W)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
      USE MOD_FiniteVolume2D_vars,ONLY: LinearWeightsOnly
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: W(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: alpha1, alpha2, alpha3
      REAL             :: beta1,  beta2,  beta3
      REAL             :: gamma1, gamma2, gamma3
      REAL             :: omega1, omega2, omega3
      REAL             :: W1, W2, W3
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
      beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
         + 11.0*Q(-2)*Q(+0) - 31.0*Q(-1)*Q(+0) + 10.0*Q(+0)*Q(+0))
      beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(+0) + 13.0*Q(+0)*Q(+0) &
         +  5.0*Q(-1)*Q(+1) - 13.0*Q(+0)*Q(+1) +  4.0*Q(+1)*Q(+1))
      beta3 = (1.0/3.0)*(10.0*Q(+0)*Q(+0) - 31.0*Q(+0)*Q(+1) + 25.0*Q(+1)*Q(+1) &
         + 11.0*Q(+0)*Q(+2) - 19.0*Q(+1)*Q(+2) +  4.0*Q(+2)*Q(+2))

!------------------------------!
! Point: x_{j-1/2*sqrt(3/5)}   !
!------------------------------!

! Linear Weights
      gamma1 = (71.0*SQRT(15.0))/5240.0 + 126.0/655.0
      gamma2 = 403.0/655.0
      gamma3 = -(71.0*SQRT(15.0))/5240.0 + 126.0/655.0

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
      W1     = (1.0/30.0 - SQRT(15.)/20.0)*Q(-2) &
         + (SQRT(15.)/5. - 1./15.)*Q(-1) &
         + (31./30. - (3.*SQRT(15.))/20.)*Q(+0)
      W2     = + (SQRT(15.)/20. + 1./30.)*Q(-1) &
         +         14./15.*Q(+0) &
         + (1./30. - SQRT(15.)/20.)*Q(+1)
      W3     = + (3.*SQRT(15.)/20. + 31./30.)*Q(+0) &
         - (SQRT(15.)/5. + 1./15.)*Q(+1) &
         + (SQRT(15.)/20. + 1./30.)*Q(+2)

      IF (LinearWeightsOnly) THEN
         W(1)   = (gamma1*W1 + gamma2*W2 + gamma3*W3)
      ELSE
         W(1)   = (omega1*W1 + omega2*W2 + omega3*W3)
      END IF

!------------------------------!
! Point: x_{j}                 !
!------------------------------!

! Linear Weights
      gamma1 = -9./80.
      gamma2 = 49./40.
      gamma3 = -9./80.

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
      W1     = -1./24. *Q(-2) &
         +1./12. *Q(-1) &
         +23./24. *Q(+0)
      W2     = - 1./24. *Q(-1) &
         + 13./12.*Q(+0) &
         - 1./24. *Q(+1)
      W3     = + 23./24.*Q(+0) &
         + 1. /12.*Q(+1) &
         - 1./24. *Q(+2)

      IF (LinearWeightsOnly) THEN
         W(2)   = (gamma1*W1 + gamma2*W2 + gamma3*W3)
      ELSE
         W(2)   = (omega1*W1 + omega2*W2 + omega3*W3)
      END IF


!------------------------------!
! Point: x_{j+1/2*sqrt(3/5)}   !
!------------------------------!

! Linear Weights
      gamma1 = -(71.0*SQRT(15.0))/5240.0 + 126.0/655.0
      gamma2 = 403.0/655.0
      gamma3 = +(71.0*SQRT(15.0))/5240.0 + 126.0/655.0

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
      W1     = (1.0/30.0 + SQRT(15.)/20.0)*Q(-2) &
         - (SQRT(15.)/5. + 1./15.)*Q(-1) &
         + (31./30. + (3.*SQRT(15.))/20.)*Q(+0)
      W2     = + (-SQRT(15.)/20. + 1./30.)*Q(-1) &
         +         14./15.*Q(+0) &
         + (1./30. + SQRT(15.)/20.)*Q(+1)
      W3     = + (-3.*SQRT(15.)/20. + 31./30.)*Q(+0) &
         + (SQRT(15.)/5. - 1./15.)*Q(+1) &
         + (-SQRT(15.)/20. + 1./30.)*Q(+2)

      IF (LinearWeightsOnly) THEN
         W(3)   = (gamma1*W1 + gamma2*W2 + gamma3*W3)
      ELSE
         W(3)   = (omega1*W1 + omega2*W2 + omega3*W3)
      END IF

!-------------------------------------------------------------------------------!
   END SUBROUTINE WENO5_SecondSweep3nGPs
!===============================================================================!
!
!
!

!===============================================================================!
SUBROUTINE WENO5_SecondSweep(Q,QCoeff,W) !4nGPS
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
REAL,INTENT(IN)  :: QCoeff(-nGhosts:nGhosts)
REAL,INTENT(OUT) :: W(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2, alpha3
REAL             :: beta1,  beta2,  beta3
REAL             :: gamma1, gamma2, gamma3
REAL             :: omega1, omega2, omega3
REAL             :: W1, W2, W3
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (1.0/3.0)*( 4.0*QCoeff(-2)*QCoeff(-2) - 19.0*QCoeff(-2)*QCoeff(-1) + 25.0*QCoeff(-1)*QCoeff(-1) &
                  + 11.0*QCoeff(-2)*QCoeff(+0) - 31.0*QCoeff(-1)*QCoeff(+0) + 10.0*QCoeff(+0)*QCoeff(+0))
beta2 = (1.0/3.0)*( 4.0*QCoeff(-1)*QCoeff(-1) - 13.0*QCoeff(-1)*QCoeff(+0) + 13.0*QCoeff(+0)*QCoeff(+0) &
                  +  5.0*QCoeff(-1)*QCoeff(+1) - 13.0*QCoeff(+0)*QCoeff(+1) +  4.0*QCoeff(+1)*QCoeff(+1))
beta3 = (1.0/3.0)*(10.0*QCoeff(+0)*QCoeff(+0) - 31.0*QCoeff(+0)*QCoeff(+1) + 25.0*QCoeff(+1)*QCoeff(+1) &
                  + 11.0*QCoeff(+0)*QCoeff(+2) - 19.0*QCoeff(+1)*QCoeff(+2) +  4.0*QCoeff(+2)*QCoeff(+2))

!--------------------------------------------!
! Point: x_{j-1/2*sqrt(3/7+2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 =0.2658420974778319
gamma2 =0.6112504900322486
gamma3 =0.1229074124899195

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1 = -0.1642562761719537 * Q(-2)+0.7590807081409336 * Q(-1)+0.4051755680310201 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)
W2 = +0.0000000000000000 * Q(-2)+0.2663118796250726 * Q(-1)+0.8979443965468811 * Q(0)-0.1642562761719537 * Q(1)+0.0000000000000000 * Q(2)
W3 = +0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+1.6968800354220990 * Q(0)-0.9631919150471715 * Q(1)+0.2663118796250726 * Q(2)
W(1)   = (omega1*W1 + omega2*W2 + omega3*W3)

!--------------------------------------------!
! Point: x_{j-1/2*sqrt(3/7-2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights

gamma1 =0.1281641584355059
gamma2 =0.5219691498503563
gamma3 =0.3498666917141379

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)


W1 = -0.1122135388132497 * Q(-2)+0.3944175994189276 * Q(-1)+0.7177959393943222 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)
W2 = +0.0000000000000000 * Q(-2)+0.0577769829791784 * Q(-1)+1.0544365558340714 * Q(0)-0.1122135388132497 * Q(1)+0.0000000000000000 * Q(2)
W3 = +0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+1.2277675047716066 * Q(0)-0.2855444877507849 * Q(1)+0.0577769829791784 * Q(2)

W(2)   = (omega1*W1 + omega2*W2 + omega3*W3)


!--------------------------------------------!
! Point: x_{j+1/2*sqrt(3/7-2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 =0.3498666917141379
gamma2 =0.5219691498503563
gamma3 =0.1281641584355059

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial

W1 = +0.0577769829791784 * Q(-2)-0.2855444877507849 * Q(-1)+1.2277675047716066 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)
W2 = +0.0000000000000000 * Q(-2)-0.1122135388132497 * Q(-1)+1.0544365558340714 * Q(0)+0.0577769829791784 * Q(1)+0.0000000000000000 * Q(2)
W3 = +0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.7177959393943222 * Q(0)+0.3944175994189276 * Q(1)-0.1122135388132497 * Q(2)

W(3)   = (omega1*W1 + omega2*W2 + omega3*W3)



!--------------------------------------------!
! Point: x_{j+1/2*sqrt(3/7+2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 =0.1229074124899195
gamma2 =0.6112504900322486
gamma3 =0.2658420974778319


alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1 = +0.2663118796250726 * Q(-2)-0.9631919150471715 * Q(-1)+1.6968800354220990 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)
W2 = +0.0000000000000000 * Q(-2)-0.1642562761719537 * Q(-1)+0.8979443965468811 * Q(0)+0.2663118796250726 * Q(1)+0.0000000000000000 * Q(2)
W3 = +0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.4051755680310201 * Q(0)+0.7590807081409336 * Q(1)-0.1642562761719537 * Q(2)
W(4)   = (omega1*W1 + omega2*W2 + omega3*W3)

!-------------------------------------------------------------------------------!
END SUBROUTINE WENO5_SecondSweep
!===============================================================================!
!
!
!
!
!
!===============================================================================!
   SUBROUTINE WENO7_FirstSweep(Q,WM,WP)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
      USE MOD_FiniteVolume2D_vars,ONLY: LinearWeightsOnly
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: WM
      REAL,INTENT(OUT) :: WP
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: alpha1, alpha2, alpha3, alpha4
      REAL             :: beta1,  beta2,  beta3,  beta4
      REAL             :: gamma1, gamma2, gamma3, gamma4
      REAL             :: omega1, omega2, omega3, omega4
      REAL             :: W1, W2, W3, W4
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
! beta1 = +(547./240.*Q(-3)*Q(-3))+(-647./80.*Q(-3)*Q(-2))+(2321./240.*Q(-3)*Q(-1))+(-309./80.*Q(-3)*Q(0))+(-647./80.*Q(-2)*Q(-3))+(7043./240.*Q(-2)*Q(-2))+(-8623./240.*Q(-2)*Q(-1))+(3521./240.*Q(-2)*Q(0))+(2321./240.*Q(-1)*Q(-3))+(-8623./240.*Q(-1)*Q(-2))+(11003/240*Q(-1)*Q(-1))+(-1567/80*Q(-1)*Q(0))+(-309/80*Q(0)*Q(-3))+(3521/240*Q(0)*Q(-2))+(-1567/80*Q(0)*Q(-1))+(2107/240*Q(0)*Q(0))
! beta2 = +(89/80*Q(-2)*Q(-2))+(-821/240*Q(-2)*Q(-1))+(267/80*Q(-2)*Q(0))+(-247/240*Q(-2)*Q(1))+(-821/240*Q(-1)*Q(-2))+(2843/240*Q(-1)*Q(-1))+(-2983/240*Q(-1)*Q(0))+(961/240*Q(-1)*Q(1))+(267/80*Q(0)*Q(-2))+(-2983/240*Q(0)*Q(-1))+(3443/240*Q(0)*Q(0))+(-1261/240*Q(0)*Q(1))+(-247/240*Q(1)*Q(-2))+(961/240*Q(1)*Q(-1))+(-1261/240*Q(1)*Q(0))+(547/240*Q(1)*Q(1))
! beta3 = +(547/240*Q(-1)*Q(-1))+(-1261/240*Q(-1)*Q(0))+(961/240*Q(-1)*Q(1))+(-247/240*Q(-1)*Q(2))+(-1261/240*Q(0)*Q(-1))+(3443/240*Q(0)*Q(0))+(-2983/240*Q(0)*Q(1))+(267/80*Q(0)*Q(2))+(961/240*Q(1)*Q(-1))+(-2983/240*Q(1)*Q(0))+(2843/240*Q(1)*Q(1))+(-821/240*Q(1)*Q(2))+(-247/240*Q(2)*Q(-1))+(267/80*Q(2)*Q(0))+(-821/240*Q(2)*Q(1))+(89/80*Q(2)*Q(2))
! beta4 = +(2107/240*Q(0)*Q(0))+(-1567/80*Q(0)*Q(1))+(3521/240*Q(0)*Q(2))+(-309/80*Q(0)*Q(3))+(-1567/80*Q(1)*Q(0))+(11003/240*Q(1)*Q(1))+(-8623/240*Q(1)*Q(2))+(2321/240*Q(1)*Q(3))+(3521/240*Q(2)*Q(0))+(-8623/240*Q(2)*Q(1))+(7043/240*Q(2)*Q(2))+(-647/80*Q(2)*Q(3))+(-309/80*Q(3)*Q(0))+(2321/240*Q(3)*Q(1))+(-647/80*Q(3)*Q(2))+(547/240*Q(3)*Q(3))

      beta1 = +(2.27916666666666679*Q(-3)*Q(-3))+(-8.08750000000000036*Q(-3)*Q(-2))+(9.67083333333333250*Q(-3)*Q(-1))+(-3.8625000000000000*Q(-3)*Q(0))+(-8.08750000000000036*Q(-2)*Q(-3))+(29.34583333333333499*Q(-2)*Q(-2))+(-35.92916666666666714*Q(-2)*Q(-1))+(14.67083333333333250*Q(-2)*Q(0))+(9.67083333333333250*Q(-1)*Q(-3))+(-35.92916666666666714*Q(-1)*Q(-2))+(45.84583333333333144*Q(-1)*Q(-1))+(-19.5875000000000000*Q(-1)*Q(0))+(-3.8625000000000000*Q(0)*Q(-3))+(14.67083333333333250*Q(0)*Q(-2))+(-19.58749999999999858*Q(0)*Q(-1))+(8.77916666666666679*Q(0)*Q(0))
      beta2 = +(1.11250000000000004*Q(-2)*Q(-2))+(-3.42083333333333339*Q(-2)*Q(-1))+(3.33749999999999991*Q(-2)*Q(0))+(-1.02916666666666656*Q(-2)*Q(1))+(-3.42083333333333339*Q(-1)*Q(-2))+(11.84583333333333321*Q(-1)*Q(-1))+(-12.42916666666666714*Q(-1)*Q(0))+(4.00416666666666643*Q(-1)*Q(1))+(3.3375000000000000*Q(0)*Q(-2))+(-12.42916666666666714*Q(0)*Q(-1))+(14.34583333333333321*Q(0)*Q(0))+(-5.25416666666666643*Q(0)*Q(1))+(-1.02916666666666656*Q(1)*Q(-2))+(4.00416666666666643*Q(1)*Q(-1))+(-5.25416666666666643*Q(1)*Q(0))+(2.27916666666666679*Q(1)*Q(1))
      beta3 = +(2.27916666666666679*Q(-1)*Q(-1))+(-5.25416666666666643*Q(-1)*Q(0))+(4.00416666666666643*Q(-1)*Q(1))+(-1.02916666666666656*Q(-1)*Q(2))+(-5.25416666666666643*Q(0)*Q(-1))+(14.34583333333333321*Q(0)*Q(0))+(-12.42916666666666714*Q(0)*Q(1))+(3.33749999999999991*Q(0)*Q(2))+(4.00416666666666643*Q(1)*Q(-1))+(-12.42916666666666714*Q(1)*Q(0))+(11.84583333333333321*Q(1)*Q(1))+(-3.42083333333333339*Q(1)*Q(2))+(-1.02916666666666656*Q(2)*Q(-1))+(3.3375000000000000*Q(2)*Q(0))+(-3.42083333333333339*Q(2)*Q(1))+(1.11250000000000004*Q(2)*Q(2))
      beta4 = +(8.77916666666666679*Q(0)*Q(0))+(-19.5875000000000000*Q(0)*Q(1))+(14.67083333333333250*Q(0)*Q(2))+(-3.8625000000000000*Q(0)*Q(3))+(-19.5875000000000000*Q(1)*Q(0))+(45.84583333333333144*Q(1)*Q(1))+(-35.92916666666666714*Q(1)*Q(2))+(9.67083333333333250*Q(1)*Q(3))+(14.67083333333333250*Q(2)*Q(0))+(-35.92916666666666714*Q(2)*Q(1))+(29.34583333333333499*Q(2)*Q(2))+(-8.08750000000000036*Q(2)*Q(3))+(-3.86249999999999982*Q(3)*Q(0))+(9.67083333333333250*Q(3)*Q(1))+(-8.08750000000000036*Q(3)*Q(2))+(2.27916666666666679*Q(3)*Q(3))


!------------------------------!
! WM: x_{i-1/2}                !
!------------------------------!

! Linear Weights
      gamma1 =0.1142857142857143
      gamma2 =0.5142857142857142
      gamma3 =0.3428571428571429
      gamma4 =0.0285714285714286

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
      alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4)
      omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4)

      W1 = +0.0833333333333333 * Q(-3)-0.4166666666666667 * Q(-2)+1.0833333333333333 * Q(-1)+0.2500000000000000 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
      W2 = +0.0000000000000000 * Q(-3)-0.0833333333333333 * Q(-2)+0.5833333333333334 * Q(-1)+0.5833333333333334 * Q(0)-0.0833333333333333 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
      W3 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.2500000000000000 * Q(-1)+1.0833333333333333 * Q(0)-0.4166666666666667 * Q(1)+0.0833333333333333 * Q(2)+0.0000000000000000 * Q(3)
      W4 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+2.0833333333333335 * Q(0)-1.9166666666666667 * Q(1)+1.0833333333333333 * Q(2)-0.2500000000000000 * Q(3)



      IF ( LinearWeightsOnly ) THEN
         WM = gamma1*W1 + gamma2*W2 + gamma3*W3 + gamma4*W4
      ELSE
         WM = omega1*W1 + omega2*W2 + omega3*W3 + omega4*W4
      END IF


!------------------------------!
! WP: x_{i+1/2}                !
!------------------------------!

! Linear Weights
      gamma1 =0.0285714285714286
      gamma2 =0.3428571428571429
      gamma3 =0.5142857142857142
      gamma4 =0.1142857142857143

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
      alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4)
      omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4)

      W1 = -0.2500000000000000 * Q(-3)+1.0833333333333333 * Q(-2)-1.9166666666666667 * Q(-1)+2.0833333333333335 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
      W2 = +0.0000000000000000 * Q(-3)+0.0833333333333333 * Q(-2)-0.4166666666666667 * Q(-1)+1.0833333333333333 * Q(0)+0.2500000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
      W3 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)-0.0833333333333333 * Q(-1)+0.5833333333333334 * Q(0)+0.5833333333333334 * Q(1)-0.0833333333333333 * Q(2)+0.0000000000000000 * Q(3)
      W4 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.2500000000000000 * Q(0)+1.0833333333333333 * Q(1)-0.4166666666666667 * Q(2)+0.0833333333333333 * Q(3)

      IF ( LinearWeightsOnly ) THEN
         WP = gamma1*W1 + gamma2*W2 + gamma3*W3 + gamma4*W4
      ELSE
         WP = omega1*W1 + omega2*W2 + omega3*W3 + omega4*W4
      END IF

!-------------------------------------------------------------------------------!
   END SUBROUTINE WENO7_FirstSweep
!===============================================================================!
!
!
!
!===============================================================================!
   SUBROUTINE WENO7_SecondSweep(Q,W) !4nGPS
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
      USE MOD_FiniteVolume2D_vars,ONLY: LinearWeightsOnly
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: W(1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: alpha1, alpha2, alpha3, alpha4
      REAL             :: beta1,  beta2,  beta3,  beta4
      REAL             :: gamma1, gamma2, gamma3, gamma4
      REAL             :: omega1, omega2, omega3, omega4
      REAL             :: W1, W2, W3, W4
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
      beta1 = +(2.27916666666666679*Q(-3)*Q(-3))+(-8.08750000000000036*Q(-3)*Q(-2))+(9.67083333333333250*Q(-3)*Q(-1))+(-3.86249999999999982*Q(-3)*Q(0))+(-8.08750000000000036*Q(-2)*Q(-3))+(29.34583333333333499*Q(-2)*Q(-2))+(-35.92916666666666714*Q(-2)*Q(-1))+(14.67083333333333250*Q(-2)*Q(0))+(9.67083333333333250*Q(-1)*Q(-3))+(-35.92916666666666714*Q(-1)*Q(-2))+(45.84583333333333144*Q(-1)*Q(-1))+(-19.58749999999999858*Q(-1)*Q(0))+(-3.86249999999999982*Q(0)*Q(-3))+(14.67083333333333250*Q(0)*Q(-2))+(-19.58749999999999858*Q(0)*Q(-1))+(8.77916666666666679*Q(0)*Q(0))
      beta2 = +(1.11250000000000004*Q(-2)*Q(-2))+(-3.42083333333333339*Q(-2)*Q(-1))+(3.33749999999999991*Q(-2)*Q(0))+(-1.02916666666666656*Q(-2)*Q(1))+(-3.42083333333333339*Q(-1)*Q(-2))+(11.84583333333333321*Q(-1)*Q(-1))+(-12.42916666666666714*Q(-1)*Q(0))+(4.00416666666666643*Q(-1)*Q(1))+(3.33749999999999991*Q(0)*Q(-2))+(-12.42916666666666714*Q(0)*Q(-1))+(14.34583333333333321*Q(0)*Q(0))+(-5.25416666666666643*Q(0)*Q(1))+(-1.02916666666666656*Q(1)*Q(-2))+(4.00416666666666643*Q(1)*Q(-1))+(-5.25416666666666643*Q(1)*Q(0))+(2.27916666666666679*Q(1)*Q(1))
      beta3 = +(2.27916666666666679*Q(-1)*Q(-1))+(-5.25416666666666643*Q(-1)*Q(0))+(4.00416666666666643*Q(-1)*Q(1))+(-1.02916666666666656*Q(-1)*Q(2))+(-5.25416666666666643*Q(0)*Q(-1))+(14.34583333333333321*Q(0)*Q(0))+(-12.42916666666666714*Q(0)*Q(1))+(3.33749999999999991*Q(0)*Q(2))+(4.00416666666666643*Q(1)*Q(-1))+(-12.42916666666666714*Q(1)*Q(0))+(11.84583333333333321*Q(1)*Q(1))+(-3.42083333333333339*Q(1)*Q(2))+(-1.02916666666666656*Q(2)*Q(-1))+(3.33749999999999991*Q(2)*Q(0))+(-3.42083333333333339*Q(2)*Q(1))+(1.11250000000000004*Q(2)*Q(2))
      beta4 = +(8.77916666666666679*Q(0)*Q(0))+(-19.58749999999999858*Q(0)*Q(1))+(14.67083333333333250*Q(0)*Q(2))+(-3.86249999999999982*Q(0)*Q(3))+(-19.58749999999999858*Q(1)*Q(0))+(45.84583333333333144*Q(1)*Q(1))+(-35.92916666666666714*Q(1)*Q(2))+(9.67083333333333250*Q(1)*Q(3))+(14.67083333333333250*Q(2)*Q(0))+(-35.92916666666666714*Q(2)*Q(1))+(29.34583333333333499*Q(2)*Q(2))+(-8.08750000000000036*Q(2)*Q(3))+(-3.86249999999999982*Q(3)*Q(0))+(9.67083333333333250*Q(3)*Q(1))+(-8.08750000000000036*Q(3)*Q(2))+(2.27916666666666679*Q(3)*Q(3))

!--------------------------------------------!
! Point: x_{j-1/2*sqrt(3/7+2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
      gamma1 =0.0978973393748262
      gamma2 =0.4945893001737974
      gamma3 =0.3706410089227926
      gamma4 =0.0368723515285839

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
      alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4)
      omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4)

! Reconstructed Polynomial
      W1 = +0.0878583391504589 * Q(-3)-0.4278312936233303 * Q(-2)+1.0226557255923103 * Q(-1)+0.3173172288805612 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
      W2 = +0.0000000000000000 * Q(-3)-0.0763979370214948 * Q(-2)+0.4955056906895569 * Q(-1)+0.6687505854823967 * Q(0)-0.0878583391504589 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
      W3 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.1899139426035779 * Q(-1)+1.1271382076113654 * Q(0)-0.3934500872364380 * Q(1)+0.0763979370214948 * Q(2)+0.0000000000000000 * Q(3)
      W4 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+1.8867939780256768 * Q(0)-1.5329337428579051 * Q(1)+0.8360537074358062 * Q(2)-0.1899139426035779 * Q(3)

      IF (LinearWeightsOnly) THEN
         W(1)   = (gamma1*W1 + gamma2*W2 + gamma3*W3 +gamma4*W4)
      ELSE
         W(1)   = (omega1*W1 + omega2*W2 + omega3*W3 +omega4*W4)
      END IF


!--------------------------------------------!
! Point: x_{j-1/2*sqrt(3/7-2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights

      gamma1 =0.0422160570229445
      gamma2 =0.3488436391325038
      gamma3 =0.4308046279496520
      gamma4 =0.1781356758948997

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
      alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4)
      omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4)


      W1 = +0.0776175431540304 * Q(-3)-0.3450661682753410 * Q(-2)+0.6272702288810189 * Q(-1)+0.6401783962402917 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
      W2 = +0.0000000000000000 * Q(-3)-0.0345959956592193 * Q(-2)+0.1615649699568364 * Q(-1)+0.9506485688564134 * Q(0)-0.0776175431540304 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
      W3 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0231809873199591 * Q(-1)+1.1582245428117293 * Q(0)-0.2160015257909077 * Q(1)+0.0345959956592193 * Q(2)+0.0000000000000000 * Q(3)
      W4 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+1.2509484920915657 * Q(0)-0.3550874497106621 * Q(1)+0.1273199449390556 * Q(2)-0.0231809873199591 * Q(3)

      IF (LinearWeightsOnly) THEN
         W(2)   = (gamma1*W1 + gamma2*W2 + gamma3*W3 +gamma4*W4)
      ELSE
         W(2)   = (omega1*W1 + omega2*W2 + omega3*W3 +omega4*W4)
      END IF



!--------------------------------------------!
! Point: x_{j+1/2*sqrt(3/7-2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
      gamma1 =0.1781356758948997
      gamma2 =0.4308046279496520
      gamma3 =0.3488436391325038
      gamma4 =0.0422160570229445

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
      alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4)
      omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4)


! Reconstructed Polynomial

      W1 = -0.0231809873199591 * Q(-3)+0.1273199449390556 * Q(-2)-0.3550874497106621 * Q(-1)+1.2509484920915657 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
      W2 = +0.0000000000000000 * Q(-3)+0.0345959956592193 * Q(-2)-0.2160015257909077 * Q(-1)+1.1582245428117293 * Q(0)+0.0231809873199591 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
      W3 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)-0.0776175431540304 * Q(-1)+0.9506485688564134 * Q(0)+0.1615649699568364 * Q(1)-0.0345959956592193 * Q(2)+0.0000000000000000 * Q(3)
      W4 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.6401783962402917 * Q(0)+0.6272702288810189 * Q(1)-0.3450661682753410 * Q(2)+0.0776175431540304 * Q(3)

      IF (LinearWeightsOnly) THEN
         W(3)   = (gamma1*W1 + gamma2*W2 + gamma3*W3 +gamma4*W4)
      ELSE
         W(3)   = (omega1*W1 + omega2*W2 + omega3*W3 +omega4*W4)
      END IF



!--------------------------------------------!
! Point: x_{j+1/2*sqrt(3/7+2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
      gamma1 =0.0368723515285839
      gamma2 =0.3706410089227926
      gamma3 =0.4945893001737974
      gamma4 =0.0978973393748262

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
      alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4)
      omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4)

! Reconstructed Polynomial
      W1 = -0.1899139426035779 * Q(-3)+0.8360537074358062 * Q(-2)-1.5329337428579051 * Q(-1)+1.8867939780256768 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
      W2 = +0.0000000000000000 * Q(-3)+0.0763979370214948 * Q(-2)-0.3934500872364380 * Q(-1)+1.1271382076113654 * Q(0)+0.1899139426035779 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)
      W3 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)-0.0878583391504589 * Q(-1)+0.6687505854823967 * Q(0)+0.4955056906895569 * Q(1)-0.0763979370214948 * Q(2)+0.0000000000000000 * Q(3)
      W4 = +0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.3173172288805612 * Q(0)+1.0226557255923103 * Q(1)-0.4278312936233303 * Q(2)+0.0878583391504589 * Q(3)

      IF (LinearWeightsOnly) THEN
         W(4)   = (gamma1*W1 + gamma2*W2 + gamma3*W3 +gamma4*W4)
      ELSE
         W(4)   = (omega1*W1 + omega2*W2 + omega3*W3 +omega4*W4)
      END IF
!-------------------------------------------------------------------------------!
   END SUBROUTINE WENO7_SecondSweep
!===============================================================================!
!
!
!
!
!===============================================================================!
   SUBROUTINE WENO9_FirstSweep(Q,WM,WP)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
      USE MOD_FiniteVolume2D_vars,ONLY: LinearWeightsOnly
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: WM
      REAL,INTENT(OUT) :: WP
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: alpha1, alpha2, alpha3, alpha4,  alpha5
      REAL             :: beta1,  beta2,  beta3,  beta4,   beta5
      REAL             :: gamma1, gamma2, gamma3, gamma4,  gamma5
      REAL             :: omega1, omega2, omega3, omega4,  omega5
      REAL             :: W1, W2, W3, W4, W5
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!

      beta1 = +(4.49563492063492021*Q(-4)*Q(-4))+(-20.68462301587301511*Q(-4)*Q(-3))+(36.19672619047619122*Q(-4)*Q(-2))+(-28.57212301587301440*Q(-4)*Q(-1))+(8.56438492063491985*Q(-4)*Q(0))+(-20.68462301587301511*Q(-3)*Q(-4))+(95.82599206349206611*Q(-3)*Q(-3))+(-169.08690476190474783*Q(-3)*Q(-2))+(134.76765873015872899*Q(-3)*Q(-1))+(-40.82212301587301795*Q(-3)*Q(0))+(36.19672619047619122*Q(-2)*Q(-4))+(-169.08690476190474783*Q(-2)*Q(-3))+(301.86369047619047024*Q(-2)*Q(-2))+(-244.25357142857143344*Q(-2)*Q(-1))+(75.28005952380952692*Q(-2)*Q(0))+(-28.57212301587301440*Q(-1)*Q(-4))+(134.76765873015872899*Q(-1)*Q(-3))+(-244.25357142857143344*Q(-1)*Q(-2))+(202.49265873015872330*Q(-1)*Q(-1))+(-64.43462301587301511*Q(-1)*Q(0))+(8.56438492063491985*Q(0)*Q(-4))+(-40.82212301587301795*Q(0)*Q(-3))+(75.28005952380952692*Q(0)*Q(-2))+(-64.43462301587301511*Q(0)*Q(-1))+(21.41230158730158806*Q(0)*Q(0))
      beta2 = +(1.37063492063492065*Q(-3)*Q(-3))+(-6.03878968253968296*Q(-3)*Q(-2))+(9.84255952380952337*Q(-3)*Q(-1))+(-6.96795634920634921*Q(-3)*Q(0))+(1.79355158730158726*Q(-3)*Q(1))+(-6.03878968253968296*Q(-2)*Q(-3))+(27.49265873015873041*Q(-2)*Q(-2))+(-46.12857142857142634*Q(-2)*Q(-1))+(33.43432539682539328*Q(-2)*Q(0))+(-8.75962301587301617*Q(-2)*Q(1))+(9.84255952380952337*Q(-1)*Q(-3))+(-46.12857142857142634*Q(-1)*Q(-2))+(80.61369047619047024*Q(-1)*Q(-1))+(-60.71190476190476204*Q(-1)*Q(0))+(16.38422619047619122*Q(-1)*Q(1))+(-6.96795634920634921*Q(0)*Q(-3))+(33.43432539682539328*Q(0)*Q(-2))+(-60.71190476190476204*Q(0)*Q(-1))+(48.15932539682539471*Q(0)*Q(0))+(-13.91378968253968296*Q(0)*Q(1))+(1.79355158730158726*Q(1)*Q(-3))+(-8.75962301587301617*Q(1)*Q(-2))+(16.38422619047619122*Q(1)*Q(-1))+(-13.91378968253968296*Q(1)*Q(0))+(4.49563492063492021*Q(1)*Q(1))
      beta3 = +(1.37063492063492065*Q(-2)*Q(-2))+(-5.05962301587301599*Q(-2)*Q(-1))+(6.73839285714285730*Q(-2)*Q(0))+(-3.86378968253968269*Q(-2)*Q(1))+(0.81438492063492063*Q(-2)*Q(2))+(-5.05962301587301599*Q(-1)*Q(-2))+(20.82599206349206256*Q(-1)*Q(-1))+(-29.67023809523809419*Q(-1)*Q(0))+(17.76765873015872899*Q(-1)*Q(1))+(-3.86378968253968269*Q(-1)*Q(2))+(6.73839285714285730*Q(0)*Q(-2))+(-29.67023809523809419*Q(0)*Q(-1))+(45.86369047619047734*Q(0)*Q(0))+(-29.67023809523809419*Q(0)*Q(1))+(6.73839285714285730*Q(0)*Q(2))+(-3.86378968253968269*Q(1)*Q(-2))+(17.76765873015872899*Q(1)*Q(-1))+(-29.67023809523809419*Q(1)*Q(0))+(20.82599206349206256*Q(1)*Q(1))+(-5.05962301587301599*Q(1)*Q(2))+(0.81438492063492063*Q(2)*Q(-2))+(-3.86378968253968269*Q(2)*Q(-1))+(6.73839285714285730*Q(2)*Q(0))+(-5.05962301587301599*Q(2)*Q(1))+(1.37063492063492065*Q(2)*Q(2))
      beta4 = +(4.49563492063492021*Q(-1)*Q(-1))+(-13.91378968253968296*Q(-1)*Q(0))+(16.38422619047619122*Q(-1)*Q(1))+(-8.75962301587301617*Q(-1)*Q(2))+(1.79355158730158726*Q(-1)*Q(3))+(-13.91378968253968296*Q(0)*Q(-1))+(48.15932539682539471*Q(0)*Q(0))+(-60.71190476190476204*Q(0)*Q(1))+(33.43432539682539328*Q(0)*Q(2))+(-6.96795634920634921*Q(0)*Q(3))+(16.38422619047619122*Q(1)*Q(-1))+(-60.71190476190476204*Q(1)*Q(0))+(80.61369047619047024*Q(1)*Q(1))+(-46.12857142857142634*Q(1)*Q(2))+(9.84255952380952337*Q(1)*Q(3))+(-8.75962301587301617*Q(2)*Q(-1))+(33.43432539682539328*Q(2)*Q(0))+(-46.12857142857142634*Q(2)*Q(1))+(27.49265873015873041*Q(2)*Q(2))+(-6.03878968253968296*Q(2)*Q(3))+(1.79355158730158726*Q(3)*Q(-1))+(-6.96795634920634921*Q(3)*Q(0))+(9.84255952380952337*Q(3)*Q(1))+(-6.03878968253968296*Q(3)*Q(2))+(1.37063492063492065*Q(3)*Q(3))
      beta5 = +(21.41230158730158806*Q(0)*Q(0))+(-64.43462301587301511*Q(0)*Q(1))+(75.28005952380952692*Q(0)*Q(2))+(-40.82212301587301795*Q(0)*Q(3))+(8.56438492063491985*Q(0)*Q(4))+(-64.43462301587301511*Q(1)*Q(0))+(202.49265873015872330*Q(1)*Q(1))+(-244.25357142857143344*Q(1)*Q(2))+(134.76765873015872899*Q(1)*Q(3))+(-28.57212301587301440*Q(1)*Q(4))+(75.28005952380952692*Q(2)*Q(0))+(-244.25357142857143344*Q(2)*Q(1))+(301.86369047619047024*Q(2)*Q(2))+(-169.08690476190474783*Q(2)*Q(3))+(36.19672619047619122*Q(2)*Q(4))+(-40.82212301587301795*Q(3)*Q(0))+(134.76765873015872899*Q(3)*Q(1))+(-169.08690476190474783*Q(3)*Q(2))+(95.82599206349206611*Q(3)*Q(3))+(-20.68462301587301511*Q(3)*Q(4))+(8.56438492063491985*Q(4)*Q(0))+(-28.57212301587301440*Q(4)*Q(1))+(36.19672619047619122*Q(4)*Q(2))+(-20.68462301587301511*Q(4)*Q(3))+(4.49563492063492021*Q(4)*Q(4))


!------------------------------!
! WM: x_{i-1/2}                !
!------------------------------!

! Linear Weights
      gamma1 =0.0396825396825397
      gamma2 =0.3174603174603174
      gamma3 =0.4761904761904762
      gamma4 =0.1587301587301587
      gamma5 =0.0079365079365079

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
      alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP
      alpha5 = gamma5/(WENOEPS + beta5)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5)
      omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5)
      omega5 = alpha5/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5)

      W1 = -0.0500000000000000 * Q(-4)+0.2833333333333333 * Q(-3)-0.7166666666666667 * Q(-2)+1.2833333333333334 * Q(-1)+0.2000000000000000 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)
      W2 = +0.0000000000000000 * Q(-4)+0.0333333333333333 * Q(-3)-0.2166666666666667 * Q(-2)+0.7833333333333333 * Q(-1)+0.4500000000000000 * Q(0)-0.0500000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)
      W3 = +0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)-0.0500000000000000 * Q(-2)+0.4500000000000000 * Q(-1)+0.7833333333333333 * Q(0)-0.2166666666666667 * Q(1)+0.0333333333333333 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)
      W4 = +0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.2000000000000000 * Q(-1)+1.2833333333333334 * Q(0)-0.7166666666666667 * Q(1)+0.2833333333333333 * Q(2)-0.0500000000000000 * Q(3)+0.0000000000000000 * Q(4)
      W5 = +0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+2.2833333333333332 * Q(0)-2.7166666666666668 * Q(1)+2.2833333333333332 * Q(2)-1.0500000000000000 * Q(3)+0.2000000000000000 * Q(4)



      IF ( LinearWeightsOnly ) THEN
         WM = gamma1*W1 + gamma2*W2 + gamma3*W3 + gamma4*W4 + gamma5*W5
      ELSE
         WM = omega1*W1 + omega2*W2 + omega3*W3 + omega4*W4 + omega5*W5
      END IF


!------------------------------!
! WP: x_{i+1/2}                !
!------------------------------!

! Linear Weights
      gamma1 =0.0079365079365079
      gamma2 =0.1587301587301587
      gamma3 =0.4761904761904762
      gamma4 =0.3174603174603174
      gamma5 =0.0396825396825397

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
      alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP
      alpha5 = gamma5/(WENOEPS + beta5)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5)
      omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5)
      omega5 = alpha5/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5)

      W1 = +0.2000000000000000 * Q(-4)-1.0500000000000000 * Q(-3)+2.2833333333333332 * Q(-2)-2.7166666666666668 * Q(-1)+2.2833333333333332 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)
      W2 = +0.0000000000000000 * Q(-4)-0.0500000000000000 * Q(-3)+0.2833333333333333 * Q(-2)-0.7166666666666667 * Q(-1)+1.2833333333333334 * Q(0)+0.2000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)
      W3 = +0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0333333333333333 * Q(-2)-0.2166666666666667 * Q(-1)+0.7833333333333333 * Q(0)+0.4500000000000000 * Q(1)-0.0500000000000000 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)
      W4 = +0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)-0.0500000000000000 * Q(-1)+0.4500000000000000 * Q(0)+0.7833333333333333 * Q(1)-0.2166666666666667 * Q(2)+0.0333333333333333 * Q(3)+0.0000000000000000 * Q(4)
      W5 = +0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.2000000000000000 * Q(0)+1.2833333333333334 * Q(1)-0.7166666666666667 * Q(2)+0.2833333333333333 * Q(3)-0.0500000000000000 * Q(4)

      IF ( LinearWeightsOnly ) THEN
         WP = gamma1*W1 + gamma2*W2 + gamma3*W3 + gamma4*W4 + gamma5*W5
      ELSE
         WP = omega1*W1 + omega2*W2 + omega3*W3 + omega4*W4 + omega5*W5
      END IF

!-------------------------------------------------------------------------------!
   END SUBROUTINE WENO9_FirstSweep
!===============================================================================!
!
!
!
!
!===============================================================================!
   SUBROUTINE WENO11_FirstSweep(Q,WM,WP)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
      USE MOD_FiniteVolume2D_vars,ONLY: LinearWeightsOnly
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: WM
      REAL,INTENT(OUT) :: WP
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: alpha1, alpha2, alpha3, alpha4, alpha5, alpha6
      REAL             :: beta1,  beta2,  beta3,  beta4,  beta5,  beta6
      REAL             :: gamma1, gamma2, gamma3, gamma4, gamma5, gamma6
      REAL             :: omega1, omega2, omega3, omega4, omega5, omega6
      REAL             :: W1, W2, W3, W4, W5, W6
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!

      beta1 = +(9.52844742063492056*Q(-5)*Q(-5))+(-53.53085317460317327*Q(-5)*Q(-4))+(121.70244708994708560*Q(-5)*Q(-3))+(-140.20669642857143344*Q(-5)*Q(-2))+(81.98722718253968367*Q(-5)*Q(-1))+(-19.48057208994708844*Q(-5)*Q(0))+(-53.53085317460317327*Q(-4)*Q(-5))+(301.59298115079366198*Q(-4)*Q(-4))+(-688.08301917989422236*Q(-4)*Q(-3))+(796.11636904761905953*Q(-4)*Q(-2))+(-467.95133928571431170*Q(-4)*Q(-1))+(111.85586144179893608*Q(-4)*Q(0))+(121.70244708994708560*Q(-3)*Q(-5))+(-688.08301917989422236*Q(-3)*Q(-4))+(1577.03019179894181434*Q(-3)*Q(-3))+(-1835.33359788359780396*Q(-3)*Q(-2))+(1086.72979497354504019*Q(-3)*Q(-1))+(-262.04581679894181434*Q(-3)*Q(0))+(-140.20669642857143344*Q(-2)*Q(-5))+(796.11636904761905953*Q(-2)*Q(-4))+(-1835.33359788359780396*Q(-2)*Q(-3))+(2153.15287698412703321*Q(-2)*Q(-2))+(-1288.73695436507932754*Q(-2)*Q(-1))+(315.00800264550264274*Q(-2)*Q(0))+(81.98722718253968367*Q(-1)*Q(-5))+(-467.95133928571431170*Q(-1)*Q(-4))+(1086.72979497354504019*Q(-1)*Q(-3))+(-1288.73695436507932754*Q(-1)*Q(-2))+(784.15374503968257613*Q(-1)*Q(-1))+(-196.18247354497353285*Q(-1)*Q(0))+(-19.48057208994708844*Q(0)*Q(-5))+(111.85586144179893608*Q(0)*Q(-4))+(-262.04581679894181434*Q(0)*Q(-3))+(315.00800264550264274*Q(0)*Q(-2))+(-196.18247354497353285*Q(0)*Q(-1))+(50.84499834656084971*Q(0)*Q(0))
      beta2 = +(2.24685019841269851*Q(-4)*Q(-4))+(-12.46580687830687850*Q(-4)*Q(-3))+(27.67281746031746081*Q(-4)*Q(-2))+(-30.62544642857142918*Q(-4)*Q(-1))+(16.81141699735449890*Q(-4)*Q(0))+(-3.63983134920634921*Q(-4)*Q(1))+(-12.46580687830687850*Q(-3)*Q(-4))+(69.85744874338624300*Q(-3)*Q(-3))+(-156.71843584656085113*Q(-3)*Q(-2))+(175.28535052910052627*Q(-3)*Q(-1))+(-97.18282076719576423*Q(-3)*Q(0))+(21.22426421957672105*Q(-3)*Q(1))+(27.67281746031746081*Q(-2)*Q(-4))+(-156.71843584656085113*Q(-2)*Q(-3))+(356.26398809523811906*Q(-2)*Q(-2))+(-404.42619047619047024*Q(-2)*Q(-1))+(227.57007275132275481*Q(-2)*Q(0))+(-50.36225198412698489*Q(-2)*Q(1))+(-30.62544642857142918*Q(-1)*Q(-4))+(175.28535052910052627*Q(-1)*Q(-3))+(-404.42619047619047024*Q(-1)*Q(-2))+(468.43759920634920491*Q(-1)*Q(-1))+(-269.61079695767193698*Q(-1)*Q(0))+(60.93948412698412653*Q(-1)*Q(1))+(16.81141699735449890*Q(0)*Q(-4))+(-97.18282076719576423*Q(0)*Q(-3))+(227.57007275132275481*Q(0)*Q(-2))+(-269.61079695767193698*Q(0)*Q(-1))+(160.10224041005290019*Q(0)*Q(0))+(-37.69011243386243137*Q(0)*Q(1))+(-3.63983134920634921*Q(1)*Q(-4))+(21.22426421957672105*Q(1)*Q(-3))+(-50.36225198412698489*Q(1)*Q(-2))+(60.93948412698412653*Q(1)*Q(-1))+(-37.69011243386243137*Q(1)*Q(0))+(9.52844742063492056*Q(1)*Q(1))
      beta3 = +(1.15437334656084656*Q(-3)*Q(-3))+(-5.91094576719576725*Q(-3)*Q(-2))+(11.83855820105820023*Q(-3)*Q(-1))+(-11.54373346560846514*Q(-3)*Q(0))+(5.47704199735449748*Q(-3)*Q(1))+(-1.01529431216931210*Q(-3)*Q(2))+(-5.91094576719576725*Q(-2)*Q(-3))+(31.62075892857142989*Q(-2)*Q(-2))+(-65.64320436507937018*Q(-2)*Q(-1))+(65.84785052910052627*Q(-2)*Q(0))+(-31.94439484126984041*Q(-2)*Q(1))+(6.02993551587301546*Q(-2)*Q(2))+(11.83855820105820023*Q(-1)*Q(-3))+(-65.64320436507937018*Q(-1)*Q(-2))+(142.15982142857143344*Q(-1)*Q(-1))+(-148.05582010582011776*Q(-1)*Q(0))+(74.01220238095237391*Q(-1)*Q(1))+(-14.31155753968253919*Q(-1)*Q(2))+(-11.54373346560846514*Q(0)*Q(-3))+(65.84785052910052627*Q(0)*Q(-2))+(-148.05582010582011776*Q(0)*Q(-1))+(161.30102513227512873*Q(0)*Q(0))+(-84.44065806878306546*Q(0)*Q(1))+(16.89133597883597915*Q(0)*Q(2))+(5.47704199735449748*Q(1)*Q(-3))+(-31.94439484126984041*Q(1)*Q(-2))+(74.01220238095237391*Q(1)*Q(-1))+(-84.44065806878306546*Q(1)*Q(0))+(46.73707837301587631*Q(1)*Q(1))+(-9.84126984126984183*Q(1)*Q(2))+(-1.01529431216931210*Q(2)*Q(-3))+(6.02993551587301546*Q(2)*Q(-2))+(-14.31155753968253919*Q(2)*Q(-1))+(16.89133597883597915*Q(2)*Q(0))+(-9.84126984126984183*Q(2)*Q(1))+(2.24685019841269851*Q(2)*Q(2))
      beta4 = +(2.24685019841269851*Q(-2)*Q(-2))+(-9.84126984126984183*Q(-2)*Q(-1))+(16.89133597883597915*Q(-2)*Q(0))+(-14.31155753968253919*Q(-2)*Q(1))+(6.02993551587301546*Q(-2)*Q(2))+(-1.01529431216931210*Q(-2)*Q(3))+(-9.84126984126984183*Q(-1)*Q(-2))+(46.73707837301587631*Q(-1)*Q(-1))+(-84.44065806878306546*Q(-1)*Q(0))+(74.01220238095237391*Q(-1)*Q(1))+(-31.94439484126984041*Q(-1)*Q(2))+(5.47704199735449748*Q(-1)*Q(3))+(16.89133597883597915*Q(0)*Q(-2))+(-84.44065806878306546*Q(0)*Q(-1))+(161.30102513227512873*Q(0)*Q(0))+(-148.05582010582011776*Q(0)*Q(1))+(65.84785052910052627*Q(0)*Q(2))+(-11.54373346560846514*Q(0)*Q(3))+(-14.31155753968253919*Q(1)*Q(-2))+(74.01220238095237391*Q(1)*Q(-1))+(-148.05582010582011776*Q(1)*Q(0))+(142.15982142857143344*Q(1)*Q(1))+(-65.64320436507937018*Q(1)*Q(2))+(11.83855820105820023*Q(1)*Q(3))+(6.02993551587301546*Q(2)*Q(-2))+(-31.94439484126984041*Q(2)*Q(-1))+(65.84785052910052627*Q(2)*Q(0))+(-65.64320436507937018*Q(2)*Q(1))+(31.62075892857142989*Q(2)*Q(2))+(-5.91094576719576725*Q(2)*Q(3))+(-1.01529431216931210*Q(3)*Q(-2))+(5.47704199735449748*Q(3)*Q(-1))+(-11.54373346560846514*Q(3)*Q(0))+(11.83855820105820023*Q(3)*Q(1))+(-5.91094576719576725*Q(3)*Q(2))+(1.15437334656084656*Q(3)*Q(3))
      beta5 = +(9.52844742063492056*Q(-1)*Q(-1))+(-37.69011243386243137*Q(-1)*Q(0))+(60.93948412698412653*Q(-1)*Q(1))+(-50.36225198412698489*Q(-1)*Q(2))+(21.22426421957672105*Q(-1)*Q(3))+(-3.63983134920634921*Q(-1)*Q(4))+(-37.69011243386243137*Q(0)*Q(-1))+(160.10224041005290019*Q(0)*Q(0))+(-269.61079695767193698*Q(0)*Q(1))+(227.57007275132275481*Q(0)*Q(2))+(-97.18282076719576423*Q(0)*Q(3))+(16.81141699735449890*Q(0)*Q(4))+(60.93948412698412653*Q(1)*Q(-1))+(-269.61079695767193698*Q(1)*Q(0))+(468.43759920634920491*Q(1)*Q(1))+(-404.42619047619047024*Q(1)*Q(2))+(175.28535052910052627*Q(1)*Q(3))+(-30.62544642857142918*Q(1)*Q(4))+(-50.36225198412698489*Q(2)*Q(-1))+(227.57007275132275481*Q(2)*Q(0))+(-404.42619047619047024*Q(2)*Q(1))+(356.26398809523811906*Q(2)*Q(2))+(-156.71843584656085113*Q(2)*Q(3))+(27.67281746031746081*Q(2)*Q(4))+(21.22426421957672105*Q(3)*Q(-1))+(-97.18282076719576423*Q(3)*Q(0))+(175.28535052910052627*Q(3)*Q(1))+(-156.71843584656085113*Q(3)*Q(2))+(69.85744874338624300*Q(3)*Q(3))+(-12.46580687830687850*Q(3)*Q(4))+(-3.63983134920634921*Q(4)*Q(-1))+(16.81141699735449890*Q(4)*Q(0))+(-30.62544642857142918*Q(4)*Q(1))+(27.67281746031746081*Q(4)*Q(2))+(-12.46580687830687850*Q(4)*Q(3))+(2.24685019841269851*Q(4)*Q(4))
      beta6 = +(50.84499834656084971*Q(0)*Q(0))+(-196.18247354497353285*Q(0)*Q(1))+(315.00800264550264274*Q(0)*Q(2))+(-262.04581679894181434*Q(0)*Q(3))+(111.85586144179893608*Q(0)*Q(4))+(-19.48057208994708844*Q(0)*Q(5))+(-196.18247354497353285*Q(1)*Q(0))+(784.15374503968257613*Q(1)*Q(1))+(-1288.73695436507932754*Q(1)*Q(2))+(1086.72979497354504019*Q(1)*Q(3))+(-467.95133928571431170*Q(1)*Q(4))+(81.98722718253968367*Q(1)*Q(5))+(315.00800264550264274*Q(2)*Q(0))+(-1288.73695436507932754*Q(2)*Q(1))+(2153.15287698412703321*Q(2)*Q(2))+(-1835.33359788359780396*Q(2)*Q(3))+(796.11636904761905953*Q(2)*Q(4))+(-140.20669642857143344*Q(2)*Q(5))+(-262.04581679894181434*Q(3)*Q(0))+(1086.72979497354504019*Q(3)*Q(1))+(-1835.33359788359780396*Q(3)*Q(2))+(1577.03019179894181434*Q(3)*Q(3))+(-688.08301917989422236*Q(3)*Q(4))+(121.70244708994708560*Q(3)*Q(5))+(111.85586144179893608*Q(4)*Q(0))+(-467.95133928571431170*Q(4)*Q(1))+(796.11636904761905953*Q(4)*Q(2))+(-688.08301917989422236*Q(4)*Q(3))+(301.59298115079366198*Q(4)*Q(4))+(-53.53085317460317327*Q(4)*Q(5))+(-19.48057208994708844*Q(5)*Q(0))+(81.98722718253968367*Q(5)*Q(1))+(-140.20669642857143344*Q(5)*Q(2))+(121.70244708994708560*Q(5)*Q(3))+(-53.53085317460317327*Q(5)*Q(4))+(9.52844742063492056*Q(5)*Q(5))



!------------------------------!
! WM: x_{i-1/2}                !
!------------------------------!

! Linear Weights
      gamma1 =0.0129870129870130
      gamma2 =0.1623376623376623
      gamma3 =0.4329004329004329
      gamma4 =0.3246753246753247
      gamma5 =0.0649350649350649
      gamma6 =0.0021645021645022

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
      alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP
      alpha5 = gamma5/(WENOEPS + beta5)**WENOEXP
      alpha6 = gamma6/(WENOEPS + beta6)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6)
      omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6)
      omega5 = alpha5/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6)
      omega6 = alpha6/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6)

      W1 = +0.0333333333333333 * Q(-5)-0.2166666666666667 * Q(-4)+0.6166666666666667 * Q(-3)-1.0500000000000000 * Q(-2)+1.4500000000000000 * Q(-1)+0.1666666666666667 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)
      W2 = +0.0000000000000000 * Q(-5)-0.0166666666666667 * Q(-4)+0.1166666666666667 * Q(-3)-0.3833333333333334 * Q(-2)+0.9500000000000000 * Q(-1)+0.3666666666666666 * Q(0)-0.0333333333333333 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)
      W3 = +0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)+0.0166666666666667 * Q(-3)-0.1333333333333333 * Q(-2)+0.6166666666666667 * Q(-1)+0.6166666666666667 * Q(0)-0.1333333333333333 * Q(1)+0.0166666666666667 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)
      W4 = +0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)-0.0333333333333333 * Q(-2)+0.3666666666666666 * Q(-1)+0.9500000000000000 * Q(0)-0.3833333333333334 * Q(1)+0.1166666666666667 * Q(2)-0.0166666666666667 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)
      W5 = +0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.1666666666666667 * Q(-1)+1.4500000000000000 * Q(0)-1.0500000000000000 * Q(1)+0.6166666666666667 * Q(2)-0.2166666666666667 * Q(3)+0.0333333333333333 * Q(4)+0.0000000000000000 * Q(5)
      W6 = +0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+2.4500000000000002 * Q(0)-3.5499999999999998 * Q(1)+3.9500000000000002 * Q(2)-2.7166666666666668 * Q(3)+1.0333333333333334 * Q(4)-0.1666666666666667 * Q(5)



      IF ( LinearWeightsOnly ) THEN
         WM = gamma1*W1 + gamma2*W2 + gamma3*W3 + gamma4*W4 + gamma5*W5 + gamma6*W6
      ELSE
         WM = omega1*W1 + omega2*W2 + omega3*W3 + omega4*W4 + omega5*W5 + omega6*W6
      END IF


!------------------------------!
! WP: x_{i+1/2}                !
!------------------------------!

! Linear Weights
      gamma1 =0.0021645021645022
      gamma2 =0.0649350649350649
      gamma3 =0.3246753246753247
      gamma4 =0.4329004329004329
      gamma5 =0.1623376623376623
      gamma6 =0.0129870129870130

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
      alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP
      alpha5 = gamma5/(WENOEPS + beta5)**WENOEXP
      alpha6 = gamma6/(WENOEPS + beta6)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6)
      omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6)
      omega5 = alpha5/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6)
      omega6 = alpha6/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6)


      W1 = -0.1666666666666667 * Q(-5)+1.0333333333333334 * Q(-4)-2.7166666666666668 * Q(-3)+3.9500000000000002 * Q(-2)-3.5499999999999998 * Q(-1)+2.4500000000000002 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)
      W2 = +0.0000000000000000 * Q(-5)+0.0333333333333333 * Q(-4)-0.2166666666666667 * Q(-3)+0.6166666666666667 * Q(-2)-1.0500000000000000 * Q(-1)+1.4500000000000000 * Q(0)+0.1666666666666667 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)
      W3 = +0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)-0.0166666666666667 * Q(-3)+0.1166666666666667 * Q(-2)-0.3833333333333334 * Q(-1)+0.9500000000000000 * Q(0)+0.3666666666666666 * Q(1)-0.0333333333333333 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)
      W4 = +0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0166666666666667 * Q(-2)-0.1333333333333333 * Q(-1)+0.6166666666666667 * Q(0)+0.6166666666666667 * Q(1)-0.1333333333333333 * Q(2)+0.0166666666666667 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)
      W5 = +0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)-0.0333333333333333 * Q(-1)+0.3666666666666666 * Q(0)+0.9500000000000000 * Q(1)-0.3833333333333334 * Q(2)+0.1166666666666667 * Q(3)-0.0166666666666667 * Q(4)+0.0000000000000000 * Q(5)
      W6 = +0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.1666666666666667 * Q(0)+1.4500000000000000 * Q(1)-1.0500000000000000 * Q(2)+0.6166666666666667 * Q(3)-0.2166666666666667 * Q(4)+0.0333333333333333 * Q(5)
      IF ( LinearWeightsOnly ) THEN
         WP = gamma1*W1 + gamma2*W2 + gamma3*W3 + gamma4*W4 + gamma5*W5 + gamma6*W6
      ELSE
         WP = omega1*W1 + omega2*W2 + omega3*W3 + omega4*W4 + omega5*W5 + omega6*W6
      END IF

!-------------------------------------------------------------------------------!
   END SUBROUTINE WENO11_FirstSweep
!===============================================================================!
!
!
!
!
!===============================================================================!
   SUBROUTINE WENO13_FirstSweep(Q,WM,WP)
!-------------------------------------------------------------------------------!
      USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
      USE MOD_FiniteVolume2D_vars,ONLY: nGPs
      USE MOD_FiniteVolume2D_vars,ONLY: WENOEPS, WENOEXP
      USE MOD_FiniteVolume2D_vars,ONLY: LinearWeightsOnly
!-------------------------------------------------------------------------------!
      IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
      REAL,INTENT(IN)  :: Q(-nGhosts:nGhosts)
      REAL,INTENT(OUT) :: WM
      REAL,INTENT(OUT) :: WP
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
      REAL             :: alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7
      REAL             :: beta1,  beta2,  beta3,  beta4,  beta5,  beta6 , beta7
      REAL             :: gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7
      REAL             :: omega1, omega2, omega3, omega4, omega5, omega6, omega7
      REAL             :: W1, W2, W3, W4, W5, W6, W7
!-------------------------------------------------------------------------------!

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!

      beta1 = +(21.01414174816952496*Q(-6)*Q(-6))+(-139.20628098945806528*Q(-6)*Q(-5))+(387.72963641875099938*Q(-6)*Q(-4))+(-582.04506049249107491*Q(-6)*Q(-3))+(497.51201489097320518*Q(-6)*Q(-2))+(-230.02750618787078452*Q(-6)*Q(-1))+(45.02305461192614189*Q(-6)*Q(0))+(-139.20628098945806528*Q(-5)*Q(-6))+(923.49471636002886044*Q(-5)*Q(-5))+(-2576.73012919372285978*Q(-5)*Q(-4))+(3876.40142005371171763*Q(-5)*Q(-3))+(-3322.10024328102463187*Q(-5)*Q(-2))+(1540.88084731240974179*Q(-5)*Q(-1))+(-302.74033026194484819*Q(-5)*Q(0))+(387.72963641875099938*Q(-4)*Q(-6))+(-2576.73012919372285978*Q(-4)*Q(-5))+(7205.30018037518038909*Q(-4)*Q(-4))+(-10869.10913049142254749*Q(-4)*Q(-3))+(9346.55924535533995368*Q(-4)*Q(-2))+(-4353.46899328102426807*Q(-4)*Q(-1))+(859.71919081689918585*Q(-4)*Q(0))+(-582.04506049249107491*Q(-3)*Q(-6))+(3876.40142005371171763*Q(-3)*Q(-5))+(-10869.10913049142254749*Q(-3)*Q(-4))+(16453.17591223077397444*Q(-3)*Q(-3))+(-14212.00727863957035879*Q(-3)*Q(-2))+(6657.85327190556381538*Q(-3)*Q(-1))+(-1324.26913456656507151*Q(-3)*Q(0))+(497.51201489097320518*Q(-2)*Q(-6))+(-3322.10024328102463187*Q(-2)*Q(-5))+(9346.55924535533995368*Q(-2)*Q(-4))+(-14212.00727863957035879*Q(-2)*Q(-3))+(12350.33143037517947960*Q(-2)*Q(-2))+(-5832.94887919372285978*Q(-2)*Q(-1))+(1172.65371049282498461*Q(-2)*Q(0))+(-230.02750618787078452*Q(-1)*Q(-6))+(1540.88084731240974179*Q(-1)*Q(-5))+(-4353.46899328102426807*Q(-1)*Q(-4))+(6657.85327190556381538*Q(-1)*Q(-3))+(-5832.94887919372285978*Q(-1)*Q(-2))+(2787.97471636002865125*Q(-1)*Q(-1))+(-570.26345691538404026*Q(-1)*Q(0))+(45.02305461192614189*Q(0)*Q(-6))+(-302.74033026194484819*Q(0)*Q(-5))+(859.71919081689918585*Q(0)*Q(-4))+(-1324.26913456656507151*Q(0)*Q(-3))+(1172.65371049282498461*Q(0)*Q(-2))+(-570.26345691538404026*Q(0)*Q(-1))+(119.87696582224360498*Q(0)*Q(0))
      beta2 = +(4.29972816792261270*Q(-5)*Q(-5))+(-28.36961046476671555*Q(-5)*Q(-4))+(78.33089027677569050*Q(-5)*Q(-3))+(-115.76103271471326650*Q(-5)*Q(-2))+(96.43524019961519400*Q(-5)*Q(-1))+(-42.82792671256213168*Q(-5)*Q(0))+(7.89271124772860855*Q(-5)*Q(1))+(-28.36961046476671555*Q(-4)*Q(-5))+(187.89196173039923110*Q(-4)*Q(-4))+(-521.01977039742666875*Q(-4)*Q(-3))+(773.66384289321786127*Q(-4)*Q(-2))+(-647.80550948472819073*Q(-4)*Q(-1))+(289.20642601611353939*Q(-4)*Q(0))+(-53.56734029280904252*Q(-4)*Q(1))+(78.33089027677569050*Q(-3)*Q(-5))+(-521.01977039742666875*Q(-3)*Q(-4))+(1452.34531926406930324*Q(-3)*Q(-3))+(-2169.83328172599021855*Q(-3)*Q(-2))+(1829.33846989237622438*Q(-3)*Q(-1))+(-822.61152800324680356*Q(-3)*Q(0))+(153.44990069344237327*Q(-3)*Q(1))+(-115.76103271471326650*Q(-2)*Q(-5))+(773.66384289321786127*Q(-2)*Q(-4))+(-2169.83328172599021855*Q(-2)*Q(-3))+(3266.81402951472409768*Q(-2)*Q(-2))+(-2779.62803481240962356*Q(-2)*Q(-1))+(1262.72742314013157738*Q(-2)*Q(0))+(-237.98294629496018615*Q(-2)*Q(1))+(96.43524019961519400*Q(-1)*Q(-5))+(-647.80550948472819073*Q(-1)*Q(-4))+(1829.33846989237622438*Q(-1)*Q(-3))+(-2779.62803481240962356*Q(-1)*Q(-2))+(2394.05596741221734192*Q(-1)*Q(-1))+(-1103.66560373076003998*Q(-1)*Q(0))+(211.26947052368927871*Q(-1)*Q(1))+(-42.82792671256213168*Q(0)*Q(-5))+(289.20642601611353939*Q(0)*Q(-4))+(-822.61152800324680356*Q(0)*Q(-3))+(1262.72742314013157738*Q(0)*Q(-2))+(-1103.66560373076003998*Q(0)*Q(-1))+(519.24714691558438062*Q(0)*Q(0))+(-102.07593762526053638*Q(0)*Q(1))+(7.89271124772860855*Q(1)*Q(-5))+(-53.56734029280904252*Q(1)*Q(-4))+(153.44990069344237327*Q(1)*Q(-3))+(-237.98294629496018615*Q(1)*Q(-2))+(211.26947052368927871*Q(1)*Q(-1))+(-102.07593762526053638*Q(1)*Q(0))+(21.01414174816952496*Q(1)*Q(1))
      beta3 = +(1.40409545187322959*Q(-4)*Q(-4))+(-9.00175938451980073*Q(-4)*Q(-3))+(23.83364876443001279*Q(-4)*Q(-2))+(-33.25586296162685329*Q(-4)*Q(-1))+(25.70915995270161858*Q(-4)*Q(0))+(-10.41776853354978272*Q(-4)*Q(1))+(1.72848671069157178*Q(-4)*Q(2))+(-9.00175938451980073*Q(-3)*Q(-4))+(58.62804969336219330*Q(-3)*Q(-3))+(-157.57063845298219462*Q(-3)*Q(-2))+(222.91746943642777978*Q(-3)*Q(-1))+(-174.45649328102453524*Q(-3)*Q(0))+(71.44677323833573723*Q(-3)*Q(1))+(-11.96340124959916551*Q(-3)*Q(2))+(23.83364876443001279*Q(-2)*Q(-4))+(-157.57063845298219462*Q(-2)*Q(-3))+(430.70874518999517022*Q(-2)*Q(-2))+(-619.92548851611354621*Q(-2)*Q(-1))+(493.06850461459833923*Q(-2)*Q(0))+(-204.84422476250600198*Q(-2)*Q(1))+(34.72945316257816017*Q(-2)*Q(2))+(-33.25586296162685329*Q(-1)*Q(-4))+(222.91746943642777978*Q(-1)*Q(-3))+(-619.92548851611354621*Q(-1)*Q(-2))+(910.75615914435354625*Q(-1)*Q(-1))+(-739.84832802228640958*Q(-1)*Q(0))+(313.41129659692160203*Q(-1)*Q(1))+(-54.05524567767623267*Q(-1)*Q(2))+(25.70915995270161858*Q(0)*Q(-4))+(-174.45649328102453524*Q(0)*Q(-3))+(493.06850461459833923*Q(0)*Q(-2))+(-739.84832802228640958*Q(0)*Q(-1))+(616.65434704184701786*Q(0)*Q(0))+(-268.59355511964889729*Q(0)*Q(1))+(47.46636481381273143*Q(0)*Q(2))+(-10.41776853354978272*Q(1)*Q(-4))+(71.44677323833573723*Q(1)*Q(-3))+(-204.84422476250600198*Q(1)*Q(-2))+(313.41129659692160203*Q(1)*Q(-1))+(-268.59355511964889729*Q(1)*Q(0))+(121.20286450817701507*Q(1)*Q(1))+(-22.20538592772967945*Q(1)*Q(2))+(1.72848671069157178*Q(2)*Q(-4))+(-11.96340124959916551*Q(2)*Q(-3))+(34.72945316257816017*Q(2)*Q(-2))+(-54.05524567767623267*Q(2)*Q(-1))+(47.46636481381273143*Q(2)*Q(0))+(-22.20538592772967945*Q(2)*Q(1))+(4.29972816792261270*Q(2)*Q(2))
      beta4 = +(1.40409545187322959*Q(-3)*Q(-3))+(-8.10018145242103493*Q(-3)*Q(-2))+(19.06823595578804031*Q(-3)*Q(-1))+(-23.43418086286141744*Q(-3)*Q(0))+(15.88747785393618805*Q(-3)*Q(1))+(-5.65235572490780847*Q(-3)*Q(2))+(0.82690877859280643*Q(-3)*Q(3))+(-8.10018145242103493*Q(-2)*Q(-3))+(48.90159136002885987*Q(-2)*Q(-2))+(-119.38481669372293936*Q(-2)*Q(-1))+(151.00859597963764713*Q(-2)*Q(0))+(-104.77055578102452671*Q(-2)*Q(1))+(37.99772231240981313*Q(-2)*Q(2))+(-5.65235572490780847*Q(-2)*Q(3))+(19.06823595578804031*Q(-1)*Q(-3))+(-119.38481669372293936*Q(-1)*Q(-2))+(302.86268037518038909*Q(-1)*Q(-1))+(-396.08945456549622577*Q(-1)*Q(0))+(282.42643285533910102*Q(-1)*Q(1))+(-104.77055578102452671*Q(-1)*Q(2))+(15.88747785393618805*Q(-1)*Q(3))+(-23.43418086286141744*Q(0)*Q(-3))+(151.00859597963764713*Q(0)*Q(-2))+(-396.08945456549622577*Q(0)*Q(-1))+(537.03007889744003478*Q(0)*Q(0))+(-396.08945456549622577*Q(0)*Q(1))+(151.00859597963764713*Q(0)*Q(2))+(-23.43418086286141744*Q(0)*Q(3))+(15.88747785393618805*Q(1)*Q(-3))+(-104.77055578102452671*Q(1)*Q(-2))+(282.42643285533910102*Q(1)*Q(-1))+(-396.08945456549622577*Q(1)*Q(0))+(302.86268037518038909*Q(1)*Q(1))+(-119.38481669372293936*Q(1)*Q(2))+(19.06823595578804031*Q(1)*Q(3))+(-5.65235572490780847*Q(2)*Q(-3))+(37.99772231240981313*Q(2)*Q(-2))+(-104.77055578102452671*Q(2)*Q(-1))+(151.00859597963764713*Q(2)*Q(0))+(-119.38481669372293936*Q(2)*Q(1))+(48.90159136002885987*Q(2)*Q(2))+(-8.10018145242103493*Q(2)*Q(3))+(0.82690877859280643*Q(3)*Q(-3))+(-5.65235572490780847*Q(3)*Q(-2))+(15.88747785393618805*Q(3)*Q(-1))+(-23.43418086286141744*Q(3)*Q(0))+(19.06823595578804031*Q(3)*Q(1))+(-8.10018145242103493*Q(3)*Q(2))+(1.40409545187322959*Q(3)*Q(3))
      beta5 = +(4.29972816792261270*Q(-2)*Q(-2))+(-22.20538592772967945*Q(-2)*Q(-1))+(47.46636481381273143*Q(-2)*Q(0))+(-54.05524567767623267*Q(-2)*Q(1))+(34.72945316257816017*Q(-2)*Q(2))+(-11.96340124959916551*Q(-2)*Q(3))+(1.72848671069157178*Q(-2)*Q(4))+(-22.20538592772967945*Q(-1)*Q(-2))+(121.20286450817701507*Q(-1)*Q(-1))+(-268.59355511964889729*Q(-1)*Q(0))+(313.41129659692160203*Q(-1)*Q(1))+(-204.84422476250600198*Q(-1)*Q(2))+(71.44677323833573723*Q(-1)*Q(3))+(-10.41776853354978272*Q(-1)*Q(4))+(47.46636481381273143*Q(0)*Q(-2))+(-268.59355511964889729*Q(0)*Q(-1))+(616.65434704184701786*Q(0)*Q(0))+(-739.84832802228640958*Q(0)*Q(1))+(493.06850461459833923*Q(0)*Q(2))+(-174.45649328102453524*Q(0)*Q(3))+(25.70915995270161858*Q(0)*Q(4))+(-54.05524567767623267*Q(1)*Q(-2))+(313.41129659692160203*Q(1)*Q(-1))+(-739.84832802228640958*Q(1)*Q(0))+(910.75615914435354625*Q(1)*Q(1))+(-619.92548851611354621*Q(1)*Q(2))+(222.91746943642777978*Q(1)*Q(3))+(-33.25586296162685329*Q(1)*Q(4))+(34.72945316257816017*Q(2)*Q(-2))+(-204.84422476250600198*Q(2)*Q(-1))+(493.06850461459833923*Q(2)*Q(0))+(-619.92548851611354621*Q(2)*Q(1))+(430.70874518999517022*Q(2)*Q(2))+(-157.57063845298219462*Q(2)*Q(3))+(23.83364876443001279*Q(2)*Q(4))+(-11.96340124959916551*Q(3)*Q(-2))+(71.44677323833573723*Q(3)*Q(-1))+(-174.45649328102453524*Q(3)*Q(0))+(222.91746943642777978*Q(3)*Q(1))+(-157.57063845298219462*Q(3)*Q(2))+(58.62804969336219330*Q(3)*Q(3))+(-9.00175938451980073*Q(3)*Q(4))+(1.72848671069157178*Q(4)*Q(-2))+(-10.41776853354978272*Q(4)*Q(-1))+(25.70915995270161858*Q(4)*Q(0))+(-33.25586296162685329*Q(4)*Q(1))+(23.83364876443001279*Q(4)*Q(2))+(-9.00175938451980073*Q(4)*Q(3))+(1.40409545187322959*Q(4)*Q(4))
      beta6 = +(21.01414174816952496*Q(-1)*Q(-1))+(-102.07593762526053638*Q(-1)*Q(0))+(211.26947052368927871*Q(-1)*Q(1))+(-237.98294629496018615*Q(-1)*Q(2))+(153.44990069344237327*Q(-1)*Q(3))+(-53.56734029280904252*Q(-1)*Q(4))+(7.89271124772860855*Q(-1)*Q(5))+(-102.07593762526053638*Q(0)*Q(-1))+(519.24714691558438062*Q(0)*Q(0))+(-1103.66560373076003998*Q(0)*Q(1))+(1262.72742314013157738*Q(0)*Q(2))+(-822.61152800324680356*Q(0)*Q(3))+(289.20642601611353939*Q(0)*Q(4))+(-42.82792671256213168*Q(0)*Q(5))+(211.26947052368927871*Q(1)*Q(-1))+(-1103.66560373076003998*Q(1)*Q(0))+(2394.05596741221734192*Q(1)*Q(1))+(-2779.62803481240962356*Q(1)*Q(2))+(1829.33846989237622438*Q(1)*Q(3))+(-647.80550948472819073*Q(1)*Q(4))+(96.43524019961519400*Q(1)*Q(5))+(-237.98294629496018615*Q(2)*Q(-1))+(1262.72742314013157738*Q(2)*Q(0))+(-2779.62803481240962356*Q(2)*Q(1))+(3266.81402951472409768*Q(2)*Q(2))+(-2169.83328172599021855*Q(2)*Q(3))+(773.66384289321786127*Q(2)*Q(4))+(-115.76103271471326650*Q(2)*Q(5))+(153.44990069344237327*Q(3)*Q(-1))+(-822.61152800324680356*Q(3)*Q(0))+(1829.33846989237622438*Q(3)*Q(1))+(-2169.83328172599021855*Q(3)*Q(2))+(1452.34531926406930324*Q(3)*Q(3))+(-521.01977039742666875*Q(3)*Q(4))+(78.33089027677569050*Q(3)*Q(5))+(-53.56734029280904252*Q(4)*Q(-1))+(289.20642601611353939*Q(4)*Q(0))+(-647.80550948472819073*Q(4)*Q(1))+(773.66384289321786127*Q(4)*Q(2))+(-521.01977039742666875*Q(4)*Q(3))+(187.89196173039923110*Q(4)*Q(4))+(-28.36961046476671555*Q(4)*Q(5))+(7.89271124772860855*Q(5)*Q(-1))+(-42.82792671256213168*Q(5)*Q(0))+(96.43524019961519400*Q(5)*Q(1))+(-115.76103271471326650*Q(5)*Q(2))+(78.33089027677569050*Q(5)*Q(3))+(-28.36961046476671555*Q(5)*Q(4))+(4.29972816792261270*Q(5)*Q(5))
      beta7 = +(119.87696582224360498*Q(0)*Q(0))+(-570.26345691538404026*Q(0)*Q(1))+(1172.65371049282498461*Q(0)*Q(2))+(-1324.26913456656507151*Q(0)*Q(3))+(859.71919081689918585*Q(0)*Q(4))+(-302.74033026194484819*Q(0)*Q(5))+(45.02305461192614189*Q(0)*Q(6))+(-570.26345691538404026*Q(1)*Q(0))+(2787.97471636002865125*Q(1)*Q(1))+(-5832.94887919372285978*Q(1)*Q(2))+(6657.85327190556381538*Q(1)*Q(3))+(-4353.46899328102426807*Q(1)*Q(4))+(1540.88084731240974179*Q(1)*Q(5))+(-230.02750618787078452*Q(1)*Q(6))+(1172.65371049282498461*Q(2)*Q(0))+(-5832.94887919372285978*Q(2)*Q(1))+(12350.33143037517947960*Q(2)*Q(2))+(-14212.00727863957035879*Q(2)*Q(3))+(9346.55924535533995368*Q(2)*Q(4))+(-3322.10024328102463187*Q(2)*Q(5))+(497.51201489097320518*Q(2)*Q(6))+(-1324.26913456656507151*Q(3)*Q(0))+(6657.85327190556381538*Q(3)*Q(1))+(-14212.00727863957035879*Q(3)*Q(2))+(16453.17591223077397444*Q(3)*Q(3))+(-10869.10913049142254749*Q(3)*Q(4))+(3876.40142005371171763*Q(3)*Q(5))+(-582.04506049249107491*Q(3)*Q(6))+(859.71919081689918585*Q(4)*Q(0))+(-4353.46899328102426807*Q(4)*Q(1))+(9346.55924535533995368*Q(4)*Q(2))+(-10869.10913049142254749*Q(4)*Q(3))+(7205.30018037518038909*Q(4)*Q(4))+(-2576.73012919372285978*Q(4)*Q(5))+(387.72963641875099938*Q(4)*Q(6))+(-302.74033026194484819*Q(5)*Q(0))+(1540.88084731240974179*Q(5)*Q(1))+(-3322.10024328102463187*Q(5)*Q(2))+(3876.40142005371171763*Q(5)*Q(3))+(-2576.73012919372285978*Q(5)*Q(4))+(923.49471636002886044*Q(5)*Q(5))+(-139.20628098945806528*Q(5)*Q(6))+(45.02305461192614189*Q(6)*Q(0))+(-230.02750618787078452*Q(6)*Q(1))+(497.51201489097320518*Q(6)*Q(2))+(-582.04506049249107491*Q(6)*Q(3))+(387.72963641875099938*Q(6)*Q(4))+(-139.20628098945806528*Q(6)*Q(5))+(21.01414174816952496*Q(6)*Q(6))



!------------------------------!
! WM: x_{i-1/2}                !
!------------------------------!

! Linear Weights
      gamma1 =0.0040792540792541
      gamma2 =0.0734265734265734
      gamma3 =0.3059440559440559
      gamma4 =0.4079254079254079
      gamma5 =0.1835664335664336
      gamma6 =0.0244755244755245
      gamma7 =0.0005827505827506

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
      alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP
      alpha5 = gamma5/(WENOEPS + beta5)**WENOEXP
      alpha6 = gamma6/(WENOEPS + beta6)**WENOEXP
      alpha7 = gamma7/(WENOEPS + beta7)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7)
      omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7)
      omega5 = alpha5/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7)
      omega6 = alpha6/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7)
      omega7 = alpha7/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7)

      W1 = -0.0238095238095238 * Q(-6)+0.1761904761904762 * Q(-5)-0.5738095238095238 * Q(-4)+1.0928571428571427 * Q(-3)-1.4071428571428573 * Q(-2)+1.5928571428571427 * Q(-1)+0.1428571428571428 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)+0.0000000000000000 * Q(6)
      W2 = +0.0000000000000000 * Q(-6)+0.0095238095238095 * Q(-5)-0.0738095238095238 * Q(-4)+0.2595238095238095 * Q(-3)-0.5738095238095238 * Q(-2)+1.0928571428571427 * Q(-1)+0.3095238095238095 * Q(0)-0.0238095238095238 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)+0.0000000000000000 * Q(6)
      W3 = +0.0000000000000000 * Q(-6)+0.0000000000000000 * Q(-5)-0.0071428571428571 * Q(-4)+0.0595238095238095 * Q(-3)-0.2404761904761905 * Q(-2)+0.7595238095238095 * Q(-1)+0.5095238095238095 * Q(0)-0.0904761904761905 * Q(1)+0.0095238095238095 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)+0.0000000000000000 * Q(6)
      W4 = +0.0000000000000000 * Q(-6)+0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)+0.0095238095238095 * Q(-3)-0.0904761904761905 * Q(-2)+0.5095238095238095 * Q(-1)+0.7595238095238095 * Q(0)-0.2404761904761905 * Q(1)+0.0595238095238095 * Q(2)-0.0071428571428571 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)+0.0000000000000000 * Q(6)
      W5 = +0.0000000000000000 * Q(-6)+0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)-0.0238095238095238 * Q(-2)+0.3095238095238095 * Q(-1)+1.0928571428571427 * Q(0)-0.5738095238095238 * Q(1)+0.2595238095238095 * Q(2)-0.0738095238095238 * Q(3)+0.0095238095238095 * Q(4)+0.0000000000000000 * Q(5)+0.0000000000000000 * Q(6)
      W6 = +0.0000000000000000 * Q(-6)+0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.1428571428571428 * Q(-1)+1.5928571428571427 * Q(0)-1.4071428571428573 * Q(1)+1.0928571428571427 * Q(2)-0.5738095238095238 * Q(3)+0.1761904761904762 * Q(4)-0.0238095238095238 * Q(5)+0.0000000000000000 * Q(6)
      W7 = +0.0000000000000000 * Q(-6)+0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+2.5928571428571430 * Q(0)-4.4071428571428575 * Q(1)+6.0928571428571425 * Q(2)-5.5738095238095235 * Q(3)+3.1761904761904760 * Q(4)-1.0238095238095237 * Q(5)+0.1428571428571428 * Q(6)



      IF ( LinearWeightsOnly ) THEN
         WM = gamma1*W1 + gamma2*W2 + gamma3*W3 + gamma4*W4 + gamma5*W5 + gamma6*W6 + gamma7*W7
      ELSE
         WM = omega1*W1 + omega2*W2 + omega3*W3 + omega4*W4 + omega5*W5 + omega6*W6 + omega7*W7
      END IF


!------------------------------!
! WP: x_{i+1/2}                !
!------------------------------!

! Linear Weights
      gamma1 =0.0005827505827506
      gamma2 =0.0244755244755245
      gamma3 =0.1835664335664336
      gamma4 =0.4079254079254079
      gamma5 =0.3059440559440559
      gamma6 =0.0734265734265734
      gamma7 =0.0040792540792541

      alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
      alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
      alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP
      alpha4 = gamma4/(WENOEPS + beta4)**WENOEXP
      alpha5 = gamma5/(WENOEPS + beta5)**WENOEXP
      alpha6 = gamma6/(WENOEPS + beta6)**WENOEXP
      alpha7 = gamma7/(WENOEPS + beta7)**WENOEXP

! Nonlinear Weights
      omega1 = alpha1/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7)
      omega2 = alpha2/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7)
      omega3 = alpha3/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7)
      omega4 = alpha4/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7)
      omega5 = alpha5/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7)
      omega6 = alpha6/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7)
      omega7 = alpha7/(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7)


      W1 = +0.1428571428571428 * Q(-6)-1.0238095238095237 * Q(-5)+3.1761904761904760 * Q(-4)-5.5738095238095235 * Q(-3)+6.0928571428571425 * Q(-2)-4.4071428571428575 * Q(-1)+2.5928571428571430 * Q(0)+0.0000000000000000 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)+0.0000000000000000 * Q(6)
      W2 = +0.0000000000000000 * Q(-6)-0.0238095238095238 * Q(-5)+0.1761904761904762 * Q(-4)-0.5738095238095238 * Q(-3)+1.0928571428571427 * Q(-2)-1.4071428571428573 * Q(-1)+1.5928571428571427 * Q(0)+0.1428571428571428 * Q(1)+0.0000000000000000 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)+0.0000000000000000 * Q(6)
      W3 = +0.0000000000000000 * Q(-6)+0.0000000000000000 * Q(-5)+0.0095238095238095 * Q(-4)-0.0738095238095238 * Q(-3)+0.2595238095238095 * Q(-2)-0.5738095238095238 * Q(-1)+1.0928571428571427 * Q(0)+0.3095238095238095 * Q(1)-0.0238095238095238 * Q(2)+0.0000000000000000 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)+0.0000000000000000 * Q(6)
      W4 = +0.0000000000000000 * Q(-6)+0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)-0.0071428571428571 * Q(-3)+0.0595238095238095 * Q(-2)-0.2404761904761905 * Q(-1)+0.7595238095238095 * Q(0)+0.5095238095238095 * Q(1)-0.0904761904761905 * Q(2)+0.0095238095238095 * Q(3)+0.0000000000000000 * Q(4)+0.0000000000000000 * Q(5)+0.0000000000000000 * Q(6)
      W5 = +0.0000000000000000 * Q(-6)+0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0095238095238095 * Q(-2)-0.0904761904761905 * Q(-1)+0.5095238095238095 * Q(0)+0.7595238095238095 * Q(1)-0.2404761904761905 * Q(2)+0.0595238095238095 * Q(3)-0.0071428571428571 * Q(4)+0.0000000000000000 * Q(5)+0.0000000000000000 * Q(6)
      W6 = +0.0000000000000000 * Q(-6)+0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)-0.0238095238095238 * Q(-1)+0.3095238095238095 * Q(0)+1.0928571428571427 * Q(1)-0.5738095238095238 * Q(2)+0.2595238095238095 * Q(3)-0.0738095238095238 * Q(4)+0.0095238095238095 * Q(5)+0.0000000000000000 * Q(6)
      W7 = +0.0000000000000000 * Q(-6)+0.0000000000000000 * Q(-5)+0.0000000000000000 * Q(-4)+0.0000000000000000 * Q(-3)+0.0000000000000000 * Q(-2)+0.0000000000000000 * Q(-1)+0.1428571428571428 * Q(0)+1.5928571428571427 * Q(1)-1.4071428571428573 * Q(2)+1.0928571428571427 * Q(3)-0.5738095238095238 * Q(4)+0.1761904761904762 * Q(5)-0.0238095238095238 * Q(6)

      IF ( LinearWeightsOnly ) THEN
         WP = gamma1*W1 + gamma2*W2 + gamma3*W3 + gamma4*W4 + gamma5*W5 + gamma6*W6 + gamma7*W7
      ELSE
         WP = omega1*W1 + omega2*W2 + omega3*W3 + omega4*W4 + omega5*W5 + omega6*W6 + omega7*W7
      END IF

!-------------------------------------------------------------------------------!
   END SUBROUTINE WENO13_FirstSweep
END MODULE MOD_Reconstruction
!-------------------------------------------------------------------------------!
