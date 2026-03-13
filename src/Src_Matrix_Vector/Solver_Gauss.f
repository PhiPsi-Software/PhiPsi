!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
      SUBROUTINE Solver_Gauss(LIN, LIMROW, LIMCOL, N, X, SINGUL) 
C Gaussian elimination method for solving systems of linear equations
       
c*GAUSS*****************************************************************
c Subroutine to find solution of a linear system of N equations in N   *
c unknowns using Gaussian elimination, provided a unique solution      *
c exists.  The coefficients and constants of the linear system are     *
c stored in the matrix LIN, which has LIMROW rows and LIMCOL columns.  *
c If the system is singular, SINGUL is returned as true, and the       *
c solution X is undefined.  Local identifiers used are:                *                        
*     I,J,K  : subscripts                                              *
c     MULT   : multiplier used to eliminate an unknown                 *
c     ASBPIV : absolute value of pivot element                         *
c     PIVROW : row containing pivot element                            *
c     EPSIL  : a small positive real value ("almost zero")             *
c     TEMP   : used to interchange rows of matrix                      *
c                                                                      *
c Accepts: Two-dimensional array LIM, integers LIMROW, LIMCOL, and N   *
c Returns: One-dimensional array X and logical value SINGUL            *
c***********************************************************************
      use Global_Float_Type
      implicit none
      integer LIMROW, LIMCOL
      real(kind=FT) LIN(LIMROW, LIMCOL), X(LIMROW), TEMP, MULT, EPSIL
      PARAMETER (EPSIL = 1D-15)
      INTEGER N, PIVROW
      real(kind=FT) ABSPIV
      LOGICAL SINGUL
      integer I,J,K
      
      SINGUL = .FALSE.
      DO 50 I = 1, N
         ABSPIV = ABS(LIN(I,I))
         PIVROW = I
         DO 10 K = I + 1, N
            IF (ABS(LIN(K,I)) .GT. ABSPIV) THEN
               ABSPIV = ABS(LIN(K,I))
               PIVROW = K
            END IF
10       CONTINUE
 
         IF (ABSPIV .LT. EPSIL) THEN
            SINGUL = .TRUE.
            RETURN
         END IF
         
         IF (PIVROW .NE. I) THEN
            DO 20 J = 1, N + 1
               TEMP = LIN(I,J)
               LIN(I,J) = LIN(PIVROW,J)
               LIN(PIVROW,J) = TEMP
20          CONTINUE
         END IF

 

         DO 40 J = I + 1, N
            MULT = -LIN(J,I) / LIN(I,I)
            DO 30 K = I, N + 1
               LIN(J,K) = LIN(J,K) +  MULT * LIN(I,K)
30          CONTINUE
40       CONTINUE
50    CONTINUE

      X(N) = LIN(N, N + 1) / LIN(N,N)
      DO 70 J = N - 1, 1, -1
         X(J) = LIN(J, N + 1)
         DO 60 K = J + 1, N
            X(J) = X(J) - LIN(J,K) * X(K)
60       CONTINUE
         X(J) = X(J) / LIN(J,J)
70    CONTINUE

      return
      END SUBROUTINE Solver_Gauss
    


