!-----------------------------------------------------------
! Brief: Solve a linear system via Gaussian elimination with pivoting.
!
! Parameters:
!   In/Out: LIN - augmented coefficient matrix (LIMROW x LIMCOL)
!   Input:  LIMROW, LIMCOL - matrix dimensions
!           N               - number of equations
!   Output: X               - solution vector of length N
!           SINGUL          - true if the system is (near-)singular
!
! Notes:   Uses partial pivoting by largest absolute pivot. The
!          augmented matrix LIN is overwritten during elimination.
!-----------------------------------------------------------

      subroutine Solver_Gauss(LIN, LIMROW, LIMCOL, N, X, SINGUL) 
! Gaussian elimination method for solving systems of linear equations
       
!*GAUSS*****************************************************************
! Subroutine to find solution of a linear system of N equations in N   *
! unknowns using Gaussian elimination, provided a unique solution      *
! exists.  The coefficients and constants of the linear system are     *
! stored in the matrix LIN, which has LIMROW rows and LIMCOL columns.  *
! If the system is singular, SINGUL is returned as true, and the       *
! solution X is undefined.  Local identifiers used are:                *                        
!     I,J,K  : subscripts                                              *
!     MULT   : multiplier used to eliminate an unknown                 *
!     ASBPIV : absolute value of pivot element                         *
!     PIVROW : row containing pivot element                            *
!     EPSIL  : a small positive real value ("almost zero")             *
!     TEMP   : used to interchange rows of matrix                      *
!                                                                      *
! Accepts: Two-dimensional array LIM, integers LIMROW, LIMCOL, and N   *
! Returns: One-dimensional array X and logical value SINGUL            *
!***********************************************************************
      use Global_Float_Type
      implicit none
      integer LIMROW, LIMCOL
      real(kind=FT) LIN(LIMROW, LIMCOL), X(LIMROW), TEMP, MULT, EPSIL
      parameter (EPSIL = 1D-15)
      integer N, PIVROW
      real(kind=FT) :: ABSPIV
      logical :: SINGUL
      integer I,J,K
      
      SINGUL = .false.
      do 50 I = 1, N
         ABSPIV = ABS(LIN(I,I))
         PIVROW = I
         do 10 K = I + 1, N
            if (ABS(LIN(K,I)) .gt. ABSPIV) then
               ABSPIV = ABS(LIN(K,I))
               PIVROW = K
            end if
10       continue
 
         if (ABSPIV .lt. EPSIL) then
            SINGUL = .true.
            return
         end if
         
         if (PIVROW .ne. I) then
            do 20 J = 1, N + 1
               TEMP = LIN(I,J)
               LIN(I,J) = LIN(PIVROW,J)
               LIN(PIVROW,J) = TEMP
20          continue
         end if

 

         do 40 J = I + 1, N
            MULT = -LIN(J,I) / LIN(I,I)
            do 30 K = I, N + 1
               LIN(J,K) = LIN(J,K) +  MULT * LIN(I,K)
30          continue
40       continue
50    continue

      X(N) = LIN(N, N + 1) / LIN(N,N)
      do 70 J = N - 1, 1, -1
         X(J) = LIN(J, N + 1)
         do 60 K = J + 1, N
            X(J) = X(J) - LIN(J,K) * X(K)
60       continue
         X(J) = X(J) / LIN(J,J)
70    continue

      return
      END subroutine Solver_Gauss
    


