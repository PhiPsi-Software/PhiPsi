 
      SUBROUTINE Solver_Gauss(LIN, LIMROW, LIMCOL, N, X, SINGUL) 
       
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
    


