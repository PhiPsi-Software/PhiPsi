 
      SUBROUTINE Vector_Sort_Dou_with_Index_v_HSL_KB07(COUNT,N,INDEX)




      INTEGER N
      DOUBLE PRECISION COUNT(*)
      INTEGER INDEX(*)
      DOUBLE PRECISION AV,X
      INTEGER I,IF,IFK,IFKA,INT,INTEST,IP,IS,IS1,IY,J,K,K1,LA,LNGTH,M,
     +        MLOOP
      INTEGER MARK(100)
      DO 10 I = 1,N
        INDEX(I) = I
   10 CONTINUE
      IF (N.EQ.1) GO TO 200
      IF (N.GE.1) GO TO 30
      WRITE (6,FMT=20)

   20 FORMAT (/,/,/,20X,' ***KB07AD***NO NUMBERS TO BE SORTED ** ',
     + 'RETURN TO CALLING PROGRAM')

      GO TO 200
   30 M = 12
      LA = 2
      IS = 1
      IF = N
      DO 190 MLOOP = 1,N
        IFKA = IF - IS
        IF ((IFKA+1).GT.M) GO TO 70
        IS1 = IS + 1
        DO 60 J = IS1,IF
          I = J
   40     IF (COUNT(I-1).LT.COUNT(I)) GO TO 60
          IF (COUNT(I-1).GT.COUNT(I)) GO TO 50
          IF (INDEX(I-1).LT.INDEX(I)) GO TO 60
   50     AV = COUNT(I-1)
          COUNT(I-1) = COUNT(I)
          COUNT(I) = AV
          INT = INDEX(I-1)
          INDEX(I-1) = INDEX(I)
          INDEX(I) = INT
          I = I - 1
          IF (I.GT.IS) GO TO 40
   60   CONTINUE
        LA = LA - 2
        GO TO 170
   70   IY = (IS+IF)/2
        X = COUNT(IY)
        INTEST = INDEX(IY)
        COUNT(IY) = COUNT(IF)
        INDEX(IY) = INDEX(IF)
        K = 1
        IFK = IF
        DO 110 I = IS,IF
          IF (X.GT.COUNT(I)) GO TO 110
          IF (X.LT.COUNT(I)) GO TO 80
          IF (INTEST.GT.INDEX(I)) GO TO 110
   80     IF (I.GE.IFK) GO TO 120
          COUNT(IFK) = COUNT(I)
          INDEX(IFK) = INDEX(I)
          K1 = K
          DO 100 K = K1,IFKA
            IFK = IF - K
            IF (COUNT(IFK).GT.X) GO TO 100
            IF (COUNT(IFK).LT.X) GO TO 90
            IF (INTEST.LE.INDEX(IFK)) GO TO 100
   90       IF (I.GE.IFK) GO TO 130
            COUNT(I) = COUNT(IFK)
            INDEX(I) = INDEX(IFK)
            GO TO 110

  100     CONTINUE
          GO TO 120

  110   CONTINUE
  120   COUNT(IFK) = X
        INDEX(IFK) = INTEST
        IP = IFK
        GO TO 140

  130   COUNT(I) = X
        INDEX(I) = INTEST
        IP = I
  140   IF ((IP-IS).GT. (IF-IP)) GO TO 150
        MARK(LA) = IF
        MARK(LA-1) = IP + 1
        IF = IP - 1
        GO TO 160

  150   MARK(LA) = IP - 1
        MARK(LA-1) = IS
        IS = IP + 1
  160   LNGTH = IF - IS
        IF (LNGTH.LE.0) GO TO 180
        LA = LA + 2
        GO TO 190

  170   IF (LA.LE.0) GO TO 200
  180   IF = MARK(LA)
        IS = MARK(LA-1)
  190 CONTINUE
  200 RETURN

      END
