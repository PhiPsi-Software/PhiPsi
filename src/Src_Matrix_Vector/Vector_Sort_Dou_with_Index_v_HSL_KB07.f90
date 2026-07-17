
      subroutine Vector_Sort_Dou_with_Index_v_HSL_KB07(COUNT,N,INDEX)




      integer :: N
      double precision COUNT(*)
      integer :: INDEX(*)
      double precision AV,X
integer I,IF,IFK,IFKA,INT,INTEST,IP,IS,IS1,IY,J,K,K1,LA,LNGTH,M, MLOOP
      integer :: MARK(100)
      do 10 I = 1,N
        INDEX(I) = I
   10 continue
      if (N.eq.1) GO TO 200
      if (N.ge.1) GO TO 30
      write (6,FMT=20)

20 format (/,/,/,20X,' ***KB07AD***NO NUMBERS TO BE SORTED ** ', 'RETURN TO CALLING PROGRAM')

      GO TO 200
   30 M = 12
      LA = 2
      IS = 1
      IF = N
      do 190 MLOOP = 1,N
        IFKA = IF - IS
        if ((IFKA+1).gt.M) GO TO 70
        IS1 = IS + 1
        do 60 J = IS1,IF
          I = J
   40     if (COUNT(I-1).lt.COUNT(I)) GO TO 60
          if (COUNT(I-1).gt.COUNT(I)) GO TO 50
          if (INDEX(I-1).lt.INDEX(I)) GO TO 60
   50     AV = COUNT(I-1)
          COUNT(I-1) = COUNT(I)
          COUNT(I) = AV
          INT = INDEX(I-1)
          INDEX(I-1) = INDEX(I)
          INDEX(I) = INT
          I = I - 1
          if (I.gt.IS) GO TO 40
   60   continue
        LA = LA - 2
        GO TO 170
   70   IY = (IS+IF)/2
        X = COUNT(IY)
        INTEST = INDEX(IY)
        COUNT(IY) = COUNT(IF)
        INDEX(IY) = INDEX(IF)
        K = 1
        IFK = IF
        do 110 I = IS,IF
          if (X.gt.COUNT(I)) GO TO 110
          if (X.lt.COUNT(I)) GO TO 80
          if (INTEST.gt.INDEX(I)) GO TO 110
   80     if (I.ge.IFK) GO TO 120
          COUNT(IFK) = COUNT(I)
          INDEX(IFK) = INDEX(I)
          K1 = K
          do 100 K = K1,IFKA
            IFK = IF - K
            if (COUNT(IFK).gt.X) GO TO 100
            if (COUNT(IFK).lt.X) GO TO 90
            if (INTEST.le.INDEX(IFK)) GO TO 100
   90       if (I.ge.IFK) GO TO 130
            COUNT(I) = COUNT(IFK)
            INDEX(I) = INDEX(IFK)
            GO TO 110

  100     continue
          GO TO 120

  110   continue
  120   COUNT(IFK) = X
        INDEX(IFK) = INTEST
        IP = IFK
        GO TO 140

  130   COUNT(I) = X
        INDEX(I) = INTEST
        IP = I
  140   if ((IP-IS).gt. (IF-IP)) GO TO 150
        MARK(LA) = IF
        MARK(LA-1) = IP + 1
        IF = IP - 1
        GO TO 160

  150   MARK(LA) = IP - 1
        MARK(LA-1) = IS
        IS = IP + 1
  160   LNGTH = IF - IS
        if (LNGTH.le.0) GO TO 180
        LA = LA + 2
        GO TO 190

  170   if (LA.le.0) GO TO 200
  180   IF = MARK(LA)
        IS = MARK(LA-1)
  190 continue
  200 return

      END
