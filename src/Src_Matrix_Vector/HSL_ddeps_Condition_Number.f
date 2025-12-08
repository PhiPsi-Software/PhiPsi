 
      SUBROUTINE MC59AD(ICNTL,NC,NR,NE,IRN,LJCN,JCN,LA,A,LIP,IP,
     &                  LIW,IW,INFO)
      INTEGER LA,LIP,LIW,LJCN,NC,NE,NR
      DOUBLE PRECISION A(LA)
      INTEGER ICNTL(10),IP(LIP),INFO(10),IRN(NE),IW(LIW),JCN(LJCN)
      INTEGER I,ICNTL1,ICNTL2,ICNTL3,ICNTL6,LAA
      INTEGER IDUP,IOUT,IUP,JOUT,LP,MP,KNE,PART
      LOGICAL LCHECK
      EXTERNAL MC59BD,MC59CD,MC59DD,MC59ED,MC59FD
      INTRINSIC MAX

      DO 10 I = 1,10
         INFO(I) = 0
   10 CONTINUE

      ICNTL1 = ICNTL(1)
      ICNTL2 = ICNTL(2)
      ICNTL3 = ICNTL(3)
      ICNTL6 = ICNTL(6)
      LCHECK = (ICNTL1.EQ.0)
      LP = ICNTL(4)
      MP = ICNTL(5)

      IF (ICNTL2.GT.2 .OR. ICNTL2.LT.0) THEN
         INFO(1) = -1
         INFO(2) = ICNTL2
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9010) ICNTL2
         END IF
         GO TO 70
      END IF

      IF (ICNTL6.GT.2 .OR. ICNTL6.LT.-2) THEN
         INFO(1) = -11
         INFO(2) = ICNTL6
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9150) ICNTL6
         END IF
         GO TO 70
      END IF

      IF (NC.LT.1) THEN
        INFO(1) = -2
        INFO(2) = NC
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9020) NC
        END IF
        GO TO 70
      END IF

      IF (NR.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9030) NR
        END IF
        GO TO 70
      END IF

      IF (ICNTL6.NE.0 .AND. NR.NE.NC) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9035) NC,NR
        END IF
        GO TO 70
      END IF

      IF (NE.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NE
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9040) NE
        END IF
        GO TO 70
      END IF

      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.1) THEN
        IF (LJCN.LT.NE) THEN
          INFO(1) = -5
          INFO(2) = NE
        END IF
      ELSE
        IF (LJCN.LT.1) THEN
          INFO(1) = -5
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-5) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9050) LJCN,INFO(2)
         END IF
         GO TO 70
      END IF

      IF (ICNTL3.EQ.0) THEN
        IF (LA.LT.NE) THEN
          INFO(1) = -6
          INFO(2) = NE
        END IF
      ELSE
        IF (LA.LT.1) THEN
          INFO(1) = -6
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-6) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9060) LA,INFO(2)
         END IF
         GO TO 70
      END IF

      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.2) THEN
        IF (LIP.LT.NC+1) THEN
          INFO(1) = -7
          INFO(2) = NC+1
        END IF
      ELSE IF (LIP.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -7
        INFO(2) = MAX(NR,NC)+1
      END IF
      IF (INFO(1).EQ.-7) THEN
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9065) LIP,INFO(2)
        END IF
        GO TO 70
      END IF

      IF (LIW.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -8
        INFO(2) = MAX(NR,NC)+1
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9070) LIW,INFO(2)
        END IF
        GO TO 70
      END IF

      LAA = NE
      IF (ICNTL3.NE.0) LAA = 1
      IOUT = 0
      JOUT = 0
      IDUP = 0
      IUP = 0

      PART = 0
      IF (ICNTL6.NE.0) PART = 1

      IF (ICNTL2.EQ.0) THEN

        CALL MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LAA,A,IP,IW,
     +              IOUT,JOUT,KNE)
        IF (KNE.EQ.0) GO TO 50

        IF (LCHECK) CALL MC59ED(NC,NR,NE,IRN,LIP,IP,LAA,A,IW,IDUP,
     &                          KNE,ICNTL6)

      ELSE IF (ICNTL2.EQ.1) THEN

        IF (ICNTL6.NE.0) PART = -1
        CALL MC59BD(LCHECK,PART,NR,NC,NE,JCN,IRN,LAA,A,IW,IP,
     +              JOUT,IOUT,KNE)
        IF (KNE.EQ.0) GO TO 50

        IF (LCHECK) CALL MC59ED(NR,NC,NE,JCN,NR+1,IW,LAA,A,IP,
     &                          IDUP,KNE,ICNTL6)

        CALL MC59CD(NC,NR,KNE,IRN,JCN,LAA,A,IP,IW)

      ELSE IF (ICNTL2.EQ.2) THEN
        IF (LCHECK) THEN
          CALL MC59FD(NC,NR,NE,IRN,NC+1,IP,LAA,A,LIW,IW,IDUP,
     +                IOUT,IUP,KNE,ICNTL6,INFO)
          IF (INFO(1).EQ.-9) GO TO 40
          IF (KNE.EQ.0) GO TO 50
        ELSE
           KNE = NE
        END IF

        CALL MC59DD(NC,KNE,IRN,IP,LAA,A)

      END IF

      INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(6) = KNE
      INFO(7) = IUP

      IF (IDUP.GT.0) INFO(1) = INFO(1) + 1
      IF (IOUT.GT.0) INFO(1) = INFO(1) + 2
      IF (JOUT.GT.0) INFO(1) = INFO(1) + 4
      IF (INFO(1).GT.0 .AND. MP.GT.0) THEN
        WRITE (MP,FMT=9080) INFO(1)
        IF (IOUT.GT.0) WRITE (MP,FMT=9090) IOUT
        IF (JOUT.GT.0) WRITE (MP,FMT=9110) JOUT
        IF (IDUP.GT.0) WRITE (MP,FMT=9100) IDUP
        IF (IUP.GT.0)  WRITE (MP,FMT=9130) IUP
      END IF
      GO TO 70

   40 INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(7) = IUP
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9140)
      END IF
      GO TO 70

   50 INFO(1) = -10
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(2) = IOUT + JOUT
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9120)
      END IF
   70 RETURN

 9000 FORMAT (/,' *** Error return from MC59AD *** INFO(1) = ',I3)
 9010 FORMAT (1X,'ICNTL(2) = ',I2,' is out of range')
 9020 FORMAT (1X,'NC = ',I6,' is out of range')
 9030 FORMAT (1X,'NR = ',I6,' is out of range')
 9035 FORMAT (1X,'Symmetric case. NC = ',I6,' but NR = ',I6)
 9040 FORMAT (1X,'NE = ',I10,' is out of range')
 9050 FORMAT (1X,'Increase LJCN from ',I10,' to at least ',I10)
 9060 FORMAT (1X,'Increase LA from ',I10,' to at least ',I10)
 9065 FORMAT (1X,'Increase LIP from ',I8,' to at least ',I10)
 9070 FORMAT (1X,'Increase LIW from ',I8,' to at least ',I10)
 9080 FORMAT (/,' *** Warning message from MC59AD *** INFO(1) = ',I3)
 9090 FORMAT (1X,I8,' entries in IRN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9100 FORMAT (1X,I8,' duplicate entries were supplied by the user')
 9110 FORMAT (1X,I8,' entries in JCN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9120 FORMAT (1X,'All entries out of range')
 9130 FORMAT (1X,I8,' of these entries were in the upper triangular ',
     +       /,'       part of matrix')
 9140 FORMAT (1X,'Entries in IP are not monotonic increasing')
 9150 FORMAT (1X,'ICNTL(6) = ',I2,' is out of range')
      END
      SUBROUTINE MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LA,A,IP,IW,IOUT,
     +                  JOUT,KNE)

      INTEGER LA,NC,NE,NR,IOUT,JOUT,KNE,PART
      LOGICAL LCHECK
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NC+1),JCN(NE)
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,JCE,JCEP,K,L,LOC

      DO 10 J = 1,NC + 1
        IW(J) = 0
   10 CONTINUE

      KNE = 0
      IOUT = 0
      JOUT = 0
      IF (LCHECK) THEN
        IF (LA.GT.1) THEN
          IF (PART.EQ.0) THEN
            DO 20 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                A(KNE) = A(K)
                IW(J) = IW(J) + 1
              END IF
   20       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 21 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   21       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 22 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   22       CONTINUE
          END IF
        ELSE
          IF (PART.EQ.0) THEN
            DO 25 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                IW(J) = IW(J) + 1
              END IF
   25       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 26 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   26       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 27 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   27       CONTINUE
          END IF
        END IF
        IF (KNE.EQ.0) GO TO 130

      ELSE

        KNE = NE
        IF (PART.EQ.0) THEN
          DO 30 K = 1,NE
            J = JCN(K)
            IW(J) = IW(J) + 1
   30     CONTINUE
        ELSE IF (PART.EQ.1) THEN
          DO 35 K = 1,NE
            I = IRN(K)
            J = JCN(K)
            IF (I.LT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   35     CONTINUE
        ELSE IF (PART.EQ.-1) THEN
          DO 36 K = 1,NE
            I = IRN(K)
            J = JCN(K)
            IF (I.GT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   36     CONTINUE
        END IF
      END IF



      IP(1) = 1
      DO 37 J = 2,NC + 1
        IP(J) = IW(J-1) + IP(J-1)
        IW(J-1) = IP(J-1)
   37 CONTINUE


      IF (LA.EQ.1) THEN
        DO 70 L = 1,NC
          DO 60 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            DO 40 J = 1,NE
              IF (JCE.EQ.L) GO TO 50
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
   40       CONTINUE
   50       JCN(K) = JCE
            IRN(K) = ICE
   60     CONTINUE
   70   CONTINUE
      ELSE

        DO 120 L = 1,NC
          DO 110 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            ACE = A(K)
            DO 90 J = 1,NE
              IF (JCE.EQ.L) GO TO 100
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
              ACEP = A(LOC)
              A(LOC) = ACE
              ACE = ACEP
   90       CONTINUE
  100       JCN(K) = JCE
            IRN(K) = ICE
            A(K) = ACE
  110     CONTINUE
  120   CONTINUE
      END IF

  130 CONTINUE

      RETURN
      END
      SUBROUTINE MC59CD(NC,NR,NE,IRN,JCN,LA,A,IP,IW)
      INTEGER LA,NC,NE,NR
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NR+1),JCN(NE)
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,J1,J2,K,L,LOC,LOCP


      DO 10 J = 1,NC
        IP(J) = 0
   10 CONTINUE

      IF (LA.GT.1) THEN

        DO 20 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
          IRN(K) = JCN(K)
   20   CONTINUE
        IP(NC+1) = NE + 1


        IP(1) = IP(1) + 1
        DO 30 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
   30   CONTINUE

        DO 50 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 40 J = J1,J2
            K = IRN(J)
            L = IP(K) - 1
            JCN(J) = L
            IRN(J) = I
            IP(K) = L
   40     CONTINUE
   50   CONTINUE
        IP(NC+1) = NE + 1
        DO 70 J = 1,NE
          LOC = JCN(J)
          IF (LOC.EQ.0) GO TO 70
          ICE = IRN(J)
          ACE = A(J)
          JCN(J) = 0
          DO 60 K = 1,NE
            LOCP = JCN(LOC)
            ICEP = IRN(LOC)
            ACEP = A(LOC)
            JCN(LOC) = 0
            IRN(LOC) = ICE
            A(LOC) = ACE
            IF (LOCP.EQ.0) GO TO 70
            ICE = ICEP
            ACE = ACEP
            LOC = LOCP
   60     CONTINUE
   70   CONTINUE
      ELSE



        DO 90 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
   90   CONTINUE
        IP(NC+1) = NE + 1


        IP(1) = IP(1) + 1
        DO 100 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
  100   CONTINUE

        DO 120 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 110 J = J1,J2
            K = JCN(J)
            L = IP(K) - 1
            IRN(L) = I
            IP(K) = L
  110     CONTINUE
  120   CONTINUE

      END IF

      RETURN
      END


      SUBROUTINE MC59DD(NC,NE,IRN,IP,LA,A)
      INTEGER LA,NC,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(NC)
      DOUBLE PRECISION ACE
      INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
      INTRINSIC ABS

      IF (LA.GT.1) THEN
        KMAX = NE
        DO 50 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 40
          KOR = KMAX
          DO 30 KDUMMY = KLO,KMAX
            ACE = A(KOR-1)
            ICE = IRN(KOR-1)
            DO 10 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 20
              IRN(K-1) = IK
              A(K-1) = A(K)
   10       CONTINUE
            K = KMAX + 1
   20       IRN(K-1) = ICE
            A(K-1) = ACE
            KOR = KOR - 1
   30     CONTINUE
   40     KMAX = KLO - 2
   50   CONTINUE
      ELSE

        KMAX = NE
        DO 150 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 140
          KOR = KMAX
          DO 130 KDUMMY = KLO,KMAX
            ICE = IRN(KOR-1)
            DO 110 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 120
              IRN(K-1) = IK
  110       CONTINUE
            K = KMAX + 1
  120       IRN(K-1) = ICE
            KOR = KOR - 1
  130     CONTINUE
  140     KMAX = KLO - 2
  150   CONTINUE
      END IF
      END

      SUBROUTINE MC59ED(NC,NR,NE,IRN,LIP,IP,LA,A,IW,IDUP,KNE,ICNTL6)

      INTEGER ICNTL6,IDUP,KNE,LIP,LA,NC,NR,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(NR)
      INTEGER I,J,K,KSTART,KSTOP,NZJ

      IDUP = 0
      KNE = 0
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE

      KSTART = IP(1)
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
              IDUP = IDUP + 1
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE

      ELSE

        DO 50 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF

      RETURN
      END

      SUBROUTINE MC59FD(NC,NR,NE,IRN,LIP,IP,LA,A,LIW,IW,IDUP,IOUT,
     +                  IUP,KNE,ICNTL6,INFO)

      INTEGER ICNTL6,IDUP,IOUT,IUP,KNE,LA,LIP,LIW,NC,NR,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(LIW),INFO(2)
      INTEGER I,J,K,KSTART,KSTOP,NZJ,LOWER

      IDUP = 0
      IOUT = 0
      IUP = 0
      KNE = 0
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE

      KSTART = IP(1)
      LOWER = 1
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.LT.J) IUP = IUP + 1
            ELSE IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
              IDUP = IDUP + 1
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE

      ELSE

        DO 50 J = 1,NC
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO  40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.GT.1) IUP = IUP + 1
            ELSE IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF

      RETURN
      END



      SUBROUTINE MA50AD(M,N,NE,LA,A,IRN,JCN,IQ,CNTL,ICNTL,IP,NP,JFIRST,
     +                  LENR,LASTR,NEXTR,IW,IFIRST,LENC,LASTC,NEXTC,
     +                  INFO,RINFO)



      INTEGER M,N,NE,LA
      DOUBLE PRECISION A(LA)
      DOUBLE PRECISION CNTL(10)
      INTEGER IRN(LA),JCN(LA),IQ(N)
      INTEGER ICNTL(20),IP(M),NP,JFIRST(M),LENR(M),LASTR(M),NEXTR(M),
     +        IW(M),IFIRST(N),LENC(N),LASTC(N),NEXTC(N),INFO(15)
      DOUBLE PRECISION RINFO(10)


      INTEGER IDAMAX
      EXTERNAL IDAMAX
      EXTERNAL MA50DD
      INTRINSIC ABS,MAX,MIN

      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0D0,ONE=1.0D0)

      DOUBLE PRECISION ALEN,AMULT,ANEW,ASW,AU,COST,CPIV
      INTEGER DISPC,DISPR,EYE,I,IDROP,IDUMMY,IEND,IFILL,IFIR,II,IJ,
     +        IJPOS,IOP,IPIV,IPOS,ISRCH,I1,I2,J,JBEG,JEND,JJ,JLAST,
     +        JMORE,JNEW,JPIV,JPOS,J1,J2,L,LC,LEN,LENPIV,LP,LR
      DOUBLE PRECISION MAXENT
      INTEGER MINC,MORD,MP,MSRCH,NC,NDROP,NEFACT,NEPR,NERED,NE1,NORD,
     +        NORD1,NR,NULLC,NULLI,NULLJ,NULLR,PIVBEG,PIVCOL,PIVEND,
     +        PIVOT
      DOUBLE PRECISION PIVR,PIVRAT,U


      LP = ICNTL(1)
      IF (ICNTL(3).LE.0) LP = 0
      MP = ICNTL(2)
      IF (ICNTL(3).LE.1) MP = 0
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = NE
      INFO(4) = NE
      INFO(5) = 0
      INFO(6) = 0
      INFO(7) = 0
      RINFO(1) = ZERO
      NP = 0

      IF (M.LT.1 .OR. N.LT.1) GO TO 690
      IF (NE.LT.1) GO TO 700
      IF (LA.LT.NE) THEN
         INFO(3) = NE
         GO TO 710
      END IF

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/2(A,I6),A,I8,A,I8/A,1P,4E10.2/A,7I4)')
     +     ' Entering MA50AD with M =',M,' N =',N,' NE =',NE,' LA =',LA,
     +     ' CNTL =',(CNTL(I),I=1,4),' ICNTL =',(ICNTL(I),I=1,7)
         IF (N.EQ.1 .OR. ICNTL(3).GT.3) THEN
            DO 10 J = 1,N - 1
               IF (IQ(J).LT.IQ(J+1)) WRITE (MP,
     +             '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J,
     +             (A(II),IRN(II),II=IQ(J),IQ(J+1)-1)
   10       CONTINUE
            IF (IQ(N).LE.NE) WRITE (MP,'(A,I5,(T13,3(1P,E12.4,I5)))')
     +          ' Column',N, (A(II),IRN(II),II=IQ(N),NE)
         ELSE
            IF (IQ(1).LT.IQ(2)) WRITE (MP,
     +          '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',1,
     +          (A(II),IRN(II),II=IQ(1),IQ(2)-1)
         END IF
         IF (ICNTL(7).EQ.2) THEN
            WRITE (MP,'(A,(T10,10(I7)))') ' IP = ',IP
            WRITE (MP,'(A,(T10,10(I7)))') ' IFIRST = ',IFIRST
         END IF
      END IF

      MINC = 1
      NERED = NE
      U = MIN(CNTL(2),ONE)
      U = MAX(U,ZERO)
      MSRCH = ICNTL(4)
      IF (MSRCH.EQ.0) MSRCH = N
      JLAST = N - ICNTL(6)
      IF (JLAST.LT.1 .OR. JLAST.GT.N) JLAST = N
      NULLI = 0
      NULLJ = 0
      MORD = 0
      NORD = 0
      NDROP = 0
      NEFACT = 0
      DO 20 I = 1,N - 1
         LENC(I) = IQ(I+1) - IQ(I)
   20 CONTINUE
      LENC(N) = NE + 1 - IQ(N)

      IF (CNTL(3).GT.ZERO) THEN
         NERED = 0
         DO 40 J = 1,N
            I = IQ(J)
            IQ(J) = NERED + 1
            DO 30 II = I,I + LENC(J) - 1
               IF (ABS(A(II)).GE.CNTL(3)) THEN
                  NERED = NERED + 1
                  A(NERED) = A(II)
                  IRN(NERED) = IRN(II)
               ELSE
                  INFO(6) = INFO(6) + 1
               END IF
   30       CONTINUE
            LENC(J) = NERED + 1 - IQ(J)
   40    CONTINUE
      END IF

      IF (ICNTL(7).EQ.2) THEN
         DO 50 I = 1,M
            NEXTR(I) = IP(I)
   50    CONTINUE
         IF (ICNTL(4).NE.1) GO TO 740
      END IF

      DISPR = NERED + 1
      DISPC = NERED + 1
      DO 60 I = 1,M
         IW(I) = 0
         LENR(I) = 0
         JFIRST(I) = 0
   60 CONTINUE
      DO 70 II = 1,NERED
         I = IRN(II)
         LENR(I) = LENR(I) + 1
   70 CONTINUE
      IP(1) = LENR(1) + 1
      DO 80 I = 2,M
         IP(I) = IP(I-1) + LENR(I)
   80 CONTINUE
      DO 100 J = 1,N
         I = IQ(J)
         DO 90 II = I,I + LENC(J) - 1
            I = IRN(II)
            IF (IW(I).EQ.J) GO TO 720
            IW(I) = J
            IPOS = IP(I) - 1
            JCN(IPOS) = J
            IP(I) = IPOS
   90    CONTINUE
  100 CONTINUE
      DO 110 I = 1,M
         IW(I) = 0
  110 CONTINUE

      IF (ICNTL(4).LE.0) THEN
         DO 120 I = 1,N
            IFIRST(I) = 0
  120    CONTINUE
         DO 130 I = M,1,-1
            NE1 = LENR(I)
            IF (NE1.GT.0) THEN
               IFIR = IFIRST(NE1)
               IFIRST(NE1) = I
               LASTR(I) = 0
               NEXTR(I) = IFIR
               IF (IFIR.GT.0) LASTR(IFIR) = I
            ELSE
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            END IF
  130    CONTINUE
      ELSE
         DO 140 I = M,1,-1
            NE1 = LENR(I)
            IF (NE1.EQ.0) THEN
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            END IF
  140    CONTINUE
      END IF
      DO 150 J = N,1,-1
         NE1 = LENC(J)
         IF (NE1.EQ.0) THEN
            IF (ICNTL(7).NE.2) THEN
               IF (J.LE.JLAST) THEN
                  NORD = NORD + 1
                  IQ(J) = -NORD
                  IF (NORD.EQ.JLAST) THEN
                     NORD = NORD + NULLJ
                     JLAST = N
                     NULLJ = 0
                  END IF
               ELSE
                  NULLJ = NULLJ + 1
                  IQ(J) = - (JLAST+NULLJ)
               END IF
               LASTC(J) = 0
               NEXTC(J) = 0
            END IF
         ELSE
            IFIR = JFIRST(NE1)
            JFIRST(NE1) = J
            NEXTC(J) = IFIR
            LASTC(J) = 0
            IF (IFIR.GT.0) LASTC(IFIR) = J
         END IF
  150 CONTINUE
      IF (INFO(6).EQ.0) THEN
         NULLC = NORD + NULLJ
         NULLR = NULLI
      END IF

      DO 630 PIVOT = 1,N
         IF (NERED.GE. (MIN(CNTL(1),ONE)*(N-NORD))*
     +       (M-MORD)) GO TO 640

         IF (ICNTL(7).EQ.2) THEN
            IPIV = 0
            J = IFIRST(PIVOT)
            IF (J.LT.1 .OR. J.GT.N) GO TO 730
            IF (IQ(J).LT.0) GO TO 730
            LEN = LENC(J)
            IF (LEN.LE.0) GO TO 320
            ALEN = LEN - 1
            I1 = IQ(J)
            I2 = I1 + LEN - 1
            II = IDAMAX(LEN,A(I1),1)
            MAXENT = ABS(A(I1+II-1))
            IF (MAXENT.LE.CNTL(4)) GO TO 320
            AU = MAX(MAXENT*U,CNTL(4))
            DO 160 II = I1,I2
               IF (ABS(A(II)).LT.AU) GO TO 160
               I = IRN(II)
               IF (IPIV.NE.0) THEN
                  IF (NEXTR(I).GE.NEXTR(IPIV)) GO TO 160
               END IF
               CPIV = ALEN*(LENR(I)-1)
               IJPOS = II
               IPIV = I
               JPIV = J
  160       CONTINUE
            GO TO 330
         END IF

         LEN = MINC
         DO 170 MINC = LEN,M - MORD
            IF (JFIRST(MINC).NE.0) GO TO 180
            IF (ICNTL(4).LE.0) THEN
               IF (IFIRST(MINC).NE.0) GO TO 180
            END IF
  170    CONTINUE

  180    CPIV = M
         CPIV = CPIV*N
         PIVRAT = ZERO
         ISRCH = 0
         DO 300 LEN = MINC,M - MORD
            ALEN = LEN - 1
            IF (CPIV.LE.ALEN**2 .AND. ICNTL(4).LE.0) GO TO 310
            IJ = JFIRST(LEN)
            DO 220 IDUMMY = 1,N
               IF (IJ.LE.0) GO TO 230
               J = IJ
               IJ = NEXTC(J)
               IF (J.GT.JLAST) GO TO 220
               MAXENT = ZERO
               I1 = IQ(J)
               I2 = I1 + LEN - 1
               II = IDAMAX(LEN,A(I1),1)
               MAXENT = ABS(A(I1+II-1))
               IF (MAXENT.LE.CNTL(4)) GO TO 320
               AU = MAX(MAXENT*U,CNTL(4))
               IF (ICNTL(7).EQ.1) THEN
                  DO 190 II = I1,I2
                     IF (IRN(II).EQ.J) GO TO 200
  190             CONTINUE
                  GO TO 220
  200             I1 = II
                  I2 = II
               END IF
               DO 210 II = I1,I2
                  IF (ABS(A(II)).LT.AU) GO TO 210
                  I = IRN(II)
                  COST = ALEN*(LENR(I)-1)
                  IF (COST.GT.CPIV) GO TO 210
                  PIVR = ABS(A(II))/MAXENT
                  IF (COST.EQ.CPIV) THEN
                     IF (PIVR.LE.PIVRAT) GO TO 210
                  END IF
                  CPIV = COST
                  IJPOS = II
                  IPIV = I
                  JPIV = J
                  IF (CPIV.LE.ALEN**2 .AND. ICNTL(4).LE.0) GO TO 330
                  PIVRAT = PIVR
  210          CONTINUE
               ISRCH = ISRCH + 1
               IF (ISRCH.GE.MSRCH) THEN
                  IF (PIVRAT.GT.ZERO) GO TO 330
               END IF
  220       CONTINUE
  230       IF (ICNTL(4).GT.0) GO TO 300
            IF (CPIV.LE.ALEN*(ALEN+1)) GO TO 310
            IF (LEN.GT.N-NORD) GO TO 300
            IJ = IFIRST(LEN)
            DO 290 IDUMMY = 1,M
               IF (IJ.EQ.0) GO TO 300
               I = IJ
               IJ = NEXTR(IJ)
               J1 = IP(I)
               J2 = J1 + LEN - 1
               IF (ICNTL(7).EQ.1) THEN
                  DO 240 JJ = J1,J2
                     IF (JCN(JJ).EQ.I) GO TO 250
  240             CONTINUE
                  GO TO 290
  250             J1 = JJ
                  J2 = JJ
               END IF
               DO 280 JJ = J1,J2
                  J = JCN(JJ)
                  IF (J.GT.JLAST) GO TO 280
                  COST = ALEN*(LENC(J)-1)
                  IF (COST.GE.CPIV) GO TO 280
                  I1 = IQ(J)
                  I2 = I1 + LENC(J) - 1
                  II = IDAMAX(LENC(J),A(I1),1)
                  MAXENT = ABS(A(I1+II-1))
                  DO 260 II = I1,I2 - 1
                     IF (IRN(II).EQ.I) GO TO 270
  260             CONTINUE
  270             JPOS = II
                  IF (MAXENT.LE.CNTL(4)) GO TO 320
                  IF (ABS(A(JPOS)).LT.MAXENT*U) GO TO 280
                  CPIV = COST
                  IPIV = I
                  JPIV = J
                  IJPOS = JPOS
                  PIVRAT = ABS(A(JPOS))/MAXENT
                  IF (CPIV.LE.ALEN*(ALEN+1)) GO TO 330
  280          CONTINUE

  290       CONTINUE
  300    CONTINUE
  310    IF (PIVRAT.GT.ZERO) GO TO 330
         INFO(1) = INFO(1) + 2
         IF (MP.GT.0) WRITE (MP,'(A/A)')
     +       ' Warning message from MA50AD: no suitable diagonal pivot',
     +       ' found, so switched to full matrix processing.'
         GO TO 640

  320    IPIV = 0
         JPIV = J

  330    NEFACT = NEFACT + LENC(JPIV)
         PIVBEG = IQ(JPIV)
         PIVEND = PIVBEG + LENC(JPIV) - 1
         NORD = NORD + 1
         NORD1 = NORD
         IF (NORD.EQ.JLAST) THEN
            NORD = NORD + NULLJ
            JLAST = N
            NULLJ = 0
         END IF
         IF (ICNTL(4).LE.0) THEN
            DO 340 II = PIVBEG,PIVEND
               I = IRN(II)
               LR = LASTR(I)
               NR = NEXTR(I)
               IF (NR.NE.0) LASTR(NR) = LR
               IF (LR.EQ.0) THEN
                  NE1 = LENR(I)
                  IFIRST(NE1) = NR
               ELSE
                  NEXTR(LR) = NR
               END IF
  340       CONTINUE
         END IF
         IF (IPIV.GT.0) THEN
            NEPR = LENR(IPIV) - 1
            NEFACT = NEFACT + NEPR
            RINFO(1) = RINFO(1) + CPIV*2 + LENR(IPIV)
            J1 = IP(IPIV)
            DO 350 JJ = J1,J1 + NEPR
               J = JCN(JJ)
               LC = LASTC(J)
               NC = NEXTC(J)
               IF (NC.NE.0) LASTC(NC) = LC
               IF (LC.EQ.0) THEN
                  NE1 = LENC(J)
                  JFIRST(NE1) = NC
               ELSE
                  NEXTC(LC) = NC
               END IF
  350       CONTINUE
            IF (PIVBEG.NE.IJPOS) THEN
               ASW = A(PIVBEG)
               A(PIVBEG) = A(IJPOS)
               A(IJPOS) = ASW
               IRN(IJPOS) = IRN(PIVBEG)
               IRN(PIVBEG) = IPIV
            END IF
         ELSE
            NEPR = 0
            NE1 = LENC(JPIV)
            IF (CNTL(3).GT.ZERO) NDROP = NDROP + NE1
            IF (NE1.GT.0) THEN
               LC = LASTC(JPIV)
               NC = NEXTC(JPIV)
               IF (NC.NE.0) LASTC(NC) = LC
               IF (LC.EQ.0) THEN
                  JFIRST(NE1) = NC
               ELSE
                  NEXTC(LC) = NC
               END IF
            END IF
         END IF
         DO 360 II = PIVBEG + 1,PIVEND
            I = IRN(II)
            IW(I) = II - PIVBEG
  360    CONTINUE
         LENPIV = PIVEND - PIVBEG
         DO 390 II = PIVBEG,PIVEND
            I = IRN(II)
            LENR(I) = LENR(I) - 1
            J1 = IP(I)
            J2 = J1 + LENR(I)
            DO 370 JJ = J1,J2 - 1
               IF (JCN(JJ).EQ.JPIV) GO TO 380
  370       CONTINUE
  380       JCN(JJ) = JCN(J2)
            JCN(J2) = 0
  390    CONTINUE

         DO 600 EYE = 1,NEPR
            J = JCN(IP(IPIV)+EYE-1)
            IDROP = 0
            JBEG = IQ(J)
            JEND = JBEG + LENC(J) - 1
            DO 400 II = JBEG,JEND - 1
               IF (IRN(II).EQ.IPIV) GO TO 410
  400       CONTINUE
  410       AMULT = -A(II)/A(IQ(JPIV))
            A(II) = A(JEND)
            IRN(II) = IRN(JEND)
            LENC(J) = LENC(J) - 1
            IRN(JEND) = 0
            JEND = JEND - 1
            IF (LENPIV.EQ.0) GO TO 600
            IOP = 0
            DO 420 II = JBEG,JEND
               I = IRN(II)
               IF (IW(I).GT.0) THEN
                  IOP = IOP + 1
                  PIVCOL = IQ(JPIV) + IW(I)
                  IW(I) = -IW(I)
                  A(II) = A(II) + AMULT*A(PIVCOL)
               END IF
  420       CONTINUE

            IF (CNTL(3).GT.ZERO) THEN
               JNEW = JBEG
               DO 450 II = JBEG,JEND
                  IF (ABS(A(II)).GE.CNTL(3)) THEN
                     A(JNEW) = A(II)
                     IRN(JNEW) = IRN(II)
                     JNEW = JNEW + 1
                  ELSE
                     I = IRN(II)
                     J1 = IP(I)
                     J2 = J1 + LENR(I) - 1
                     DO 430 JJ = J1,J2 - 1
                        IF (JCN(JJ).EQ.J) GO TO 440
  430                CONTINUE
  440                JCN(JJ) = JCN(J2)
                     JCN(J2) = 0
                     LENR(I) = LENR(I) - 1
                  END IF
  450          CONTINUE
               DO 460 II = JNEW,JEND
                  IRN(II) = 0
  460          CONTINUE
               IDROP = JEND + 1 - JNEW
               JEND = JNEW - 1
               LENC(J) = LENC(J) - IDROP
               NERED = NERED - IDROP
               INFO(6) = INFO(6) + IDROP
            END IF

            IFILL = LENPIV - IOP
            NERED = NERED + IFILL
            INFO(3) = MAX(INFO(3),NERED+LENC(J))

            IF (IFILL.EQ.0) THEN
               DO 470 II = PIVBEG + 1,PIVEND
                  I = IRN(II)
                  IW(I) = -IW(I)
  470          CONTINUE
               GO TO 600
            END IF

            DO 480 IPOS = JEND + 1,MIN(JEND+IFILL,DISPC-1)
               IF (IRN(IPOS).NE.0) GO TO 490
  480       CONTINUE
            IF (IPOS.EQ.JEND+IFILL+1) GO TO 540
            IF (JEND+IFILL+1.LE.LA+1) THEN
               DISPC = JEND + IFILL + 1
               GO TO 540
            END IF
            IPOS = LA
            DISPC = LA + 1
  490       JMORE = JEND + IFILL - IPOS + 1
            DO 500 IPOS = JBEG - 1,MAX(JBEG-JMORE,1),-1
               IF (IRN(IPOS).NE.0) GO TO 510
  500       CONTINUE
            IPOS = IPOS + 1
            IF (IPOS.EQ.JBEG-JMORE) GO TO 520
  510       IF (DISPC+LENC(J)+IFILL.GT.LA+1) THEN
               INFO(2) = INFO(2) + 1
               CALL MA50DD(LA,A,IRN,IQ,N,DISPC,.TRUE.)
               JBEG = IQ(J)
               JEND = JBEG + LENC(J) - 1
               PIVBEG = IQ(JPIV)
               PIVEND = PIVBEG + LENC(JPIV) - 1
               IF (DISPC+LENC(J)+IFILL.GT.LA+1) GO TO 705
            END IF
            IPOS = DISPC
            DISPC = DISPC + LENC(J) + IFILL
  520       IQ(J) = IPOS
            DO 530 II = JBEG,JEND
               A(IPOS) = A(II)
               IRN(IPOS) = IRN(II)
               IPOS = IPOS + 1
               IRN(II) = 0
  530       CONTINUE
            JBEG = IQ(J)
            JEND = IPOS - 1
  540       IDROP = 0
            DO 580 II = PIVBEG + 1,PIVEND
               I = IRN(II)
               INFO(3) = MAX(INFO(3),NERED+LENR(I)+1)
               IF (IW(I).LT.0) THEN
                  IW(I) = -IW(I)
                  GO TO 580
               END IF
               ANEW = AMULT*A(II)
               IF (ABS(ANEW).LT.CNTL(3)) THEN
                  IDROP = IDROP + 1
               ELSE
                  JEND = JEND + 1
                  A(JEND) = ANEW
                  IRN(JEND) = I

                  IEND = IP(I) + LENR(I)
                  IF (IEND.LT.DISPR) THEN
                     IF (JCN(IEND).EQ.0) GO TO 560
                  ELSE
                     IF (DISPR.LE.LA) THEN
                        DISPR = DISPR + 1
                        GO TO 560
                     END IF
                  END IF
                  IF (IP(I).GT.1) THEN
                     IF (JCN(IP(I)-1).EQ.0) THEN
                        IEND = IEND - 1
                        DO 545 JJ = IP(I),IEND
                           JCN(JJ-1) = JCN(JJ)
  545                   CONTINUE
                        IP(I) = IP(I) - 1
                        GO TO 560
                     END IF
                  END IF
                  IF (DISPR+LENR(I).GT.LA) THEN
                     INFO(2) = INFO(2) + 1
                     CALL MA50DD(LA,A,JCN,IP,M,DISPR,.FALSE.)
                     IF (DISPR+LENR(I).GT.LA) GO TO 705
                  END IF
                  J1 = IP(I)
                  J2 = IP(I) + LENR(I) - 1
                  IP(I) = DISPR
                  DO 550 JJ = J1,J2
                     JCN(DISPR) = JCN(JJ)
                     JCN(JJ) = 0
                     DISPR = DISPR + 1
  550             CONTINUE
                  IEND = DISPR
                  DISPR = IEND + 1
  560             JCN(IEND) = J
                  LENR(I) = LENR(I) + 1
               END IF
  580       CONTINUE
            INFO(6) = INFO(6) + IDROP
            NERED = NERED - IDROP
            DO 590 II = 1,IDROP
               IRN(JEND+II) = 0
  590       CONTINUE
            LENC(J) = LENC(J) + IFILL - IDROP
  600    CONTINUE


         DO 610 EYE = 1,NEPR
            JJ = IP(IPIV) + EYE - 1
            J = JCN(JJ)
            JCN(JJ) = 0
            NE1 = LENC(J)
            LASTC(J) = 0
            IF (NE1.GT.0) THEN
               IFIR = JFIRST(NE1)
               JFIRST(NE1) = J
               NEXTC(J) = IFIR
               IF (IFIR.NE.0) LASTC(IFIR) = J
               MINC = MIN(MINC,NE1)
            ELSE IF (ICNTL(7).NE.2) THEN
               IF (INFO(6).EQ.0) NULLC = NULLC + 1
               IF (J.LE.JLAST) THEN
                  NORD = NORD + 1
                  IQ(J) = -NORD
                  IF (NORD.EQ.JLAST) THEN
                     NORD = NORD + NULLJ
                     JLAST = N
                     NULLJ = 0
                  END IF
               ELSE
                  NULLJ = NULLJ + 1
                  IQ(J) = - (JLAST+NULLJ)
               END IF
            END IF
  610    CONTINUE
         NERED = NERED - NEPR

         IF (IPIV.NE.0) THEN
            LENR(IPIV) = 0
            IW(IPIV) = 0
            IRN(PIVBEG) = 0
            MORD = MORD + 1
            PIVBEG = PIVBEG + 1
            IP(IPIV) = -MORD
         END IF
         NERED = NERED - LENPIV - 1
         DO 620 II = PIVBEG,PIVEND
            I = IRN(II)
            IW(I) = 0
            IRN(II) = 0
            NE1 = LENR(I)
            IF (NE1.EQ.0) THEN
               IF (INFO(6).EQ.0) NULLR = NULLR + 1
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            ELSE IF (ICNTL(4).LE.0) THEN
               IFIR = IFIRST(NE1)
               LASTR(I) = 0
               NEXTR(I) = IFIR
               IFIRST(NE1) = I
               IF (IFIR.NE.0) LASTR(IFIR) = I
               MINC = MIN(MINC,NE1)
            END IF
  620    CONTINUE
         IQ(JPIV) = -NORD1
  630 CONTINUE


  640 INFO(5) = MORD + MIN(M-MORD-NULLI,N-NORD-NULLJ)
      DO 650 L = 1,MIN(M-MORD,N-NORD)
         RINFO(1) = RINFO(1) + M - MORD - L + 1 +
     +                   REAL(M-MORD-L)*(N-NORD-L)*2
  650 CONTINUE
      NP = NORD
      INFO(4) = 2 + NEFACT + M*2 + MAX(N-NORD+M-MORD,
     +          (N-NORD)*(M-MORD))
      INFO(6) = INFO(6) + NDROP
      INFO(7) = M - MORD
      DO 660 L = 1,M
         IF (IP(L).LT.0) THEN
            IP(L) = -IP(L)
         ELSE
            MORD = MORD + 1
            IP(L) = MORD
         END IF
  660 CONTINUE
      DO 670 L = 1,N
         IF (IQ(L).LT.0) THEN
            LASTC(L) = -IQ(L)
         ELSE
            IF (NORD.EQ.JLAST) NORD = NORD + NULLJ
            NORD = NORD + 1
            LASTC(L) = NORD
         END IF
  670 CONTINUE
      DO 680 L = 1,N
         IQ(LASTC(L)) = L
  680 CONTINUE

      IF (INFO(5).LT.MIN(M,N)) INFO(1) = INFO(1) + 1

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(A,I6,A,F12.1/A,7I8)') ' Leaving MA50AD with NP =',
     +     NP,' RINFO(1) =',RINFO(1),' INFO =',(INFO(I),I=1,7)
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            WRITE (MP,'(A,(T6,10(I7)))') ' IQ = ',IQ
         END IF
      END IF

      GO TO 750

  690 INFO(1) = -1
      IF (LP.GT.0) WRITE (LP,'(/A/(2(A,I8)))')
     +    ' **** Error return from MA50AD ****',' M =',M,' N =',N
      GO TO 750
  700 INFO(1) = -2
      IF (LP.GT.0) WRITE (LP,'(/A/(A,I10))')
     +    ' **** Error return from MA50AD ****',' NE =',NE
      GO TO 750
  705 INFO(4) =  NEFACT + NERED
      INFO(6) = INFO(6) + NDROP
  710 INFO(1) = -3
      IF (LP.GT.0) WRITE (LP,'(/A/A,I9,A,I9)')
     +    ' **** Error return from MA50AD ****',
     +    ' LA  must be increased from',LA,' to at least',INFO(3)
      GO TO 750
  720 INFO(1) = -4
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50AD ****',' Entry in row',I,
     +    ' and column',J,' duplicated'
      GO TO 750
  730 INFO(1) = -5
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50AD ****',' Fault in component ',
     +    PIVOT,' of column permutation given in IFIRST'
      GO TO 750
  740 INFO(1) = -6
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50AD ****',' ICNTL(4) = ',ICNTL(4),
     +    ' when ICNTL(6) = 2'
  750 END


      SUBROUTINE MA50BD(M,N,NE,JOB,AA,IRNA,IPTRA,CNTL,ICNTL,IP,IQ,NP,
     +                  LFACT,FACT,IRNF,IPTRL,IPTRU,W,IW,INFO,RINFO)
      INTEGER M,N,NE,JOB
      DOUBLE PRECISION AA(NE)
      INTEGER IRNA(NE),IPTRA(N)
      DOUBLE PRECISION CNTL(10)
      INTEGER ICNTL(20),IP(M),IQ(*),NP,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER IRNF(LFACT),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION W(M)
      INTEGER IW(M+2*N),INFO(15)
      DOUBLE PRECISION RINFO(10)

      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0D0,ONE=1.0D0)
      DOUBLE PRECISION AMULT,ASW
      INTEGER BEGCOL
      LOGICAL DROP
      INTEGER ENDCOL,EYE,EYE1,I,IA1,IA2,IF1,IF2,II,IL1,IL2,IPIV,IQPIV,
     +        IU1,IU2,ISW,J,JDUMMY,JJ,JLAST,K,LP
      DOUBLE PRECISION MAXENT
      INTEGER MF,MORD,MP,NEU,NF,NULLC
      DOUBLE PRECISION PIVLIM
      INTEGER RANK
      DOUBLE PRECISION U
      EXTERNAL MA50ED,MA50FD,MA50GD
      INTRINSIC ABS,MAX,MIN

      INFO(1) = 0
      INFO(4) = 0
      INFO(5) = 0
      INFO(6) = 0
      INFO(7) = 0
      RINFO(1) = ZERO
      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (ICNTL(3).LE.0) LP = 0
      IF (ICNTL(3).LE.1) MP = 0

      IF (M.LT.1 .OR. N.LT.1) THEN
         INFO(1) = -1
         IF (LP.GT.0) WRITE (LP,'(/A/A,I8,A,I8)')
     +       ' **** Error return from MA50BD ****',' M =',M,' N =',N
         GO TO 550
      END IF
      IF (NE.LE.0) THEN
         INFO(1) = -2
         IF (LP.GT.0) WRITE (LP,'(/A/A,I6)')
     +       ' **** Error return from MA50BD ****',' NE =',NE
         GO TO 550
      END IF
      IF (NP.LT.0 .OR. NP.GT.N) THEN
         INFO(1) = -7
         IF (LP.GT.0) WRITE (LP,'(/A/A,I8,A,I8)')
     +       ' **** Error return from MA50BD ****',' NP =',NP,' N =',N
         GO TO 550
      END IF
      IF (LFACT.LT.MAX(M,NE+2)) THEN
         INFO(4) = MAX(M,NE+2)
         GO TO 520
      END IF
      IF (JOB.EQ.1) THEN
      ELSE IF (JOB.EQ.2 .OR. JOB.EQ.3) THEN
         IF (IRNF(1).NE.0) THEN
            INFO(1) = -6
            IF (LP.GT.0) WRITE (LP,'(/A/A,I1,A)')
     +          ' **** Error return from MA50BD ***',' Call with JOB=',
     +          JOB,' follows JOB=1 call in which entries were dropped'
            GO TO 550
         END IF
      ELSE
         INFO(1) = -5
         IF (LP.GT.0) WRITE (LP,'(/A/A,I2)')
     +       ' **** Error return from MA50BD ****',' JOB =',JOB
         GO TO 550
      END IF

      IF (MP.GT.0) THEN
         IF (ICNTL(3).GT.2) WRITE (MP,
     +       '(/2(A,I6),A,I8,A,I3/A,I8,A,I7/A,1P,4E10.2/A,7I8)')
     +       ' Entering MA50BD with M =',M,' N =',N,' NE =',NE,' JOB =',
     +       JOB,' LFACT =',LFACT,' NP =',NP,' CNTL =',(CNTL(I),I=1,4),
     +       ' ICNTL =',(ICNTL(I),I=1,7)
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            IF (IQ(1).GT.0) THEN
               WRITE (MP,'(A,(T6,10(I7)))') ' IQ = ',(IQ(J),J=1,N)
            ELSE
               WRITE (MP,'(A,(T6,I7))') ' IQ = ',IQ(1)
            END IF
            DO 10 J = 1,N - 1
               IF (IPTRA(J).LT.IPTRA(J+1)) WRITE (MP,
     +             '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J,
     +             (AA(II),IRNA(II),II=IPTRA(J),IPTRA(J+1)-1)
   10       CONTINUE
            IF (IPTRA(N).LE.NE) WRITE (MP,
     +          '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',N,
     +          (AA(II),IRNA(II),II=IPTRA(N),NE)
         END IF
      END IF

      JLAST = 0
      NULLC = 0
      IF (JOB.GT.1 .AND. ICNTL(6).GT.0 .AND.
     +    ICNTL(6).LT.N) JLAST = MIN(NP,N-ICNTL(6))

      U = MIN(CNTL(2),ONE)
      U = MAX(U,ZERO)
      DO 20 I = 1,M
         IW(I+N) = 0
         W(I) = ZERO
   20 CONTINUE
      MORD = 0
      IF1 = LFACT + 1
      IF2 = 0
      NF = N - NP
      MF = 0
      IL2 = 2
      IF (JLAST.GT.0) IL2 = IPTRL(JLAST)
      NEU = 0

      IF (JOB.EQ.2) GO TO 370

      IF (JOB.EQ.3) THEN
         DO 30 J = 1,NP
            IA1 = IPTRU(J) + 1
            IF (IA1.GT.IPTRL(J)) GO TO 30
            IF (J.LE.JLAST) THEN
               MORD = MORD + 1
               IP(IRNF(IA1)) = -J
            ELSE
               IP(IRNF(IA1)) = J
            END IF
   30    CONTINUE
         MF = IRNF(2)
         IA1 = IPTRL(N)
         DO 40 J = 1,MF
            IP(IRNF(IA1+J)) = NP + J
   40    CONTINUE
      END IF

      DO 50 K = 1,JLAST
         IW(M+N+K) = IPTRL(K)
   50 CONTINUE

      DO 310 K = JLAST + 1,N
         DROP = .FALSE.
         IF (K.EQ.NP+1) THEN
            MF = M - MORD
            IF1 = LFACT + 1 - MF
            II = 0
            DO 60 I = 1,M
               IF (IP(I).GT.0) THEN
                  IW(I+N) = N
                  IRNF(IF1+II) = I
                  II = II + 1
                  IP(I) = NP + II
               END IF
   60       CONTINUE
            IF1 = LFACT + 1 - MAX(MF*NF,MF+NF)
            IF2 = IF1 - 1 + MF*MAX(0,JLAST-NP)
         END IF
         J = K
         IF (IQ(1).GT.0) J = IQ(K)
         IA1 = IPTRA(J)
         IA2 = NE
         IF (J.NE.N) IA2 = IPTRA(J+1) - 1
         IU1 = IL2 + 1
         IU2 = IU1 - 1
         IL1 = IF1 - 1 + IA1 - IA2
         IL2 = IL1 - 1
         INFO(4) = MAX(INFO(4),NEU+LFACT-IL1+IU2+M+1)
         IF (IL1-IU2.LE.M) THEN
            IF (INFO(1).NE.-3) THEN
               INFO(1) = -3
               IF (ICNTL(8).NE.0) GO TO 480
               NEU = IL2 + LFACT + 1 - MF - IF1
               IF1 = LFACT + 1 - MF
               IF2 = IF1 - 1
               IL2 = 0
               EYE = 0
               DO 80 J = 1,MIN(K-1,NP)
                  IU2 = IPTRU(J)
                  IPTRU(J) = EYE
                  IL2 = IPTRL(J)
                  NEU = NEU + IU2 - IL2
                  DO 70 II = IU2 + 1,IL2
                     EYE = EYE + 1
                     IRNF(EYE) = IRNF(II)
                     FACT(EYE) = FACT(II)
   70             CONTINUE
                  IPTRL(J) = EYE
                  IW(M+N+J) = EYE
   80          CONTINUE
               IU1 = EYE + 1
               IU2 = EYE
               IL1 = IF1 - 1 + IA1 - IA2
               IL2 = IL1 - 1
            END IF
            IF (IL1-IU2.LE.M) GO TO 480
         END IF
         EYE = IL1
         DO 90 II = IA1,IA2
            I = IRNA(II)
            IF (IW(I+N).EQ.-1) GO TO 540
            IW(I+N) = -1
            W(I) = AA(II)
            IRNF(EYE) = I
            EYE = EYE + 1
   90    CONTINUE
         IPTRL(K) = EYE - 1
         IW(M+N+K) = EYE - 1
         IW(K) = IL1
         J = K
         DO 120 JDUMMY = 1,2*K
            DO 100 II = IW(J),ABS(IW(M+N+J))
               I = IRNF(II)
               IF (IW(I+N).GE.K) GO TO 100
               IF (IP(I).LE.0) GO TO 110
               IW(I+N) = K
               IL1 = IL1 - 1
               IRNF(IL1) = I
  100       CONTINUE
            IF (J.EQ.K) GO TO 130
            IU2 = IU2 + 1
            I = IRNF(IPTRU(J)+1)
            IRNF(IU2) = I
            J = -IW(I+N)
            IW(I+N) = K
            GO TO 120
  110       IW(I+N) = -J
            IW(J) = II + 1
            J = -IP(I)
            IW(J) = IPTRU(J) + 2
  120    CONTINUE
  130    DO 150 II = IU2,IU1,-1
            I = IRNF(II)
            J = -IP(I)
            EYE1 = IPTRU(J) + 1
            IF (ABS(W(I)).LT.CNTL(3)) GO TO 150
            AMULT = -W(I)*FACT(EYE1)
            W(I) = AMULT
            DO 140 EYE = EYE1 + 1,IPTRL(J)
               I = IRNF(EYE)
               W(I) = W(I) + AMULT*FACT(EYE)
  140       CONTINUE
            RINFO(1) = RINFO(1) + ONE + 2*(IPTRL(J)-EYE1)
  150    CONTINUE

         IF (CNTL(3).GT.ZERO) THEN
            EYE = IU1
            DO 160 II = IU1,IU2
               I = IRNF(II)
               IF (ABS(W(I)).LT.CNTL(3)) THEN
                  INFO(6) = INFO(6) + 1
               ELSE
                  IRNF(EYE) = -IP(I)
                  FACT(EYE) = W(I)
                  EYE = EYE + 1
               END IF
               W(I) = ZERO
  160       CONTINUE
            IU2 = EYE - 1
         ELSE
            DO 170 II = IU1,IU2
               I = IRNF(II)
               IRNF(II) = -IP(I)
               FACT(II) = W(I)
               W(I) = ZERO
  170       CONTINUE
         END IF
         IF (INFO(1).EQ.-3) THEN
            NEU = NEU + IU2 - IU1 + 1
            IU2 = IU1 - 1
         END IF
         IPTRU(K) = IU2
         IF (K.LE.NP) THEN
            MAXENT = ZERO
            IF (CNTL(3).GT.ZERO) THEN
               EYE = IL1
               DO 180 II = IL1,IL2
                  I = IRNF(II)
                  IF (ABS(W(I)).LT.CNTL(3)) THEN
                     INFO(6) = INFO(6) + 1
                     W(I) = ZERO
                     DROP = .TRUE.
                  ELSE
                     IRNF(EYE) = I
                     EYE = EYE + 1
                     MAXENT = MAX(ABS(W(I)),MAXENT)
                  END IF
  180          CONTINUE
               IL2 = EYE - 1
            ELSE
               DO 190 II = IL1,IL2
                  MAXENT = MAX(ABS(W(IRNF(II))),MAXENT)
  190          CONTINUE
            END IF
            PIVLIM = U*MAXENT
            EYE = IU2
            IQPIV = M + N
            IF (IL1.GT.IL2) NULLC = NULLC + 1
            DO 200 II = IL1,IL2
               I = IRNF(II)
               EYE = EYE + 1
               IRNF(EYE) = I
               FACT(EYE) = W(I)
               W(I) = ZERO
               IF (ABS(FACT(EYE)).GE.PIVLIM) THEN
                  IF (ABS(FACT(EYE)).GT.CNTL(4)) THEN
                     IF (IP(I).LT.IQPIV) THEN
                        IQPIV = IP(I)
                        IPIV = EYE
                     END IF
                  END IF
               END IF
  200       CONTINUE
            IL1 = IU2 + 1
            IL2 = EYE
            IF (IL1.LE.IL2) THEN
               IF (IQPIV.EQ.M+N) THEN
                  IF (CNTL(3).GT.ZERO) INFO(6) = INFO(6) + EYE - IU2
                  IL2 = IU2
               ELSE
                  IF (IL1.NE.IPIV) THEN
                     ASW = FACT(IPIV)
                     FACT(IPIV) = FACT(IL1)
                     FACT(IL1) = ASW
                     ISW = IRNF(IL1)
                     IRNF(IL1) = IRNF(IPIV)
                     IRNF(IPIV) = ISW
                  END IF
                  INFO(5) = INFO(5) + 1
                  FACT(IL1) = ONE/FACT(IL1)
                  RINFO(1) = RINFO(1) + ONE
                  MORD = MORD + 1
                  IP(IRNF(IL1)) = -K
               END IF
            END IF
         ELSE
            IL2 = IPTRU(K)
            DO 210 II = LFACT - MF + 1,LFACT
               I = IRNF(II)
               IF2 = IF2 + 1
               FACT(IF2) = W(I)
               W(I) = ZERO
  210       CONTINUE
            IF (INFO(1).EQ.-3) IF2 = IF2 - MF
         END IF
         IW(M+N+K) = IL2
         IPTRL(K) = IL2
         IF (DROP) GO TO 310
         DO 300 II = IU1,IU2
            I = IRNF(II)
            IF (IW(M+N+I).LT.0) GO TO 300
            BEGCOL = IPTRU(I) + 2
            ENDCOL = IPTRL(I)
            IF (K.LE.NP) THEN
               DO 220 JJ = BEGCOL,ENDCOL
                  IF (IP(IRNF(JJ)).EQ.-K) GO TO 230
  220          CONTINUE
               GO TO 300
            END IF
  230       DO 280 JDUMMY = BEGCOL,ENDCOL
               JJ = BEGCOL
               DO 240 BEGCOL = JJ,ENDCOL
                  IF (IP(IRNF(BEGCOL)).GT.0) GO TO 250
  240          CONTINUE
               GO TO 290
  250          JJ = ENDCOL
               DO 260 ENDCOL = JJ,BEGCOL,-1
                  IF (IP(IRNF(ENDCOL)).LT.0) GO TO 270
  260          CONTINUE
               GO TO 290
  270          ASW = FACT(BEGCOL)
               FACT(BEGCOL) = FACT(ENDCOL)
               FACT(ENDCOL) = ASW
               J = IRNF(BEGCOL)
               IRNF(BEGCOL) = IRNF(ENDCOL)
               IRNF(ENDCOL) = J
               BEGCOL = BEGCOL + 1
               ENDCOL = ENDCOL - 1
  280       CONTINUE
  290       IW(M+N+I) = -ENDCOL
  300    CONTINUE
  310 CONTINUE
      IF (N.EQ.NP) THEN
         MF = M - MORD
         IF1 = LFACT + 1 - MF
         II = 0
         DO 320 I = 1,M
            IF (IP(I).GT.0) THEN
               IW(I+N) = N
               IRNF(IF1+II) = I
               II = II + 1
               IP(I) = NP + II
            END IF
  320    CONTINUE
         IF1 = LFACT + 1 - MAX(MF*NF,MF+NF)
         IF2 = IF1 - 1 + MF*MAX(0,JLAST-NP)
      END IF
      IF (INFO(5).EQ.MIN(M,N)) THEN
         DO 330 I = 1,M
            IP(I) = ABS(IP(I))
  330    CONTINUE
      ELSE
         MORD = NP
         DO 340 I = 1,M
            IF (IP(I).LT.0) THEN
               IP(I) = -IP(I)
            ELSE
               MORD = MORD + 1
               IP(I) = MORD
            END IF
  340    CONTINUE
      END IF
      IRNF(1) = INFO(6)
      IRNF(2) = MF
      INFO(7) = MF
      FACT(1) = CNTL(3)
      FACT(2) = CNTL(4)
      IF (INFO(1).EQ.-3) GO TO 520
      IF2 = IF2 - MF*NF
      DO 350 II = 1,MF*NF
         FACT(IL2+II) = FACT(IF1-1+II)
  350 CONTINUE
      DO 360 II = 1,MF
         IRNF(IL2+II) = IRNF(LFACT-MF+II)
  360 CONTINUE
      IF1 = IL2 + 1
      GO TO 440
  370 MF = IRNF(2)
      IF1 = IPTRL(N) + 1
      IF2 = IF1 - 1
      DO 430 K = JLAST + 1,N
         J = K
         IF (IQ(1).GT.0) J = IQ(K)
         IA1 = IPTRA(J)
         IA2 = NE
         IF (J.NE.N) IA2 = IPTRA(J+1) - 1
         IU1 = IL2 + 1
         IU2 = IPTRU(K)
         IL1 = IU2 + 1
         IL2 = IPTRL(K)
         DO 380 II = IA1,IA2
            W(IRNA(II)) = AA(II)
  380    CONTINUE
         DO 400 II = IU2,IU1,-1
            J = IRNF(II)
            I = IRNF(IPTRU(J)+1)
            EYE1 = IPTRU(J) + 1
            AMULT = -W(I)*FACT(EYE1)
            FACT(II) = AMULT
            W(I) = ZERO
            DO 390 EYE = EYE1 + 1,IPTRL(J)
               I = IRNF(EYE)
               W(I) = W(I) + AMULT*FACT(EYE)
  390       CONTINUE
            RINFO(1) = RINFO(1) + ONE + 2*(IPTRL(J)-EYE1)
  400    CONTINUE
         IF (K.LE.NP) THEN
            IF (IL1.LE.IL2) THEN
               DO 410 II = IL1,IL2
                  I = IRNF(II)
                  FACT(II) = W(I)
                  W(I) = ZERO
  410          CONTINUE
               IF (ABS(FACT(IL1)).LE.CNTL(4)) THEN
                  GO TO 530
               ELSE
                  INFO(5) = INFO(5) + 1
                  FACT(IL1) = ONE/FACT(IL1)
                  RINFO(1) = RINFO(1) + ONE
               END IF
            END IF
         ELSE
            DO 420 II = IF1,IF1 + MF - 1
               I = IRNF(II)
               IF2 = IF2 + 1
               FACT(IF2) = W(I)
               W(I) = ZERO
  420       CONTINUE
         END IF
  430 CONTINUE
      INFO(4) = MAX(IF1+MF+NF-1,IF2)

  440 IF (MF.GT.0 .AND. NF.GT.0) THEN
         IF (ICNTL(5).GT.1) CALL MA50GD(MF,NF,FACT(IF1),MF,ICNTL(5),
     +                                  CNTL(4),IRNF(IF1+MF),RANK)
         IF (ICNTL(5).EQ.1) CALL MA50FD(MF,NF,FACT(IF1),MF,CNTL(4),
     +                                  IRNF(IF1+MF),RANK)
         IF (ICNTL(5).LE.0) CALL MA50ED(MF,NF,FACT(IF1),MF,CNTL(4),
     +                                  IRNF(IF1+MF),RANK)
         INFO(5) = INFO(5) + RANK
         DO 450 I = 1,MIN(MF,NF)
            RINFO(1) = RINFO(1) + MF - I + 1 + REAL(MF-I)*(NF-I)*2
  450    CONTINUE
      END IF
      IF (INFO(5).LT.MIN(M,N)) INFO(1) = 1
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(A,I6,A,F12.1/A,I3,A,4I8)')
     +     ' Leaving MA50BD with IRNF(2) =',IRNF(2),
     +     ' RINFO(1) =',RINFO(1),
     +     ' INFO(1) =',INFO(1),' INFO(4:7) =', (INFO(J),J=4,7)
         IF (ICNTL(3).GT.3) THEN
            IF (JOB.NE.2) WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            DO 460 J = 1,N
               IF (J.GT.1) THEN
                  IF (IPTRL(J-1).LT.IPTRU(J)) WRITE (MP,
     +                '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,
     +                ' of U', (FACT(II),IRNF(II),II=IPTRL(J-1)+1,
     +                IPTRU(J))
               END IF
               IF (IPTRU(J).LT.IPTRL(J)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L',
     +              (FACT(II),IRNF(II),II=IPTRU(J)+1,IPTRL(J))
  460       CONTINUE
            WRITE (MP,'(A)') ' Full part'
            WRITE (MP,'((6I12))') (IRNF(IF1+MF+J),J=0,NF-1)
            DO 470 I = 0,MF - 1
               WRITE (MP,'(I4,1P,6E12.4:/(4X,1P,6E12.4))') IRNF(IF1+I),
     +            (FACT(IF1+I+J*MF),J=0,NF-1)
  470       CONTINUE
         END IF
      END IF
      GO TO 550

  480 DO 490 I = 1,M
         IW(I) = 0
  490 CONTINUE
      DO 500 I = 1,M
         IF (IP(I).GT.0) THEN
            IW(IP(I)) = I
         ELSE
            IP(I) = -IP(I)
         END IF
  500 CONTINUE
      DO 510 I = 1,M
         IF (IW(I).GT.0) THEN
            IP(IW(I)) = K
            K = K + 1
         END IF
  510 CONTINUE
  520 INFO(1) = -3
      IF (LP.GT.0) THEN
         WRITE (LP,'(/A)')' **** Error return from MA50BD **** '
         IF (ICNTL(8).EQ.0) THEN
            WRITE (LP,'(A,I7,A,I7)')' LFACT must be increased from',
     +         LFACT,' to at least',INFO(4)
         ELSE
            WRITE (LP,'(A,I7)')' LFACT must be increased from',LFACT
         END IF
      END IF
      GO TO 550
  530 INFO(1) = - (7+K)
      IF (LP.GT.0) WRITE (LP,'(/A/A,I6,A)')
     +    ' **** Error return from MA50BD **** ',
     +    ' Small pivot found in column',K,' of the permuted matrix.'
      GO TO 550
  540 INFO(1) = -4
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50BD ****',' Entry in row',I,
     +    ' and column',J,' duplicated'
  550 END

      SUBROUTINE MA50CD(M,N,ICNTL,IQ,NP,TRANS,LFACT,FACT,IRNF,IPTRL,
     +                  IPTRU,B,X,W,INFO)
      INTEGER M,N,ICNTL(20),IQ(*),NP
      LOGICAL TRANS
      INTEGER LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER IRNF(LFACT),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION B(*),X(*),W(*)
      INTEGER INFO(15)

      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0D0)
      INTEGER I,II,IA1,IF1,J,LP,MF,MP,NF
      DOUBLE PRECISION PROD

      EXTERNAL MA50HD

      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (ICNTL(3).LE.0) LP = 0
      IF (ICNTL(3).LE.1) MP = 0

      IF (M.LT.1 .OR. N.LT.1) GO TO 250

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) WRITE (MP,
     +    '(/2(A,I6),A,I4,A,L2)') ' Entering MA50CD with M=',M,' N =',N,
     +    ' NP =',NP,' TRANS =',TRANS
      IF1 = IPTRL(N) + 1
      MF = IRNF(2)
      NF = N - NP
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) WRITE (MP,
     +    '(A,I5,A,I5)') ' Size of full submatrix',MF,' by',NF
      IF (MP.GT.0 .AND. ICNTL(3).GT.3) THEN
         DO 10 J = 1,N
            IF (J.GT.1) THEN
               IF (IPTRL(J-1).LT.IPTRU(J)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of U',
     +              (FACT(II),IRNF(II),II=IPTRL(J-1)+1,IPTRU(J))
            END IF
            IF (IPTRU(J).LT.IPTRL(J)) WRITE (MP,
     +          '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L',
     +          (FACT(II),IRNF(II),II=IPTRU(J)+1,IPTRL(J))
   10    CONTINUE
         WRITE (MP,'(A)') ' Full part'
         WRITE (MP,'((6I12))') (IRNF(IF1+MF+J),J=0,NF-1)
         DO 20 I = 0,MF - 1
            WRITE (MP,'(I4,1P,6E12.4:/(4X,1P,6E12.4))') IRNF(IF1+I),
     +        (FACT(IF1+I+J*MF),J=0,NF-1)
   20    CONTINUE
      END IF

      IF (TRANS) THEN
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A4,5F10.4:/(4X,5F10.4))') ' B =', (B(I),I=1,N)
         IF (IQ(1).GT.0) THEN
            DO 30 I = 1,N
               W(I) = B(IQ(I))
   30       CONTINUE
         ELSE
            DO 40 I = 1,N
               W(I) = B(I)
   40       CONTINUE
         END IF
         DO 50 I = 1,M
            X(I) = ZERO
   50    CONTINUE
         DO 70 I = 2,N
            PROD = ZERO
            DO 60 II = IPTRL(I-1) + 1,IPTRU(I)
               PROD = PROD + FACT(II)*W(IRNF(II))
   60       CONTINUE
            W(I) = W(I) + PROD
   70    CONTINUE
         DO 80 I = 1,NF
            X(I) = W(NP+I)
   80    CONTINUE
         IF (MF.GT.0 .AND. NF.GT.0) THEN
            CALL MA50HD(TRANS,MF,NF,FACT(IF1),MF,IRNF(IF1+MF),X,
     +                  ICNTL(5))
         ELSE
            DO 90 I = 1,MF
               X(I) = ZERO
   90       CONTINUE
         END IF
         DO 100 I = MF,1,-1
            J = IRNF(IF1+I-1)
            IF (J.NE.I) X(J) = X(I)
  100    CONTINUE
         DO 120 I = NP,1,-1
            IA1 = IPTRU(I) + 1
            IF (IA1.GT.IPTRL(I)) GO TO 120
            PROD = ZERO
            DO 110 II = IA1 + 1,IPTRL(I)
               PROD = PROD + FACT(II)*X(IRNF(II))
  110       CONTINUE
            X(IRNF(IA1)) = (W(I)-PROD)*FACT(IA1)
  120    CONTINUE
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A/(4X,5F10.4))') ' Leaving MA50CD with X =', (X(I),I=1,M)
      ELSE
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A4,5F10.4:/(4X,5F10.4))') ' B =', (B(I),I=1,M)
         DO 130 I = 1,M
            W(I) = B(I)
  130    CONTINUE
         DO 150 I = 1,NP
            IA1 = IPTRU(I) + 1
            IF (IA1.LE.IPTRL(I)) THEN
               X(I) = W(IRNF(IA1))*FACT(IA1)
               IF (X(I).NE.ZERO) THEN
                  DO 140 II = IA1 + 1,IPTRL(I)
                     W(IRNF(II)) = W(IRNF(II)) - FACT(II)*X(I)
  140             CONTINUE
               END IF
            END IF
  150    CONTINUE
         IF (MF.GT.0 .AND. NF.GT.0) THEN
            DO 160 I = 1,MF
               W(I) = W(IRNF(IF1+I-1))
  160       CONTINUE
            CALL MA50HD(TRANS,MF,NF,FACT(IF1),MF,IRNF(IF1+MF),W,
     +                  ICNTL(5))
            DO 170 I = 1,NF
               X(NP+I) = W(I)
  170       CONTINUE
         ELSE
            DO 180 I = 1,NF
               X(NP+I) = ZERO
  180       CONTINUE
         END IF
         DO 200 J = N,MAX(2,NP+1),-1
            PROD = X(J)
            DO 190 II = IPTRL(J-1) + 1,IPTRU(J)
               X(IRNF(II)) = X(IRNF(II)) + FACT(II)*PROD
  190       CONTINUE
  200    CONTINUE
         DO 220 J = NP,2,-1
            IA1 = IPTRU(J)
            IF (IA1.GE.IPTRL(J)) THEN
               X(J) = ZERO
            ELSE
               PROD = X(J)
               DO 210 II = IPTRL(J-1) + 1,IA1
                  X(IRNF(II)) = X(IRNF(II)) + FACT(II)*PROD
  210          CONTINUE
            END IF
  220    CONTINUE
         IF (NP.GE.1 .AND. IPTRU(1).GE.IPTRL(1)) X(1) = ZERO
         IF (IQ(1).GT.0) THEN
            DO 230 I = 1,N
               W(I) = X(I)
  230       CONTINUE
            DO 240 I = 1,N
               X(IQ(I)) = W(I)
  240       CONTINUE
         END IF
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A/(4X,5F10.4))') ' Leaving MA50CD with X =', (X(I),I=1,N)
      END IF
      RETURN
  250 INFO(1) = -1
      IF (LP.GT.0) WRITE (LP,'(/A/2(A,I8))')
     +    ' **** Error return from MA50CD ****',' M =',M,' N =',N
      END

      SUBROUTINE MA50DD(LA,A,IND,IPTR,N,DISP,REALS)
      INTEGER LA,N,DISP
      DOUBLE PRECISION A(LA)
      INTEGER IPTR(N)
      LOGICAL REALS
      INTEGER IND(LA)
      INTEGER J,K,KN
      DO 10 J = 1,N
         K = IPTR(J)
         IF (K.GT.0) THEN
            IPTR(J) = IND(K)
            IND(K) = -J
         END IF
   10 CONTINUE
      KN = 0
      DO 20 K = 1,DISP - 1
         IF (IND(K).EQ.0) GO TO 20
         KN = KN + 1
         IF (REALS) A(KN) = A(K)
         IF (IND(K).LE.0) THEN
            J = -IND(K)
            IND(K) = IPTR(J)
            IPTR(J) = KN
         END IF
         IND(KN) = IND(K)
   20 CONTINUE
      DISP = KN + 1
      END


      SUBROUTINE MA50ED(M,N,A,LDA,PIVTOL,IPIV,RANK)
      INTEGER LDA,M,N,RANK
      DOUBLE PRECISION PIVTOL

      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N)



      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)

      INTEGER I,J,JP,K
      LOGICAL PIVOT

      INTEGER IDAMAX
      EXTERNAL IDAMAX

      EXTERNAL DAXPY,DSCAL,DSWAP
      INTRINSIC ABS

      J = 1
      DO 30 K = 1,N

         DO 10 I = 1,J - 1
            IF (M.GT.I) CALL DAXPY(M-I,-A(I,J),A(I+1,I),1,A(I+1,J),1)
   10    CONTINUE

         IF (J.LE.M) THEN
            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .GT. PIVTOL
         ELSE
            PIVOT = .FALSE.
         END IF
         IF (PIVOT) THEN

            IF (JP.NE.J) CALL DSWAP(N+J-K,A(J,1),LDA,A(JP,1),LDA)
            IF (J.LT.M) CALL DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)

            J = J + 1
         ELSE
            DO 20 I = J,M
               A(I,J) = ZERO
   20       CONTINUE
            IF (K.LT.N) CALL DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            IPIV(N-K+J) = -J
         END IF
   30 CONTINUE

      RANK = J - 1
      END


      SUBROUTINE MA50FD(M,N,A,LDA,PIVTOL,IPIV,RANK)
      INTEGER LDA,M,N,RANK
      DOUBLE PRECISION PIVTOL

      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N)



      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)

      INTEGER I,J,JP,K
      LOGICAL PIVOT

      INTEGER IDAMAX
      EXTERNAL IDAMAX

      EXTERNAL DGEMV,DSCAL,DSWAP

      INTRINSIC ABS

      J = 1
      DO 20 K = 1,N

         IF (J.LE.M) THEN
            CALL DGEMV('No transpose',M-J+1,J-1,-ONE,A(J,1),LDA,A(1,J),
     +                 1,ONE,A(J,J),1)
            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .GT. PIVTOL
         ELSE
            PIVOT = .FALSE.
         END IF

         IF (PIVOT) THEN

            IF (JP.NE.J) CALL DSWAP(N+J-K,A(J,1),LDA,A(JP,1),LDA)
            IF (J.LT.M) CALL DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)

            IF (J.LT.N) THEN
               CALL DGEMV('Transpose',J-1,N-J,-ONE,A(1,J+1),LDA,A(J,1),
     +                    LDA,ONE,A(J,J+1),LDA)
            END IF

            J = J + 1
         ELSE
            DO 10 I = J,M
               A(I,J) = ZERO
   10       CONTINUE
            IF (K.LT.N) CALL DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            IPIV(N-K+J) = -J
         END IF
   20 CONTINUE

      RANK = J - 1
      END


      SUBROUTINE MA50GD(M,N,A,LDA,NB,PIVTOL,IPIV,RANK)
      INTEGER LDA,M,N,NB,RANK
      DOUBLE PRECISION PIVTOL

      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N)

      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)

      INTEGER I,J,JJ,JP,J1,J2,K
      LOGICAL PIVOT
      DOUBLE PRECISION TEMP


      EXTERNAL DGEMM,DGEMV,DSWAP,DSCAL,DTRSM,DTRSV

      INTEGER IDAMAX
      EXTERNAL IDAMAX

      INTRINSIC ABS,MIN

      J = 1
      J1 = 1
      J2 = MIN(N,NB)
      DO 70 K = 1,N

         IF (J.LE.M) THEN

            CALL DGEMV('No transpose',M-J+1,J-J1,-ONE,A(J,J1),LDA,
     +                 A(J1,J),1,ONE,A(J,J),1)

            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .GT. PIVTOL
         ELSE
            PIVOT = .FALSE.
         END IF

         IF (PIVOT) THEN

            IF (JP.NE.J) CALL DSWAP(J2-J1+1,A(J,J1),LDA,A(JP,J1),LDA)
            IF (J.LT.M) CALL DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)

            IF (J+1.LE.J2) THEN
               CALL DGEMV('Transpose',J-J1,J2-J,-ONE,A(J1,J+1),LDA,
     +                    A(J,J1),LDA,ONE,A(J,J+1),LDA)
            END IF

            J = J + 1
         ELSE

            DO 10 I = J,M
               A(I,J) = ZERO
   10       CONTINUE
            IPIV(N-K+J) = -J
            IF (K.NE.N) CALL DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            IF (N-K+J.GT.J2) THEN
               DO 20 I = J1,J - 1
                  JP = IPIV(I)
                  TEMP = A(I,J)
                  A(I,J) = A(JP,J)
                  A(JP,J) = TEMP
   20          CONTINUE
               IF(J.GT.J1) CALL DTRSV('Lower','No transpose','Unit',
     +                                J-J1,A(J1,J1),LDA,A(J1,J),1)
            ELSE
               J2 = J2 - 1
            END IF
         END IF

         IF (J.GT.J2) THEN
            DO 40 JJ = 1,J1 - 1
               DO 30 I = J1,J2
                  JP = IPIV(I)
                  TEMP = A(I,JJ)
                  A(I,JJ) = A(JP,JJ)
                  A(JP,JJ) = TEMP
   30          CONTINUE
   40       CONTINUE
            DO 60 JJ = J2 + 1,N - K + J - 1
               DO 50 I = J1,J2
                  JP = IPIV(I)
                  TEMP = A(I,JJ)
                  A(I,JJ) = A(JP,JJ)
                  A(JP,JJ) = TEMP
   50          CONTINUE
   60       CONTINUE

            IF (K.NE.N) THEN
               CALL DTRSM('Left','Lower','No transpose','Unit',J2-J1+1,
     +                    N-K,ONE,A(J1,J1),LDA,A(J1,J2+1),LDA)
               IF (M.GT.J2) CALL DGEMM('No transpose','No transpose',
     +                                 M-J2,N-K,J2-J1+1,-ONE,A(J2+1,J1),
     +                                 LDA,A(J1,J2+1),LDA,ONE,
     +                                 A(J2+1,J2+1),LDA)
            END IF

            J1 = J2 + 1
            J2 = MIN(J2+NB,N-K+J-1)

         END IF

   70 CONTINUE
      RANK = J - 1
      END

      SUBROUTINE MA50HD(TRANS,M,N,A,LDA,IPIV,B,ICNTL5)
      LOGICAL TRANS
      INTEGER LDA,M,N,ICNTL5
      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N),B(*)


      INTEGER I,K,RANK

      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0D0)
      DOUBLE PRECISION TEMP
      INTRINSIC MIN
      EXTERNAL DAXPY,DDOT,DTRSV
      DOUBLE PRECISION DDOT

      RANK = 0
      DO 10 RANK = MIN(M,N),1,-1
         IF (IPIV(RANK).GT.0) GO TO 20
   10 CONTINUE

   20 IF (.NOT.TRANS) THEN
         DO 30 I = 1,RANK
            K = IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   30    CONTINUE
         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('L','NoTrans','Unit',RANK,A,LDA,B,
     +                                1)
         ELSE
            DO 40 K = 1,RANK - 1
               IF (B(K).NE.ZERO) CALL DAXPY(RANK-K,-B(K),A(K+1,K),1,
     +                                B(K+1),1)
   40       CONTINUE
         END IF

         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('U','NoTrans','NonUnit',RANK,A,
     +                                LDA,B,1)
         ELSE
            DO 50 K = RANK,2,-1
               IF (B(K).NE.ZERO) THEN
                  B(K) = B(K)/A(K,K)
                  CALL DAXPY(K-1,-B(K),A(1,K),1,B(1),1)
               END IF
   50       CONTINUE
            IF (RANK.GT.0) B(1) = B(1)/A(1,1)
         END IF

         DO 60 K = RANK + 1,N
            B(K) = ZERO
   60    CONTINUE
         DO 70 I = RANK + 1,N
            K = -IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   70    CONTINUE

      ELSE

         DO 80 I = N,RANK + 1,-1
            K = -IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   80    CONTINUE
         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('U','Trans','NonUnit',RANK,A,LDA,
     +                                B,1)
         ELSE
            IF (RANK.GT.0) B(1) = B(1)/A(1,1)
            DO 90 I = 2,RANK
               TEMP = B(I) - DDOT(I-1,A(1,I),1,B(1),1)
               B(I) = TEMP/A(I,I)
   90       CONTINUE
         END IF

         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('L','Trans','Unit',RANK,A,LDA,B,1)
         ELSE
            DO 100 I = RANK - 1,1,-1
               B(I) = B(I) - DDOT(RANK-I,A(I+1,I),1,B(I+1),1)
  100       CONTINUE
         END IF

         DO 110 I = RANK + 1,M
            B(I) = ZERO
  110    CONTINUE
         DO 120 I = RANK,1,-1
            K = IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
  120    CONTINUE
      END IF

      END

      SUBROUTINE MA50ID(CNTL,ICNTL)

      DOUBLE PRECISION CNTL(10)
      INTEGER I,ICNTL(20)

      CNTL(1) = 0.5D0
      CNTL(2) = 0.1D0
      DO 10 I = 3,10
         CNTL(I) = 0.0D0
   10 CONTINUE

      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 1
      ICNTL(4) = 3
      ICNTL(5) = 32
      DO 20 I = 6,20
        ICNTL(I) = 0
   20 CONTINUE

      END

      SUBROUTINE MC13DD(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW)
      INTEGER LICN,N,NUM
      INTEGER IB(N),ICN(LICN),IOR(N),IP(N),IW(N,3),LENR(N)
      EXTERNAL MC13ED
      CALL MC13ED(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW(1,1),IW(1,2),IW(1,3))
      RETURN

      END
      SUBROUTINE MC13ED(N,ICN,LICN,IP,LENR,ARP,IB,NUM,LOWL,NUMB,PREV)
      INTEGER LICN,N,NUM
      INTEGER ARP(N),IB(N),ICN(LICN),IP(N),LENR(N),LOWL(N),NUMB(N),
     +        PREV(N)
      INTEGER DUMMY,I,I1,I2,ICNT,II,ISN,IST,IST1,IV,IW,J,K,LCNT,NNM1,STP
      INTRINSIC MIN
      ICNT = 0
      NUM = 0
      NNM1 = N + N - 1
      DO 20 J = 1,N
        NUMB(J) = 0
        ARP(J) = LENR(J) - 1
   20 CONTINUE
      DO 120 ISN = 1,N
        IF (NUMB(ISN).NE.0) GO TO 120
        IV = ISN
        IST = 1
        LOWL(IV) = 1
        NUMB(IV) = 1
        IB(N) = IV
        DO 110 DUMMY = 1,NNM1
          I1 = ARP(IV)
          IF (I1.LT.0) GO TO 60
          I2 = IP(IV) + LENR(IV) - 1
          I1 = I2 - I1
          DO 50 II = I1,I2
            IW = ICN(II)
            IF (NUMB(IW).EQ.0) GO TO 100
            LOWL(IV) = MIN(LOWL(IV),LOWL(IW))
   50     CONTINUE
          ARP(IV) = -1
   60     IF (LOWL(IV).LT.NUMB(IV)) GO TO 90
          NUM = NUM + 1
          IST1 = N + 1 - IST
          LCNT = ICNT + 1
          DO 70 STP = IST1,N
            IW = IB(STP)
            LOWL(IW) = N + 1
            ICNT = ICNT + 1
            NUMB(IW) = ICNT
            IF (IW.EQ.IV) GO TO 80
   70     CONTINUE
   80     IST = N - STP
          IB(NUM) = LCNT
          IF (IST.NE.0) GO TO 90
          IF (ICNT.LT.N) GO TO 120
          GO TO 130
   90     IW = IV
          IV = PREV(IV)
          LOWL(IV) = MIN(LOWL(IV),LOWL(IW))
          GO TO 110
  100     ARP(IV) = I2 - II - 1
          PREV(IW) = IV
          IV = IW
          IST = IST + 1
          LOWL(IV) = IST
          NUMB(IV) = IST
          K = N + 1 - IST
          IB(K) = IV
  110   CONTINUE
  120 CONTINUE
  130 DO 140 I = 1,N
        II = NUMB(I)
        ARP(II) = I
  140 CONTINUE
      RETURN

      END

      SUBROUTINE MC21AD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
      INTEGER LICN,N,NUMNZ
      INTEGER ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
      EXTERNAL MC21BD
      CALL MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),
     +            IW(1,3),IW(1,4))
      RETURN
      END
      SUBROUTINE MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
      INTEGER LICN,N,NUMNZ
      INTEGER ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
      INTEGER I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
      DO 10 I = 1,N
        ARP(I) = LENR(I) - 1
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      NUMNZ = 0
      DO 100 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENR(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = ICN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
          ARP(J) = -1
   30     CONTINUE
          OUT(J) = LENR(J) - 1
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
            DO 40 II = IN1,IN2
              I = ICN(II)
              IF (CV(I).EQ.JORD) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 100
   60     CONTINUE
   70   CONTINUE
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUMNZ = NUMNZ + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 100
          II = IP(J) + LENR(J) - OUT(J) - 2
          I = ICN(II)
          IPERM(I) = J
   90   CONTINUE
  100 CONTINUE
      IF (NUMNZ.EQ.N) RETURN
      DO 110 I = 1,N
        ARP(I) = 0
  110 CONTINUE
      K = 0
      DO 130 I = 1,N
        IF (IPERM(I).NE.0) GO TO 120
        K = K + 1
        OUT(K) = I
        GO TO 130
  120   CONTINUE
        J = IPERM(I)
        ARP(J) = I
  130 CONTINUE
      K = 0
      DO 140 I = 1,N
        IF (ARP(I).NE.0) GO TO 140
        K = K + 1
        IOUTK = OUT(K)
        IPERM(IOUTK) = I
  140 CONTINUE
      RETURN
      END

      SUBROUTINE MC29AD(M,N,NE,A,IRN,ICN,R,C,W,LP,IFAIL)
      INTEGER M,N,NE
      DOUBLE PRECISION A(NE)
      INTEGER IRN(NE),ICN(NE)
      DOUBLE PRECISION R(M),C(N),W(M*2+N*3)
      INTEGER LP,IFAIL

      INTRINSIC LOG,ABS,MIN

      INTEGER MAXIT
      PARAMETER (MAXIT=100)
      DOUBLE PRECISION ONE,SMIN,ZERO
      PARAMETER (ONE=1D0,SMIN=0.1,ZERO=0D0)

      INTEGER I,I1,I2,I3,I4,I5,ITER,J,K
      DOUBLE PRECISION E,E1,EM,Q,Q1,QM,S,S1,SM,U,V

      IFAIL = 0
      IF (M.LT.1 .OR. N.LT.1) THEN
         IFAIL = -1
         GO TO 220
      ELSE IF (NE.LE.0) THEN
         IFAIL = -2
         GO TO 220
      END IF

      I1 = 0
      I2 = M
      I3 = M + N
      I4 = M + N*2
      I5 = M + N*3

      DO 10 I = 1,M
         R(I) = ZERO
         W(I1+I) = ZERO
   10 CONTINUE
      DO 20 J = 1,N
         C(J) = ZERO
         W(I2+J) = ZERO
         W(I3+J) = ZERO
         W(I4+J) = ZERO
   20 CONTINUE

      DO 30 K = 1,NE
         U = ABS(A(K))
         IF (U.EQ.ZERO) GO TO 30
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 30
         U = LOG(U)
         W(I1+I) = W(I1+I) + ONE
         W(I2+J) = W(I2+J) + ONE
         R(I) = R(I) + U
         W(I3+J) = W(I3+J) + U
   30 CONTINUE
      DO 40 I = 1,M
         IF (W(I1+I).EQ.ZERO) W(I1+I) = ONE
         R(I) = R(I)/W(I1+I)
         W(I5+I) = R(I)
   40 CONTINUE
      DO 50 J = 1,N
         IF (W(I2+J).EQ.ZERO) W(I2+J) = ONE
         W(I3+J) = W(I3+J)/W(I2+J)
   50 CONTINUE
      SM = SMIN*NE

      DO 60 K = 1,NE
         IF (A(K).EQ.ZERO) GO TO 60
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 60
         R(I) = R(I) - W(I3+J)/W(I1+I)
   60 CONTINUE
      E = ZERO
      Q = ONE
      S = ZERO
      DO 70 I = 1,M
         S = S + W(I1+I)*R(I)**2
   70 CONTINUE
      IF (S.LE.SM) GO TO 160

      DO 150 ITER = 1,MAXIT
         DO 80 K = 1,NE
            IF (A(K).EQ.ZERO) GO TO 80
            J = ICN(K)
            I = IRN(K)
            IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 80
            C(J) = C(J) + R(I)
   80    CONTINUE
         S1 = S
         S = ZERO
         DO 90 J = 1,N
            V = -C(J)/Q
            C(J) = V/W(I2+J)
            S = S + V*C(J)
   90    CONTINUE
         E1 = E
         E = Q*S/S1
         Q = ONE - E
         IF (S.LE.SM) E = ZERO
         DO 100 I = 1,M
            R(I) = R(I)*E*W(I1+I)
  100    CONTINUE
         IF (S.LE.SM) GO TO 180
         EM = E*E1
         DO 110 K = 1,NE
            IF (A(K).EQ.ZERO) GO TO 110
            I = IRN(K)
            J = ICN(K)
            IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 110
            R(I) = R(I) + C(J)
  110    CONTINUE
         S1 = S
         S = ZERO
         DO 120 I = 1,M
            V = -R(I)/Q
            R(I) = V/W(I1+I)
            S = S + V*R(I)
  120    CONTINUE
         E1 = E
         E = Q*S/S1
         Q1 = Q
         Q = ONE - E
         IF (S.LE.SM) Q = ONE
         QM = Q*Q1
         DO 130 J = 1,N
            W(I4+J) = (EM*W(I4+J)+C(J))/QM
            W(I3+J) = W(I3+J) + W(I4+J)
  130    CONTINUE
         IF (S.LE.SM) GO TO 160
         DO 140 J = 1,N
            C(J) = C(J)*E*W(I2+J)
  140    CONTINUE
  150 CONTINUE
  160 DO 170 I = 1,M
         R(I) = R(I)*W(I1+I)
  170 CONTINUE
  180 DO 190 K = 1,NE
         IF (A(K).EQ.ZERO) GO TO 190
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 190
         R(I) = R(I) + W(I3+J)
  190 CONTINUE
      DO 200 I = 1,M
         R(I) = R(I)/W(I1+I) - W(I5+I)
  200 CONTINUE
      DO 210 J = 1,N
         C(J) = -W(I3+J)
  210 CONTINUE
      RETURN

  220 IF (LP.GT.0) WRITE (LP,'(/A/A,I3)')
     +    ' **** Error return from MC29AD ****',' IFAIL =',IFAIL

      END


      SUBROUTINE MC71AD(N,KASE,X,EST,W,IW,KEEP)
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION EST
      INTEGER KASE,N
      DOUBLE PRECISION W(*),X(*)
      INTEGER IW(*),KEEP(5)
      DOUBLE PRECISION ALTSGN,TEMP
      INTEGER I,ITER,J,JLAST,JUMP
      INTEGER IDAMAX
      EXTERNAL IDAMAX
      INTRINSIC ABS,SIGN,NINT,DBLE
      IF (N.LE.0) THEN
        KASE = -1
        RETURN

      END IF

      IF (KASE.EQ.0) THEN
        DO 10 I = 1,N
          X(I) = ONE/DBLE(N)
   10   CONTINUE
        KASE = 1
        JUMP = 1
        KEEP(1) = JUMP
        KEEP(2) = 0
        KEEP(3) = 0
        KEEP(4) = 0
        RETURN

      END IF
      JUMP  = KEEP(1)
      ITER  = KEEP(2)
      J     = KEEP(3)
      JLAST = KEEP(4)
      GO TO (100,200,300,400,500) JUMP
  100 CONTINUE
      IF (N.EQ.1) THEN
        W(1) = X(1)
        EST = ABS(W(1))
        GO TO 510

      END IF
      DO 110 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  110 CONTINUE
      KASE = 2
      JUMP = 2
      GO TO 1010
  200 CONTINUE
      J = IDAMAX(N,X,1)
      ITER = 2
  220 CONTINUE
      DO 230 I = 1,N
        X(I) = ZERO
  230 CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      GO TO 1010
  300 CONTINUE
      DO 310 I = 1,N
        W(I) = X(I)
  310 CONTINUE
      DO 320 I = 1,N
        IF (NINT(SIGN(ONE,X(I))).NE.IW(I)) GO TO 330
  320 CONTINUE
      GO TO 410
  330 CONTINUE
      DO 340 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  340 CONTINUE
      KASE = 2
      JUMP = 4
      GO TO 1010
  400 CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF ((ABS(X(JLAST)).NE.ABS(X(J))) .AND. (ITER.LT.ITMAX)) THEN
        ITER = ITER + 1
        GO TO 220

      END IF
  410 CONTINUE
      EST = ZERO
      DO 420 I = 1,N
        EST = EST + ABS(W(I))
  420 CONTINUE
      ALTSGN = ONE
      DO 430 I = 1,N
        X(I) = ALTSGN* (ONE+DBLE(I-1)/DBLE(N-1))
        ALTSGN = -ALTSGN
  430 CONTINUE
      KASE = 1
      JUMP = 5
      GO TO 1010
  500 CONTINUE
      TEMP = ZERO
      DO 520 I = 1,N
        TEMP = TEMP + ABS(X(I))
  520 CONTINUE
      TEMP = 2.0*TEMP/DBLE(3*N)
      IF (TEMP.GT.EST) THEN
        DO 530 I = 1,N
          W(I) = X(I)
  530   CONTINUE
        EST = TEMP
      END IF
  510 KASE = 0
 1010 CONTINUE
      KEEP(1) = JUMP
      KEEP(2) = ITER
      KEEP(3) = J
      KEEP(4) = JLAST
      RETURN
      END
      SUBROUTINE MA48AD(M,N,NE,JOB,LA,A,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,
     +                  RINFO)
      INTEGER M,N,NE,JOB,LA
      DOUBLE PRECISION A(LA)
      INTEGER IRN(LA),JCN(LA),KEEP(*)
      DOUBLE PRECISION CNTL(10)
      DOUBLE PRECISION RINFO(10)
      INTEGER ICNTL(20),IW(6*M+3*N),INFO(20)


      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.D0)
      DOUBLE PRECISION CNTL5(10)
      INTEGER EYE,HEADC,I,IB,ICNTL5(20),IDUMMY,INFO5(15),IP,IPTRA,IPTRD,
     +        IPTRO,IPTRP,IQ,ISW,IW13,IW21,IW50,J,JAY,JB,JFIRST,J1,J2,
     +        J3,K,KB,KBLOCK,KD,KK,KL,KO,L,LASTR,LASTC,LBLOCK
      LOGICAL LDUP
      INTEGER LENC,LENP,LENR,LP,MBLOCK,MINBLK,MP,NB,NBLOCK,NC,NDIAG,
     +        NEXTC,NEXTR,NP,NR,NXTEYE,NZA,NZB,NZD,PTRD,PTRO
      DOUBLE PRECISION RINFO5(10),TOL
      EXTERNAL MC13DD,MC21AD,MA50AD
      INTRINSIC ABS,MAX
      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (M.LE.0 .OR. N.LE.0) THEN
         INFO(1) = -1
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9100) M,N
         GO TO 530
      END IF
      IF (NE.LE.0) THEN
         INFO(1) = -2
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9110) NE
         GO TO 530
      END IF
      IF (LA.LT.2*NE) THEN
         INFO(1) = -3
         INFO(3) = 2*NE
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9120) LA,INFO(3)
         GO TO 530
      END IF
      IF (JOB.LT.1 .OR. JOB.GT.3) THEN
         INFO(1) = -6
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9130) JOB
         GO TO 530
      END IF
      IF (JOB.EQ.2) THEN
         DO 10 I = 1,MAX(M,N)
            IW(N+I) = 0
   10    CONTINUE
         DO 20 I = 1,M
            J = KEEP(I)
            IF (J.LT.1 .OR. J.GT.M) GO TO 40
            IF (IW(N+J).EQ.1) GO TO 40
            IW(N+J) = 1
   20    CONTINUE
         DO 30 I = 1,N
            J = KEEP(M+I)
            IF (J.LT.1 .OR. J.GT.N) GO TO 40
            IF (IW(N+J).EQ.2) GO TO 40
            IW(N+J) = 2
   30    CONTINUE
         GO TO 50
   40    INFO(1) = -5
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9140)
         GO TO 530
      END IF

   50 IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,
     +     '(/A/A,I7,A,I6,A,I7,A,I2,A,I7/A,1P,4D12.4/A,4I8/A,3I8)')
     +     ' Entering MA48A/AD with',' M =',M,'     N =',N,'     NE =',
     +     NE,'     JOB =',JOB,'     LA =',LA,' CNTL (1:4) =',
     +     (CNTL(I),I=1,4),' ICNTL(1:4) = ', (ICNTL(I),I=1,4),
     +     ' ICNTL(6:8) = ', (ICNTL(I),I=6,8)
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9000) (A(K),IRN(K),JCN(K),K=1,NE)
 9000       FORMAT (' Entries:'/3 (1P,D12.4,2I6))
         ELSE
            WRITE (MP,9000) (A(K),IRN(K),JCN(K),K=1,MIN(9,NE))
         END IF
         IF (JOB.EQ.2) THEN
            WRITE (MP,'(A)') ' Permutations input (JOB=2)'
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9010) (KEEP(I),I=1,M)
 9010          FORMAT (' Positions of original rows in the permuted ma',
     +                'trix'/ (10I6))
               WRITE (MP,9020) (KEEP(M+I),I=1,N)
 9020          FORMAT (' Positions of columns of permuted matrix ','in',
     +                ' or','iginal matrix '/ (10I6))
            ELSE
               WRITE (MP,9010) (KEEP(I),I=1,MIN(10,M))
               WRITE (MP,9020) (KEEP(M+I),I=1,MIN(10,N))
            END IF
         END IF
         IF (ICNTL(8).NE.0) THEN
            WRITE (MP,'(A,I6)')
     +        ' Value of IW entries on call with ICNTL(8) =',ICNTL(8)
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9030) (IW(I),I=1,N)
 9030          FORMAT (10I6)
            ELSE
               WRITE (MP,9030) (IW(I),I=1,MIN(10,N))
            END IF
         END IF
      END IF
      DO 53 I = 1,20
         INFO(I) = 0
   53 CONTINUE
      INFO(3) = NE*2
      INFO(4) = NE
      INFO(10) = MIN(M,N)
      DO 56 I = 1,5
        RINFO(I) = ZERO
   56 CONTINUE
      DO 60 I = 1,4
         CNTL5(I) = CNTL(I)
         ICNTL5(I) = ICNTL(I)
   60 CONTINUE
      ICNTL5(3) = 0
      ICNTL5(5) = ICNTL(5)
      ICNTL5(6) = 0
      ICNTL5(7) = 0
      IF (JOB.EQ.2) THEN
         ICNTL5(7) = 2
         ICNTL5(4) = 1
      END IF
      IF (JOB.EQ.3) ICNTL5(7) = 1

      TOL = MAX(ZERO,CNTL(4))
      MINBLK = MAX(1,ICNTL(6))
      LDUP = .FALSE.

      IPTRD = M + 3*N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1

      HEADC = N + 1
      LASTC = HEADC + N
      DO 70 J = 1,N
         IW(HEADC+J) = 0
   70 CONTINUE
      DO 80 K = 1,NE
         I = IRN(K)
         J = JCN(K)
         IF (I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N) THEN
            INFO(12) = INFO(12) + 1
            IF (MP.GT.0 .AND. INFO(12).LE.10 .AND.
     +          ICNTL(3).GE.2) WRITE (MP,'(A,I7,A,2I6)')
     +          ' Message from MA48A/AD .. indices for entry ',K,' are',
     +          I,J
            JCN(K) = 0
         ELSE
            JCN(K) = IW(HEADC+J)
            IW(HEADC+J) = K
         END IF
   80 CONTINUE

      IF (MINBLK.GE.N .OR. M.NE.N .OR. JOB.GT.1) GO TO 190

      IW21 = 2*N
      IW13 = IW21
      IPTRA = LASTC + N
      LENC = IPTRA + N
      IB = LENC
      IP = LENC + N
      IPTRP = IP + N
      LENP = IPTRP + N

      DO 90 I = 1,N
         IW(LASTC+I) = 0
   90 CONTINUE

      LDUP = .TRUE.
      K = 1
      DO 120 J = 1,N
         EYE = IW(HEADC+J)
         IW(IPTRA+J) = K
         DO 100 IDUMMY = 1,NE
            IF (EYE.EQ.0) GO TO 110
            I = IRN(EYE)
            IF (IW(LASTC+I).NE.J) THEN
               IW(LASTC+I) = J
               IRN(NE+K) = I
               K = K + 1
            ELSE
               INFO(11) = INFO(11) + 1
               IF (MP.GT.0 .AND. INFO(11).LE.10 .AND.
     +             ICNTL(3).GE.2) WRITE (MP,'(A,I7,A,2I6)')
     +             ' Message from MA48A/AD .. duplicate in position ',K,
     +             ' with indices',I,J
            END IF
            EYE = JCN(EYE)
  100    CONTINUE
  110    IW(LENC+J) = K - IW(IPTRA+J)
  120 CONTINUE

      NZA = K - 1

      CALL MC21AD(N,IRN(NE+1),NZA,IW(IPTRA+1),IW(LENC+1),IW(IP+1),NDIAG,
     +           KEEP(IW21+1))
      INFO(10) = NDIAG
      IF (NDIAG.LT.N) THEN
         IF (ICNTL(7).NE.0) THEN
            INFO(1) = -4
            IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,
     +          '(A,A/A,I7,A,I7)')
     +        ' Error return from MA48A/AD because matrix structurally '
     +          ,' singular',' order is ',N,' and structural rank',NDIAG
            GO TO 530
         END IF
         GO TO 190
      END IF
      DO 130 J = 1,N
         JAY = IW(IP+J)
         IW(IPTRP+J) = IW(IPTRA+JAY)
         IW(LENP+J) = IW(LENC+JAY)
  130 CONTINUE

      CALL MC13DD(N,IRN(NE+1),NZA,IW(IPTRP+1),IW(LENP+1),KEEP(M+1),
     +           IW(IB+1),NB,KEEP(IW13+1))

      DO 140 JB = 2,NB
         IW(IB+JB-1) = IW(IB+JB) - IW(IB+JB-1)
  140 CONTINUE
      IW(IB+NB) = N + 1 - IW(IB+NB)
      IF (IW(IB+1).EQ.1) IW(IB+1) = -1
      KB = 1
      DO 150 JB = 2,NB
         L = IW(IB+JB)
         IF (L.EQ.1 .AND. IW(IB+KB).LE.0) THEN
            IW(IB+KB) = IW(IB+KB) - 1
         ELSE
            KB = KB + 1
            IF (L.EQ.1) THEN
               IW(IB+KB) = -1
            ELSE
               IW(IB+KB) = L
            END IF
         END IF
  150 CONTINUE
      NB = KB
      KB = 1
      DO 160 JB = 2,NB
         IF (ABS(IW(IB+KB)).LT.MINBLK) THEN
            IW(IB+KB) = ABS(IW(IB+KB)) + ABS(IW(IB+JB))
         ELSE
            KB = KB + 1
            IW(IB+KB) = IW(IB+JB)
         END IF
  160 CONTINUE
      NB = KB
      DO 170 JB = 1,NB
         KEEP(NBLOCK+3*JB) = ABS(IW(IB+JB))
         KEEP(MBLOCK+3*JB) = IW(IB+JB)
  170 CONTINUE
      DO 180 J = 1,N
         KEEP(KEEP(M+J)) = J
         KEEP(M+J) = IW(IP+KEEP(M+J))
  180 CONTINUE
      GO TO 220
  190 NB = 1
      IF (JOB.EQ.1 .OR. JOB.EQ.3) THEN
         DO 200 I = 1,M
            KEEP(I) = I
  200    CONTINUE
         DO 210 I = 1,N
            KEEP(M+I) = I
  210    CONTINUE
      END IF
      KEEP(NBLOCK+3) = N
      KEEP(MBLOCK+3) = 0

  220 IF (ICNTL(8).NE.0) THEN
         LBLOCK = KBLOCK + 3*NB
         IF (JOB.EQ.2) THEN
            DO 230 I = 1,N
               IF (IW(I).EQ.0) GO TO 240
  230       CONTINUE
  240       KEEP(LBLOCK+1) = N - I + 1
         ELSE
            J = 1
            DO 270 JB = 1,NB
               KEEP(LBLOCK+JB) = 0
               J2 = J + KEEP(NBLOCK+3*JB) - 1
               J1 = J2
               IF (KEEP(MBLOCK+3*JB).LT.0) GO TO 260
  250          IF (J.EQ.J2) GO TO 260
               IF (IW(KEEP(M+J)).EQ.0) THEN
                  KEEP(LBLOCK+JB) = KEEP(LBLOCK+JB) + 1
                  ISW = KEEP(M+J2)
                  KEEP(M+J2) = KEEP(M+J)
                  KEEP(M+J) = ISW
                  J2 = J2 - 1
               ELSE
                  J = J + 1
               END IF
               GO TO 250
  260          J = J1 + 1
  270       CONTINUE
         END IF
      END IF

      LASTR = LASTC + M
      DO 280 I = 1,M
         IW(LASTC+I) = 0
  280 CONTINUE
      KEEP(KBLOCK+3) = NB
      K = 1
      KK = NE
      J2 = 0
      DO 310 JB = 1,NB
         J1 = J2 + 1
         J2 = J1 + KEEP(NBLOCK+3*JB) - 1
         DO 300 JAY = J1,J2
            J = KEEP(M+JAY)
            EYE = IW(HEADC+J)
            KEEP(IPTRD+JAY) = K
            KEEP(IPTRO+JAY) = KK
            IF (KEEP(MBLOCK+3*JB).LT.0) THEN
               IW(LASTC+JAY) = JAY
               A(NE+K) = ZERO
               IRN(NE+K) = JAY
               IW(LASTR+JAY) = K
               K = K + 1
            END IF
            DO 290 IDUMMY = 1,NE
               IF (EYE.EQ.0) GO TO 300
               NXTEYE = JCN(EYE)
               I = KEEP(IRN(EYE))
               IF (IW(LASTC+I).NE.JAY) THEN
                  IW(LASTC+I) = JAY
                  IF ((I.GE.J1.AND.I.LE.J2) .OR. (M.NE.N)) THEN
                     A(NE+K) = A(EYE)
                     IRN(NE+K) = I
                     IW(LASTR+I) = K
                     JCN(EYE) = K
                     K = K + 1
                  ELSE
                     A(NE+KK) = A(EYE)
                     IRN(NE+KK) = I
                     IW(LASTR+I) = KK
                     JCN(EYE) = KK
                     KK = KK - 1
                  END IF
               ELSE
                  IF (.NOT.LDUP) THEN
                     INFO(11) = INFO(11) + 1
                     IF (MP.GT.0 .AND. INFO(11).LE.10 .AND.
     +                   ICNTL(3).GE.2) WRITE (MP,'(A,I7,A,2I6)')
     +                ' Message from MA48A/AD .. duplicate in position '
     +                   ,EYE,' with indices',IRN(EYE),J
                  END IF
                  KL = IW(LASTR+I)
                  JCN(EYE) = KL
                  A(NE+KL) = A(NE+KL) + A(EYE)
               END IF
               EYE = NXTEYE
  290       CONTINUE
  300    CONTINUE
  310 CONTINUE
      KEEP(IPTRD+N+1) = K
      KEEP(IPTRO+N+1) = KK
      NZD = K - 1
      DO 320 I = 1,K-1
         IRN(I) = IRN(NE+I)
  320 CONTINUE
      DO 325 I = KK+1,NE
         IRN(I) = IRN(NE+I)
  325 CONTINUE
      DO 326 I = K,KK
         IRN(I) = 0
  326 CONTINUE
      IP = 0
      IQ = M
      PTRD = M + N
      JFIRST = M + 2*N
      LENR = 2*M + 2*N
      LASTR = 3*M + 2*N
      NEXTR = 4*M + 2*N
      IW50 = 5*M + 2*N
      NEXTC = 6*M + 2*N
      PTRO = NEXTC
      LASTC = M + N
      LENC = LASTC + N
      J1 = N + 1
      KB = 0
      DO 390 JB = NB,1,-1
         NC = KEEP(NBLOCK+3*JB)
         J2 = J1 - 1
         J1 = J2 + 1 - NC
         IF (KEEP(MBLOCK+3*JB).LT.0) THEN
            DO 330 J = J1,J2
               IF (ABS(A(NE+KEEP(IPTRD+J))).GT.TOL)
     +             INFO(5) = INFO(5) + 1
               IW(IP+J) = J
               IW(IQ+J) = J
  330       CONTINUE
         ELSE
            NZB = KEEP(IPTRD+J2+1) - KEEP(IPTRD+J1)
            DO 340 K = KEEP(IPTRD+J1),KEEP(IPTRD+J2+1) - 1
               IRN(NE+K) = IRN(K) - J1 + 1
  340       CONTINUE
            K = KEEP(IPTRD+J1) - 1
            DO 350 J = J1,J2
               IW(IQ+J) = KEEP(IPTRD+J) - K
  350       CONTINUE
            NR = NC
            IF (NB.EQ.1) NR = M
            IF (JOB.EQ.2) THEN
               DO 360 J = J1,J1 + NR - 1
                  IW(IP+J) = J - J1 + 1
  360          CONTINUE
               DO 370 J = J1,J2
                  IW(PTRD+J) = J - J1 + 1
  370          CONTINUE
            END IF
            INFO(7) = MAX(INFO(7),NR)
            INFO(8) = INFO(8) + NC
            INFO(9) = INFO(9) + NZB
            IF (ICNTL(8).NE.0) ICNTL5(6) = KEEP(LBLOCK+JB)
            CALL MA50AD(NR,NC,NZB,LA-NE-K,A(NE+K+1),IRN(NE+K+1),
     +                  JCN(NE+1),IW(IQ+J1),CNTL5,ICNTL5,IW(IP+J1),
     +                  NP,IW(JFIRST+1),IW(LENR+1),
     +                  IW(LASTR+1),IW(NEXTR+1),IW(IW50+1),IW(PTRD+J1),
     +                  KEEP(LENC+KB+1),KEEP(LASTC+KB+1),IW(NEXTC+1),
     +                  INFO5,RINFO5)
            KEEP(MBLOCK+3*JB) = NP
            DO 380 J = J1,J1+NR-1
               IW(IP+J) = IW(IP+J) + J1 - 1
  380       CONTINUE
            DO 385 J = J1,J2
               IW(IQ+J) = IW(IQ+J) + J1 - 1
  385       CONTINUE
            IF (INFO5(1).EQ.1) THEN
               IF (INFO(1).EQ.0 .OR. INFO(1).EQ.4) INFO(1) = INFO(1) + 2
            END IF
            IF (INFO5(1).EQ.2) THEN
               IF (INFO(1).EQ.0 .OR. INFO(1).EQ.2) INFO(1) = INFO(1) + 4
            END IF
            IF (INFO5(1).EQ.3 .AND. INFO(1).GE.0) INFO(1) = 6
            IF (INFO5(1).EQ.-3) INFO(1) = -3
            INFO(2) = INFO(2) + INFO5(2)
            INFO(3) = MAX(INFO(3),NE+K+INFO5(3))

            INFO(5) = INFO(5) + INFO5(5)
            INFO(6) = INFO(6) + INFO5(6)
            RINFO(1) = RINFO(1) + RINFO5(1)
            KB = KB + 1
            KEEP(LENC+KB) = INFO5(4) - 2*NR
            KEEP(LASTC+KB) = INFO5(4)
         END IF
  390 CONTINUE

      INFO(4) = NE*2
      K = NE
      DO 400 JB = KB,1,-1
         INFO(4) = MAX(INFO(4),K+KEEP(LASTC+JB))
         K = K + KEEP(LENC+JB)
  400 CONTINUE

      IF (INFO(1).EQ.-3) THEN
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9120) LA,INFO(3)
         GO TO 530
      END IF

      DO 410 K = 1,NE
         IRN(NE+K) = IRN(K)
  410 CONTINUE
      DO 420 J = 1,N
         IW(PTRD+J) = KEEP(IPTRD+J)
         IW(PTRO+J) = KEEP(IPTRO+J+1) + 1
  420 CONTINUE
      KD = 1
      KO = NZD + 1
      DO 450 J = 1,N
         KEEP(IPTRD+J) = KD
         JAY = IW(IQ+J)
         KL = NZD
         IF (JAY.NE.N) KL = IW(PTRD+JAY+1) - 1
         DO 430 KK = IW(PTRD+JAY),KL
            IRN(KD) = IW(IP+IRN(NE+KK))
            IRN(NE+KK) = KD
            KD = KD + 1
  430    CONTINUE
         KEEP(IPTRO+J) = KO
         KL = NE
         IF (JAY.NE.1) KL = IW(PTRO+JAY-1) - 1
         DO 440 KK = IW(PTRO+JAY),KL
            IRN(KO) = IW(IP+IRN(NE+KK))
            IRN(NE+KK) = KO
            KO = KO + 1
  440    CONTINUE
  450 CONTINUE
      KEEP(IPTRO+N+1) = KO
      DO 460 I = 1,M
         KEEP(I) = IW(IP+KEEP(I))
  460 CONTINUE
      DO 465 I = 1,N
         IW(IQ+I) = KEEP(M+IW(IQ+I))
  465 CONTINUE
      DO 470 I = 1,N
         KEEP(M+I) = IW(IQ+I)
  470 CONTINUE
      IW(1) = IRN(NE)
      IRN(NE) = NE
      DO 480 K = 1,NE
         JCN(K) = IRN(NE+JCN(K))
  480 CONTINUE
      IRN(NE) = IW(1)

      IF (INFO(11).GT.0 .OR. INFO(12).GT.0) THEN
         INFO(1) = INFO(1) + 1
         IF (MP.GT.0 .AND. ICNTL(3).GE.2) THEN
            IF (INFO(11).GT.0) WRITE (MP,9150) INFO(11)
            IF (INFO(12).GT.0) WRITE (MP,9160) INFO(12)
         END IF
         IF (INFO(11).GT.0) JCN(1) = -JCN(1)
      END IF

      IF (INFO(10).LT.INFO(5)) THEN
         INFO(5) = INFO(10)
         IF (INFO(1).NE.2 .AND. INFO(1).NE.3 .AND.
     +       INFO(1).LT.6) INFO(1) = INFO(1) + 2
      END IF

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/A/A/12I6/A,F12.1)') ' Leaving MA48A/AD with',
     +     ' INFO(1:12)  =',(INFO(I),I=1,12),' RINFO(1) =',RINFO(1)
         WRITE (MP,'(A)') ' Permuted matrix by blocks (IRN)'
         KB = NB
         IF (ICNTL(3).EQ.3) THEN
            WRITE (MP,'(A)')
     +        ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         END IF
         WRITE (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         DO 500 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            IF (J1.LE.J2) WRITE (MP,'(A,I6)') ' Block',JB
            J3 = J2
            IF (ICNTL(3).EQ.3) J3 = J1
            DO 490 J = J1,J3
               WRITE (MP,'(A,I5,(T13,10I6))') ' Column',J,
     +           (IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
  490       CONTINUE
            J1 = J2 + 1
  500    CONTINUE
         IF (KEEP(IPTRO+N+1).GT.KEEP(IPTRD+N+1)) THEN
            WRITE (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            DO 520 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               IF (ICNTL(3).EQ.3) J3 = J1
               DO 510 J = J1,J3
                  IF (KEEP(IPTRO+J+1).GT.KEEP(IPTRO+J)) WRITE (MP,
     +                '(A,I5,(T13,10I6))') ' Column',J,
     +                (IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
  510          CONTINUE
               J1 = J2 + 1
  520       CONTINUE
         END IF
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9040) (JCN(K),K=1,NE)
 9040       FORMAT (' JCN (MAP) ='/ (6X,10I6))
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9010) (KEEP(I),I=1,M)
            WRITE (MP,9020) (KEEP(M+I),I=1,N)
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9050) (KEEP(IPTRD+J),J=1,N+1)
 9050       FORMAT (' IPTRD ='/ (8X,10I6))
            WRITE (MP,9060) (KEEP(IPTRO+J),J=1,N+1)
 9060       FORMAT (' IPTRO ='/ (8X,10I6))
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9070) (KEEP(NBLOCK+3*JB),JB=1,NB)
            WRITE (MP,9080) (KEEP(MBLOCK+3*JB),JB=1,NB)
 9070       FORMAT (' NBLOCK (order blocks) ='/ (8X,10I6))
 9080       FORMAT (' MBLOCK (triangular flag and number packed rows) ='
     +             / (8X,10I6))
 9090       FORMAT (' LBLOCK (number of changed columns) ='/ (8X,10I6))
            IF (ICNTL(8).NE.0) WRITE (MP,9090) (KEEP(LBLOCK+JB),JB=1,NB)
         ELSE
            WRITE (MP,9040) (JCN(K),K=1,MIN(10,NE))
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9010) (KEEP(I),I=1,MIN(10,M))
            WRITE (MP,9020) (KEEP(M+I),I=1,MIN(10,N))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9050) (KEEP(IPTRD+J),J=1,MIN(10,N+1))
            WRITE (MP,9060) (KEEP(IPTRO+J),J=1,MIN(10,N+1))
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9070) (KEEP(NBLOCK+3*JB),JB=1,MIN(10,NB))
            WRITE (MP,9080) (KEEP(MBLOCK+3*JB),JB=1,MIN(10,NB))
            IF (ICNTL(8).NE.0) WRITE (MP,9090) (KEEP(LBLOCK+JB),JB=1,
     +          MIN(10,NB))
         END IF
      END IF
  530 RETURN

 9100 FORMAT (' Error return from MA48A/AD because M =',I10,' and N =',
     +       I10)
 9110 FORMAT (' Error return from MA48A/AD because NE =',I10)
 9120 FORMAT (' Error return from MA48A/AD because LA is',I10/' and ',
     +       'must be at least',I10)
 9130 FORMAT (' Error return from MA48A/AD because ','JOB = ',I10)
 9140 FORMAT (' Error return from MA48A/AD because ','faulty permutati',
     +       'ons input when JOB = 2')
 9150 FORMAT (' Message from MA48A/AD ..',I8,' duplicates found')
 9160 FORMAT (' Message from MA48A/AD ..',I8,' out-of-range indices fo',
     +       'und')
      END


      SUBROUTINE MA48BD(M,N,NE,JOB,LA,A,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,
     +                  INFO,RINFO)
      INTEGER M,N,NE,JOB,LA
      DOUBLE PRECISION A(LA)
      INTEGER IRN(LA),JCN(NE),KEEP(*)
      DOUBLE PRECISION CNTL(10)
      INTEGER ICNTL(20)
      DOUBLE PRECISION W(M)
      INTEGER IW(2*M+2*N),INFO(20)
      DOUBLE PRECISION RINFO(10)


      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.D0)
      DOUBLE PRECISION CNTL5(10)
      INTEGER I,ICNTL5(20),INFO5(15),IPTRD,IPTRL,IPTRO,IPTRU,IQB(1),
     +        ITRY,J,JB,JOB5,J1,J2,J3,K,KB,KBLOCK,KK,LBLOCK,LP,
     +        MBLOCK,MP,NB,NBLOCK,NEWNE,NC,NP,NR,NRF,NZB
      DOUBLE PRECISION RINFO5(10),TOL
      LOGICAL TRISNG

      EXTERNAL MA50BD
      INTRINSIC MAX

      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (M.LE.0 .OR. N.LE.0) THEN
         INFO(1) = -1
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9160) M,N
         GO TO 240
      END IF
      IF (NE.LE.0) THEN
         INFO(1) = -2
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9170) NE
         GO TO 240
      END IF
      IF (LA.LT.2*NE) THEN
         INFO(1) = -3
         INFO(4) = 2*NE
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9180) LA,INFO(4)
         GO TO 240
      END IF
      IF (JOB.LT.1 .OR. JOB.GT.3) THEN
         INFO(1) = -6
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9190) JOB
         GO TO 240
      END IF
      INFO(1) = 0
      DO 10 I = 1,4
         CNTL5(I) = CNTL(I)
         ICNTL5(I) = ICNTL(I)
   10 CONTINUE
      ICNTL5(3) = 0
      ICNTL5(5) = ICNTL(5)
      ICNTL5(6) = 0
      ICNTL5(7) = 0
      ICNTL5(8) = ICNTL(10)

      IPTRL = M + N
      IPTRU = IPTRL + N
      IPTRD = IPTRU + N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1
      NB = KEEP(KBLOCK+3)
      LBLOCK = KBLOCK + 3*NB
      NEWNE = KEEP(IPTRO+N+1) - 1

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/A/3(A,I8),A,I2,A,I8/A,1P,3D12.4/A,3I8/A,I8/A,I8)')
     +     ' Entering MA48B/BD with',' M =',M,'     N =',N,'     NE ='
     +     ,NE,'     JOB =',JOB,'     LA =',LA,' CNTL (2:4) =',
     +     (CNTL(I),I=2,4),' ICNTL(1:3) =', (ICNTL(I),I=1,3),
     +     ' ICNTL(5)   =',ICNTL(5),' ICNTL(8)   =',ICNTL(8)
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9000) (A(K),K=1,NE)
         ELSE
            WRITE (MP,9000) (A(K),K=1,MIN(10,NE))
         END IF
 9000    FORMAT (' A ='/ (4X,1P,5D12.4))
         WRITE (MP,'(A)') ' Indices for permuted matrix by blocks'
         KB = NB
         IF (ICNTL(3).EQ.3) THEN
            WRITE (MP,'(A)')
     +        ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         END IF
         WRITE (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         DO 30 JB = 1,KB
            WRITE (MP,'(A,I6)') ' Block',JB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            J3 = J2
            IF(ICNTL(3).EQ.3) J3 = J1
            DO 20 J = J1,J3
               WRITE (MP,'(A,I5,(T13,10I6))') ' Column',J,
     +           (IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
   20       CONTINUE
            J1 = J2 + 1
   30    CONTINUE
         IF (KEEP(IPTRO+N+1).GT.KEEP(IPTRD+N+1)) THEN
            WRITE (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            DO 50 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               IF(ICNTL(3).EQ.3) J3 = J1
               DO 40 J = J1,J3
                  IF (KEEP(IPTRO+J+1).GT.KEEP(IPTRO+J)) WRITE (MP,
     +                '(A,I5,(T13,10I6))') ' Column',J,
     +                (IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
   40          CONTINUE
               J1 = J2 + 1
   50       CONTINUE
         END IF
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9010) (JCN(K),K=1,NE)
 9010       FORMAT (' JCN (MAP) ='/ (6X,10I6))
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9020) (KEEP(I),I=1,M)
            WRITE (MP,9030) (KEEP(M+I),I=1,N)
 9020       FORMAT (' Positions of original rows in the permuted matrix'
     +             / (10I6))
 9030       FORMAT (' Positions of columns of permuted matrix ','in ',
     +             'original matrix '/ (10I6))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9040) (KEEP(IPTRD+J),J=1,N+1)
 9040       FORMAT (' IPTRD ='/ (8X,10I6))
            WRITE (MP,9050) (KEEP(IPTRO+J),J=1,N+1)
 9050       FORMAT (' IPTRO ='/ (8X,10I6))
            IF (JOB.GT.1) THEN
               WRITE (MP,9060) (KEEP(IPTRL+J),J=1,N)
 9060          FORMAT (' IPTRL ='/ (8X,10I6))
               WRITE (MP,9070) (KEEP(IPTRU+J),J=1,N)
 9070          FORMAT (' IPTRU ='/ (8X,10I6))
            END IF
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,NB)
            WRITE (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,NB)
 9080       FORMAT (' NBLOCK (order blocks) ='/ (8X,10I6))
 9090       FORMAT (' MBLOCK (triangular flag and number packed rows) ='
     +             / (8X,10I6))
 9100       FORMAT (' KBLOCK (position of beginning of block) ='/
     +             (8X,10I6))
 9110       FORMAT (' LBLOCK (number of changed columns) ='/ (8X,10I6))
            IF (JOB.GT.1) THEN
               WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,NB)
               IF (ICNTL(8).NE.0) WRITE (MP,9110) (KEEP(LBLOCK+JB),JB=1,
     +             NB)
            END IF
         ELSE
            WRITE (MP,9010) (JCN(K),K=1,MIN(10,NE))
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9020) (KEEP(I),I=1,MIN(10,M))
            WRITE (MP,9030) (KEEP(M+I),I=1,MIN(10,N))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9040) (KEEP(IPTRD+J),J=1,MIN(10,N+1))
            WRITE (MP,9050) (KEEP(IPTRO+J),J=1,MIN(10,N+1))
            IF (JOB.GT.1) THEN
               WRITE (MP,9060) (KEEP(IPTRL+J),J=1,MIN(10,N))
               WRITE (MP,9070) (KEEP(IPTRU+J),J=1,MIN(10,N))
            END IF
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,MIN(10,NB))
            WRITE (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,MIN(10,NB))
            IF (JOB.GT.1) THEN
               WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,MIN(10,NB))
               IF (ICNTL(8).NE.0) WRITE (MP,9110) (KEEP(LBLOCK+JB),JB=1,
     +             MIN(10,NB))
            END IF
         END IF
      END IF
      INFO(4) = NE
      INFO(5) = 0
      INFO(6) = 0
      RINFO(1) = ZERO
      TOL = MAX(ZERO,CNTL(4))
      IF (JCN(1).GT.0) THEN
         DO 60 K = 1,NE
            A(NE+K) = A(K)
   60    CONTINUE
         DO 70 K = 1,NE
            A(JCN(K)) = A(NE+K)
   70    CONTINUE
      ELSE
         DO 80 K = 1,NE
            A(NE+K) = A(K)
            A(K) = ZERO
   80    CONTINUE
         A(-JCN(1)) = A(NE+1)
         DO 90 K = 2,NE
            KK = JCN(K)
            A(KK) = A(KK) + A(NE+K)
   90    CONTINUE
      END IF
      IQB(1) = 0
      KK = 0
      J2 = 0
      JOB5 = JOB
      DO 150 JB = 1,NB
         NC = KEEP(NBLOCK+3*JB)
         J1 = J2 + 1
         J2 = J1 + NC - 1
         KEEP(KBLOCK+3*JB) = 0
         IF (KEEP(MBLOCK+3*JB).LT.0) THEN
            TRISNG = .FALSE.
            DO 100 J = J1,J2
               IF (ABS(A(KEEP(IPTRD+J))).LE.TOL) TRISNG = .TRUE.
               KEEP(IPTRL+J) = 0
               KEEP(IPTRU+J) = 0
  100       CONTINUE
            IF (.NOT.TRISNG)  THEN
               INFO(5) = INFO(5) + NC
               GO TO 150
            ENDIF
            IF (JOB.EQ.2) THEN
               INFO(1) = -7
               IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,'(A)')
     +              ' Error return from MA48B/BD with JOB=2 because ',
     +              ' the matrix is incompatible with expectations'
               GO TO 240
            ELSE
               KEEP(MBLOCK+3*JB) = NC
            ENDIF
         END IF
       DO 145 ITRY = 1,2
         NR = NC
         IF (NB.EQ.1) NR = M
         NZB = KEEP(IPTRD+J2+1) - KEEP(IPTRD+J1)
         IF (ICNTL(8).NE.0) ICNTL5(6) = KEEP(LBLOCK+JB)
         DO 110 K = KEEP(IPTRD+J1),KEEP(IPTRD+J2+1) - 1
            IRN(K) = IRN(K) - J1 + 1
  110    CONTINUE
         K = KEEP(IPTRD+J1) - 1
         DO 115 J = J1,J1+NC-1
            KEEP(IPTRD+J) = KEEP(IPTRD+J) - K
  115    CONTINUE
         DO 120 J = J1,J1+NR-1
            IW(J) = J - J1 + 1
  120    CONTINUE
         NP = KEEP(MBLOCK+3*JB)
         CALL MA50BD(NR,NC,NZB,JOB5,A(K+1),IRN(K+1),KEEP(IPTRD+J1),
     +               CNTL5,ICNTL5,IW(J1),IQB,NP,LA-NEWNE-KK,
     +               A(NEWNE+KK+1),IRN(NEWNE+KK+1),KEEP(IPTRL+J1),
     +               KEEP(IPTRU+J1),W,IW(M+1),INFO5,RINFO5)
         DO 130 J = J1,J2
            KEEP(IPTRD+J) = KEEP(IPTRD+J) + K
  130    CONTINUE
         DO 140 K = KEEP(IPTRD+J1),KEEP(IPTRD+J2+1) - 1
            IRN(K) = IRN(K) + J1 - 1
  140    CONTINUE
         IF (INFO5(1).EQ.-6) THEN
            INFO(1) = -6
            IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,'(A)')
     +          ' Error return from MA48B/BD with JOB greater than 1',
     +          ' and entries dropped during previous factorization'
            GO TO 240
         END IF
         IF (INFO5(1).LT.-7) THEN
            IF (ICNTL(11).EQ.1 .AND. JOB.EQ.2) THEN
               JOB5 = 1
               IF (LP.GT.0 .AND. ICNTL(3).GE.2) WRITE(LP,'(A,2(A,I4))')
     +          ' Warning from MA48B/BD. Switched from JOB=2 to JOB=1',
     +          ' in block', JB, ' of ',NB
               INFO(1) = INFO(1) + 1
               GO TO 145
            END IF
            INFO(1) = -7
            IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,'(A)')
     +          ' Error return from MA48B/BD with JOB=2 because ',
     +          ' the matrix is incompatible with expectations'
            GO TO 240
         ELSE
            GO TO 147
         END IF
  145  CONTINUE
  147    IF (INFO5(1).EQ.-3) THEN
            INFO(1) = -3
            IF (ICNTL(10).EQ.1) THEN
               KEEP(KBLOCK+3) = NB
               IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE(LP,'(A,2(A,I4))')
     +          ' Error return from MA48B/BD because LA is too small.',
     +          ' In block', JB, ' of ',NB
               GO TO 240
            END IF
         END IF
         IF (INFO(1).EQ.-3) THEN
            INFO(4) = INFO(4) + INFO5(4)
            KK = 0
         ELSE
            INFO(4) = MAX(INFO(4),KK+NEWNE+INFO5(4))
            NRF = IRN(NEWNE+KK+2)
            KEEP(KBLOCK+3*JB) = KK + 1
            KK = KK + KEEP(IPTRL+J2) + MAX((NC-KEEP(MBLOCK+3*JB))*
     +           (NRF), (NC-KEEP(MBLOCK+3*JB))+(NRF))
         END IF
         IF (INFO5(1).EQ.1) THEN
            IF (INFO(1).NE.-3) INFO(1) = MIN(INFO(1)+2,3)
         END IF
         RINFO(1) = RINFO(1) + RINFO5(1)
         INFO(5) = INFO(5) + INFO5(5)
         INFO(6) = INFO(6) + INFO5(6)
  150 CONTINUE

      INFO(4) = MAX(NE*2,INFO(4))
      KEEP(KBLOCK+3) = NB

      IF (INFO(1).EQ.-3) THEN
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9180) LA,INFO(4)
         GO TO 240
      END IF

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/A/A,I6/A,3I6/A,F12.1)') ' Leaving MA48B/BD with',
     +     ' INFO(1)   = ',INFO(1),' INFO(4:6) = ', (INFO(I),I=4,6),
     +     ' RINFO(1)     =',RINFO(1)
         WRITE (MP,'(A)') ' Permuted matrix by blocks'
         KB = NB
         IF (ICNTL(3).EQ.3) THEN
            WRITE (MP,'(A)')
     +        ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         END IF
         WRITE (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         DO 170 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            IF (J1.LE.J2) WRITE (MP,'(A,I6)') ' Block',JB
            J3 = J2
            IF (ICNTL(3).EQ.3) J3 = MIN(J1,J2)
            DO 160 J = J1,J3
               WRITE (MP,'(A,I5,(T13,3(1PD12.4,I5)))') ' Column',J,
     +           (A(I),IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
  160       CONTINUE
            J1 = J2 + 1
  170    CONTINUE
         IF (KEEP(IPTRO+N+1).GT.KEEP(IPTRD+N+1)) THEN
            WRITE (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            DO 190 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               IF (ICNTL(3).EQ.3) J3 = MIN(J1,J2)
               DO 180 J = J1,J3
                  IF (KEEP(IPTRO+J+1).GT.KEEP(IPTRO+J)) WRITE (MP,
     +                '(A,I5,(T13,3(1P,D12.4,I5)))') ' Column',J,
     +                (A(I),IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
  180          CONTINUE
               J1 = J2 + 1
  190       CONTINUE
         END IF
         WRITE (MP,'(A)') ' Factorized matrix by blocks'
         J1 = 1
         DO 230 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            IF (KEEP(MBLOCK+3*JB).LT.0) GO TO 220
            NC = J2 - J1 + 1
            NR = NC
            IF (KB.EQ.1) NR = M
            WRITE (MP,'(A,I6)') ' Block',JB
            K = NEWNE
            IF (JB.GT.1) K = KEEP(KBLOCK+3*JB) + NEWNE - 1
            IF (KEEP(IPTRL+J1).GT.KEEP(IPTRU+J1)) WRITE (MP,
     +          '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J1,' of L',
     +          (A(K+I),IRN(K+I),I=KEEP(IPTRU+J1)+1,KEEP(IPTRL+J1))
            IF (ICNTL(3).EQ.3) GO TO 210
            DO 200 J = J1 + 1,J2
               IF (KEEP(IPTRU+J).GT.KEEP(IPTRL+J-1)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J,' of U',
     +              (A(K+I),IRN(K+I),I=KEEP(IPTRL+J-1)+1,KEEP(IPTRU+J))
               IF (KEEP(IPTRU+J).LT.KEEP(IPTRL+J)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J,' of L',
     +              (A(K+I),IRN(K+I),I=KEEP(IPTRU+J)+1,KEEP(IPTRL+J))
  200       CONTINUE
  210       WRITE (MP,'(A)') ' Full block'
            WRITE (MP,'(A)') ' Row indices'
            NRF = IRN(K+2)
            K = K + KEEP(IPTRL+J2)
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9120) (IRN(K+I),I=1,NRF)
               WRITE (MP,'(A)') ' Column pivoting information'
               WRITE (MP,9120) (IRN(K+I),I=NRF+1,
     +           NRF+NC-KEEP(MBLOCK+3*JB))
               WRITE (MP,'(A)') ' Reals by columns'
               WRITE (MP,9130) (A(K+I),I=1,
     +           (NRF)* (NC-KEEP(MBLOCK+3*JB)))
 9120          FORMAT (10I6)
 9130          FORMAT (1P,5D12.4)
            ELSE
               WRITE (MP,9120) (IRN(K+I),I=1,MIN(10,NRF))
               WRITE (MP,'(A)') ' Column pivoting information'
               WRITE (MP,9120) (IRN(K+I),I=NRF+1,
     +           NRF+MIN(10,NC-KEEP(MBLOCK+3*JB)))
               WRITE (MP,'(A)') ' Reals by columns'
               WRITE (MP,9130) (A(K+I),I=1,
     +           MIN(10, (NRF)* (NC-KEEP(MBLOCK+3*JB))))
            END IF
  220       J1 = J2 + 1
  230    CONTINUE
         IF (JOB.EQ.1 .OR. JOB.EQ.3) THEN
            WRITE (MP,'(A)') ' Contents of KEEP array'
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9020) (KEEP(I),I=1,M)
               WRITE (MP,9030) (KEEP(M+I),I=1,N)
               WRITE (MP,'(A)') ' Pointer information from KEEP'
               WRITE (MP,9140) (KEEP(IPTRL+J),J=1,N)
 9140          FORMAT (' IPTRL ='/ (8X,10I6))
               WRITE (MP,9150) (KEEP(IPTRU+J),J=1,N)
 9150          FORMAT (' IPTRU ='/ (8X,10I6))
            ELSE
               WRITE (MP,9020) (KEEP(I),I=1,MIN(10,M))
               WRITE (MP,9030) (KEEP(M+I),I=1,MIN(10,N))
               WRITE (MP,'(A)') ' Pointer information from KEEP'
               WRITE (MP,9140) (KEEP(IPTRL+J),J=1,MIN(10,N))
               WRITE (MP,9150) (KEEP(IPTRU+J),J=1,MIN(10,N))
            END IF
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,KB)
         END IF
      END IF

  240 RETURN

 9160 FORMAT (' Error return from MA48B/BD because M =',I10,' and N =',
     +       I10)
 9170 FORMAT (' Error return from MA48B/BD because NE =',I10)
 9180 FORMAT (' Error return from MA48B/BD because LA is',I10/' and mu',
     +       'st be at least',I10)
 9190 FORMAT (' Error return from MA48B/BD because ','JOB = ',I10)
      END


      SUBROUTINE MA48CD(M,N,TRANS,JOB,LA,A,IRN,KEEP,CNTL,ICNTL,RHS,X,
     +                  ERROR,W,IW,INFO)

      INTEGER M,N
      LOGICAL TRANS
      INTEGER JOB,LA
      DOUBLE PRECISION A(LA)
      INTEGER IRN(LA),KEEP(*)
      DOUBLE PRECISION CNTL(10)
      INTEGER ICNTL(20)
      DOUBLE PRECISION RHS(*),X(*),ERROR(3),W(*)
      INTEGER IW(*),INFO(20)

      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.D0,ONE=1.0D0)
      DOUBLE PRECISION COND(2),CTAU,DXMAX
      INTEGER I,ICNTL5(20),IPTRD,IPTRL,IPTRO,IPTRU,IQB(1),J,JB,JJ,J1,J2,
     +        J3,K,KASE,KB,KBLOCK,KK,KEEP71(5)
      LOGICAL LCOND(2)
      INTEGER LP,MBLOCK,MP,NB,NBLOCK,NC,NE,NEQ,NRF,NVAR
      DOUBLE PRECISION OLDOMG(2),OMEGA(2),OM1,OM2,TAU
      EXTERNAL MA48DD,MA50CD,MC71AD
      DOUBLE PRECISION EPS
      INTRINSIC ABS,MAX

      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (N.LE.0 .OR. M.LE.0) THEN
         INFO(1) = -1
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9140) M,N
         GO TO 380
      END IF
      IF (JOB.GT.4 .OR. JOB.LT.1) THEN
         INFO(1) = -6
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9150) JOB
         GO TO 380
      END IF
      INFO(1) = 0

      IPTRL = M + N
      IPTRU = IPTRL + N
      IPTRD = IPTRU + N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1
      NB = KEEP(KBLOCK+3)

      NE = KEEP(IPTRO+N+1) - 1

      OMEGA(1) = ZERO
      OMEGA(2) = ZERO
      ERROR(3) = ZERO
      EPS = EPSILON(EPS)
      CTAU = 1000.*EPS

      IQB(1) = 0

      DO 10 I = 1,7
         ICNTL5(I) = 0
   10 CONTINUE
      ICNTL5(5) = ICNTL(5)

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,
     +     '(/A/3(A,I8),A,I2/A,L2,A,I7/A,1P,E12.4/A,3I6,2(/A,I6))')
     +     ' Entering MA48C/CD with',' M =',M,'     N =',N,'     LA =',
     +     LA,'      JOB =',JOB,'   TRANS =',TRANS,
     +     '      No. of blocks =',
     +     NB,'   CNTL(5)    = ',CNTL(5),'   ICNTL(1:3) = ',
     +     (ICNTL(I),I=1,3),'   ICNTL(5)   = ',ICNTL(5),
     +     '   ICNTL(9)   = ',ICNTL(9)
         WRITE (MP,'(A)') ' Permuted matrix by blocks'
         KB = NB
         IF (ICNTL(3).EQ.3) THEN
            WRITE (MP,'(A)')
     +        ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         END IF
         WRITE (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         DO 30 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            IF (J1.LE.J2) WRITE (MP,'(A,I6)') ' Block',JB
            J3 = J2
            IF (ICNTL(3).EQ.3) J3 = MIN(J1,J2)
            DO 20 J = J1,J3
               WRITE (MP,'(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J,
     +           (A(I),IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
   20       CONTINUE
            J1 = J2 + 1
   30    CONTINUE
         IF (KEEP(IPTRO+N+1).GT.KEEP(IPTRD+N+1)) THEN
            WRITE (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            DO 50 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               IF (ICNTL(3).EQ.3) J3 = MIN(J1,J2)
               DO 40 J = J1,J3
                  IF (KEEP(IPTRO+J+1).GT.KEEP(IPTRO+J)) WRITE (MP,
     +                '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J,
     +                (A(I),IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
   40          CONTINUE
               J1 = J2 + 1
   50       CONTINUE
         END IF
         WRITE (MP,'(A)') ' Factorized matrix by blocks'
         J1 = 1
         DO 90 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            IF (KEEP(MBLOCK+3*JB).LT.0) GO TO 80
            NC = J2 - J1 + 1
            WRITE (MP,'(A,I6)') ' Block',JB
            K = NE
            IF (JB.GT.1) K = KEEP(KBLOCK+3*JB) + NE - 1
            IF (KEEP(IPTRL+J1).GT.KEEP(IPTRU+J1)) WRITE (MP,
     +          '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J1,' of L',
     +          (A(K+I),IRN(K+I),I=KEEP(IPTRU+J1)+1,KEEP(IPTRL+J1))
            IF (ICNTL(3).EQ.3) GO TO 70
            DO 60 J = J1 + 1,J2
               IF (KEEP(IPTRU+J).GT.KEEP(IPTRL+J-1)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of U',
     +              (A(K+I),IRN(K+I),I=KEEP(IPTRL+J-1)+1,KEEP(IPTRU+J))
               IF (KEEP(IPTRU+J).LT.KEEP(IPTRL+J)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L',
     +              (A(K+I),IRN(K+I),I=KEEP(IPTRU+J)+1,KEEP(IPTRL+J))
   60       CONTINUE
   70       WRITE (MP,'(A)') ' Full block'
            WRITE (MP,'(A)') ' Row indices'
            NRF = IRN(K+2)
            K = K + KEEP(IPTRL+J2)
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9000) (IRN(K+I),I=1,NRF)
               WRITE (MP,'(A)') ' Column pivoting information'
               WRITE (MP,9000) (IRN(K+I),I=NRF+1,
     +           NRF+NC-KEEP(MBLOCK+3*JB))
               WRITE (MP,'(A)') ' Reals by columns'
               WRITE (MP,9010) (A(K+I),I=1,
     +           (NRF)* (NC-KEEP(MBLOCK+3*JB)))
 9000          FORMAT (10I6)
 9010          FORMAT (1P,5D12.4)
            ELSE
               WRITE (MP,9000) (IRN(K+I),I=1,MIN(10,NRF))
               WRITE (MP,'(A)') ' Column pivoting information'
               WRITE (MP,9000) (IRN(K+I),I=NRF+1,
     +           NRF+MIN(10,NC-KEEP(MBLOCK+3*JB)))
               WRITE (MP,'(A)') ' Reals by columns'
               WRITE (MP,9010) (A(K+I),I=1,
     +           MIN(10, (NRF)* (NC-KEEP(MBLOCK+3*JB))))
            END IF
   80       J1 = J2 + 1
   90    CONTINUE
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9020) (KEEP(I),I=1,M)
            WRITE (MP,9030) (KEEP(M+I),I=1,N)
 9020       FORMAT (' Positions of original rows in the permuted matrix'
     +             / (10I6))
 9030       FORMAT (' Positions of columns of permuted matrix ','in or',
     +             'ig','inal matrix '/ (10I6))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9040) (KEEP(IPTRD+J),J=1,N+1)
 9040       FORMAT (' IPTRD ='/ (8X,10I6))
            WRITE (MP,9050) (KEEP(IPTRO+J),J=1,N+1)
 9050       FORMAT (' IPTRO ='/ (8X,10I6))
            WRITE (MP,9060) (KEEP(IPTRL+J),J=1,N)
 9060       FORMAT (' IPTRL ='/ (8X,10I6))
            WRITE (MP,9070) (KEEP(IPTRU+J),J=1,N)
 9070       FORMAT (' IPTRU ='/ (8X,10I6))
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,NB)
            WRITE (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,NB)
            WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,NB)
 9080       FORMAT (' NBLOCK (order blocks) ='/ (8X,10I6))
 9090       FORMAT (' MBLOCK (triangular flag and number packed rows) ='
     +             / (8X,10I6))
 9100       FORMAT (' KBLOCK (position of beginning of block) ='/
     +             (8X,10I6))
            IF (TRANS) THEN
               WRITE (MP,9110) (RHS(I),I=1,N)
 9110          FORMAT (' RHS =  ',1P,5D12.4/ (8X,5D12.4))
            ELSE
               WRITE (MP,9110) (RHS(I),I=1,M)
            END IF
         ELSE
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9020) (KEEP(I),I=1,MIN(10,M))
            WRITE (MP,9030) (KEEP(M+I),I=1,MIN(10,N))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9040) (KEEP(IPTRD+K),K=1,MIN(10,N+1))
            WRITE (MP,9050) (KEEP(IPTRO+K),K=1,MIN(10,N+1))
            WRITE (MP,9060) (KEEP(IPTRL+J),J=1,MIN(10,N))
            WRITE (MP,9070) (KEEP(IPTRU+J),J=1,MIN(10,N))
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,KB)
            WRITE (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,KB)
            WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,KB)
            IF (TRANS) THEN
               WRITE (MP,9110) (RHS(I),I=1,MIN(10,N))
            ELSE
               WRITE (MP,9110) (RHS(I),I=1,MIN(10,M))
            END IF
         END IF
      END IF

      IF (TRANS) THEN
         NEQ = N
         NVAR = M
         DO 100 I = 1,NEQ
            W(I) = RHS(KEEP(M+I))
  100    CONTINUE
      ELSE
         NEQ = M
         NVAR = N
         DO 110 I = 1,NEQ
            W(KEEP(I)) = RHS(I)
  110    CONTINUE
      END IF
      IF (JOB.EQ.1) THEN
         IF (NB.EQ.1 .AND. KEEP(MBLOCK+3).GE.0) THEN
            CALL MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),TRANS,LA-NE,
     +                  A(NE+1),IRN(NE+1),KEEP(IPTRL+1),KEEP(IPTRU+1),W,
     +                  X,W(NEQ+1),INFO)
         ELSE
            CALL MA48DD
     +           (N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN,KEEP(IPTRD+1),
     +            KEEP(IPTRO+1),NB,KEEP(NBLOCK+3),KEEP(IPTRL+1),
     +            KEEP(IPTRU+1),W,X,TRANS,ICNTL5,W(NEQ+1))
         END IF
         GO TO 340
      END IF

      DO 120 I = 1,NVAR
         X(I) = ZERO
  120 CONTINUE
      DO 130 I = 1,NEQ
         RHS(I) = W(I)
  130 CONTINUE

      OM1 = ZERO
      DO 260 K = 1,ICNTL(9)
         IF (NB.EQ.1 .AND. KEEP(MBLOCK+3).GE.0) THEN
            CALL MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),TRANS,LA-NE,
     +                A(NE+1),IRN(NE+1),KEEP(IPTRL+1),KEEP(IPTRU+1),
     +                W,W(NEQ+1),W(M+N+1),INFO)
         ELSE
            CALL MA48DD
     +           (N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN,KEEP(IPTRD+1),
     +            KEEP(IPTRO+1),NB,KEEP(NBLOCK+3),KEEP(IPTRL+1),
     +            KEEP(IPTRU+1),W,W(NEQ+1),TRANS,ICNTL5,W(M+N+1))
         END IF
         DO 140 I = 1,NVAR
            X(I) = X(I) + W(NEQ+I)
  140    CONTINUE
         DO 150 I = 1,NEQ
            W(I) = RHS(I)
            W(NEQ+I) = ZERO
            W(2*NEQ+I) = ZERO
  150    CONTINUE
         IF (TRANS) THEN
            DO 180 J = 1,N
               DO 160 JJ = KEEP(IPTRD+J),KEEP(IPTRD+J+1) - 1
                  I = IRN(JJ)
                  W(J) = W(J) - A(JJ)*X(I)
                  W(NEQ+J) = W(NEQ+J) + ABS(A(JJ)*X(I))
                  W(2*NEQ+J) = W(2*NEQ+J) + ABS(A(JJ))
  160          CONTINUE
               DO 170 JJ = KEEP(IPTRO+J),KEEP(IPTRO+J+1) - 1
                  I = IRN(JJ)
                  W(J) = W(J) - A(JJ)*X(I)
                  W(NEQ+J) = W(NEQ+J) + ABS(A(JJ)*X(I))
                  W(2*NEQ+J) = W(2*NEQ+J) + ABS(A(JJ))
  170          CONTINUE
  180       CONTINUE
         ELSE
            DO 210 J = 1,N
               DO 190 JJ = KEEP(IPTRD+J),KEEP(IPTRD+J+1) - 1
                  I = IRN(JJ)
                  W(I) = W(I) - A(JJ)*X(J)
                  W(NEQ+I) = W(NEQ+I) + ABS(A(JJ)*X(J))
                  W(2*NEQ+I) = W(2*NEQ+I) + ABS(A(JJ))
  190          CONTINUE
               DO 200 JJ = KEEP(IPTRO+J),KEEP(IPTRO+J+1) - 1
                  I = IRN(JJ)
                  W(I) = W(I) - A(JJ)*X(J)
                  W(NEQ+I) = W(NEQ+I) + ABS(A(JJ)*X(J))
                  W(2*NEQ+I) = W(2*NEQ+I) + ABS(A(JJ))
  200          CONTINUE
  210       CONTINUE
         END IF
         DXMAX = ZERO
         DO 220 I = 1,NVAR
            DXMAX = MAX(DXMAX,ABS(X(I)))
  220    CONTINUE
         OMEGA(1) = ZERO
         OMEGA(2) = ZERO
         DO 230 I = 1,NEQ
            TAU = (W(2*NEQ+I)*DXMAX+ABS(RHS(I)))*NVAR*CTAU
            IF ((W(NEQ+I)+ABS(RHS(I))).GT.TAU) THEN
               OMEGA(1) = MAX(OMEGA(1),ABS(W(I))/
     +                    (W(NEQ+I)+ABS(RHS(I))))
               IW(I) = 1
            ELSE
               IF (TAU.GT.ZERO) THEN
                  OMEGA(2) = MAX(OMEGA(2),ABS(W(I))/
     +                       (W(NEQ+I)+W(2*NEQ+I)*DXMAX))
               END IF
               IW(I) = 2
            END IF
  230    CONTINUE
         IF (JOB.EQ.2) GO TO 340
         OM2 = OMEGA(1) + OMEGA(2)
         IF (OM2.LE.EPS) GO TO 270
         IF (K.GT.1 .AND. OM2.GT.OM1*CNTL(5)) THEN
            IF (OM2.GT.OM1) THEN
               OMEGA(1) = OLDOMG(1)
               OMEGA(2) = OLDOMG(2)
               DO 240 I = 1,NVAR
                  X(I) = W(3*NEQ+I)
  240          CONTINUE
            END IF
            GO TO 270
         END IF
         DO 250 I = 1,NVAR
            W(3*NEQ+I) = X(I)
  250    CONTINUE
         OLDOMG(1) = OMEGA(1)
         OLDOMG(2) = OMEGA(2)
         OM1 = OM2
  260 CONTINUE
      INFO(1) = -8
      IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9170) INFO(1),ICNTL(9)
      GO TO 340

  270 IF (JOB.LE.3) GO TO 340
      IF (M.NE.N) GO TO 340
      LCOND(1) = .FALSE.
      LCOND(2) = .FALSE.
      DO 280 I = 1,NEQ
         IF (IW(I).EQ.1) THEN
            W(I) = W(NEQ+I) + ABS(RHS(I))
            W(NEQ+I) = ZERO
            LCOND(1) = .TRUE.
         ELSE

            W(NEQ+I) = W(NEQ+I) + W(2*NEQ+I)*DXMAX
            W(I) = ZERO
            LCOND(2) = .TRUE.
         END IF
  280 CONTINUE
      KASE = 0
      DO 330 K = 1,2
         IF (LCOND(K)) THEN
            DO 310 KK = 1,40
               CALL MC71AD(N,KASE,W(3*NEQ+1),COND(K),RHS,IW,KEEP71)
               IF (KASE.EQ.0) GO TO 320
               IF (KASE.EQ.1) THEN
                  IF (NB.EQ.1 .AND. KEEP(MBLOCK+3).GE.0) THEN
                     CALL MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),
     +                           .NOT.TRANS,LA-NE,A(NE+1),IRN(NE+1),
     +                           KEEP(IPTRL+1),KEEP(IPTRU+1),W(3*NEQ+1),
     +                           W(2*NEQ+1),RHS,INFO)
                  ELSE
                     CALL MA48DD(N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN,
     +                           KEEP(IPTRD+1),KEEP(IPTRO+1),NB,
     +                           KEEP(NBLOCK+3),KEEP(IPTRL+1),
     +                           KEEP(IPTRU+1),W(3*NEQ+1),W(2*NEQ+1),
     +                           .NOT.TRANS,ICNTL5,RHS)
                  END IF

                  DO 290 I = 1,M
                     W(3*NEQ+I) = W((K-1)*NEQ+I)*W(2*NEQ+I)
  290             CONTINUE
               END IF
               IF (KASE.EQ.2) THEN
                  DO 300 I = 1,N
                     W(2*NEQ+I) = W((K-1)*NEQ+I)*W(3*NEQ+I)
  300             CONTINUE
                  IF (NB.EQ.1 .AND. KEEP(MBLOCK+3).GE.0) THEN
                     CALL MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),TRANS,
     +                           LA-NE,A(NE+1),IRN(NE+1),KEEP(IPTRL+1),
     +                           KEEP(IPTRU+1),W(2*NEQ+1),W(3*NEQ+1),
     +                           RHS,INFO)
                  ELSE
                     CALL MA48DD(N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN,
     +                           KEEP(IPTRD+1),KEEP(IPTRO+1),NB,
     +                           KEEP(NBLOCK+3),KEEP(IPTRL+1),
     +                           KEEP(IPTRU+1),W(2*NEQ+1),W(3*NEQ+1),
     +                           TRANS,ICNTL5,RHS)
                  END IF
               END IF
  310       CONTINUE
            INFO(1) = -9
            IF (LP.NE.0 .AND. ICNTL(3).GE.1) WRITE (LP,9160)
            GO TO 340
  320       IF (DXMAX.GT.ZERO) COND(K) = COND(K)/DXMAX
            ERROR(3) = ERROR(3) + OMEGA(K)*COND(K)
         END IF
  330 CONTINUE

  340 DO 350 I = 1,NVAR
         W(I) = X(I)
  350 CONTINUE
      IF (.NOT.TRANS) THEN
         DO 360 I = 1,NVAR
            X(KEEP(M+I)) = W(I)
  360    CONTINUE
      ELSE
         DO 370 I = 1,NVAR
            X(I) = W(KEEP(I))
  370    CONTINUE
      END IF
      IF (JOB.GE.2) THEN
         ERROR(1) = OMEGA(1)
         ERROR(2) = OMEGA(2)
      END IF

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/A,I6)') ' Leaving MA48C/CD with INFO(1) =',INFO(1)
         IF (JOB.GT.1) THEN
            K = 2
            IF (JOB.EQ.4 .AND. INFO(1).NE.-9) K = 3
            WRITE (MP,9120) (ERROR(I),I=1,K)
 9120       FORMAT (' ERROR =',1P,3D12.4)
         END IF
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9130) (X(I),I=1,NVAR)
 9130       FORMAT (' X =    ',1P,5D12.4:/ (8X,5D12.4))
         ELSE
            WRITE (MP,9130) (X(I),I=1,MIN(10,NVAR))
         END IF
      END IF
  380 RETURN

 9140 FORMAT (' Error return from MA48C/CD because M =',I10,' and N =',
     +       I10)
 9150 FORMAT (' Error return from MA48C/CD because ','JOB = ',I10)
 9160 FORMAT (' Error return from MA48C/CD because of ','error in MC71',
     +       'A/AD'/' ERROR(3) not calculated')
 9170 FORMAT (' Error return from MA48C/CD because of ','nonconvergenc',
     +       'e of iterative refinement'/' Error INFO(1) = ',I2,'  wit',
     +       'h ICNTL','(9) = ',I10)
      END


      SUBROUTINE MA48DD(N,NE,LA,A,AA,IRN,IRNA,IPTRD,IPTRO,NB,IBLOCK,
     +                  IPTRL,IPTRU,RHS,X,TRANS,ICNTL5,W)
      INTEGER N,NE,LA
      DOUBLE PRECISION A(LA),AA(NE)
      INTEGER IRN(LA),IRNA(NE),IPTRD(N+1),IPTRO(N+1),NB
      INTEGER IBLOCK(3,NB),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION RHS(N),X(N)
      LOGICAL TRANS
      INTEGER ICNTL5(20)
      DOUBLE PRECISION W(N)
      INTEGER I,IFLAG(15),IQB(1),J,JB,JJ,J1,K1,K2,NC,NUMB
      EXTERNAL MA50CD
      IQB(1) = 0
      NUMB = IBLOCK(3,1)
      IF (.NOT.TRANS) THEN
         K1 = N + 1
         DO 50 JB = NUMB,1,-1
            NC = IBLOCK(1,JB)
            K2 = K1 - 1
            K1 = K1 - NC
            IF (IBLOCK(2,JB).LT.0) THEN
               DO 20 J = K2,K1,-1
                  X(J) = RHS(J)/AA(IPTRD(J))
                  DO 10 JJ = IPTRD(J) + 1,IPTRD(J+1) - 1
                     I = IRNA(JJ)
                     RHS(I) = RHS(I) - AA(JJ)*X(J)
   10             CONTINUE
   20          CONTINUE
            ELSE
               J1 = 1
               IF (JB.GT.1) J1 = IBLOCK(3,JB)
               CALL MA50CD(NC,NC,ICNTL5,IQB,IBLOCK(2,JB),TRANS,LA+1-J1,
     +                     A(J1),IRN(J1),IPTRL(K1),IPTRU(K1),RHS(K1),
     +                     X(K1),W,IFLAG)
            END IF
            IF (JB.EQ.1) GO TO 50
            DO 40 J = K1,K2
               DO 30 JJ = IPTRO(J),IPTRO(J+1) - 1
                  I = IRNA(JJ)
                  RHS(I) = RHS(I) - AA(JJ)*X(J)
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      ELSE
         K2 = 0
         DO 100 JB = 1,NUMB
            NC = IBLOCK(1,JB)
            K1 = K2 + 1
            K2 = K2 + NC
            IF (JB.GT.1) THEN
               DO 70 J = K1,K2
                  DO 60 JJ = IPTRO(J),IPTRO(J+1) - 1
                     I = IRNA(JJ)
                     RHS(J) = RHS(J) - AA(JJ)*X(I)
   60             CONTINUE
   70          CONTINUE
            END IF
            IF (IBLOCK(2,JB).LT.0) THEN
               DO 90 J = K1,K2
                  DO 80 JJ = IPTRD(J) + 1,IPTRD(J+1) - 1
                     I = IRNA(JJ)
                     RHS(J) = RHS(J) - AA(JJ)*X(I)
   80             CONTINUE
                  X(J) = RHS(J)/AA(IPTRD(J))
   90          CONTINUE
            ELSE
               J1 = 1
               IF (JB.GT.1) J1 = IBLOCK(3,JB)
               CALL MA50CD(NC,NC,ICNTL5,IQB,IBLOCK(2,JB),TRANS,LA+1-J1,
     +                     A(J1),IRN(J1),IPTRL(K1),IPTRU(K1),RHS(K1),
     +                     X(K1),W,IFLAG)
            END IF
  100    CONTINUE
      END IF

      RETURN
      END

      SUBROUTINE MA48ID(CNTL,ICNTL)

      DOUBLE PRECISION CNTL(10)
      INTEGER ICNTL(20)


      INTEGER I

      DO 10 I = 3, 10
         CNTL(I) = 0.0D0
   10 CONTINUE
      CNTL(1) = 0.5D0
      CNTL(2) = 0.1D0
      CNTL(5) = 0.5D0
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 2
      ICNTL(4) = 3
      ICNTL(5) = 32
      ICNTL(6) = 1
      ICNTL(7) = 1
      ICNTL(8) = 0
      ICNTL(9) = 10
      DO 20 I = 10, 20
         ICNTL(I) = 0
   20 CONTINUE

      END

