 
      SUBROUTINE MC75ID(ICNTL)
      INTEGER ICNTL(5)

      ICNTL(1) = 6
      ICNTL(2) = 0
      ICNTL(3) = 0
      ICNTL(4) = 0
      ICNTL(5) = 0

      RETURN
      END
      SUBROUTINE MC75AD(N,NZ,LA,A,IRN,JCN,COND,LIW,IW,LW,W,ICNTL,INFO)
      DOUBLE PRECISION ZERO,U
      PARAMETER (ZERO=0.0D0,U=0.1D0)
      INTEGER LA,N,NZ,LIW,LW
      DOUBLE PRECISION A(LA),COND(2),W(LW)
      INTEGER IRN(LA),IW(LIW),JCN(LA),ICNTL(5),INFO(5)
      DOUBLE PRECISION ERROR(3)
      INTEGER KEEP71(5)
      DOUBLE PRECISION CNTL48(10),RINF48(10)
      INTEGER ICNT48(20),INFO48(20)
      DOUBLE PRECISION DANRM
      INTEGER I,IDUMMY,K,KASE,IW71
      INTEGER IIW,LIIW,IKEEP,LIKEEP,MINIW,IW75B,IWRHS,IWX,IW48,MINW
      INTEGER IDAMAX
      EXTERNAL IDAMAX
      EXTERNAL MA48AD,MA48BD,MA48CD,MC71AD,MC75BD
      INTRINSIC ABS
      INFO(1) = 0
      IF (N.LE.0) THEN
        INFO(1) = -7
        IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9700) INFO(1)
 9700   FORMAT (1X,'Error in MC75: no.',I5,' (N LE 0)')
        GO TO 1010
      END IF

      IF (NZ.LE.0) THEN
        INFO(1) = -6
        IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9600) INFO(1)
 9600   FORMAT (1X,'Error in MC75: no.',I5,' (NZ LE 0)')
        GO TO 1010
      END IF

      IF (LA.LT.NZ) THEN
        INFO(1) = -5
        IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9500) INFO(1)
 9500   FORMAT (1X,'Error in MC75: no.',I5,' (LA LT NZ)')
        GO TO 1010
      END IF

      COND(1) = ZERO
      COND(2) = ZERO

      IIW = 0
      LIIW = 9*N
      IKEEP = IIW + LIIW
      LIKEEP = 10*N + 7
      MINIW = IKEEP + LIKEEP

      IF (MINIW.GT.LIW) THEN
        INFO(1) = -4
        INFO(2) = MINIW
        IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9400) INFO(1)
 9400   FORMAT (1X,'Error in MC75: no.',I5,' (IW array too short)')
        GO TO 1010
      END IF
      IW75B = 0
      IW71  = IW75B + N
      IWRHS = IW71 + N
      IWX   = IWRHS + N
      IW48  = IWX + N
      MINW  = IW48 + 4*N

      IF (MINW.GT.LW) THEN
        INFO(1) = -3
        INFO(2) = MINW
        IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9300) INFO(1)
 9300   FORMAT (1X,'Error in MC75: no.',I5,' (W array too short)')
        GO TO 1010
      END IF
      DO 140 I = 1,N
        W(IW75B+I) = ZERO
  140 CONTINUE
      DO 150 K = 1,NZ
        I = IRN(K)
        W(IW75B+I) = W(IW75B+I) + ABS(A(K))
  150 CONTINUE
      DANRM = ABS(W(IW75B+IDAMAX(N,W(IW75B+1),1)))
      CALL MA48ID(CNTL48,ICNT48)
      ICNT48(1) = ICNTL(2)
      ICNT48(2) = ICNTL(2)
      ICNT48(7) = 0
      CALL MA48AD(N,N,NZ,1,LA,A,IRN,JCN,IW(IKEEP+1),CNTL48,
     +           ICNT48,IW(IIW+1),INFO48,RINF48)
      INFO(3) = INFO48(3)
      INFO(4) = INFO48(4)
      INFO(5) = INFO48(5)
      IF (INFO48(1).EQ.1) THEN
        IF (INFO48(12).GT.0) THEN
          INFO(1) = -8
        ELSE
          INFO(1) = 1
        ENDIF
      ELSE IF (INFO48(1).EQ.-3) THEN
        INFO(1) = -1
      ELSE IF (INFO48(1).LT.0) THEN
        INFO(1) = -9
        INFO(2) = INFO48(1)
        IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9009) INFO(1)
 9009   FORMAT (1X,'Error in MC75: no.',I5,' (MA48A/AD error)')
        GO TO 1010
      END IF

      IF (INFO(1).LT.0) GO TO 9999
      CALL MA48BD(N,N,NZ,1,LA,A,IRN,JCN,IW(IKEEP+1),CNTL48,
     +           ICNT48,W(IW48+1),IW(IIW+1),INFO48,RINF48)
      INFO(3) = INFO48(3)
      INFO(4) = INFO48(4)
      INFO(5) = INFO48(5)
      IF (INFO48(1).EQ.-3) THEN
        INFO(1) = -2
      ELSE IF (INFO48(1).EQ.2) THEN
        INFO(1) = 2
      ELSE IF (INFO48(1).LT.0) THEN
        INFO(1) = -10
        INFO(2) = INFO48(1)
        IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9010) INFO(1)
 9010   FORMAT (1X,'Error in MC75: no.',I5,' (MA48B/BD error)')
        GO TO 1010
      END IF

      IF (INFO(1).LT.0) GO TO 9999
      KASE = 0
      IF (INFO(1).LE.1) THEN
        DO 170 IDUMMY = 1,100
          CALL MC71AD(N,KASE,W(IWX+1),COND(2),W(IW71+1),IW(IIW+1),
     +                KEEP71)
          IF (KASE.LE.0) GO TO 160
          IF (KASE.EQ.1) CALL MC75BD(N,W(IWX+1),W(IW75B+1))
          DO 165 I = 1,N
            W(IWRHS+I) = W(IWX+I)
  165     CONTINUE
          CALL MA48CD(N,N,KASE.EQ.1,1,LA,A,IRN,IW(IKEEP+1),CNTL48,
     +               ICNT48,W(IWRHS+1),W(IWX+1),ERROR,
     +               W(IW48+1),IW(IIW+1),INFO48)
          IF (INFO48(1).LT.0) THEN
             INFO(1) = -11
             INFO(2) = INFO48(1)
             IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9011) INFO(1)
 9011        FORMAT (1X,'Error in MC75: no.',I5,' (MA48C/CD error)')
             GO TO 1010
          END IF
          IF (KASE.EQ.2) CALL MC75BD(N,W(IWX+1),W(IW75B+1))
  170   CONTINUE
      END IF

  160 KASE = 0
      DO 190 IDUMMY = 1,100
        CALL MC71AD(N,KASE,W(IWX+1),COND(1),W(IW71+1),IW(IIW+1),
     +              KEEP71)
        IF (KASE.LE.0) GO TO 180
        DO 175 I = 1,N
          W(IWRHS+I) = W(IWX+I)
  175   CONTINUE
        CALL MA48CD(N,N,KASE.EQ.1,1,LA,A,IRN,IW(IKEEP+1),CNTL48,
     +             ICNT48,W(IWRHS+1),W(IWX+1),ERROR,
     +             W(IW48+1),IW(IIW+1),INFO48)
        IF (INFO48(1).LT.0) THEN
           INFO(1) = -11
           INFO(2) = INFO48(1)
           IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9011) INFO(1)
           GO TO 1010
        END IF
  190 CONTINUE
  180 COND(1) = COND(1)*DANRM
      GO TO 1010

 9999 IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9000) INFO(1)
 9000 FORMAT (1X,'Error in MC75: no.',I5)

 1010 RETURN
      END
      SUBROUTINE MC75BD(N,R,W)
      INTEGER N
      DOUBLE PRECISION R(*),W(*)
      INTEGER I
      DO 100 I = 1,N
        R(I) = R(I)*W(I)
  100 CONTINUE
      RETURN
      END
