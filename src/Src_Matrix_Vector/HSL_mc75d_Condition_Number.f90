
      subroutine MC75ID(ICNTL)
      integer :: ICNTL(5)

      ICNTL(1) = 6
      ICNTL(2) = 0
      ICNTL(3) = 0
      ICNTL(4) = 0
      ICNTL(5) = 0

      return
      END
      subroutine MC75AD(N,NZ,LA,A,IRN,JCN,COND,LIW,IW,LW,W,ICNTL,INFO)
      double precision ZERO,U
      parameter (ZERO=0.0D0,U=0.1D0)
      integer LA,N,NZ,LIW,LW
      double precision A(LA),COND(2),W(LW)
      integer IRN(LA),IW(LIW),JCN(LA),ICNTL(5),INFO(5)
      double precision ERROR(3)
      integer :: KEEP71(5)
      double precision CNTL48(10),RINF48(10)
      integer ICNT48(20),INFO48(20)
      double precision DANRM
      integer I,IDUMMY,K,KASE,IW71
      integer IIW,LIIW,IKEEP,LIKEEP,MINIW,IW75B,IWRHS,IWX,IW48,MINW
      integer :: IDAMAX
      external IDAMAX
      external MA48AD,MA48BD,MA48CD,MC71AD,MC75BD
      intrinsic ABS
      INFO(1) = 0
      if (N.le.0) then
        INFO(1) = -7
        if (ICNTL(1).gt.0) write (ICNTL(1),FMT=9700) INFO(1)
 9700   format (1X,'Error in MC75: no.',I5,' (N LE 0)')
        GO TO 1010
      end if

      if (NZ.le.0) then
        INFO(1) = -6
        if (ICNTL(1).gt.0) write (ICNTL(1),FMT=9600) INFO(1)
 9600   format (1X,'Error in MC75: no.',I5,' (NZ LE 0)')
        GO TO 1010
      end if

      if (LA.lt.NZ) then
        INFO(1) = -5
        if (ICNTL(1).gt.0) write (ICNTL(1),FMT=9500) INFO(1)
 9500   format (1X,'Error in MC75: no.',I5,' (LA LT NZ)')
        GO TO 1010
      end if

      COND(1) = ZERO
      COND(2) = ZERO

      IIW = 0
      LIIW = 9*N
      IKEEP = IIW + LIIW
      LIKEEP = 10*N + 7
      MINIW = IKEEP + LIKEEP

      if (MINIW.gt.LIW) then
        INFO(1) = -4
        INFO(2) = MINIW
        if (ICNTL(1).gt.0) write (ICNTL(1),FMT=9400) INFO(1)
 9400   format (1X,'Error in MC75: no.',I5,' (IW array too short)')
        GO TO 1010
      end if
      IW75B = 0
      IW71  = IW75B + N
      IWRHS = IW71 + N
      IWX   = IWRHS + N
      IW48  = IWX + N
      MINW  = IW48 + 4*N

      if (MINW.gt.LW) then
        INFO(1) = -3
        INFO(2) = MINW
        if (ICNTL(1).gt.0) write (ICNTL(1),FMT=9300) INFO(1)
 9300   format (1X,'Error in MC75: no.',I5,' (W array too short)')
        GO TO 1010
      end if
      do 140 I = 1,N
        W(IW75B+I) = ZERO
  140 continue
      do 150 K = 1,NZ
        I = IRN(K)
        W(IW75B+I) = W(IW75B+I) + ABS(A(K))
  150 continue
      DANRM = ABS(W(IW75B+IDAMAX(N,W(IW75B+1),1)))
      call MA48ID(CNTL48,ICNT48)
      ICNT48(1) = ICNTL(2)
      ICNT48(2) = ICNTL(2)
      ICNT48(7) = 0
call MA48AD(N,N,NZ,1,LA,A,IRN,JCN,IW(IKEEP+1),CNTL48, ICNT48,IW(IIW+1),INFO48,RINF48)
      INFO(3) = INFO48(3)
      INFO(4) = INFO48(4)
      INFO(5) = INFO48(5)
      if (INFO48(1).eq.1) then
        if (INFO48(12).gt.0) then
          INFO(1) = -8
        else
          INFO(1) = 1
        endif
      else if (INFO48(1).eq.-3) then
        INFO(1) = -1
      else if (INFO48(1).lt.0) then
        INFO(1) = -9
        INFO(2) = INFO48(1)
        if (ICNTL(1).gt.0) write (ICNTL(1),FMT=9009) INFO(1)
 9009   format (1X,'Error in MC75: no.',I5,' (MA48A/AD error)')
        GO TO 1010
      end if

      if (INFO(1).lt.0) GO TO 9999
call MA48BD(N,N,NZ,1,LA,A,IRN,JCN,IW(IKEEP+1),CNTL48, ICNT48,W(IW48+1),IW(IIW+1),INFO48,RINF48)
      INFO(3) = INFO48(3)
      INFO(4) = INFO48(4)
      INFO(5) = INFO48(5)
      if (INFO48(1).eq.-3) then
        INFO(1) = -2
      else if (INFO48(1).eq.2) then
        INFO(1) = 2
      else if (INFO48(1).lt.0) then
        INFO(1) = -10
        INFO(2) = INFO48(1)
        if (ICNTL(1).gt.0) write (ICNTL(1),FMT=9010) INFO(1)
 9010   format (1X,'Error in MC75: no.',I5,' (MA48B/BD error)')
        GO TO 1010
      end if

      if (INFO(1).lt.0) GO TO 9999
      KASE = 0
      if (INFO(1).le.1) then
        do 170 IDUMMY = 1,100
call MC71AD(N,KASE,W(IWX+1),COND(2),W(IW71+1),IW(IIW+1), KEEP71)
          if (KASE.le.0) GO TO 160
          if (KASE.eq.1) call MC75BD(N,W(IWX+1),W(IW75B+1))
          do 165 I = 1,N
            W(IWRHS+I) = W(IWX+I)
  165     continue
call MA48CD(N,N,KASE.eq.1,1,LA,A,IRN,IW(IKEEP+1),CNTL48, ICNT48,W(IWRHS+1),W(IWX+1),ERROR, W(IW48+1),IW(IIW+1),INFO48)
          if (INFO48(1).lt.0) then
             INFO(1) = -11
             INFO(2) = INFO48(1)
             if (ICNTL(1).gt.0) write (ICNTL(1),FMT=9011) INFO(1)
 9011        format (1X,'Error in MC75: no.',I5,' (MA48C/CD error)')
             GO TO 1010
          end if
          if (KASE.eq.2) call MC75BD(N,W(IWX+1),W(IW75B+1))
  170   continue
      end if

  160 KASE = 0
      do 190 IDUMMY = 1,100
call MC71AD(N,KASE,W(IWX+1),COND(1),W(IW71+1),IW(IIW+1), KEEP71)
        if (KASE.le.0) GO TO 180
        do 175 I = 1,N
          W(IWRHS+I) = W(IWX+I)
  175   continue
call MA48CD(N,N,KASE.eq.1,1,LA,A,IRN,IW(IKEEP+1),CNTL48, ICNT48,W(IWRHS+1),W(IWX+1),ERROR, W(IW48+1),IW(IIW+1),INFO48)
        if (INFO48(1).lt.0) then
           INFO(1) = -11
           INFO(2) = INFO48(1)
           if (ICNTL(1).gt.0) write (ICNTL(1),FMT=9011) INFO(1)
           GO TO 1010
        end if
  190 continue
  180 COND(1) = COND(1)*DANRM
      GO TO 1010

 9999 if (ICNTL(1).gt.0) write (ICNTL(1),FMT=9000) INFO(1)
 9000 format (1X,'Error in MC75: no.',I5)

 1010 return
      END
      subroutine MC75BD(N,R,W)
      integer :: N
      double precision R(*),W(*)
      integer :: I
      do 100 I = 1,N
        R(I) = R(I)*W(I)
  100 continue
      return
      END
