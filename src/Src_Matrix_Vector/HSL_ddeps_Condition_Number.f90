
subroutine MC59AD(ICNTL,NC,NR,NE,IRN,LJCN,JCN,LA,A,LIP,IP, LIW,IW,INFO)
      integer LA,LIP,LIW,LJCN,NC,NE,NR
      double precision A(LA)
      integer ICNTL(10),IP(LIP),INFO(10),IRN(NE),IW(LIW),JCN(LJCN)
      integer I,ICNTL1,ICNTL2,ICNTL3,ICNTL6,LAA
      integer IDUP,IOUT,IUP,JOUT,LP,MP,KNE,PART
      logical :: LCHECK
      external MC59BD,MC59CD,MC59DD,MC59ED,MC59FD
      intrinsic MAX

      do 10 I = 1,10
         INFO(I) = 0
   10 continue

      ICNTL1 = ICNTL(1)
      ICNTL2 = ICNTL(2)
      ICNTL3 = ICNTL(3)
      ICNTL6 = ICNTL(6)
      LCHECK = (ICNTL1.eq.0)
      LP = ICNTL(4)
      MP = ICNTL(5)

      if (ICNTL2.gt.2 .or. ICNTL2.lt.0) then
         INFO(1) = -1
         INFO(2) = ICNTL2
         if (LP.gt.0) then
            write (LP,FMT=9000) INFO(1)
            write (LP,FMT=9010) ICNTL2
         end if
         GO TO 70
      end if

      if (ICNTL6.gt.2 .or. ICNTL6.lt.-2) then
         INFO(1) = -11
         INFO(2) = ICNTL6
         if (LP.gt.0) then
            write (LP,FMT=9000) INFO(1)
            write (LP,FMT=9150) ICNTL6
         end if
         GO TO 70
      end if

      if (NC.lt.1) then
        INFO(1) = -2
        INFO(2) = NC
        if (LP.gt.0) then
          write (LP,FMT=9000) INFO(1)
          write (LP,FMT=9020) NC
        end if
        GO TO 70
      end if

      if (NR.lt.1) then
        INFO(1) = -3
        INFO(2) = NR
        if (LP.gt.0) then
          write (LP,FMT=9000) INFO(1)
          write (LP,FMT=9030) NR
        end if
        GO TO 70
      end if

      if (ICNTL6.ne.0 .and. NR.ne.NC) then
        INFO(1) = -3
        INFO(2) = NR
        if (LP.gt.0) then
          write (LP,FMT=9000) INFO(1)
          write (LP,FMT=9035) NC,NR
        end if
        GO TO 70
      end if

      if (NE.lt.1) then
        INFO(1) = -4
        INFO(2) = NE
        if (LP.gt.0) then
          write (LP,FMT=9000) INFO(1)
          write (LP,FMT=9040) NE
        end if
        GO TO 70
      end if

      if (ICNTL2.eq.0 .or. ICNTL2.eq.1) then
        if (LJCN.lt.NE) then
          INFO(1) = -5
          INFO(2) = NE
        end if
      else
        if (LJCN.lt.1) then
          INFO(1) = -5
          INFO(2) = 1
        end if
      end if
      if (INFO(1).eq.-5) then
         if (LP.gt.0) then
            write (LP,FMT=9000) INFO(1)
            write (LP,FMT=9050) LJCN,INFO(2)
         end if
         GO TO 70
      end if

      if (ICNTL3.eq.0) then
        if (LA.lt.NE) then
          INFO(1) = -6
          INFO(2) = NE
        end if
      else
        if (LA.lt.1) then
          INFO(1) = -6
          INFO(2) = 1
        end if
      end if
      if (INFO(1).eq.-6) then
         if (LP.gt.0) then
            write (LP,FMT=9000) INFO(1)
            write (LP,FMT=9060) LA,INFO(2)
         end if
         GO TO 70
      end if

      if (ICNTL2.eq.0 .or. ICNTL2.eq.2) then
        if (LIP.lt.NC+1) then
          INFO(1) = -7
          INFO(2) = NC+1
        end if
      else if (LIP.lt.MAX(NR,NC)+1) then
        INFO(1) = -7
        INFO(2) = MAX(NR,NC)+1
      end if
      if (INFO(1).eq.-7) then
        if (LP.gt.0) then
          write (LP,FMT=9000) INFO(1)
          write (LP,FMT=9065) LIP,INFO(2)
        end if
        GO TO 70
      end if

      if (LIW.lt.MAX(NR,NC)+1) then
        INFO(1) = -8
        INFO(2) = MAX(NR,NC)+1
        if (LP.gt.0) then
          write (LP,FMT=9000) INFO(1)
          write (LP,FMT=9070) LIW,INFO(2)
        end if
        GO TO 70
      end if

      LAA = NE
      if (ICNTL3.ne.0) LAA = 1
      IOUT = 0
      JOUT = 0
      IDUP = 0
      IUP = 0

      PART = 0
      if (ICNTL6.ne.0) PART = 1

      if (ICNTL2.eq.0) then

call MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LAA,A,IP,IW, IOUT,JOUT,KNE)
        if (KNE.eq.0) GO TO 50

if (LCHECK) call MC59ED(NC,NR,NE,IRN,LIP,IP,LAA,A,IW,IDUP, KNE,ICNTL6)

      else if (ICNTL2.eq.1) then

        if (ICNTL6.ne.0) PART = -1
call MC59BD(LCHECK,PART,NR,NC,NE,JCN,IRN,LAA,A,IW,IP, JOUT,IOUT,KNE)
        if (KNE.eq.0) GO TO 50

if (LCHECK) call MC59ED(NR,NC,NE,JCN,NR+1,IW,LAA,A,IP, IDUP,KNE,ICNTL6)

        call MC59CD(NC,NR,KNE,IRN,JCN,LAA,A,IP,IW)

      else if (ICNTL2.eq.2) then
        if (LCHECK) then
call MC59FD(NC,NR,NE,IRN,NC+1,IP,LAA,A,LIW,IW,IDUP, IOUT,IUP,KNE,ICNTL6,INFO)
          if (INFO(1).eq.-9) GO TO 40
          if (KNE.eq.0) GO TO 50
        else
           KNE = NE
        end if

        call MC59DD(NC,KNE,IRN,IP,LAA,A)

      end if

      INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(6) = KNE
      INFO(7) = IUP

      if (IDUP.gt.0) INFO(1) = INFO(1) + 1
      if (IOUT.gt.0) INFO(1) = INFO(1) + 2
      if (JOUT.gt.0) INFO(1) = INFO(1) + 4
      if (INFO(1).gt.0 .and. MP.gt.0) then
        write (MP,FMT=9080) INFO(1)
        if (IOUT.gt.0) write (MP,FMT=9090) IOUT
        if (JOUT.gt.0) write (MP,FMT=9110) JOUT
        if (IDUP.gt.0) write (MP,FMT=9100) IDUP
        if (IUP.gt.0)  write (MP,FMT=9130) IUP
      end if
      GO TO 70

   40 INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(7) = IUP
      if (LP.gt.0) then
        write (LP,FMT=9000) INFO(1)
        write (LP,FMT=9140)
      end if
      GO TO 70

   50 INFO(1) = -10
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(2) = IOUT + JOUT
      if (LP.gt.0) then
        write (LP,FMT=9000) INFO(1)
        write (LP,FMT=9120)
      end if
   70 return

 9000 format (/,' *** Error return from MC59AD *** INFO(1) = ',I3)
 9010 format (1X,'ICNTL(2) = ',I2,' is out of range')
 9020 format (1X,'NC = ',I6,' is out of range')
 9030 format (1X,'NR = ',I6,' is out of range')
 9035 format (1X,'Symmetric case. NC = ',I6,' but NR = ',I6)
 9040 format (1X,'NE = ',I10,' is out of range')
 9050 format (1X,'Increase LJCN from ',I10,' to at least ',I10)
 9060 format (1X,'Increase LA from ',I10,' to at least ',I10)
 9065 format (1X,'Increase LIP from ',I8,' to at least ',I10)
 9070 format (1X,'Increase LIW from ',I8,' to at least ',I10)
 9080 format (/,' *** Warning message from MC59AD *** INFO(1) = ',I3)
9090 format (1X,I8,' entries in IRN supplied by the user were ', &
/,'       out of range and were ignored by the routine')
 9100 format (1X,I8,' duplicate entries were supplied by the user')
9110 format (1X,I8,' entries in JCN supplied by the user were ', &
/,'       out of range and were ignored by the routine')
 9120 format (1X,'All entries out of range')
9130 format (1X,I8,' of these entries were in the upper triangular ', /,'       part of matrix')
 9140 format (1X,'Entries in IP are not monotonic increasing')
 9150 format (1X,'ICNTL(6) = ',I2,' is out of range')
      END
subroutine MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LA,A,IP,IW,IOUT, JOUT,KNE)

      integer LA,NC,NE,NR,IOUT,JOUT,KNE,PART
      logical :: LCHECK
      double precision A(LA)
      integer IP(NC+1),IRN(NE),IW(NC+1),JCN(NE)
      double precision ACE,ACEP
      integer I,ICE,ICEP,J,JCE,JCEP,K,L,LOC

      do 10 J = 1,NC + 1
        IW(J) = 0
   10 continue

      KNE = 0
      IOUT = 0
      JOUT = 0
      if (LCHECK) then
        if (LA.gt.1) then
          if (PART.eq.0) then
            do 20 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              if (I.gt.NR .or. I.lt.1) then
                IOUT = IOUT + 1
                if (J.gt.NC .or. J.lt.1)  JOUT = JOUT + 1
              else if (J.gt.NC .or. J.lt.1) then
                JOUT = JOUT + 1
              else
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                A(KNE) = A(K)
                IW(J) = IW(J) + 1
              end if
   20       continue
          else if (PART.eq.1) then
            do 21 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              if (I.gt.NR .or. I.lt.1) then
                IOUT = IOUT + 1
                if (J.gt.NC .or. J.lt.1)  JOUT = JOUT + 1
              else if (J.gt.NC .or. J.lt.1) then
                JOUT = JOUT + 1
              else
                KNE = KNE + 1
                if (I.lt.J) then
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                else
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                end if
                A(KNE) = A(K)
              end if
   21       continue
          else if (PART.eq.-1) then
            do 22 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              if (I.gt.NR .or. I.lt.1) then
                IOUT = IOUT + 1
                if (J.gt.NC .or. J.lt.1)  JOUT = JOUT + 1
              else if (J.gt.NC .or. J.lt.1) then
                JOUT = JOUT + 1
              else
                KNE = KNE + 1
                if (I.gt.J) then
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                else
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                end if
                A(KNE) = A(K)
              end if
   22       continue
          end if
        else
          if (PART.eq.0) then
            do 25 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              if (I.gt.NR .or. I.lt.1) then
                IOUT = IOUT + 1
                if (J.gt.NC .or. J.lt.1)  JOUT = JOUT + 1
              else if (J.gt.NC .or. J.lt.1) then
                JOUT = JOUT + 1
              else
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                IW(J) = IW(J) + 1
              end if
   25       continue
          else if (PART.eq.1) then
            do 26 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              if (I.gt.NR .or. I.lt.1) then
                IOUT = IOUT + 1
                if (J.gt.NC .or. J.lt.1)  JOUT = JOUT + 1
              else if (J.gt.NC .or. J.lt.1) then
                JOUT = JOUT + 1
              else
                KNE = KNE + 1
                if (I.lt.J) then
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                else
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                end if
              end if
   26       continue
          else if (PART.eq.-1) then
            do 27 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              if (I.gt.NR .or. I.lt.1) then
                IOUT = IOUT + 1
                if (J.gt.NC .or. J.lt.1)  JOUT = JOUT + 1
              else if (J.gt.NC .or. J.lt.1) then
                JOUT = JOUT + 1
              else
                KNE = KNE + 1
                if (I.gt.J) then
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                else
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                end if
              end if
   27       continue
          end if
        end if
        if (KNE.eq.0) GO TO 130

      else

        KNE = NE
        if (PART.eq.0) then
          do 30 K = 1,NE
            J = JCN(K)
            IW(J) = IW(J) + 1
   30     continue
        else if (PART.eq.1) then
          do 35 K = 1,NE
            I = IRN(K)
            J = JCN(K)
            if (I.lt.J) then
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            else
              IW(J) = IW(J) + 1
            end if
   35     continue
        else if (PART.eq.-1) then
          do 36 K = 1,NE
            I = IRN(K)
            J = JCN(K)
            if (I.gt.J) then
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            else
              IW(J) = IW(J) + 1
            end if
   36     continue
        end if
      end if



      IP(1) = 1
      do 37 J = 2,NC + 1
        IP(J) = IW(J-1) + IP(J-1)
        IW(J-1) = IP(J-1)
   37 continue


      if (LA.eq.1) then
        do 70 L = 1,NC
          do 60 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            do 40 J = 1,NE
              if (JCE.eq.L) GO TO 50
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
   40       continue
   50       JCN(K) = JCE
            IRN(K) = ICE
   60     continue
   70   continue
      else

        do 120 L = 1,NC
          do 110 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            ACE = A(K)
            do 90 J = 1,NE
              if (JCE.eq.L) GO TO 100
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
   90       continue
  100       JCN(K) = JCE
            IRN(K) = ICE
            A(K) = ACE
  110     continue
  120   continue
      end if

  130 continue

      return
      END
      subroutine MC59CD(NC,NR,NE,IRN,JCN,LA,A,IP,IW)
      integer LA,NC,NE,NR
      double precision A(LA)
      integer IP(NC+1),IRN(NE),IW(NR+1),JCN(NE)
      double precision ACE,ACEP
      integer I,ICE,ICEP,J,J1,J2,K,L,LOC,LOCP


      do 10 J = 1,NC
        IP(J) = 0
   10 continue

      if (LA.gt.1) then

        do 20 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
          IRN(K) = JCN(K)
   20   continue
        IP(NC+1) = NE + 1


        IP(1) = IP(1) + 1
        do 30 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
   30   continue

        do 50 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          do 40 J = J1,J2
            K = IRN(J)
            L = IP(K) - 1
            JCN(J) = L
            IRN(J) = I
            IP(K) = L
   40     continue
   50   continue
        IP(NC+1) = NE + 1
        do 70 J = 1,NE
          LOC = JCN(J)
          if (LOC.eq.0) GO TO 70
          ICE = IRN(J)
          ACE = A(J)
          JCN(J) = 0
          do 60 K = 1,NE
            LOCP = JCN(LOC)
            ICEP = IRN(LOC)
            ACEP = A(LOC)
            JCN(LOC) = 0
            IRN(LOC) = ICE
            A(LOC) = ACE
            if (LOCP.eq.0) GO TO 70
            ICE = ICEP
            ACE = ACEP
            LOC = LOCP
   60     continue
   70   continue
      else



        do 90 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
   90   continue
        IP(NC+1) = NE + 1


        IP(1) = IP(1) + 1
        do 100 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
  100   continue

        do 120 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          do 110 J = J1,J2
            K = JCN(J)
            L = IP(K) - 1
            IRN(L) = I
            IP(K) = L
  110     continue
  120   continue

      end if

      return
      END


      subroutine MC59DD(NC,NE,IRN,IP,LA,A)
      integer LA,NC,NE
      double precision A(LA)
      integer IRN(NE),IP(NC)
      double precision ACE
      integer ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
      intrinsic ABS

      if (LA.gt.1) then
        KMAX = NE
        do 50 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          if (KLO.gt.KMAX) GO TO 40
          KOR = KMAX
          do 30 KDUMMY = KLO,KMAX
            ACE = A(KOR-1)
            ICE = IRN(KOR-1)
            do 10 K = KOR,KMAX
              IK = IRN(K)
              if (ABS(ICE).le.ABS(IK)) GO TO 20
              IRN(K-1) = IK
              A(K-1) = A(K)
   10       continue
            K = KMAX + 1
   20       IRN(K-1) = ICE
            A(K-1) = ACE
            KOR = KOR - 1
   30     continue
   40     KMAX = KLO - 2
   50   continue
      else

        KMAX = NE
        do 150 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          if (KLO.gt.KMAX) GO TO 140
          KOR = KMAX
          do 130 KDUMMY = KLO,KMAX
            ICE = IRN(KOR-1)
            do 110 K = KOR,KMAX
              IK = IRN(K)
              if (ABS(ICE).le.ABS(IK)) GO TO 120
              IRN(K-1) = IK
  110       continue
            K = KMAX + 1
  120       IRN(K-1) = ICE
            KOR = KOR - 1
  130     continue
  140     KMAX = KLO - 2
  150   continue
      end if
      END

      subroutine MC59ED(NC,NR,NE,IRN,LIP,IP,LA,A,IW,IDUP,KNE,ICNTL6)

      integer ICNTL6,IDUP,KNE,LIP,LA,NC,NR,NE
      double precision A(LA)
      integer IRN(NE),IP(LIP),IW(NR)
      integer I,J,K,KSTART,KSTOP,NZJ

      IDUP = 0
      KNE = 0
      do 10 I = 1,NR
        IW(I) = 0
   10 continue

      KSTART = IP(1)
      if (LA.gt.1) then
        NZJ = 0
        do 30 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          do 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            if (IW(I).le.NZJ) then
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            else
              IDUP = IDUP + 1
              if (ICNTL6.ge.0) A(IW(I)) = A(IW(I)) + A(K)
            end if
   20     continue
          KSTART = KSTOP
          NZJ = KNE
   30   continue

      else

        do 50 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          do 40 K = KSTART,KSTOP - 1
            I = IRN(K)
            if (IW(I).lt.J) then
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            else
              IDUP = IDUP + 1
            end if
   40     continue
          KSTART = KSTOP
   50   continue
      end if

      return
      END

subroutine MC59FD(NC,NR,NE,IRN,LIP,IP,LA,A,LIW,IW,IDUP,IOUT, IUP,KNE,ICNTL6,INFO)

      integer ICNTL6,IDUP,IOUT,IUP,KNE,LA,LIP,LIW,NC,NR,NE
      double precision A(LA)
      integer IRN(NE),IP(LIP),IW(LIW),INFO(2)
      integer I,J,K,KSTART,KSTOP,NZJ,LOWER

      IDUP = 0
      IOUT = 0
      IUP = 0
      KNE = 0
      do 10 I = 1,NR
        IW(I) = 0
   10 continue

      KSTART = IP(1)
      LOWER = 1
      if (LA.gt.1) then
        NZJ = 0
        do 30 J = 1,NC
          if (ICNTL6.ne.0) LOWER = J
          KSTOP = IP(J+1)
          if (KSTART.gt.KSTOP) then
            INFO(1) = -9
            INFO(2) = J
            return
          end if
          IP(J+1) = IP(J)
          do 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            if (I.gt.NR .or. I.lt.LOWER) then
              IOUT = IOUT + 1
              if (ICNTL6.ne.0 .and. I.lt.J) IUP = IUP + 1
            else if (IW(I).le.NZJ) then
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            else
              IDUP = IDUP + 1
              if (ICNTL6.ge.0) A(IW(I)) = A(IW(I)) + A(K)
            end if
   20     continue
          KSTART = KSTOP
          NZJ = KNE
   30   continue

      else

        do 50 J = 1,NC
          if (ICNTL6.ne.0) LOWER = J
          KSTOP = IP(J+1)
          if (KSTART.gt.KSTOP) then
            INFO(1) = -9
            INFO(2) = J
            return
          end if
          IP(J+1) = IP(J)
          do  40 K = KSTART,KSTOP - 1
            I = IRN(K)
            if (I.gt.NR .or. I.lt.LOWER) then
              IOUT = IOUT + 1
              if (ICNTL6.ne.0 .and. I.gt.1) IUP = IUP + 1
            else if (IW(I).lt.J) then
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            else
              IDUP = IDUP + 1
            end if
   40     continue
          KSTART = KSTOP
   50   continue
      end if

      return
      END



subroutine MA50AD(M,N,NE,LA,A,IRN,JCN,IQ,CNTL,ICNTL,IP,NP,JFIRST, LENR,LASTR,NEXTR,IW,IFIRST,LENC,LASTC,NEXTC, &
INFO,RINFO)



      integer M,N,NE,LA
      double precision A(LA)
      double precision CNTL(10)
      integer IRN(LA),JCN(LA),IQ(N)
integer ICNTL(20),IP(M),NP,JFIRST(M),LENR(M),LASTR(M),NEXTR(M), IW(M),IFIRST(N),LENC(N),LASTC(N),NEXTC(N),INFO(15)
      double precision RINFO(10)


      integer :: IDAMAX
      external IDAMAX
      external MA50DD
      intrinsic ABS,MAX,MIN

      double precision ZERO,ONE
      parameter (ZERO=0D0,ONE=1.0D0)

      double precision ALEN,AMULT,ANEW,ASW,AU,COST,CPIV
integer DISPC,DISPR,EYE,I,IDROP,IDUMMY,IEND,IFILL,IFIR,II,IJ, IJPOS,IOP,IPIV,IPOS,ISRCH,I1,I2,J,JBEG,JEND,JJ,JLAST, &
JMORE,JNEW,JPIV,JPOS,J1,J2,L,LC,LEN,LENPIV,LP,LR
      double precision MAXENT
integer MINC,MORD,MP,MSRCH,NC,NDROP,NEFACT,NEPR,NERED,NE1,NORD, NORD1,NR,NULLC,NULLI,NULLJ,NULLR,PIVBEG,PIVCOL,PIVEND, &
PIVOT
      double precision PIVR,PIVRAT,U


      LP = ICNTL(1)
      if (ICNTL(3).le.0) LP = 0
      MP = ICNTL(2)
      if (ICNTL(3).le.1) MP = 0
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = NE
      INFO(4) = NE
      INFO(5) = 0
      INFO(6) = 0
      INFO(7) = 0
      RINFO(1) = ZERO
      NP = 0

      if (M.lt.1 .or. N.lt.1) GO TO 690
      if (NE.lt.1) GO TO 700
      if (LA.lt.NE) then
         INFO(3) = NE
         GO TO 710
      end if

      if (MP.gt.0 .and. ICNTL(3).gt.2) then
write (MP,'(/2(A,I6),A,I8,A,I8/A,1P,4E10.2/A,7I4)') ' Entering MA50AD with M =',M,' N =',N,' NE =',NE,' LA =',LA, &
' CNTL =',(CNTL(I),I=1,4),' ICNTL =',(ICNTL(I),I=1,7)
         if (N.eq.1 .or. ICNTL(3).gt.3) then
            do 10 J = 1,N - 1
if (IQ(J).lt.IQ(J+1)) write (MP, '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J, (A(II),IRN(II),II=IQ(J),IQ(J+1)-1)
   10       continue
if (IQ(N).le.NE) write (MP,'(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',N, (A(II),IRN(II),II=IQ(N),NE)
         else
if (IQ(1).lt.IQ(2)) write (MP, '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',1, (A(II),IRN(II),II=IQ(1),IQ(2)-1)
         end if
         if (ICNTL(7).eq.2) then
            write (MP,'(A,(T10,10(I7)))') ' IP = ',IP
            write (MP,'(A,(T10,10(I7)))') ' IFIRST = ',IFIRST
         end if
      end if

      MINC = 1
      NERED = NE
      U = MIN(CNTL(2),ONE)
      U = MAX(U,ZERO)
      MSRCH = ICNTL(4)
      if (MSRCH.eq.0) MSRCH = N
      JLAST = N - ICNTL(6)
      if (JLAST.lt.1 .or. JLAST.gt.N) JLAST = N
      NULLI = 0
      NULLJ = 0
      MORD = 0
      NORD = 0
      NDROP = 0
      NEFACT = 0
      do 20 I = 1,N - 1
         LENC(I) = IQ(I+1) - IQ(I)
   20 continue
      LENC(N) = NE + 1 - IQ(N)

      if (CNTL(3).gt.ZERO) then
         NERED = 0
         do 40 J = 1,N
            I = IQ(J)
            IQ(J) = NERED + 1
            do 30 II = I,I + LENC(J) - 1
               if (ABS(A(II)).ge.CNTL(3)) then
                  NERED = NERED + 1
                  A(NERED) = A(II)
                  IRN(NERED) = IRN(II)
               else
                  INFO(6) = INFO(6) + 1
               end if
   30       continue
            LENC(J) = NERED + 1 - IQ(J)
   40    continue
      end if

      if (ICNTL(7).eq.2) then
         do 50 I = 1,M
            NEXTR(I) = IP(I)
   50    continue
         if (ICNTL(4).ne.1) GO TO 740
      end if

      DISPR = NERED + 1
      DISPC = NERED + 1
      do 60 I = 1,M
         IW(I) = 0
         LENR(I) = 0
         JFIRST(I) = 0
   60 continue
      do 70 II = 1,NERED
         I = IRN(II)
         LENR(I) = LENR(I) + 1
   70 continue
      IP(1) = LENR(1) + 1
      do 80 I = 2,M
         IP(I) = IP(I-1) + LENR(I)
   80 continue
      do 100 J = 1,N
         I = IQ(J)
         do 90 II = I,I + LENC(J) - 1
            I = IRN(II)
            if (IW(I).eq.J) GO TO 720
            IW(I) = J
            IPOS = IP(I) - 1
            JCN(IPOS) = J
            IP(I) = IPOS
   90    continue
  100 continue
      do 110 I = 1,M
         IW(I) = 0
  110 continue

      if (ICNTL(4).le.0) then
         do 120 I = 1,N
            IFIRST(I) = 0
  120    continue
         do 130 I = M,1,-1
            NE1 = LENR(I)
            if (NE1.gt.0) then
               IFIR = IFIRST(NE1)
               IFIRST(NE1) = I
               LASTR(I) = 0
               NEXTR(I) = IFIR
               if (IFIR.gt.0) LASTR(IFIR) = I
            else
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            end if
  130    continue
      else
         do 140 I = M,1,-1
            NE1 = LENR(I)
            if (NE1.eq.0) then
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            end if
  140    continue
      end if
      do 150 J = N,1,-1
         NE1 = LENC(J)
         if (NE1.eq.0) then
            if (ICNTL(7).ne.2) then
               if (J.le.JLAST) then
                  NORD = NORD + 1
                  IQ(J) = -NORD
                  if (NORD.eq.JLAST) then
                     NORD = NORD + NULLJ
                     JLAST = N
                     NULLJ = 0
                  end if
               else
                  NULLJ = NULLJ + 1
                  IQ(J) = - (JLAST+NULLJ)
               end if
               LASTC(J) = 0
               NEXTC(J) = 0
            end if
         else
            IFIR = JFIRST(NE1)
            JFIRST(NE1) = J
            NEXTC(J) = IFIR
            LASTC(J) = 0
            if (IFIR.gt.0) LASTC(IFIR) = J
         end if
  150 continue
      if (INFO(6).eq.0) then
         NULLC = NORD + NULLJ
         NULLR = NULLI
      end if

      do 630 PIVOT = 1,N
if (NERED.ge. (MIN(CNTL(1),ONE)*(N-NORD))* (M-MORD)) GO TO 640

         if (ICNTL(7).eq.2) then
            IPIV = 0
            J = IFIRST(PIVOT)
            if (J.lt.1 .or. J.gt.N) GO TO 730
            if (IQ(J).lt.0) GO TO 730
            len = LENC(J)
            if (len.le.0) GO TO 320
            ALEN = len - 1
            I1 = IQ(J)
            I2 = I1 + len - 1
            II = IDAMAX(len,A(I1),1)
            MAXENT = ABS(A(I1+II-1))
            if (MAXENT.le.CNTL(4)) GO TO 320
            AU = MAX(MAXENT*U,CNTL(4))
            do 160 II = I1,I2
               if (ABS(A(II)).lt.AU) GO TO 160
               I = IRN(II)
               if (IPIV.ne.0) then
                  if (NEXTR(I).ge.NEXTR(IPIV)) GO TO 160
               end if
               CPIV = ALEN*(LENR(I)-1)
               IJPOS = II
               IPIV = I
               JPIV = J
  160       continue
            GO TO 330
         end if

         len = MINC
         do 170 MINC = len,M - MORD
            if (JFIRST(MINC).ne.0) GO TO 180
            if (ICNTL(4).le.0) then
               if (IFIRST(MINC).ne.0) GO TO 180
            end if
  170    continue

  180    CPIV = M
         CPIV = CPIV*N
         PIVRAT = ZERO
         ISRCH = 0
         do 300 len = MINC,M - MORD
            ALEN = len - 1
            if (CPIV.le.ALEN**2 .and. ICNTL(4).le.0) GO TO 310
            IJ = JFIRST(len)
            do 220 IDUMMY = 1,N
               if (IJ.le.0) GO TO 230
               J = IJ
               IJ = NEXTC(J)
               if (J.gt.JLAST) GO TO 220
               MAXENT = ZERO
               I1 = IQ(J)
               I2 = I1 + len - 1
               II = IDAMAX(len,A(I1),1)
               MAXENT = ABS(A(I1+II-1))
               if (MAXENT.le.CNTL(4)) GO TO 320
               AU = MAX(MAXENT*U,CNTL(4))
               if (ICNTL(7).eq.1) then
                  do 190 II = I1,I2
                     if (IRN(II).eq.J) GO TO 200
  190             continue
                  GO TO 220
  200             I1 = II
                  I2 = II
               end if
               do 210 II = I1,I2
                  if (ABS(A(II)).lt.AU) GO TO 210
                  I = IRN(II)
                  COST = ALEN*(LENR(I)-1)
                  if (COST.gt.CPIV) GO TO 210
                  PIVR = ABS(A(II))/MAXENT
                  if (COST.eq.CPIV) then
                     if (PIVR.le.PIVRAT) GO TO 210
                  end if
                  CPIV = COST
                  IJPOS = II
                  IPIV = I
                  JPIV = J
                  if (CPIV.le.ALEN**2 .and. ICNTL(4).le.0) GO TO 330
                  PIVRAT = PIVR
  210          continue
               ISRCH = ISRCH + 1
               if (ISRCH.ge.MSRCH) then
                  if (PIVRAT.gt.ZERO) GO TO 330
               end if
  220       continue
  230       if (ICNTL(4).gt.0) GO TO 300
            if (CPIV.le.ALEN*(ALEN+1)) GO TO 310
            if (len.gt.N-NORD) GO TO 300
            IJ = IFIRST(len)
            do 290 IDUMMY = 1,M
               if (IJ.eq.0) GO TO 300
               I = IJ
               IJ = NEXTR(IJ)
               J1 = IP(I)
               J2 = J1 + len - 1
               if (ICNTL(7).eq.1) then
                  do 240 JJ = J1,J2
                     if (JCN(JJ).eq.I) GO TO 250
  240             continue
                  GO TO 290
  250             J1 = JJ
                  J2 = JJ
               end if
               do 280 JJ = J1,J2
                  J = JCN(JJ)
                  if (J.gt.JLAST) GO TO 280
                  COST = ALEN*(LENC(J)-1)
                  if (COST.ge.CPIV) GO TO 280
                  I1 = IQ(J)
                  I2 = I1 + LENC(J) - 1
                  II = IDAMAX(LENC(J),A(I1),1)
                  MAXENT = ABS(A(I1+II-1))
                  do 260 II = I1,I2 - 1
                     if (IRN(II).eq.I) GO TO 270
  260             continue
  270             JPOS = II
                  if (MAXENT.le.CNTL(4)) GO TO 320
                  if (ABS(A(JPOS)).lt.MAXENT*U) GO TO 280
                  CPIV = COST
                  IPIV = I
                  JPIV = J
                  IJPOS = JPOS
                  PIVRAT = ABS(A(JPOS))/MAXENT
                  if (CPIV.le.ALEN*(ALEN+1)) GO TO 330
  280          continue

  290       continue
  300    continue
  310    if (PIVRAT.gt.ZERO) GO TO 330
         INFO(1) = INFO(1) + 2
if (MP.gt.0) write (MP,'(A/A)') ' Warning message from MA50AD: no suitable diagonal pivot', &
' found, so switched to full matrix processing.'
         GO TO 640

  320    IPIV = 0
         JPIV = J

  330    NEFACT = NEFACT + LENC(JPIV)
         PIVBEG = IQ(JPIV)
         PIVEND = PIVBEG + LENC(JPIV) - 1
         NORD = NORD + 1
         NORD1 = NORD
         if (NORD.eq.JLAST) then
            NORD = NORD + NULLJ
            JLAST = N
            NULLJ = 0
         end if
         if (ICNTL(4).le.0) then
            do 340 II = PIVBEG,PIVEND
               I = IRN(II)
               LR = LASTR(I)
               NR = NEXTR(I)
               if (NR.ne.0) LASTR(NR) = LR
               if (LR.eq.0) then
                  NE1 = LENR(I)
                  IFIRST(NE1) = NR
               else
                  NEXTR(LR) = NR
               end if
  340       continue
         end if
         if (IPIV.gt.0) then
            NEPR = LENR(IPIV) - 1
            NEFACT = NEFACT + NEPR
            RINFO(1) = RINFO(1) + CPIV*2 + LENR(IPIV)
            J1 = IP(IPIV)
            do 350 JJ = J1,J1 + NEPR
               J = JCN(JJ)
               LC = LASTC(J)
               NC = NEXTC(J)
               if (NC.ne.0) LASTC(NC) = LC
               if (LC.eq.0) then
                  NE1 = LENC(J)
                  JFIRST(NE1) = NC
               else
                  NEXTC(LC) = NC
               end if
  350       continue
            if (PIVBEG.ne.IJPOS) then
               ASW = A(PIVBEG)
               A(PIVBEG) = A(IJPOS)
               A(IJPOS) = ASW
               IRN(IJPOS) = IRN(PIVBEG)
               IRN(PIVBEG) = IPIV
            end if
         else
            NEPR = 0
            NE1 = LENC(JPIV)
            if (CNTL(3).gt.ZERO) NDROP = NDROP + NE1
            if (NE1.gt.0) then
               LC = LASTC(JPIV)
               NC = NEXTC(JPIV)
               if (NC.ne.0) LASTC(NC) = LC
               if (LC.eq.0) then
                  JFIRST(NE1) = NC
               else
                  NEXTC(LC) = NC
               end if
            end if
         end if
         do 360 II = PIVBEG + 1,PIVEND
            I = IRN(II)
            IW(I) = II - PIVBEG
  360    continue
         LENPIV = PIVEND - PIVBEG
         do 390 II = PIVBEG,PIVEND
            I = IRN(II)
            LENR(I) = LENR(I) - 1
            J1 = IP(I)
            J2 = J1 + LENR(I)
            do 370 JJ = J1,J2 - 1
               if (JCN(JJ).eq.JPIV) GO TO 380
  370       continue
  380       JCN(JJ) = JCN(J2)
            JCN(J2) = 0
  390    continue

         do 600 EYE = 1,NEPR
            J = JCN(IP(IPIV)+EYE-1)
            IDROP = 0
            JBEG = IQ(J)
            JEND = JBEG + LENC(J) - 1
            do 400 II = JBEG,JEND - 1
               if (IRN(II).eq.IPIV) GO TO 410
  400       continue
  410       AMULT = -A(II)/A(IQ(JPIV))
            A(II) = A(JEND)
            IRN(II) = IRN(JEND)
            LENC(J) = LENC(J) - 1
            IRN(JEND) = 0
            JEND = JEND - 1
            if (LENPIV.eq.0) GO TO 600
            IOP = 0
!DIR$ IVDEP
            do 420 II = JBEG,JEND
               I = IRN(II)
               if (IW(I).gt.0) then
                  IOP = IOP + 1
                  PIVCOL = IQ(JPIV) + IW(I)
                  IW(I) = -IW(I)
                  A(II) = A(II) + AMULT*A(PIVCOL)
               end if
  420       continue

            if (CNTL(3).gt.ZERO) then
               JNEW = JBEG
               do 450 II = JBEG,JEND
                  if (ABS(A(II)).ge.CNTL(3)) then
                     A(JNEW) = A(II)
                     IRN(JNEW) = IRN(II)
                     JNEW = JNEW + 1
                  else
                     I = IRN(II)
                     J1 = IP(I)
                     J2 = J1 + LENR(I) - 1
                     do 430 JJ = J1,J2 - 1
                        if (JCN(JJ).eq.J) GO TO 440
  430                continue
  440                JCN(JJ) = JCN(J2)
                     JCN(J2) = 0
                     LENR(I) = LENR(I) - 1
                  end if
  450          continue
               do 460 II = JNEW,JEND
                  IRN(II) = 0
  460          continue
               IDROP = JEND + 1 - JNEW
               JEND = JNEW - 1
               LENC(J) = LENC(J) - IDROP
               NERED = NERED - IDROP
               INFO(6) = INFO(6) + IDROP
            end if

            IFILL = LENPIV - IOP
            NERED = NERED + IFILL
            INFO(3) = MAX(INFO(3),NERED+LENC(J))

            if (IFILL.eq.0) then
!DIR$ IVDEP
               do 470 II = PIVBEG + 1,PIVEND
                  I = IRN(II)
                  IW(I) = -IW(I)
  470          continue
               GO TO 600
            end if

            do 480 IPOS = JEND + 1,MIN(JEND+IFILL,DISPC-1)
               if (IRN(IPOS).ne.0) GO TO 490
  480       continue
            if (IPOS.eq.JEND+IFILL+1) GO TO 540
            if (JEND+IFILL+1.le.LA+1) then
               DISPC = JEND + IFILL + 1
               GO TO 540
            end if
            IPOS = LA
            DISPC = LA + 1
  490       JMORE = JEND + IFILL - IPOS + 1
            do 500 IPOS = JBEG - 1,MAX(JBEG-JMORE,1),-1
               if (IRN(IPOS).ne.0) GO TO 510
  500       continue
            IPOS = IPOS + 1
            if (IPOS.eq.JBEG-JMORE) GO TO 520
  510       if (DISPC+LENC(J)+IFILL.gt.LA+1) then
               INFO(2) = INFO(2) + 1
               call MA50DD(LA,A,IRN,IQ,N,DISPC,.true.)
               JBEG = IQ(J)
               JEND = JBEG + LENC(J) - 1
               PIVBEG = IQ(JPIV)
               PIVEND = PIVBEG + LENC(JPIV) - 1
               if (DISPC+LENC(J)+IFILL.gt.LA+1) GO TO 705
            end if
            IPOS = DISPC
            DISPC = DISPC + LENC(J) + IFILL
  520       IQ(J) = IPOS
            do 530 II = JBEG,JEND
               A(IPOS) = A(II)
               IRN(IPOS) = IRN(II)
               IPOS = IPOS + 1
               IRN(II) = 0
  530       continue
            JBEG = IQ(J)
            JEND = IPOS - 1
  540       IDROP = 0
            do 580 II = PIVBEG + 1,PIVEND
               I = IRN(II)
               INFO(3) = MAX(INFO(3),NERED+LENR(I)+1)
               if (IW(I).lt.0) then
                  IW(I) = -IW(I)
                  GO TO 580
               end if
               ANEW = AMULT*A(II)
               if (ABS(ANEW).lt.CNTL(3)) then
                  IDROP = IDROP + 1
               else
                  JEND = JEND + 1
                  A(JEND) = ANEW
                  IRN(JEND) = I

                  IEND = IP(I) + LENR(I)
                  if (IEND.lt.DISPR) then
                     if (JCN(IEND).eq.0) GO TO 560
                  else
                     if (DISPR.le.LA) then
                        DISPR = DISPR + 1
                        GO TO 560
                     end if
                  end if
                  if (IP(I).gt.1) then
                     if (JCN(IP(I)-1).eq.0) then
                        IEND = IEND - 1
                        do 545 JJ = IP(I),IEND
                           JCN(JJ-1) = JCN(JJ)
  545                   continue
                        IP(I) = IP(I) - 1
                        GO TO 560
                     end if
                  end if
                  if (DISPR+LENR(I).gt.LA) then
                     INFO(2) = INFO(2) + 1
                     call MA50DD(LA,A,JCN,IP,M,DISPR,.false.)
                     if (DISPR+LENR(I).gt.LA) GO TO 705
                  end if
                  J1 = IP(I)
                  J2 = IP(I) + LENR(I) - 1
                  IP(I) = DISPR
                  do 550 JJ = J1,J2
                     JCN(DISPR) = JCN(JJ)
                     JCN(JJ) = 0
                     DISPR = DISPR + 1
  550             continue
                  IEND = DISPR
                  DISPR = IEND + 1
  560             JCN(IEND) = J
                  LENR(I) = LENR(I) + 1
               end if
  580       continue
            INFO(6) = INFO(6) + IDROP
            NERED = NERED - IDROP
            do 590 II = 1,IDROP
               IRN(JEND+II) = 0
  590       continue
            LENC(J) = LENC(J) + IFILL - IDROP
  600    continue


         do 610 EYE = 1,NEPR
            JJ = IP(IPIV) + EYE - 1
            J = JCN(JJ)
            JCN(JJ) = 0
            NE1 = LENC(J)
            LASTC(J) = 0
            if (NE1.gt.0) then
               IFIR = JFIRST(NE1)
               JFIRST(NE1) = J
               NEXTC(J) = IFIR
               if (IFIR.ne.0) LASTC(IFIR) = J
               MINC = MIN(MINC,NE1)
            else if (ICNTL(7).ne.2) then
               if (INFO(6).eq.0) NULLC = NULLC + 1
               if (J.le.JLAST) then
                  NORD = NORD + 1
                  IQ(J) = -NORD
                  if (NORD.eq.JLAST) then
                     NORD = NORD + NULLJ
                     JLAST = N
                     NULLJ = 0
                  end if
               else
                  NULLJ = NULLJ + 1
                  IQ(J) = - (JLAST+NULLJ)
               end if
            end if
  610    continue
         NERED = NERED - NEPR

         if (IPIV.ne.0) then
            LENR(IPIV) = 0
            IW(IPIV) = 0
            IRN(PIVBEG) = 0
            MORD = MORD + 1
            PIVBEG = PIVBEG + 1
            IP(IPIV) = -MORD
         end if
         NERED = NERED - LENPIV - 1
         do 620 II = PIVBEG,PIVEND
            I = IRN(II)
            IW(I) = 0
            IRN(II) = 0
            NE1 = LENR(I)
            if (NE1.eq.0) then
               if (INFO(6).eq.0) NULLR = NULLR + 1
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            else if (ICNTL(4).le.0) then
               IFIR = IFIRST(NE1)
               LASTR(I) = 0
               NEXTR(I) = IFIR
               IFIRST(NE1) = I
               if (IFIR.ne.0) LASTR(IFIR) = I
               MINC = MIN(MINC,NE1)
            end if
  620    continue
         IQ(JPIV) = -NORD1
  630 continue


  640 INFO(5) = MORD + MIN(M-MORD-NULLI,N-NORD-NULLJ)
      do 650 L = 1,MIN(M-MORD,N-NORD)
RINFO(1) = RINFO(1) + M - MORD - L + 1 + REAL(M-MORD-L)*(N-NORD-L)*2
  650 continue
      NP = NORD
INFO(4) = 2 + NEFACT + M*2 + MAX(N-NORD+M-MORD, (N-NORD)*(M-MORD))
      INFO(6) = INFO(6) + NDROP
      INFO(7) = M - MORD
      do 660 L = 1,M
         if (IP(L).lt.0) then
            IP(L) = -IP(L)
         else
            MORD = MORD + 1
            IP(L) = MORD
         end if
  660 continue
      do 670 L = 1,N
         if (IQ(L).lt.0) then
            LASTC(L) = -IQ(L)
         else
            if (NORD.eq.JLAST) NORD = NORD + NULLJ
            NORD = NORD + 1
            LASTC(L) = NORD
         end if
  670 continue
      do 680 L = 1,N
         IQ(LASTC(L)) = L
  680 continue

      if (INFO(5).lt.MIN(M,N)) INFO(1) = INFO(1) + 1

      if (MP.gt.0 .and. ICNTL(3).gt.2) then
write (MP,'(A,I6,A,F12.1/A,7I8)') ' Leaving MA50AD with NP =', NP,' RINFO(1) =',RINFO(1),' INFO =',(INFO(I),I=1,7)
         if (ICNTL(3).gt.3) then
            write (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            write (MP,'(A,(T6,10(I7)))') ' IQ = ',IQ
         end if
      end if

      GO TO 750

  690 INFO(1) = -1
if (LP.gt.0) write (LP,'(/A/(2(A,I8)))') ' **** Error return from MA50AD ****',' M =',M,' N =',N
      GO TO 750
  700 INFO(1) = -2
if (LP.gt.0) write (LP,'(/A/(A,I10))') ' **** Error return from MA50AD ****',' NE =',NE
      GO TO 750
  705 INFO(4) =  NEFACT + NERED
      INFO(6) = INFO(6) + NDROP
  710 INFO(1) = -3
if (LP.gt.0) write (LP,'(/A/A,I9,A,I9)') ' **** Error return from MA50AD ****', &
' LA  must be increased from',LA,' to at least',INFO(3)
      GO TO 750
  720 INFO(1) = -4
if (LP.gt.0) write (LP,'(/A/(3(A,I9)))') ' **** Error return from MA50AD ****',' Entry in row',I, &
' and column',J,' duplicated'
      GO TO 750
  730 INFO(1) = -5
if (LP.gt.0) write (LP,'(/A/(3(A,I9)))') ' **** Error return from MA50AD ****',' Fault in component ', &
PIVOT,' of column permutation given in IFIRST'
      GO TO 750
  740 INFO(1) = -6
if (LP.gt.0) write (LP,'(/A/(3(A,I9)))') ' **** Error return from MA50AD ****',' ICNTL(4) = ',ICNTL(4), &
' when ICNTL(6) = 2'
  750 END


subroutine MA50BD(M,N,NE,JOB,AA,IRNA,IPTRA,CNTL,ICNTL,IP,IQ,NP, LFACT,FACT,IRNF,IPTRL,IPTRU,W,IW,INFO,RINFO)
      integer M,N,NE,JOB
      double precision AA(NE)
      integer IRNA(NE),IPTRA(N)
      double precision CNTL(10)
      integer ICNTL(20),IP(M),IQ(*),NP,LFACT
      double precision FACT(LFACT)
      integer IRNF(LFACT),IPTRL(N),IPTRU(N)
      double precision W(M)
      integer IW(M+2*N),INFO(15)
      double precision RINFO(10)

      double precision ZERO,ONE
      parameter (ZERO=0D0,ONE=1.0D0)
      double precision AMULT,ASW
      integer :: BEGCOL
      logical :: DROP
integer ENDCOL,EYE,EYE1,I,IA1,IA2,IF1,IF2,II,IL1,IL2,IPIV,IQPIV, IU1,IU2,ISW,J,JDUMMY,JJ,JLAST,K,LP
      double precision MAXENT
      integer MF,MORD,MP,NEU,NF,NULLC
      double precision PIVLIM
      integer :: RANK
      double precision U
      external MA50ED,MA50FD,MA50GD
      intrinsic ABS,MAX,MIN

      INFO(1) = 0
      INFO(4) = 0
      INFO(5) = 0
      INFO(6) = 0
      INFO(7) = 0
      RINFO(1) = ZERO
      LP = ICNTL(1)
      MP = ICNTL(2)
      if (ICNTL(3).le.0) LP = 0
      if (ICNTL(3).le.1) MP = 0

      if (M.lt.1 .or. N.lt.1) then
         INFO(1) = -1
if (LP.gt.0) write (LP,'(/A/A,I8,A,I8)') ' **** Error return from MA50BD ****',' M =',M,' N =',N
         GO TO 550
      end if
      if (NE.le.0) then
         INFO(1) = -2
if (LP.gt.0) write (LP,'(/A/A,I6)') ' **** Error return from MA50BD ****',' NE =',NE
         GO TO 550
      end if
      if (NP.lt.0 .or. NP.gt.N) then
         INFO(1) = -7
if (LP.gt.0) write (LP,'(/A/A,I8,A,I8)') ' **** Error return from MA50BD ****',' NP =',NP,' N =',N
         GO TO 550
      end if
      if (LFACT.lt.MAX(M,NE+2)) then
         INFO(4) = MAX(M,NE+2)
         GO TO 520
      end if
      if (JOB.eq.1) then
      else if (JOB.eq.2 .or. JOB.eq.3) then
         if (IRNF(1).ne.0) then
            INFO(1) = -6
if (LP.gt.0) write (LP,'(/A/A,I1,A)') ' **** Error return from MA50BD ***',' Call with JOB=', &
JOB,' follows JOB=1 call in which entries were dropped'
            GO TO 550
         end if
      else
         INFO(1) = -5
if (LP.gt.0) write (LP,'(/A/A,I2)') ' **** Error return from MA50BD ****',' JOB =',JOB
         GO TO 550
      end if

      if (MP.gt.0) then
if (ICNTL(3).gt.2) write (MP, '(/2(A,I6),A,I8,A,I3/A,I8,A,I7/A,1P,4E10.2/A,7I8)') &
' Entering MA50BD with M =',M,' N =',N,' NE =',NE,' JOB =', JOB,' LFACT =',LFACT,' NP =',NP,' CNTL =',(CNTL(I),I=1,4), &
' ICNTL =',(ICNTL(I),I=1,7)
         if (ICNTL(3).gt.3) then
            write (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            if (IQ(1).gt.0) then
               write (MP,'(A,(T6,10(I7)))') ' IQ = ',(IQ(J),J=1,N)
            else
               write (MP,'(A,(T6,I7))') ' IQ = ',IQ(1)
            end if
            do 10 J = 1,N - 1
if (IPTRA(J).lt.IPTRA(J+1)) write (MP, '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J, &
(AA(II),IRNA(II),II=IPTRA(J),IPTRA(J+1)-1)
   10       continue
if (IPTRA(N).le.NE) write (MP, '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',N, (AA(II),IRNA(II),II=IPTRA(N),NE)
         end if
      end if

      JLAST = 0
      NULLC = 0
if (JOB.gt.1 .and. ICNTL(6).gt.0 .and. ICNTL(6).LT.N) JLAST = MIN(NP,N-ICNTL(6))

      U = MIN(CNTL(2),ONE)
      U = MAX(U,ZERO)
      do 20 I = 1,M
         IW(I+N) = 0
         W(I) = ZERO
   20 continue
      MORD = 0
      IF1 = LFACT + 1
      IF2 = 0
      NF = N - NP
      MF = 0
      IL2 = 2
      if (JLAST.gt.0) IL2 = IPTRL(JLAST)
      NEU = 0

      if (JOB.eq.2) GO TO 370

      if (JOB.eq.3) then
         do 30 J = 1,NP
            IA1 = IPTRU(J) + 1
            if (IA1.gt.IPTRL(J)) GO TO 30
            if (J.le.JLAST) then
               MORD = MORD + 1
               IP(IRNF(IA1)) = -J
            else
               IP(IRNF(IA1)) = J
            end if
   30    continue
         MF = IRNF(2)
         IA1 = IPTRL(N)
         do 40 J = 1,MF
            IP(IRNF(IA1+J)) = NP + J
   40    continue
      end if

      do 50 K = 1,JLAST
         IW(M+N+K) = IPTRL(K)
   50 continue

      do 310 K = JLAST + 1,N
         DROP = .false.
         if (K.eq.NP+1) then
            MF = M - MORD
            IF1 = LFACT + 1 - MF
            II = 0
            do 60 I = 1,M
               if (IP(I).gt.0) then
                  IW(I+N) = N
                  IRNF(IF1+II) = I
                  II = II + 1
                  IP(I) = NP + II
               end if
   60       continue
            IF1 = LFACT + 1 - MAX(MF*NF,MF+NF)
            IF2 = IF1 - 1 + MF*MAX(0,JLAST-NP)
         end if
         J = K
         if (IQ(1).gt.0) J = IQ(K)
         IA1 = IPTRA(J)
         IA2 = NE
         if (J.ne.N) IA2 = IPTRA(J+1) - 1
         IU1 = IL2 + 1
         IU2 = IU1 - 1
         IL1 = IF1 - 1 + IA1 - IA2
         IL2 = IL1 - 1
         INFO(4) = MAX(INFO(4),NEU+LFACT-IL1+IU2+M+1)
         if (IL1-IU2.le.M) then
            if (INFO(1).ne.-3) then
               INFO(1) = -3
               if (ICNTL(8).ne.0) GO TO 480
               NEU = IL2 + LFACT + 1 - MF - IF1
               IF1 = LFACT + 1 - MF
               IF2 = IF1 - 1
               IL2 = 0
               EYE = 0
               do 80 J = 1,MIN(K-1,NP)
                  IU2 = IPTRU(J)
                  IPTRU(J) = EYE
                  IL2 = IPTRL(J)
                  NEU = NEU + IU2 - IL2
                  do 70 II = IU2 + 1,IL2
                     EYE = EYE + 1
                     IRNF(EYE) = IRNF(II)
                     FACT(EYE) = FACT(II)
   70             continue
                  IPTRL(J) = EYE
                  IW(M+N+J) = EYE
   80          continue
               IU1 = EYE + 1
               IU2 = EYE
               IL1 = IF1 - 1 + IA1 - IA2
               IL2 = IL1 - 1
            end if
            if (IL1-IU2.le.M) GO TO 480
         end if
         EYE = IL1
         do 90 II = IA1,IA2
            I = IRNA(II)
            if (IW(I+N).eq.-1) GO TO 540
            IW(I+N) = -1
            W(I) = AA(II)
            IRNF(EYE) = I
            EYE = EYE + 1
   90    continue
         IPTRL(K) = EYE - 1
         IW(M+N+K) = EYE - 1
         IW(K) = IL1
         J = K
         do 120 JDUMMY = 1,2*K
            do 100 II = IW(J),ABS(IW(M+N+J))
               I = IRNF(II)
               if (IW(I+N).ge.K) GO TO 100
               if (IP(I).le.0) GO TO 110
               IW(I+N) = K
               IL1 = IL1 - 1
               IRNF(IL1) = I
  100       continue
            if (J.eq.K) GO TO 130
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
  120    continue
  130    do 150 II = IU2,IU1,-1
            I = IRNF(II)
            J = -IP(I)
            EYE1 = IPTRU(J) + 1
            if (ABS(W(I)).lt.CNTL(3)) GO TO 150
            AMULT = -W(I)*FACT(EYE1)
            W(I) = AMULT
            do 140 EYE = EYE1 + 1,IPTRL(J)
               I = IRNF(EYE)
               W(I) = W(I) + AMULT*FACT(EYE)
  140       continue
            RINFO(1) = RINFO(1) + ONE + 2*(IPTRL(J)-EYE1)
  150    continue

         if (CNTL(3).gt.ZERO) then
            EYE = IU1
            do 160 II = IU1,IU2
               I = IRNF(II)
               if (ABS(W(I)).lt.CNTL(3)) then
                  INFO(6) = INFO(6) + 1
               else
                  IRNF(EYE) = -IP(I)
                  FACT(EYE) = W(I)
                  EYE = EYE + 1
               end if
               W(I) = ZERO
  160       continue
            IU2 = EYE - 1
         else
            do 170 II = IU1,IU2
               I = IRNF(II)
               IRNF(II) = -IP(I)
               FACT(II) = W(I)
               W(I) = ZERO
  170       continue
         end if
         if (INFO(1).eq.-3) then
            NEU = NEU + IU2 - IU1 + 1
            IU2 = IU1 - 1
         end if
         IPTRU(K) = IU2
         if (K.le.NP) then
            MAXENT = ZERO
            if (CNTL(3).gt.ZERO) then
               EYE = IL1
               do 180 II = IL1,IL2
                  I = IRNF(II)
                  if (ABS(W(I)).lt.CNTL(3)) then
                     INFO(6) = INFO(6) + 1
                     W(I) = ZERO
                     DROP = .true.
                  else
                     IRNF(EYE) = I
                     EYE = EYE + 1
                     MAXENT = MAX(ABS(W(I)),MAXENT)
                  end if
  180          continue
               IL2 = EYE - 1
            else
               do 190 II = IL1,IL2
                  MAXENT = MAX(ABS(W(IRNF(II))),MAXENT)
  190          continue
            end if
            PIVLIM = U*MAXENT
            EYE = IU2
            IQPIV = M + N
            if (IL1.gt.IL2) NULLC = NULLC + 1
            do 200 II = IL1,IL2
               I = IRNF(II)
               EYE = EYE + 1
               IRNF(EYE) = I
               FACT(EYE) = W(I)
               W(I) = ZERO
               if (ABS(FACT(EYE)).ge.PIVLIM) then
                  if (ABS(FACT(EYE)).gt.CNTL(4)) then
                     if (IP(I).lt.IQPIV) then
                        IQPIV = IP(I)
                        IPIV = EYE
                     end if
                  end if
               end if
  200       continue
            IL1 = IU2 + 1
            IL2 = EYE
            if (IL1.le.IL2) then
               if (IQPIV.eq.M+N) then
                  if (CNTL(3).gt.ZERO) INFO(6) = INFO(6) + EYE - IU2
                  IL2 = IU2
               else
                  if (IL1.ne.IPIV) then
                     ASW = FACT(IPIV)
                     FACT(IPIV) = FACT(IL1)
                     FACT(IL1) = ASW
                     ISW = IRNF(IL1)
                     IRNF(IL1) = IRNF(IPIV)
                     IRNF(IPIV) = ISW
                  end if
                  INFO(5) = INFO(5) + 1
                  FACT(IL1) = ONE/FACT(IL1)
                  RINFO(1) = RINFO(1) + ONE
                  MORD = MORD + 1
                  IP(IRNF(IL1)) = -K
               end if
            end if
         else
            IL2 = IPTRU(K)
!DIR$ IVDEP
            do 210 II = LFACT - MF + 1,LFACT
               I = IRNF(II)
               IF2 = IF2 + 1
               FACT(IF2) = W(I)
               W(I) = ZERO
  210       continue
            if (INFO(1).eq.-3) IF2 = IF2 - MF
         end if
         IW(M+N+K) = IL2
         IPTRL(K) = IL2
         if (DROP) GO TO 310
         do 300 II = IU1,IU2
            I = IRNF(II)
            if (IW(M+N+I).lt.0) GO TO 300
            BEGCOL = IPTRU(I) + 2
            ENDCOL = IPTRL(I)
            if (K.le.NP) then
               do 220 JJ = BEGCOL,ENDCOL
                  if (IP(IRNF(JJ)).eq.-K) GO TO 230
  220          continue
               GO TO 300
            end if
  230       do 280 JDUMMY = BEGCOL,ENDCOL
               JJ = BEGCOL
               do 240 BEGCOL = JJ,ENDCOL
                  if (IP(IRNF(BEGCOL)).gt.0) GO TO 250
  240          continue
               GO TO 290
  250          JJ = ENDCOL
               do 260 ENDCOL = JJ,BEGCOL,-1
                  if (IP(IRNF(ENDCOL)).lt.0) GO TO 270
  260          continue
               GO TO 290
  270          ASW = FACT(BEGCOL)
               FACT(BEGCOL) = FACT(ENDCOL)
               FACT(ENDCOL) = ASW
               J = IRNF(BEGCOL)
               IRNF(BEGCOL) = IRNF(ENDCOL)
               IRNF(ENDCOL) = J
               BEGCOL = BEGCOL + 1
               ENDCOL = ENDCOL - 1
  280       continue
  290       IW(M+N+I) = -ENDCOL
  300    continue
  310 continue
      if (N.eq.NP) then
         MF = M - MORD
         IF1 = LFACT + 1 - MF
         II = 0
         do 320 I = 1,M
            if (IP(I).gt.0) then
               IW(I+N) = N
               IRNF(IF1+II) = I
               II = II + 1
               IP(I) = NP + II
            end if
  320    continue
         IF1 = LFACT + 1 - MAX(MF*NF,MF+NF)
         IF2 = IF1 - 1 + MF*MAX(0,JLAST-NP)
      end if
      if (INFO(5).eq.MIN(M,N)) then
         do 330 I = 1,M
            IP(I) = ABS(IP(I))
  330    continue
      else
         MORD = NP
         do 340 I = 1,M
            if (IP(I).lt.0) then
               IP(I) = -IP(I)
            else
               MORD = MORD + 1
               IP(I) = MORD
            end if
  340    continue
      end if
      IRNF(1) = INFO(6)
      IRNF(2) = MF
      INFO(7) = MF
      FACT(1) = CNTL(3)
      FACT(2) = CNTL(4)
      if (INFO(1).eq.-3) GO TO 520
      IF2 = IF2 - MF*NF
      do 350 II = 1,MF*NF
         FACT(IL2+II) = FACT(IF1-1+II)
  350 continue
      do 360 II = 1,MF
         IRNF(IL2+II) = IRNF(LFACT-MF+II)
  360 continue
      IF1 = IL2 + 1
      GO TO 440
  370 MF = IRNF(2)
      IF1 = IPTRL(N) + 1
      IF2 = IF1 - 1
      do 430 K = JLAST + 1,N
         J = K
         if (IQ(1).gt.0) J = IQ(K)
         IA1 = IPTRA(J)
         IA2 = NE
         if (J.ne.N) IA2 = IPTRA(J+1) - 1
         IU1 = IL2 + 1
         IU2 = IPTRU(K)
         IL1 = IU2 + 1
         IL2 = IPTRL(K)
         do 380 II = IA1,IA2
            W(IRNA(II)) = AA(II)
  380    continue
         do 400 II = IU2,IU1,-1
            J = IRNF(II)
            I = IRNF(IPTRU(J)+1)
            EYE1 = IPTRU(J) + 1
            AMULT = -W(I)*FACT(EYE1)
            FACT(II) = AMULT
            W(I) = ZERO
            do 390 EYE = EYE1 + 1,IPTRL(J)
               I = IRNF(EYE)
               W(I) = W(I) + AMULT*FACT(EYE)
  390       continue
            RINFO(1) = RINFO(1) + ONE + 2*(IPTRL(J)-EYE1)
  400    continue
         if (K.le.NP) then
            if (IL1.le.IL2) then
!DIR$ IVDEP
               do 410 II = IL1,IL2
                  I = IRNF(II)
                  FACT(II) = W(I)
                  W(I) = ZERO
  410          continue
               if (ABS(FACT(IL1)).le.CNTL(4)) then
                  GO TO 530
               else
                  INFO(5) = INFO(5) + 1
                  FACT(IL1) = ONE/FACT(IL1)
                  RINFO(1) = RINFO(1) + ONE
               end if
            end if
         else
            do 420 II = IF1,IF1 + MF - 1
               I = IRNF(II)
               IF2 = IF2 + 1
               FACT(IF2) = W(I)
               W(I) = ZERO
  420       continue
         end if
  430 continue
      INFO(4) = MAX(IF1+MF+NF-1,IF2)

  440 if (MF.gt.0 .and. NF.gt.0) then
if (ICNTL(5).gt.1) call MA50GD(MF,NF,FACT(IF1),MF,ICNTL(5), CNTL(4),IRNF(IF1+MF),RANK)
if (ICNTL(5).eq.1) call MA50FD(MF,NF,FACT(IF1),MF,CNTL(4), IRNF(IF1+MF),RANK)
if (ICNTL(5).le.0) call MA50ED(MF,NF,FACT(IF1),MF,CNTL(4), IRNF(IF1+MF),RANK)
         INFO(5) = INFO(5) + RANK
         do 450 I = 1,MIN(MF,NF)
            RINFO(1) = RINFO(1) + MF - I + 1 + real(MF-I)*(NF-I)*2
  450    continue
      end if
      if (INFO(5).lt.MIN(M,N)) INFO(1) = 1
      if (MP.gt.0 .and. ICNTL(3).gt.2) then
write (MP,'(A,I6,A,F12.1/A,I3,A,4I8)') ' Leaving MA50BD with IRNF(2) =',IRNF(2), ' RINFO(1) =',RINFO(1), &
' INFO(1) =',INFO(1),' INFO(4:7) =', (INFO(J),J=4,7)
         if (ICNTL(3).gt.3) then
            if (JOB.ne.2) write (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            do 460 J = 1,N
               if (J.gt.1) then
if (IPTRL(J-1).lt.IPTRU(J)) write (MP, '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J, &
' of U', (FACT(II),IRNF(II),II=IPTRL(J-1)+1, IPTRU(J))
               end if
if (IPTRU(J).lt.IPTRL(J)) write (MP, '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L', &
(FACT(II),IRNF(II),II=IPTRU(J)+1,IPTRL(J))
  460       continue
            write (MP,'(A)') ' Full part'
            write (MP,'((6I12))') (IRNF(IF1+MF+J),J=0,NF-1)
            do 470 I = 0,MF - 1
write (MP,'(I4,1P,6E12.4:/(4X,1P,6E12.4))') IRNF(IF1+I), (FACT(IF1+I+J*MF),J=0,NF-1)
  470       continue
         end if
      end if
      GO TO 550

  480 do 490 I = 1,M
         IW(I) = 0
  490 continue
      do 500 I = 1,M
         if (IP(I).gt.0) then
            IW(IP(I)) = I
         else
            IP(I) = -IP(I)
         end if
  500 continue
      do 510 I = 1,M
         if (IW(I).gt.0) then
            IP(IW(I)) = K
            K = K + 1
         end if
  510 continue
  520 INFO(1) = -3
      if (LP.gt.0) then
         write (LP,'(/A)')' **** Error return from MA50BD **** '
         if (ICNTL(8).eq.0) then
write (LP,'(A,I7,A,I7)')' LFACT must be increased from', LFACT,' to at least',INFO(4)
         else
            write (LP,'(A,I7)')' LFACT must be increased from',LFACT
         end if
      end if
      GO TO 550
  530 INFO(1) = - (7+K)
if (LP.gt.0) write (LP,'(/A/A,I6,A)') ' **** Error return from MA50BD **** ', &
' Small pivot found in column',K,' of the permuted matrix.'
      GO TO 550
  540 INFO(1) = -4
if (LP.gt.0) write (LP,'(/A/(3(A,I9)))') ' **** Error return from MA50BD ****',' Entry in row',I, &
' and column',J,' duplicated'
  550 END

subroutine MA50CD(M,N,ICNTL,IQ,NP,TRANS,LFACT,FACT,IRNF,IPTRL, IPTRU,B,X,W,INFO)
      integer M,N,ICNTL(20),IQ(*),NP
      logical :: TRANS
      integer :: LFACT
      double precision FACT(LFACT)
      integer IRNF(LFACT),IPTRL(N),IPTRU(N)
      double precision B(*),X(*),W(*)
      integer :: INFO(15)

      double precision ZERO
      parameter (ZERO=0D0)
      integer I,II,IA1,IF1,J,LP,MF,MP,NF
      double precision PROD

      external MA50HD

      LP = ICNTL(1)
      MP = ICNTL(2)
      if (ICNTL(3).le.0) LP = 0
      if (ICNTL(3).le.1) MP = 0

      if (M.lt.1 .or. N.lt.1) GO TO 250

if (MP.gt.0 .and. ICNTL(3).gt.2) write (MP, '(/2(A,I6),A,I4,A,L2)') ' Entering MA50CD with M=',M,' N =',N, &
' NP =',NP,' TRANS =',TRANS
      IF1 = IPTRL(N) + 1
      MF = IRNF(2)
      NF = N - NP
if (MP.gt.0 .and. ICNTL(3).gt.2) write (MP, '(A,I5,A,I5)') ' Size of full submatrix',MF,' by',NF
      if (MP.gt.0 .and. ICNTL(3).gt.3) then
         do 10 J = 1,N
            if (J.gt.1) then
if (IPTRL(J-1).lt.IPTRU(J)) write (MP, '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of U', &
(FACT(II),IRNF(II),II=IPTRL(J-1)+1,IPTRU(J))
            end if
if (IPTRU(J).lt.IPTRL(J)) write (MP, '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L', &
(FACT(II),IRNF(II),II=IPTRU(J)+1,IPTRL(J))
   10    continue
         write (MP,'(A)') ' Full part'
         write (MP,'((6I12))') (IRNF(IF1+MF+J),J=0,NF-1)
         do 20 I = 0,MF - 1
write (MP,'(I4,1P,6E12.4:/(4X,1P,6E12.4))') IRNF(IF1+I), (FACT(IF1+I+J*MF),J=0,NF-1)
   20    continue
      end if

      if (TRANS) then
if (MP.gt.0 .and. ICNTL(3).gt.3) write (MP, '(A4,5F10.4:/(4X,5F10.4))') ' B =', (B(I),I=1,N)
         if (IQ(1).gt.0) then
            do 30 I = 1,N
               W(I) = B(IQ(I))
   30       continue
         else
            do 40 I = 1,N
               W(I) = B(I)
   40       continue
         end if
         do 50 I = 1,M
            X(I) = ZERO
   50    continue
         do 70 I = 2,N
            PROD = ZERO
            do 60 II = IPTRL(I-1) + 1,IPTRU(I)
               PROD = PROD + FACT(II)*W(IRNF(II))
   60       continue
            W(I) = W(I) + PROD
   70    continue
         do 80 I = 1,NF
            X(I) = W(NP+I)
   80    continue
         if (MF.gt.0 .and. NF.gt.0) then
call MA50HD(TRANS,MF,NF,FACT(IF1),MF,IRNF(IF1+MF),X, ICNTL(5))
         else
            do 90 I = 1,MF
               X(I) = ZERO
   90       continue
         end if
         do 100 I = MF,1,-1
            J = IRNF(IF1+I-1)
            if (J.ne.I) X(J) = X(I)
  100    continue
         do 120 I = NP,1,-1
            IA1 = IPTRU(I) + 1
            if (IA1.gt.IPTRL(I)) GO TO 120
            PROD = ZERO
            do 110 II = IA1 + 1,IPTRL(I)
               PROD = PROD + FACT(II)*X(IRNF(II))
  110       continue
            X(IRNF(IA1)) = (W(I)-PROD)*FACT(IA1)
  120    continue
if (MP.gt.0 .and. ICNTL(3).gt.3) write (MP, '(A/(4X,5F10.4))') ' Leaving MA50CD with X =', (X(I),I=1,M)
      else
if (MP.gt.0 .and. ICNTL(3).gt.3) write (MP, '(A4,5F10.4:/(4X,5F10.4))') ' B =', (B(I),I=1,M)
         do 130 I = 1,M
            W(I) = B(I)
  130    continue
         do 150 I = 1,NP
            IA1 = IPTRU(I) + 1
            if (IA1.le.IPTRL(I)) then
               X(I) = W(IRNF(IA1))*FACT(IA1)
               if (X(I).ne.ZERO) then
!DIR$ IVDEP
                  do 140 II = IA1 + 1,IPTRL(I)
                     W(IRNF(II)) = W(IRNF(II)) - FACT(II)*X(I)
  140             continue
               end if
            end if
  150    continue
         if (MF.gt.0 .and. NF.gt.0) then
            do 160 I = 1,MF
               W(I) = W(IRNF(IF1+I-1))
  160       continue
call MA50HD(TRANS,MF,NF,FACT(IF1),MF,IRNF(IF1+MF),W, ICNTL(5))
            do 170 I = 1,NF
               X(NP+I) = W(I)
  170       continue
         else
            do 180 I = 1,NF
               X(NP+I) = ZERO
  180       continue
         end if
         do 200 J = N,MAX(2,NP+1),-1
            PROD = X(J)
!DIR$ IVDEP
            do 190 II = IPTRL(J-1) + 1,IPTRU(J)
               X(IRNF(II)) = X(IRNF(II)) + FACT(II)*PROD
  190       continue
  200    continue
         do 220 J = NP,2,-1
            IA1 = IPTRU(J)
            if (IA1.ge.IPTRL(J)) then
               X(J) = ZERO
            else
               PROD = X(J)
!DIR$ IVDEP
               do 210 II = IPTRL(J-1) + 1,IA1
                  X(IRNF(II)) = X(IRNF(II)) + FACT(II)*PROD
  210          continue
            end if
  220    continue
         if (NP.ge.1 .and. IPTRU(1).ge.IPTRL(1)) X(1) = ZERO
         if (IQ(1).gt.0) then
            do 230 I = 1,N
               W(I) = X(I)
  230       continue
            do 240 I = 1,N
               X(IQ(I)) = W(I)
  240       continue
         end if
if (MP.gt.0 .and. ICNTL(3).gt.3) write (MP, '(A/(4X,5F10.4))') ' Leaving MA50CD with X =', (X(I),I=1,N)
      end if
      return
  250 INFO(1) = -1
if (LP.gt.0) write (LP,'(/A/2(A,I8))') ' **** Error return from MA50CD ****',' M =',M,' N =',N
      END

      subroutine MA50DD(LA,A,IND,IPTR,N,DISP,REALS)
      integer LA,N,DISP
      double precision A(LA)
      integer :: IPTR(N)
      logical :: REALS
      integer :: IND(LA)
      integer J,K,KN
      do 10 J = 1,N
         K = IPTR(J)
         if (K.gt.0) then
            IPTR(J) = IND(K)
            IND(K) = -J
         end if
   10 continue
      KN = 0
      do 20 K = 1,DISP - 1
         if (IND(K).eq.0) GO TO 20
         KN = KN + 1
         if (REALS) A(KN) = A(K)
         if (IND(K).le.0) then
            J = -IND(K)
            IND(K) = IPTR(J)
            IPTR(J) = KN
         end if
         IND(KN) = IND(K)
   20 continue
      DISP = KN + 1
      END


      subroutine MA50ED(M,N,A,LDA,PIVTOL,IPIV,RANK)
      integer LDA,M,N,RANK
      double precision PIVTOL

      integer :: IPIV(N)
      double precision A(LDA,N)



      double precision ONE,ZERO
      parameter (ONE=1.0D+0,ZERO=0.0D+0)

      integer I,J,JP,K
      logical :: PIVOT

      integer :: IDAMAX
      external IDAMAX

      external DAXPY,DSCAL,DSWAP
      intrinsic ABS

      J = 1
      do 30 K = 1,N

         do 10 I = 1,J - 1
            if (M.gt.I) call DAXPY(M-I,-A(I,J),A(I+1,I),1,A(I+1,J),1)
   10    continue

         if (J.le.M) then
            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .gt. PIVTOL
         else
            PIVOT = .false.
         end if
         if (PIVOT) then

            if (JP.ne.J) call DSWAP(N+J-K,A(J,1),LDA,A(JP,1),LDA)
            if (J.lt.M) call DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)

            J = J + 1
         else
            do 20 I = J,M
               A(I,J) = ZERO
   20       continue
            if (K.lt.N) call DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            IPIV(N-K+J) = -J
         end if
   30 continue

      RANK = J - 1
      END


      subroutine MA50FD(M,N,A,LDA,PIVTOL,IPIV,RANK)
      integer LDA,M,N,RANK
      double precision PIVTOL

      integer :: IPIV(N)
      double precision A(LDA,N)



      double precision ONE,ZERO
      parameter (ONE=1.0D+0,ZERO=0.0D+0)

      integer I,J,JP,K
      logical :: PIVOT

      integer :: IDAMAX
      external IDAMAX

      external DGEMV,DSCAL,DSWAP

      intrinsic ABS

      J = 1
      do 20 K = 1,N

         if (J.le.M) then
call DGEMV('No transpose',M-J+1,J-1,-ONE,A(J,1),LDA,A(1,J), 1,ONE,A(J,J),1)
            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .gt. PIVTOL
         else
            PIVOT = .false.
         end if

         if (PIVOT) then

            if (JP.ne.J) call DSWAP(N+J-K,A(J,1),LDA,A(JP,1),LDA)
            if (J.lt.M) call DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)

            if (J.lt.N) then
call DGEMV('Transpose',J-1,N-J,-ONE,A(1,J+1),LDA,A(J,1), LDA,ONE,A(J,J+1),LDA)
            end if

            J = J + 1
         else
            do 10 I = J,M
               A(I,J) = ZERO
   10       continue
            if (K.lt.N) call DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            IPIV(N-K+J) = -J
         end if
   20 continue

      RANK = J - 1
      END


      subroutine MA50GD(M,N,A,LDA,NB,PIVTOL,IPIV,RANK)
      integer LDA,M,N,NB,RANK
      double precision PIVTOL

      integer :: IPIV(N)
      double precision A(LDA,N)

      double precision ONE,ZERO
      parameter (ONE=1.0D+0,ZERO=0.0D+0)

      integer I,J,JJ,JP,J1,J2,K
      logical :: PIVOT
      double precision TEMP


      external DGEMM,DGEMV,DSWAP,DSCAL,DTRSM,DTRSV

      integer :: IDAMAX
      external IDAMAX

      intrinsic ABS,MIN

      J = 1
      J1 = 1
      J2 = MIN(N,NB)
      do 70 K = 1,N

         if (J.le.M) then

call DGEMV('No transpose',M-J+1,J-J1,-ONE,A(J,J1),LDA, A(J1,J),1,ONE,A(J,J),1)

            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .gt. PIVTOL
         else
            PIVOT = .false.
         end if

         if (PIVOT) then

            if (JP.ne.J) call DSWAP(J2-J1+1,A(J,J1),LDA,A(JP,J1),LDA)
            if (J.lt.M) call DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)

            if (J+1.le.J2) then
call DGEMV('Transpose',J-J1,J2-J,-ONE,A(J1,J+1),LDA, A(J,J1),LDA,ONE,A(J,J+1),LDA)
            end if

            J = J + 1
         else

            do 10 I = J,M
               A(I,J) = ZERO
   10       continue
            IPIV(N-K+J) = -J
            if (K.ne.N) call DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            if (N-K+J.gt.J2) then
               do 20 I = J1,J - 1
                  JP = IPIV(I)
                  TEMP = A(I,J)
                  A(I,J) = A(JP,J)
                  A(JP,J) = TEMP
   20          continue
if (J.gt.J1) call DTRSV('Lower','No transpose','Unit', J-J1,A(J1,J1),LDA,A(J1,J),1)
            else
               J2 = J2 - 1
            end if
         end if

         if (J.gt.J2) then
            do 40 JJ = 1,J1 - 1
               do 30 I = J1,J2
                  JP = IPIV(I)
                  TEMP = A(I,JJ)
                  A(I,JJ) = A(JP,JJ)
                  A(JP,JJ) = TEMP
   30          continue
   40       continue
            do 60 JJ = J2 + 1,N - K + J - 1
               do 50 I = J1,J2
                  JP = IPIV(I)
                  TEMP = A(I,JJ)
                  A(I,JJ) = A(JP,JJ)
                  A(JP,JJ) = TEMP
   50          continue
   60       continue

            if (K.ne.N) then
call DTRSM('Left','Lower','No transpose','Unit',J2-J1+1, N-K,ONE,A(J1,J1),LDA,A(J1,J2+1),LDA)
if (M.gt.J2) call DGEMM('No transpose','No transpose', M-J2,N-K,J2-J1+1,-ONE,A(J2+1,J1), LDA,A(J1,J2+1),LDA,ONE, &
A(J2+1,J2+1),LDA)
            end if

            J1 = J2 + 1
            J2 = MIN(J2+NB,N-K+J-1)

         end if

   70 continue
      RANK = J - 1
      END

      subroutine MA50HD(TRANS,M,N,A,LDA,IPIV,B,ICNTL5)
      logical :: TRANS
      integer LDA,M,N,ICNTL5
      integer :: IPIV(N)
      double precision A(LDA,N),B(*)


      integer I,K,RANK

      double precision ZERO
      parameter (ZERO=0D0)
      double precision TEMP
      intrinsic MIN
      external DAXPY,DDOT,DTRSV
      double precision DDOT

      RANK = 0
      do 10 RANK = MIN(M,N),1,-1
         if (IPIV(RANK).gt.0) GO TO 20
   10 continue

   20 if (.not.TRANS) then
         do 30 I = 1,RANK
            K = IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   30    continue
         if (ICNTL5.gt.0) then
if (RANK.gt.0) call DTRSV('L','NoTrans','Unit',RANK,A,LDA,B, 1)
         else
            do 40 K = 1,RANK - 1
if (B(K).ne.ZERO) call DAXPY(RANK-K,-B(K),A(K+1,K),1, B(K+1),1)
   40       continue
         end if

         if (ICNTL5.gt.0) then
if (RANK.gt.0) call DTRSV('U','NoTrans','NonUnit',RANK,A, LDA,B,1)
         else
            do 50 K = RANK,2,-1
               if (B(K).ne.ZERO) then
                  B(K) = B(K)/A(K,K)
                  call DAXPY(K-1,-B(K),A(1,K),1,B(1),1)
               end if
   50       continue
            if (RANK.gt.0) B(1) = B(1)/A(1,1)
         end if

         do 60 K = RANK + 1,N
            B(K) = ZERO
   60    continue
         do 70 I = RANK + 1,N
            K = -IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   70    continue

      else

         do 80 I = N,RANK + 1,-1
            K = -IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   80    continue
         if (ICNTL5.gt.0) then
if (RANK.gt.0) call DTRSV('U','Trans','NonUnit',RANK,A,LDA, B,1)
         else
            if (RANK.gt.0) B(1) = B(1)/A(1,1)
            do 90 I = 2,RANK
               TEMP = B(I) - DDOT(I-1,A(1,I),1,B(1),1)
               B(I) = TEMP/A(I,I)
   90       continue
         end if

         if (ICNTL5.gt.0) then
            if (RANK.gt.0) call DTRSV('L','Trans','Unit',RANK,A,LDA,B,1)
         else
            do 100 I = RANK - 1,1,-1
               B(I) = B(I) - DDOT(RANK-I,A(I+1,I),1,B(I+1),1)
  100       continue
         end if

         do 110 I = RANK + 1,M
            B(I) = ZERO
  110    continue
         do 120 I = RANK,1,-1
            K = IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
  120    continue
      end if

      END

      subroutine MA50ID(CNTL,ICNTL)

      double precision CNTL(10)
      integer I,ICNTL(20)

      CNTL(1) = 0.5D0
      CNTL(2) = 0.1D0
      do 10 I = 3,10
         CNTL(I) = 0.0D0
   10 continue

      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 1
      ICNTL(4) = 3
      ICNTL(5) = 32
      do 20 I = 6,20
        ICNTL(I) = 0
   20 continue

      END

      subroutine MC13DD(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW)
      integer LICN,N,NUM
      integer IB(N),ICN(LICN),IOR(N),IP(N),IW(N,3),LENR(N)
      external MC13ED
      call MC13ED(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW(1,1),IW(1,2),IW(1,3))
      return

      END
      subroutine MC13ED(N,ICN,LICN,IP,LENR,ARP,IB,NUM,LOWL,NUMB,PREV)
      integer LICN,N,NUM
integer ARP(N),IB(N),ICN(LICN),IP(N),LENR(N),LOWL(N),NUMB(N), PREV(N)
      integer DUMMY,I,I1,I2,ICNT,II,ISN,IST,IST1,IV,IW,J,K,LCNT,NNM1,STP
      intrinsic MIN
      ICNT = 0
      NUM = 0
      NNM1 = N + N - 1
      do 20 J = 1,N
        NUMB(J) = 0
        ARP(J) = LENR(J) - 1
   20 continue
      do 120 ISN = 1,N
        if (NUMB(ISN).ne.0) GO TO 120
        IV = ISN
        IST = 1
        LOWL(IV) = 1
        NUMB(IV) = 1
        IB(N) = IV
        do 110 DUMMY = 1,NNM1
          I1 = ARP(IV)
          if (I1.lt.0) GO TO 60
          I2 = IP(IV) + LENR(IV) - 1
          I1 = I2 - I1
          do 50 II = I1,I2
            IW = ICN(II)
            if (NUMB(IW).eq.0) GO TO 100
            LOWL(IV) = MIN(LOWL(IV),LOWL(IW))
   50     continue
          ARP(IV) = -1
   60     if (LOWL(IV).lt.NUMB(IV)) GO TO 90
          NUM = NUM + 1
          IST1 = N + 1 - IST
          LCNT = ICNT + 1
          do 70 STP = IST1,N
            IW = IB(STP)
            LOWL(IW) = N + 1
            ICNT = ICNT + 1
            NUMB(IW) = ICNT
            if (IW.eq.IV) GO TO 80
   70     continue
   80     IST = N - STP
          IB(NUM) = LCNT
          if (IST.ne.0) GO TO 90
          if (ICNT.lt.N) GO TO 120
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
  110   continue
  120 continue
  130 do 140 I = 1,N
        II = NUMB(I)
        ARP(II) = I
  140 continue
      return

      END

      subroutine MC21AD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
      integer LICN,N,NUMNZ
      integer ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
      external MC21BD
call MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2), IW(1,3),IW(1,4))
      return
      END
      subroutine MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
      integer LICN,N,NUMNZ
      integer ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
      integer I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
      do 10 I = 1,N
        ARP(I) = LENR(I) - 1
        CV(I) = 0
        IPERM(I) = 0
   10 continue
      NUMNZ = 0
      do 100 JORD = 1,N
        J = JORD
        PR(J) = -1
        do 70 K = 1,JORD
          IN1 = ARP(J)
          if (IN1.lt.0) GO TO 30
          IN2 = IP(J) + LENR(J) - 1
          IN1 = IN2 - IN1
          do 20 II = IN1,IN2
            I = ICN(II)
            if (IPERM(I).eq.0) GO TO 80
   20     continue
          ARP(J) = -1
   30     continue
          OUT(J) = LENR(J) - 1
          do 60 KK = 1,JORD
            IN1 = OUT(J)
            if (IN1.lt.0) GO TO 50
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
            do 40 II = IN1,IN2
              I = ICN(II)
              if (CV(I).eq.JORD) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       continue
   50       continue
            J = PR(J)
            if (J.eq.-1) GO TO 100
   60     continue
   70   continue
   80   continue
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUMNZ = NUMNZ + 1
        do 90 K = 1,JORD
          J = PR(J)
          if (J.eq.-1) GO TO 100
          II = IP(J) + LENR(J) - OUT(J) - 2
          I = ICN(II)
          IPERM(I) = J
   90   continue
  100 continue
      if (NUMNZ.eq.N) return
      do 110 I = 1,N
        ARP(I) = 0
  110 continue
      K = 0
      do 130 I = 1,N
        if (IPERM(I).ne.0) GO TO 120
        K = K + 1
        OUT(K) = I
        GO TO 130
  120   continue
        J = IPERM(I)
        ARP(J) = I
  130 continue
      K = 0
      do 140 I = 1,N
        if (ARP(I).ne.0) GO TO 140
        K = K + 1
        IOUTK = OUT(K)
        IPERM(IOUTK) = I
  140 continue
      return
      END

      subroutine MC29AD(M,N,NE,A,IRN,ICN,R,C,W,LP,IFAIL)
      integer M,N,NE
      double precision A(NE)
      integer IRN(NE),ICN(NE)
      double precision R(M),C(N),W(M*2+N*3)
      integer LP,IFAIL

      intrinsic LOG,ABS,MIN

      integer :: MAXIT
      parameter (MAXIT=100)
      double precision ONE,SMIN,ZERO
      parameter (ONE=1D0,SMIN=0.1,ZERO=0D0)

      integer I,I1,I2,I3,I4,I5,ITER,J,K
      double precision E,E1,EM,Q,Q1,QM,S,S1,SM,U,V

      IFAIL = 0
      if (M.lt.1 .or. N.lt.1) then
         IFAIL = -1
         GO TO 220
      else if (NE.le.0) then
         IFAIL = -2
         GO TO 220
      end if

      I1 = 0
      I2 = M
      I3 = M + N
      I4 = M + N*2
      I5 = M + N*3

      do 10 I = 1,M
         R(I) = ZERO
         W(I1+I) = ZERO
   10 continue
      do 20 J = 1,N
         C(J) = ZERO
         W(I2+J) = ZERO
         W(I3+J) = ZERO
         W(I4+J) = ZERO
   20 continue

      do 30 K = 1,NE
         U = ABS(A(K))
         if (U.eq.ZERO) GO TO 30
         I = IRN(K)
         J = ICN(K)
         if (MIN(I,J).lt.1 .or. I.gt.M .or. J.gt.N) GO TO 30
         U = LOG(U)
         W(I1+I) = W(I1+I) + ONE
         W(I2+J) = W(I2+J) + ONE
         R(I) = R(I) + U
         W(I3+J) = W(I3+J) + U
   30 continue
      do 40 I = 1,M
         if (W(I1+I).eq.ZERO) W(I1+I) = ONE
         R(I) = R(I)/W(I1+I)
         W(I5+I) = R(I)
   40 continue
      do 50 J = 1,N
         if (W(I2+J).eq.ZERO) W(I2+J) = ONE
         W(I3+J) = W(I3+J)/W(I2+J)
   50 continue
      SM = SMIN*NE

      do 60 K = 1,NE
         if (A(K).eq.ZERO) GO TO 60
         I = IRN(K)
         J = ICN(K)
         if (MIN(I,J).lt.1 .or. I.gt.M .or. J.gt.N) GO TO 60
         R(I) = R(I) - W(I3+J)/W(I1+I)
   60 continue
      E = ZERO
      Q = ONE
      S = ZERO
      do 70 I = 1,M
         S = S + W(I1+I)*R(I)**2
   70 continue
      if (S.le.SM) GO TO 160

      do 150 ITER = 1,MAXIT
         do 80 K = 1,NE
            if (A(K).eq.ZERO) GO TO 80
            J = ICN(K)
            I = IRN(K)
            if (MIN(I,J).lt.1 .or. I.gt.M .or. J.gt.N) GO TO 80
            C(J) = C(J) + R(I)
   80    continue
         S1 = S
         S = ZERO
         do 90 J = 1,N
            V = -C(J)/Q
            C(J) = V/W(I2+J)
            S = S + V*C(J)
   90    continue
         E1 = E
         E = Q*S/S1
         Q = ONE - E
         if (S.le.SM) E = ZERO
         do 100 I = 1,M
            R(I) = R(I)*E*W(I1+I)
  100    continue
         if (S.le.SM) GO TO 180
         EM = E*E1
         do 110 K = 1,NE
            if (A(K).eq.ZERO) GO TO 110
            I = IRN(K)
            J = ICN(K)
            if (MIN(I,J).lt.1 .or. I.gt.M .or. J.gt.N) GO TO 110
            R(I) = R(I) + C(J)
  110    continue
         S1 = S
         S = ZERO
         do 120 I = 1,M
            V = -R(I)/Q
            R(I) = V/W(I1+I)
            S = S + V*R(I)
  120    continue
         E1 = E
         E = Q*S/S1
         Q1 = Q
         Q = ONE - E
         if (S.le.SM) Q = ONE
         QM = Q*Q1
         do 130 J = 1,N
            W(I4+J) = (EM*W(I4+J)+C(J))/QM
            W(I3+J) = W(I3+J) + W(I4+J)
  130    continue
         if (S.le.SM) GO TO 160
         do 140 J = 1,N
            C(J) = C(J)*E*W(I2+J)
  140    continue
  150 continue
  160 do 170 I = 1,M
         R(I) = R(I)*W(I1+I)
  170 continue
  180 do 190 K = 1,NE
         if (A(K).eq.ZERO) GO TO 190
         I = IRN(K)
         J = ICN(K)
         if (MIN(I,J).lt.1 .or. I.gt.M .or. J.gt.N) GO TO 190
         R(I) = R(I) + W(I3+J)
  190 continue
      do 200 I = 1,M
         R(I) = R(I)/W(I1+I) - W(I5+I)
  200 continue
      do 210 J = 1,N
         C(J) = -W(I3+J)
  210 continue
      return

220 if (LP.gt.0) write (LP,'(/A/A,I3)') ' **** Error return from MC29AD ****',' IFAIL =',IFAIL

      END


      subroutine MC71AD(N,KASE,X,EST,W,IW,KEEP)
      integer :: ITMAX
      parameter (ITMAX=5)
      double precision ZERO,ONE
      parameter (ZERO=0.0D0,ONE=1.0D0)
      double precision EST
      integer KASE,N
      double precision W(*),X(*)
      integer IW(*),KEEP(5)
      double precision ALTSGN,TEMP
      integer I,ITER,J,JLAST,JUMP
      integer :: IDAMAX
      external IDAMAX
      intrinsic ABS,SIGN,NINT,DBLE
      if (N.le.0) then
        KASE = -1
        return

      end if

      if (KASE.eq.0) then
        do 10 I = 1,N
          X(I) = ONE/DBLE(N)
   10   continue
        KASE = 1
        JUMP = 1
        KEEP(1) = JUMP
        KEEP(2) = 0
        KEEP(3) = 0
        KEEP(4) = 0
        return

      end if
      JUMP  = KEEP(1)
      ITER  = KEEP(2)
      J     = KEEP(3)
      JLAST = KEEP(4)
      GO TO (100,200,300,400,500) JUMP
  100 continue
      if (N.eq.1) then
        W(1) = X(1)
        EST = ABS(W(1))
        GO TO 510

      end if
      do 110 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  110 continue
      KASE = 2
      JUMP = 2
      GO TO 1010
  200 continue
      J = IDAMAX(N,X,1)
      ITER = 2
  220 continue
      do 230 I = 1,N
        X(I) = ZERO
  230 continue
      X(J) = ONE
      KASE = 1
      JUMP = 3
      GO TO 1010
  300 continue
      do 310 I = 1,N
        W(I) = X(I)
  310 continue
      do 320 I = 1,N
        if (NINT(SIGN(ONE,X(I))).ne.IW(I)) GO TO 330
  320 continue
      GO TO 410
  330 continue
      do 340 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  340 continue
      KASE = 2
      JUMP = 4
      GO TO 1010
  400 continue
      JLAST = J
      J = IDAMAX(N,X,1)
      if ((ABS(X(JLAST)).ne.ABS(X(J))) .and. (ITER.lt.ITMAX)) then
        ITER = ITER + 1
        GO TO 220

      end if
  410 continue
      EST = ZERO
      do 420 I = 1,N
        EST = EST + ABS(W(I))
  420 continue
      ALTSGN = ONE
      do 430 I = 1,N
        X(I) = ALTSGN* (ONE+DBLE(I-1)/DBLE(N-1))
        ALTSGN = -ALTSGN
  430 continue
      KASE = 1
      JUMP = 5
      GO TO 1010
  500 continue
      TEMP = ZERO
      do 520 I = 1,N
        TEMP = TEMP + ABS(X(I))
  520 continue
      TEMP = 2.0*TEMP/DBLE(3*N)
      if (TEMP.gt.EST) then
        do 530 I = 1,N
          W(I) = X(I)
  530   continue
        EST = TEMP
      end if
  510 KASE = 0
 1010 continue
      KEEP(1) = JUMP
      KEEP(2) = ITER
      KEEP(3) = J
      KEEP(4) = JLAST
      return
      END
subroutine MA48AD(M,N,NE,JOB,LA,A,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO, RINFO)
      integer M,N,NE,JOB,LA
      double precision A(LA)
      integer IRN(LA),JCN(LA),KEEP(*)
      double precision CNTL(10)
      double precision RINFO(10)
      integer ICNTL(20),IW(6*M+3*N),INFO(20)


      double precision ZERO
      parameter (ZERO=0.D0)
      double precision CNTL5(10)
integer EYE,HEADC,I,IB,ICNTL5(20),IDUMMY,INFO5(15),IP,IPTRA,IPTRD, &
IPTRO,IPTRP,IQ,ISW,IW13,IW21,IW50,J,JAY,JB,JFIRST,J1,J2, J3,K,KB,KBLOCK,KD,KK,KL,KO,L,LASTR,LASTC,LBLOCK
      logical :: LDUP
integer LENC,LENP,LENR,LP,MBLOCK,MINBLK,MP,NB,NBLOCK,NC,NDIAG, NEXTC,NEXTR,NP,NR,NXTEYE,NZA,NZB,NZD,PTRD,PTRO
      double precision RINFO5(10),TOL
      external MC13DD,MC21AD,MA50AD
      intrinsic ABS,MAX
      LP = ICNTL(1)
      MP = ICNTL(2)
      if (M.le.0 .or. N.le.0) then
         INFO(1) = -1
         if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,9100) M,N
         GO TO 530
      end if
      if (NE.le.0) then
         INFO(1) = -2
         if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,9110) NE
         GO TO 530
      end if
      if (LA.lt.2*NE) then
         INFO(1) = -3
         INFO(3) = 2*NE
         if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,9120) LA,INFO(3)
         GO TO 530
      end if
      if (JOB.lt.1 .or. JOB.gt.3) then
         INFO(1) = -6
         if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,9130) JOB
         GO TO 530
      end if
      if (JOB.eq.2) then
         do 10 I = 1,MAX(M,N)
            IW(N+I) = 0
   10    continue
         do 20 I = 1,M
            J = KEEP(I)
            if (J.lt.1 .or. J.gt.M) GO TO 40
            if (IW(N+J).eq.1) GO TO 40
            IW(N+J) = 1
   20    continue
         do 30 I = 1,N
            J = KEEP(M+I)
            if (J.lt.1 .or. J.gt.N) GO TO 40
            if (IW(N+J).eq.2) GO TO 40
            IW(N+J) = 2
   30    continue
         GO TO 50
   40    INFO(1) = -5
         if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,9140)
         GO TO 530
      end if

   50 if (MP.gt.0 .and. ICNTL(3).gt.2) then
write (MP, '(/A/A,I7,A,I6,A,I7,A,I2,A,I7/A,1P,4D12.4/A,4I8/A,3I8)') &
' Entering MA48A/AD with',' M =',M,'     N =',N,'     NE =', NE,'     JOB =',JOB,'     LA =',LA,' CNTL (1:4) =', &
(CNTL(I),I=1,4),' ICNTL(1:4) = ', (ICNTL(I),I=1,4), ' ICNTL(6:8) = ', (ICNTL(I),I=6,8)
         if (ICNTL(3).gt.3) then
            write (MP,9000) (A(K),IRN(K),JCN(K),K=1,NE)
 9000       format (' Entries:'/3 (1P,D12.4,2I6))
         else
            write (MP,9000) (A(K),IRN(K),JCN(K),K=1,MIN(9,NE))
         end if
         if (JOB.eq.2) then
            write (MP,'(A)') ' Permutations input (JOB=2)'
            if (ICNTL(3).gt.3) then
               write (MP,9010) (KEEP(I),I=1,M)
9010          format (' Positions of original rows in the permuted ma', 'trix'/ (10I6))
               write (MP,9020) (KEEP(M+I),I=1,N)
9020          format (' Positions of columns of permuted matrix ','in', ' or','iginal matrix '/ (10I6))
            else
               write (MP,9010) (KEEP(I),I=1,MIN(10,M))
               write (MP,9020) (KEEP(M+I),I=1,MIN(10,N))
            end if
         end if
         if (ICNTL(8).ne.0) then
write (MP,'(A,I6)') ' Value of IW entries on call with ICNTL(8) =',ICNTL(8)
            if (ICNTL(3).gt.3) then
               write (MP,9030) (IW(I),I=1,N)
 9030          format (10I6)
            else
               write (MP,9030) (IW(I),I=1,MIN(10,N))
            end if
         end if
      end if
      do 53 I = 1,20
         INFO(I) = 0
   53 continue
      INFO(3) = NE*2
      INFO(4) = NE
      INFO(10) = MIN(M,N)
      do 56 I = 1,5
        RINFO(I) = ZERO
   56 continue
      do 60 I = 1,4
         CNTL5(I) = CNTL(I)
         ICNTL5(I) = ICNTL(I)
   60 continue
      ICNTL5(3) = 0
      ICNTL5(5) = ICNTL(5)
      ICNTL5(6) = 0
      ICNTL5(7) = 0
      if (JOB.eq.2) then
         ICNTL5(7) = 2
         ICNTL5(4) = 1
      end if
      if (JOB.eq.3) ICNTL5(7) = 1

      TOL = MAX(ZERO,CNTL(4))
      MINBLK = MAX(1,ICNTL(6))
      LDUP = .false.

      IPTRD = M + 3*N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1

      HEADC = N + 1
      LASTC = HEADC + N
      do 70 J = 1,N
         IW(HEADC+J) = 0
   70 continue
      do 80 K = 1,NE
         I = IRN(K)
         J = JCN(K)
         if (I.lt.1 .or. I.gt.M .or. J.lt.1 .or. J.gt.N) then
            INFO(12) = INFO(12) + 1
if (MP.gt.0 .and. INFO(12).le.10 .and. ICNTL(3).GE.2) WRITE (MP,'(A,I7,A,2I6)') &
' Message from MA48A/AD .. indices for entry ',K,' are', I,J
            JCN(K) = 0
         else
            JCN(K) = IW(HEADC+J)
            IW(HEADC+J) = K
         end if
   80 continue

      if (MINBLK.ge.N .or. M.ne.N .or. JOB.gt.1) GO TO 190

      IW21 = 2*N
      IW13 = IW21
      IPTRA = LASTC + N
      LENC = IPTRA + N
      IB = LENC
      IP = LENC + N
      IPTRP = IP + N
      LENP = IPTRP + N

      do 90 I = 1,N
         IW(LASTC+I) = 0
   90 continue

      LDUP = .true.
      K = 1
      do 120 J = 1,N
         EYE = IW(HEADC+J)
         IW(IPTRA+J) = K
         do 100 IDUMMY = 1,NE
            if (EYE.eq.0) GO TO 110
            I = IRN(EYE)
            if (IW(LASTC+I).ne.J) then
               IW(LASTC+I) = J
               IRN(NE+K) = I
               K = K + 1
            else
               INFO(11) = INFO(11) + 1
if (MP.gt.0 .and. INFO(11).le.10 .and. ICNTL(3).GE.2) WRITE (MP,'(A,I7,A,2I6)') &
' Message from MA48A/AD .. duplicate in position ',K, ' with indices',I,J
            end if
            EYE = JCN(EYE)
  100    continue
  110    IW(LENC+J) = K - IW(IPTRA+J)
  120 continue

      NZA = K - 1

call MC21AD(N,IRN(NE+1),NZA,IW(IPTRA+1),IW(LENC+1),IW(IP+1),NDIAG, KEEP(IW21+1))
      INFO(10) = NDIAG
      if (NDIAG.lt.N) then
         if (ICNTL(7).ne.0) then
            INFO(1) = -4
if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP, '(A,A/A,I7,A,I7)') &
' Error return from MA48A/AD because matrix structurally ' ,' singular',' order is ',N,' and structural rank',NDIAG
            GO TO 530
         end if
         GO TO 190
      end if
!DIR$ IVDEP
      do 130 J = 1,N
         JAY = IW(IP+J)
         IW(IPTRP+J) = IW(IPTRA+JAY)
         IW(LENP+J) = IW(LENC+JAY)
  130 continue

call MC13DD(N,IRN(NE+1),NZA,IW(IPTRP+1),IW(LENP+1),KEEP(M+1), IW(IB+1),NB,KEEP(IW13+1))

      do 140 JB = 2,NB
         IW(IB+JB-1) = IW(IB+JB) - IW(IB+JB-1)
  140 continue
      IW(IB+NB) = N + 1 - IW(IB+NB)
      if (IW(IB+1).eq.1) IW(IB+1) = -1
      KB = 1
      do 150 JB = 2,NB
         L = IW(IB+JB)
         if (L.eq.1 .and. IW(IB+KB).le.0) then
            IW(IB+KB) = IW(IB+KB) - 1
         else
            KB = KB + 1
            if (L.eq.1) then
               IW(IB+KB) = -1
            else
               IW(IB+KB) = L
            end if
         end if
  150 continue
      NB = KB
      KB = 1
      do 160 JB = 2,NB
         if (ABS(IW(IB+KB)).lt.MINBLK) then
            IW(IB+KB) = ABS(IW(IB+KB)) + ABS(IW(IB+JB))
         else
            KB = KB + 1
            IW(IB+KB) = IW(IB+JB)
         end if
  160 continue
      NB = KB
      do 170 JB = 1,NB
         KEEP(NBLOCK+3*JB) = ABS(IW(IB+JB))
         KEEP(MBLOCK+3*JB) = IW(IB+JB)
  170 continue
      do 180 J = 1,N
         KEEP(KEEP(M+J)) = J
         KEEP(M+J) = IW(IP+KEEP(M+J))
  180 continue
      GO TO 220
  190 NB = 1
      if (JOB.eq.1 .or. JOB.eq.3) then
         do 200 I = 1,M
            KEEP(I) = I
  200    continue
         do 210 I = 1,N
            KEEP(M+I) = I
  210    continue
      end if
      KEEP(NBLOCK+3) = N
      KEEP(MBLOCK+3) = 0

  220 if (ICNTL(8).ne.0) then
         LBLOCK = KBLOCK + 3*NB
         if (JOB.eq.2) then
            do 230 I = 1,N
               if (IW(I).eq.0) GO TO 240
  230       continue
  240       KEEP(LBLOCK+1) = N - I + 1
         else
            J = 1
            do 270 JB = 1,NB
               KEEP(LBLOCK+JB) = 0
               J2 = J + KEEP(NBLOCK+3*JB) - 1
               J1 = J2
               if (KEEP(MBLOCK+3*JB).lt.0) GO TO 260
  250          if (J.eq.J2) GO TO 260
               if (IW(KEEP(M+J)).eq.0) then
                  KEEP(LBLOCK+JB) = KEEP(LBLOCK+JB) + 1
                  ISW = KEEP(M+J2)
                  KEEP(M+J2) = KEEP(M+J)
                  KEEP(M+J) = ISW
                  J2 = J2 - 1
               else
                  J = J + 1
               end if
               GO TO 250
  260          J = J1 + 1
  270       continue
         end if
      end if

      LASTR = LASTC + M
      do 280 I = 1,M
         IW(LASTC+I) = 0
  280 continue
      KEEP(KBLOCK+3) = NB
      K = 1
      KK = NE
      J2 = 0
      do 310 JB = 1,NB
         J1 = J2 + 1
         J2 = J1 + KEEP(NBLOCK+3*JB) - 1
         do 300 JAY = J1,J2
            J = KEEP(M+JAY)
            EYE = IW(HEADC+J)
            KEEP(IPTRD+JAY) = K
            KEEP(IPTRO+JAY) = KK
            if (KEEP(MBLOCK+3*JB).lt.0) then
               IW(LASTC+JAY) = JAY
               A(NE+K) = ZERO
               IRN(NE+K) = JAY
               IW(LASTR+JAY) = K
               K = K + 1
            end if
            do 290 IDUMMY = 1,NE
               if (EYE.eq.0) GO TO 300
               NXTEYE = JCN(EYE)
               I = KEEP(IRN(EYE))
               if (IW(LASTC+I).ne.JAY) then
                  IW(LASTC+I) = JAY
                  if ((I.ge.J1.and.I.le.J2) .or. (M.ne.N)) then
                     A(NE+K) = A(EYE)
                     IRN(NE+K) = I
                     IW(LASTR+I) = K
                     JCN(EYE) = K
                     K = K + 1
                  else
                     A(NE+KK) = A(EYE)
                     IRN(NE+KK) = I
                     IW(LASTR+I) = KK
                     JCN(EYE) = KK
                     KK = KK - 1
                  end if
               else
                  if (.not.LDUP) then
                     INFO(11) = INFO(11) + 1
if (MP.gt.0 .and. INFO(11).le.10 .and. ICNTL(3).GE.2) WRITE (MP,'(A,I7,A,2I6)') &
' Message from MA48A/AD .. duplicate in position ' ,EYE,' with indices',IRN(EYE),J
                  end if
                  KL = IW(LASTR+I)
                  JCN(EYE) = KL
                  A(NE+KL) = A(NE+KL) + A(EYE)
               end if
               EYE = NXTEYE
  290       continue
  300    continue
  310 continue
      KEEP(IPTRD+N+1) = K
      KEEP(IPTRO+N+1) = KK
      NZD = K - 1
      do 320 I = 1,K-1
         IRN(I) = IRN(NE+I)
  320 continue
      do 325 I = KK+1,NE
         IRN(I) = IRN(NE+I)
  325 continue
      do 326 I = K,KK
         IRN(I) = 0
  326 continue
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
      do 390 JB = NB,1,-1
         NC = KEEP(NBLOCK+3*JB)
         J2 = J1 - 1
         J1 = J2 + 1 - NC
         if (KEEP(MBLOCK+3*JB).lt.0) then
            do 330 J = J1,J2
if (ABS(A(NE+KEEP(IPTRD+J))).gt.TOL) INFO(5) = INFO(5) + 1
               IW(IP+J) = J
               IW(IQ+J) = J
  330       continue
         else
            NZB = KEEP(IPTRD+J2+1) - KEEP(IPTRD+J1)
            do 340 K = KEEP(IPTRD+J1),KEEP(IPTRD+J2+1) - 1
               IRN(NE+K) = IRN(K) - J1 + 1
  340       continue
            K = KEEP(IPTRD+J1) - 1
            do 350 J = J1,J2
               IW(IQ+J) = KEEP(IPTRD+J) - K
  350       continue
            NR = NC
            if (NB.eq.1) NR = M
            if (JOB.eq.2) then
               do 360 J = J1,J1 + NR - 1
                  IW(IP+J) = J - J1 + 1
  360          continue
               do 370 J = J1,J2
                  IW(PTRD+J) = J - J1 + 1
  370          continue
            end if
            INFO(7) = MAX(INFO(7),NR)
            INFO(8) = INFO(8) + NC
            INFO(9) = INFO(9) + NZB
            if (ICNTL(8).ne.0) ICNTL5(6) = KEEP(LBLOCK+JB)
call MA50AD(NR,NC,NZB,LA-NE-K,A(NE+K+1),IRN(NE+K+1), JCN(NE+1),IW(IQ+J1),CNTL5,ICNTL5,IW(IP+J1), &
NP,IW(JFIRST+1),IW(LENR+1), IW(LASTR+1),IW(NEXTR+1),IW(IW50+1),IW(PTRD+J1), &
KEEP(LENC+KB+1),KEEP(LASTC+KB+1),IW(NEXTC+1), INFO5,RINFO5)
            KEEP(MBLOCK+3*JB) = NP
            do 380 J = J1,J1+NR-1
               IW(IP+J) = IW(IP+J) + J1 - 1
  380       continue
            do 385 J = J1,J2
               IW(IQ+J) = IW(IQ+J) + J1 - 1
  385       continue
            if (INFO5(1).eq.1) then
               if (INFO(1).eq.0 .or. INFO(1).eq.4) INFO(1) = INFO(1) + 2
            end if
            if (INFO5(1).eq.2) then
               if (INFO(1).eq.0 .or. INFO(1).eq.2) INFO(1) = INFO(1) + 4
            end if
            if (INFO5(1).eq.3 .and. INFO(1).ge.0) INFO(1) = 6
            if (INFO5(1).eq.-3) INFO(1) = -3
            INFO(2) = INFO(2) + INFO5(2)
            INFO(3) = MAX(INFO(3),NE+K+INFO5(3))

            INFO(5) = INFO(5) + INFO5(5)
            INFO(6) = INFO(6) + INFO5(6)
            RINFO(1) = RINFO(1) + RINFO5(1)
            KB = KB + 1
            KEEP(LENC+KB) = INFO5(4) - 2*NR
            KEEP(LASTC+KB) = INFO5(4)
         end if
  390 continue

      INFO(4) = NE*2
      K = NE
      do 400 JB = KB,1,-1
         INFO(4) = MAX(INFO(4),K+KEEP(LASTC+JB))
         K = K + KEEP(LENC+JB)
  400 continue

      if (INFO(1).eq.-3) then
         if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,9120) LA,INFO(3)
         GO TO 530
      end if

      do 410 K = 1,NE
         IRN(NE+K) = IRN(K)
  410 continue
      do 420 J = 1,N
         IW(PTRD+J) = KEEP(IPTRD+J)
         IW(PTRO+J) = KEEP(IPTRO+J+1) + 1
  420 continue
      KD = 1
      KO = NZD + 1
      do 450 J = 1,N
         KEEP(IPTRD+J) = KD
         JAY = IW(IQ+J)
         KL = NZD
         if (JAY.ne.N) KL = IW(PTRD+JAY+1) - 1
         do 430 KK = IW(PTRD+JAY),KL
            IRN(KD) = IW(IP+IRN(NE+KK))
            IRN(NE+KK) = KD
            KD = KD + 1
  430    continue
         KEEP(IPTRO+J) = KO
         KL = NE
         if (JAY.ne.1) KL = IW(PTRO+JAY-1) - 1
         do 440 KK = IW(PTRO+JAY),KL
            IRN(KO) = IW(IP+IRN(NE+KK))
            IRN(NE+KK) = KO
            KO = KO + 1
  440    continue
  450 continue
      KEEP(IPTRO+N+1) = KO
      do 460 I = 1,M
         KEEP(I) = IW(IP+KEEP(I))
  460 continue
      do 465 I = 1,N
         IW(IQ+I) = KEEP(M+IW(IQ+I))
  465 continue
      do 470 I = 1,N
         KEEP(M+I) = IW(IQ+I)
  470 continue
      IW(1) = IRN(NE)
      IRN(NE) = NE
      do 480 K = 1,NE
         JCN(K) = IRN(NE+JCN(K))
  480 continue
      IRN(NE) = IW(1)

      if (INFO(11).gt.0 .or. INFO(12).gt.0) then
         INFO(1) = INFO(1) + 1
         if (MP.gt.0 .and. ICNTL(3).ge.2) then
            if (INFO(11).gt.0) write (MP,9150) INFO(11)
            if (INFO(12).gt.0) write (MP,9160) INFO(12)
         end if
         if (INFO(11).gt.0) JCN(1) = -JCN(1)
      end if

      if (INFO(10).lt.INFO(5)) then
         INFO(5) = INFO(10)
if (INFO(1).ne.2 .and. INFO(1).ne.3 .and. INFO(1).LT.6) INFO(1) = INFO(1) + 2
      end if

      if (MP.gt.0 .and. ICNTL(3).gt.2) then
write (MP,'(/A/A/12I6/A,F12.1)') ' Leaving MA48A/AD with', ' INFO(1:12)  =',(INFO(I),I=1,12),' RINFO(1) =',RINFO(1)
         write (MP,'(A)') ' Permuted matrix by blocks (IRN)'
         KB = NB
         if (ICNTL(3).eq.3) then
write (MP,'(A)') ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         end if
         write (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         do 500 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            if (J1.le.J2) write (MP,'(A,I6)') ' Block',JB
            J3 = J2
            if (ICNTL(3).eq.3) J3 = J1
            do 490 J = J1,J3
write (MP,'(A,I5,(T13,10I6))') ' Column',J, (IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
  490       continue
            J1 = J2 + 1
  500    continue
         if (KEEP(IPTRO+N+1).gt.KEEP(IPTRD+N+1)) then
            write (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            do 520 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               if (ICNTL(3).eq.3) J3 = J1
               do 510 J = J1,J3
if (KEEP(IPTRO+J+1).gt.KEEP(IPTRO+J)) write (MP, '(A,I5,(T13,10I6))') ' Column',J, &
(IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
  510          continue
               J1 = J2 + 1
  520       continue
         end if
         if (ICNTL(3).gt.3) then
            write (MP,9040) (JCN(K),K=1,NE)
 9040       format (' JCN (MAP) ='/ (6X,10I6))
            write (MP,'(A)') ' Contents of KEEP array'
            write (MP,9010) (KEEP(I),I=1,M)
            write (MP,9020) (KEEP(M+I),I=1,N)
            write (MP,'(A)') ' Pointer information from KEEP'
            write (MP,9050) (KEEP(IPTRD+J),J=1,N+1)
 9050       format (' IPTRD ='/ (8X,10I6))
            write (MP,9060) (KEEP(IPTRO+J),J=1,N+1)
 9060       format (' IPTRO ='/ (8X,10I6))
            write (MP,'(A)') ' Block structure information from KEEP'
            write (MP,9070) (KEEP(NBLOCK+3*JB),JB=1,NB)
            write (MP,9080) (KEEP(MBLOCK+3*JB),JB=1,NB)
 9070       format (' NBLOCK (order blocks) ='/ (8X,10I6))
9080       format (' MBLOCK (triangular flag and number packed rows) =' / (8X,10I6))
 9090       format (' LBLOCK (number of changed columns) ='/ (8X,10I6))
            if (ICNTL(8).ne.0) write (MP,9090) (KEEP(LBLOCK+JB),JB=1,NB)
         else
            write (MP,9040) (JCN(K),K=1,MIN(10,NE))
            write (MP,'(A)') ' Contents of KEEP array'
            write (MP,9010) (KEEP(I),I=1,MIN(10,M))
            write (MP,9020) (KEEP(M+I),I=1,MIN(10,N))
            write (MP,'(A)') ' Pointer information from KEEP'
            write (MP,9050) (KEEP(IPTRD+J),J=1,MIN(10,N+1))
            write (MP,9060) (KEEP(IPTRO+J),J=1,MIN(10,N+1))
            write (MP,'(A)') ' Block structure information from KEEP'
            write (MP,9070) (KEEP(NBLOCK+3*JB),JB=1,MIN(10,NB))
            write (MP,9080) (KEEP(MBLOCK+3*JB),JB=1,MIN(10,NB))
if (ICNTL(8).ne.0) write (MP,9090) (KEEP(LBLOCK+JB),JB=1, MIN(10,NB))
         end if
      end if
  530 return

9100 format (' Error return from MA48A/AD because M =',I10,' and N =', I10)
 9110 format (' Error return from MA48A/AD because NE =',I10)
9120 format (' Error return from MA48A/AD because LA is',I10/' and ', 'must be at least',I10)
 9130 format (' Error return from MA48A/AD because ','JOB = ',I10)
9140 format (' Error return from MA48A/AD because ','faulty permutati', 'ons input when JOB = 2')
 9150 format (' Message from MA48A/AD ..',I8,' duplicates found')
9160 format (' Message from MA48A/AD ..',I8,' out-of-range indices fo', 'und')
      END


subroutine MA48BD(M,N,NE,JOB,LA,A,IRN,JCN,KEEP,CNTL,ICNTL,W,IW, INFO,RINFO)
      integer M,N,NE,JOB,LA
      double precision A(LA)
      integer IRN(LA),JCN(NE),KEEP(*)
      double precision CNTL(10)
      integer :: ICNTL(20)
      double precision W(M)
      integer IW(2*M+2*N),INFO(20)
      double precision RINFO(10)


      double precision ZERO
      parameter (ZERO=0.D0)
      double precision CNTL5(10)
integer I,ICNTL5(20),INFO5(15),IPTRD,IPTRL,IPTRO,IPTRU,IQB(1), ITRY,J,JB,JOB5,J1,J2,J3,K,KB,KBLOCK,KK,LBLOCK,LP, &
MBLOCK,MP,NB,NBLOCK,NEWNE,NC,NP,NR,NRF,NZB
      double precision RINFO5(10),TOL
      logical :: TRISNG

      external MA50BD
      intrinsic MAX

      LP = ICNTL(1)
      MP = ICNTL(2)
      if (M.le.0 .or. N.le.0) then
         INFO(1) = -1
         if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,9160) M,N
         GO TO 240
      end if
      if (NE.le.0) then
         INFO(1) = -2
         if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,9170) NE
         GO TO 240
      end if
      if (LA.lt.2*NE) then
         INFO(1) = -3
         INFO(4) = 2*NE
         if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,9180) LA,INFO(4)
         GO TO 240
      end if
      if (JOB.lt.1 .or. JOB.gt.3) then
         INFO(1) = -6
         if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,9190) JOB
         GO TO 240
      end if
      INFO(1) = 0
      do 10 I = 1,4
         CNTL5(I) = CNTL(I)
         ICNTL5(I) = ICNTL(I)
   10 continue
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

      if (MP.gt.0 .and. ICNTL(3).gt.2) then
write (MP,'(/A/3(A,I8),A,I2,A,I8/A,1P,3D12.4/A,3I8/A,I8/A,I8)') &
' Entering MA48B/BD with',' M =',M,'     N =',N,'     NE =' ,NE,'     JOB =',JOB,'     LA =',LA,' CNTL (2:4) =', &
(CNTL(I),I=2,4),' ICNTL(1:3) =', (ICNTL(I),I=1,3), ' ICNTL(5)   =',ICNTL(5),' ICNTL(8)   =',ICNTL(8)
         if (ICNTL(3).gt.3) then
            write (MP,9000) (A(K),K=1,NE)
         else
            write (MP,9000) (A(K),K=1,MIN(10,NE))
         end if
 9000    format (' A ='/ (4X,1P,5D12.4))
         write (MP,'(A)') ' Indices for permuted matrix by blocks'
         KB = NB
         if (ICNTL(3).eq.3) then
write (MP,'(A)') ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         end if
         write (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         do 30 JB = 1,KB
            write (MP,'(A,I6)') ' Block',JB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            J3 = J2
            if (ICNTL(3).eq.3) J3 = J1
            do 20 J = J1,J3
write (MP,'(A,I5,(T13,10I6))') ' Column',J, (IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
   20       continue
            J1 = J2 + 1
   30    continue
         if (KEEP(IPTRO+N+1).gt.KEEP(IPTRD+N+1)) then
            write (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            do 50 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               if (ICNTL(3).eq.3) J3 = J1
               do 40 J = J1,J3
if (KEEP(IPTRO+J+1).gt.KEEP(IPTRO+J)) write (MP, '(A,I5,(T13,10I6))') ' Column',J, &
(IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
   40          continue
               J1 = J2 + 1
   50       continue
         end if
         if (ICNTL(3).gt.3) then
            write (MP,9010) (JCN(K),K=1,NE)
 9010       format (' JCN (MAP) ='/ (6X,10I6))
            write (MP,'(A)') ' Contents of KEEP array'
            write (MP,9020) (KEEP(I),I=1,M)
            write (MP,9030) (KEEP(M+I),I=1,N)
9020       format (' Positions of original rows in the permuted matrix' / (10I6))
9030       format (' Positions of columns of permuted matrix ','in ', 'original matrix '/ (10I6))
            write (MP,'(A)') ' Pointer information from KEEP'
            write (MP,9040) (KEEP(IPTRD+J),J=1,N+1)
 9040       format (' IPTRD ='/ (8X,10I6))
            write (MP,9050) (KEEP(IPTRO+J),J=1,N+1)
 9050       format (' IPTRO ='/ (8X,10I6))
            if (JOB.gt.1) then
               write (MP,9060) (KEEP(IPTRL+J),J=1,N)
 9060          format (' IPTRL ='/ (8X,10I6))
               write (MP,9070) (KEEP(IPTRU+J),J=1,N)
 9070          format (' IPTRU ='/ (8X,10I6))
            end if
            write (MP,'(A)') ' Block structure information from KEEP'
            write (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,NB)
            write (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,NB)
 9080       format (' NBLOCK (order blocks) ='/ (8X,10I6))
9090       format (' MBLOCK (triangular flag and number packed rows) =' / (8X,10I6))
9100       format (' KBLOCK (position of beginning of block) ='/ (8X,10I6))
 9110       format (' LBLOCK (number of changed columns) ='/ (8X,10I6))
            if (JOB.gt.1) then
               write (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,NB)
if (ICNTL(8).ne.0) write (MP,9110) (KEEP(LBLOCK+JB),JB=1, NB)
            end if
         else
            write (MP,9010) (JCN(K),K=1,MIN(10,NE))
            write (MP,'(A)') ' Contents of KEEP array'
            write (MP,9020) (KEEP(I),I=1,MIN(10,M))
            write (MP,9030) (KEEP(M+I),I=1,MIN(10,N))
            write (MP,'(A)') ' Pointer information from KEEP'
            write (MP,9040) (KEEP(IPTRD+J),J=1,MIN(10,N+1))
            write (MP,9050) (KEEP(IPTRO+J),J=1,MIN(10,N+1))
            if (JOB.gt.1) then
               write (MP,9060) (KEEP(IPTRL+J),J=1,MIN(10,N))
               write (MP,9070) (KEEP(IPTRU+J),J=1,MIN(10,N))
            end if
            write (MP,'(A)') ' Block structure information from KEEP'
            write (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,MIN(10,NB))
            write (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,MIN(10,NB))
            if (JOB.gt.1) then
               write (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,MIN(10,NB))
if (ICNTL(8).ne.0) write (MP,9110) (KEEP(LBLOCK+JB),JB=1, MIN(10,NB))
            end if
         end if
      end if
      INFO(4) = NE
      INFO(5) = 0
      INFO(6) = 0
      RINFO(1) = ZERO
      TOL = MAX(ZERO,CNTL(4))
      if (JCN(1).gt.0) then
         do 60 K = 1,NE
            A(NE+K) = A(K)
   60    continue
!DIR$ IVDEP
         do 70 K = 1,NE
            A(JCN(K)) = A(NE+K)
   70    continue
      else
         do 80 K = 1,NE
            A(NE+K) = A(K)
            A(K) = ZERO
   80    continue
         A(-JCN(1)) = A(NE+1)
         do 90 K = 2,NE
            KK = JCN(K)
            A(KK) = A(KK) + A(NE+K)
   90    continue
      end if
      IQB(1) = 0
      KK = 0
      J2 = 0
      JOB5 = JOB
      do 150 JB = 1,NB
         NC = KEEP(NBLOCK+3*JB)
         J1 = J2 + 1
         J2 = J1 + NC - 1
         KEEP(KBLOCK+3*JB) = 0
         if (KEEP(MBLOCK+3*JB).lt.0) then
            TRISNG = .false.
            do 100 J = J1,J2
               if (ABS(A(KEEP(IPTRD+J))).le.TOL) TRISNG = .true.
               KEEP(IPTRL+J) = 0
               KEEP(IPTRU+J) = 0
  100       continue
            if (.not.TRISNG)  then
               INFO(5) = INFO(5) + NC
               GO TO 150
            endif
            if (JOB.eq.2) then
               INFO(1) = -7
if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,'(A)') ' Error return from MA48B/BD with JOB=2 because ', &
' the matrix is incompatible with expectations'
               GO TO 240
            else
               KEEP(MBLOCK+3*JB) = NC
            endif
         end if
       do 145 ITRY = 1,2
         NR = NC
         if (NB.eq.1) NR = M
         NZB = KEEP(IPTRD+J2+1) - KEEP(IPTRD+J1)
         if (ICNTL(8).ne.0) ICNTL5(6) = KEEP(LBLOCK+JB)
         do 110 K = KEEP(IPTRD+J1),KEEP(IPTRD+J2+1) - 1
            IRN(K) = IRN(K) - J1 + 1
  110    continue
         K = KEEP(IPTRD+J1) - 1
         do 115 J = J1,J1+NC-1
            KEEP(IPTRD+J) = KEEP(IPTRD+J) - K
  115    continue
         do 120 J = J1,J1+NR-1
            IW(J) = J - J1 + 1
  120    continue
         NP = KEEP(MBLOCK+3*JB)
call MA50BD(NR,NC,NZB,JOB5,A(K+1),IRN(K+1),KEEP(IPTRD+J1), CNTL5,ICNTL5,IW(J1),IQB,NP,LA-NEWNE-KK, &
A(NEWNE+KK+1),IRN(NEWNE+KK+1),KEEP(IPTRL+J1), KEEP(IPTRU+J1),W,IW(M+1),INFO5,RINFO5)
         do 130 J = J1,J2
            KEEP(IPTRD+J) = KEEP(IPTRD+J) + K
  130    continue
         do 140 K = KEEP(IPTRD+J1),KEEP(IPTRD+J2+1) - 1
            IRN(K) = IRN(K) + J1 - 1
  140    continue
         if (INFO5(1).eq.-6) then
            INFO(1) = -6
if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,'(A)') ' Error return from MA48B/BD with JOB greater than 1', &
' and entries dropped during previous factorization'
            GO TO 240
         end if
         if (INFO5(1).lt.-7) then
            if (ICNTL(11).eq.1 .and. JOB.eq.2) then
               JOB5 = 1
if (LP.gt.0 .and. ICNTL(3).ge.2) write(LP,'(A,2(A,I4))') ' Warning from MA48B/BD. Switched from JOB=2 to JOB=1', &
' in block', JB, ' of ',NB
               INFO(1) = INFO(1) + 1
               GO TO 145
            end if
            INFO(1) = -7
if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,'(A)') ' Error return from MA48B/BD with JOB=2 because ', &
' the matrix is incompatible with expectations'
            GO TO 240
         else
            GO TO 147
         end if
  145  continue
  147    if (INFO5(1).eq.-3) then
            INFO(1) = -3
            if (ICNTL(10).eq.1) then
               KEEP(KBLOCK+3) = NB
if (LP.gt.0 .and. ICNTL(3).ge.1) write(LP,'(A,2(A,I4))') ' Error return from MA48B/BD because LA is too small.', &
' In block', JB, ' of ',NB
               GO TO 240
            end if
         end if
         if (INFO(1).eq.-3) then
            INFO(4) = INFO(4) + INFO5(4)
            KK = 0
         else
            INFO(4) = MAX(INFO(4),KK+NEWNE+INFO5(4))
            NRF = IRN(NEWNE+KK+2)
            KEEP(KBLOCK+3*JB) = KK + 1
KK = KK + KEEP(IPTRL+J2) + MAX((NC-KEEP(MBLOCK+3*JB))* (NRF), (NC-KEEP(MBLOCK+3*JB))+(NRF))
         end if
         if (INFO5(1).eq.1) then
            if (INFO(1).ne.-3) INFO(1) = MIN(INFO(1)+2,3)
         end if
         RINFO(1) = RINFO(1) + RINFO5(1)
         INFO(5) = INFO(5) + INFO5(5)
         INFO(6) = INFO(6) + INFO5(6)
  150 continue

      INFO(4) = MAX(NE*2,INFO(4))
      KEEP(KBLOCK+3) = NB

      if (INFO(1).eq.-3) then
         if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,9180) LA,INFO(4)
         GO TO 240
      end if

      if (MP.gt.0 .and. ICNTL(3).gt.2) then
write (MP,'(/A/A,I6/A,3I6/A,F12.1)') ' Leaving MA48B/BD with', ' INFO(1)   = ',INFO(1),' INFO(4:6) = ', (INFO(I),I=4,6), &
' RINFO(1)     =',RINFO(1)
         write (MP,'(A)') ' Permuted matrix by blocks'
         KB = NB
         if (ICNTL(3).eq.3) then
write (MP,'(A)') ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         end if
         write (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         do 170 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            if (J1.le.J2) write (MP,'(A,I6)') ' Block',JB
            J3 = J2
            if (ICNTL(3).eq.3) J3 = MIN(J1,J2)
            do 160 J = J1,J3
write (MP,'(A,I5,(T13,3(1PD12.4,I5)))') ' Column',J, (A(I),IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
  160       continue
            J1 = J2 + 1
  170    continue
         if (KEEP(IPTRO+N+1).gt.KEEP(IPTRD+N+1)) then
            write (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            do 190 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               if (ICNTL(3).eq.3) J3 = MIN(J1,J2)
               do 180 J = J1,J3
if (KEEP(IPTRO+J+1).gt.KEEP(IPTRO+J)) write (MP, '(A,I5,(T13,3(1P,D12.4,I5)))') ' Column',J, &
(A(I),IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
  180          continue
               J1 = J2 + 1
  190       continue
         end if
         write (MP,'(A)') ' Factorized matrix by blocks'
         J1 = 1
         do 230 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            if (KEEP(MBLOCK+3*JB).lt.0) GO TO 220
            NC = J2 - J1 + 1
            NR = NC
            if (KB.eq.1) NR = M
            write (MP,'(A,I6)') ' Block',JB
            K = NEWNE
            if (JB.gt.1) K = KEEP(KBLOCK+3*JB) + NEWNE - 1
if (KEEP(IPTRL+J1).gt.KEEP(IPTRU+J1)) write (MP, '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J1,' of L', &
(A(K+I),IRN(K+I),I=KEEP(IPTRU+J1)+1,KEEP(IPTRL+J1))
            if (ICNTL(3).eq.3) GO TO 210
            do 200 J = J1 + 1,J2
if (KEEP(IPTRU+J).gt.KEEP(IPTRL+J-1)) write (MP, '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J,' of U', &
(A(K+I),IRN(K+I),I=KEEP(IPTRL+J-1)+1,KEEP(IPTRU+J))
if (KEEP(IPTRU+J).lt.KEEP(IPTRL+J)) write (MP, '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J,' of L', &
(A(K+I),IRN(K+I),I=KEEP(IPTRU+J)+1,KEEP(IPTRL+J))
  200       continue
  210       write (MP,'(A)') ' Full block'
            write (MP,'(A)') ' Row indices'
            NRF = IRN(K+2)
            K = K + KEEP(IPTRL+J2)
            if (ICNTL(3).gt.3) then
               write (MP,9120) (IRN(K+I),I=1,NRF)
               write (MP,'(A)') ' Column pivoting information'
write (MP,9120) (IRN(K+I),I=NRF+1, NRF+NC-KEEP(MBLOCK+3*JB))
               write (MP,'(A)') ' Reals by columns'
write (MP,9130) (A(K+I),I=1, (NRF)* (NC-KEEP(MBLOCK+3*JB)))
 9120          format (10I6)
 9130          format (1P,5D12.4)
            else
               write (MP,9120) (IRN(K+I),I=1,MIN(10,NRF))
               write (MP,'(A)') ' Column pivoting information'
write (MP,9120) (IRN(K+I),I=NRF+1, NRF+MIN(10,NC-KEEP(MBLOCK+3*JB)))
               write (MP,'(A)') ' Reals by columns'
write (MP,9130) (A(K+I),I=1, MIN(10, (NRF)* (NC-KEEP(MBLOCK+3*JB))))
            end if
  220       J1 = J2 + 1
  230    continue
         if (JOB.eq.1 .or. JOB.eq.3) then
            write (MP,'(A)') ' Contents of KEEP array'
            if (ICNTL(3).gt.3) then
               write (MP,9020) (KEEP(I),I=1,M)
               write (MP,9030) (KEEP(M+I),I=1,N)
               write (MP,'(A)') ' Pointer information from KEEP'
               write (MP,9140) (KEEP(IPTRL+J),J=1,N)
 9140          format (' IPTRL ='/ (8X,10I6))
               write (MP,9150) (KEEP(IPTRU+J),J=1,N)
 9150          format (' IPTRU ='/ (8X,10I6))
            else
               write (MP,9020) (KEEP(I),I=1,MIN(10,M))
               write (MP,9030) (KEEP(M+I),I=1,MIN(10,N))
               write (MP,'(A)') ' Pointer information from KEEP'
               write (MP,9140) (KEEP(IPTRL+J),J=1,MIN(10,N))
               write (MP,9150) (KEEP(IPTRU+J),J=1,MIN(10,N))
            end if
            write (MP,'(A)') ' Block structure information from KEEP'
            write (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,KB)
         end if
      end if

  240 return

9160 format (' Error return from MA48B/BD because M =',I10,' and N =', I10)
 9170 format (' Error return from MA48B/BD because NE =',I10)
9180 format (' Error return from MA48B/BD because LA is',I10/' and mu', 'st be at least',I10)
 9190 format (' Error return from MA48B/BD because ','JOB = ',I10)
      END


subroutine MA48CD(M,N,TRANS,JOB,LA,A,IRN,KEEP,CNTL,ICNTL,RHS,X, ERROR,W,IW,INFO)

      integer M,N
      logical :: TRANS
      integer JOB,LA
      double precision A(LA)
      integer IRN(LA),KEEP(*)
      double precision CNTL(10)
      integer :: ICNTL(20)
      double precision RHS(*),X(*),ERROR(3),W(*)
      integer IW(*),INFO(20)

      double precision ZERO,ONE
      parameter (ZERO=0.D0,ONE=1.0D0)
      double precision COND(2),CTAU,DXMAX
integer I,ICNTL5(20),IPTRD,IPTRL,IPTRO,IPTRU,IQB(1),J,JB,JJ,J1,J2, J3,K,KASE,KB,KBLOCK,KK,KEEP71(5)
      logical :: LCOND(2)
      integer LP,MBLOCK,MP,NB,NBLOCK,NC,NE,NEQ,NRF,NVAR
      double precision OLDOMG(2),OMEGA(2),OM1,OM2,TAU
      external MA48DD,MA50CD,MC71AD
      double precision EPS
      intrinsic ABS,MAX

      LP = ICNTL(1)
      MP = ICNTL(2)
      if (N.le.0 .or. M.le.0) then
         INFO(1) = -1
         if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,9140) M,N
         GO TO 380
      end if
      if (JOB.gt.4 .or. JOB.lt.1) then
         INFO(1) = -6
         if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,9150) JOB
         GO TO 380
      end if
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

      do 10 I = 1,7
         ICNTL5(I) = 0
   10 continue
      ICNTL5(5) = ICNTL(5)

      if (MP.gt.0 .and. ICNTL(3).gt.2) then
write (MP, '(/A/3(A,I8),A,I2/A,L2,A,I7/A,1P,E12.4/A,3I6,2(/A,I6))') &
' Entering MA48C/CD with',' M =',M,'     N =',N,'     LA =', LA,'      JOB =',JOB,'   TRANS =',TRANS, &
'      No. of blocks =', NB,'   CNTL(5)    = ',CNTL(5),'   ICNTL(1:3) = ', (ICNTL(I),I=1,3),'   ICNTL(5)   = ',ICNTL(5), &
'   ICNTL(9)   = ',ICNTL(9)
         write (MP,'(A)') ' Permuted matrix by blocks'
         KB = NB
         if (ICNTL(3).eq.3) then
write (MP,'(A)') ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         end if
         write (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         do 30 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            if (J1.le.J2) write (MP,'(A,I6)') ' Block',JB
            J3 = J2
            if (ICNTL(3).eq.3) J3 = MIN(J1,J2)
            do 20 J = J1,J3
write (MP,'(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J, (A(I),IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
   20       continue
            J1 = J2 + 1
   30    continue
         if (KEEP(IPTRO+N+1).gt.KEEP(IPTRD+N+1)) then
            write (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            do 50 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               if (ICNTL(3).eq.3) J3 = MIN(J1,J2)
               do 40 J = J1,J3
if (KEEP(IPTRO+J+1).gt.KEEP(IPTRO+J)) write (MP, '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J, &
(A(I),IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
   40          continue
               J1 = J2 + 1
   50       continue
         end if
         write (MP,'(A)') ' Factorized matrix by blocks'
         J1 = 1
         do 90 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            if (KEEP(MBLOCK+3*JB).lt.0) GO TO 80
            NC = J2 - J1 + 1
            write (MP,'(A,I6)') ' Block',JB
            K = NE
            if (JB.gt.1) K = KEEP(KBLOCK+3*JB) + NE - 1
if (KEEP(IPTRL+J1).gt.KEEP(IPTRU+J1)) write (MP, '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J1,' of L', &
(A(K+I),IRN(K+I),I=KEEP(IPTRU+J1)+1,KEEP(IPTRL+J1))
            if (ICNTL(3).eq.3) GO TO 70
            do 60 J = J1 + 1,J2
if (KEEP(IPTRU+J).gt.KEEP(IPTRL+J-1)) write (MP, '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of U', &
(A(K+I),IRN(K+I),I=KEEP(IPTRL+J-1)+1,KEEP(IPTRU+J))
if (KEEP(IPTRU+J).lt.KEEP(IPTRL+J)) write (MP, '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L', &
(A(K+I),IRN(K+I),I=KEEP(IPTRU+J)+1,KEEP(IPTRL+J))
   60       continue
   70       write (MP,'(A)') ' Full block'
            write (MP,'(A)') ' Row indices'
            NRF = IRN(K+2)
            K = K + KEEP(IPTRL+J2)
            if (ICNTL(3).gt.3) then
               write (MP,9000) (IRN(K+I),I=1,NRF)
               write (MP,'(A)') ' Column pivoting information'
write (MP,9000) (IRN(K+I),I=NRF+1, NRF+NC-KEEP(MBLOCK+3*JB))
               write (MP,'(A)') ' Reals by columns'
write (MP,9010) (A(K+I),I=1, (NRF)* (NC-KEEP(MBLOCK+3*JB)))
 9000          format (10I6)
 9010          format (1P,5D12.4)
            else
               write (MP,9000) (IRN(K+I),I=1,MIN(10,NRF))
               write (MP,'(A)') ' Column pivoting information'
write (MP,9000) (IRN(K+I),I=NRF+1, NRF+MIN(10,NC-KEEP(MBLOCK+3*JB)))
               write (MP,'(A)') ' Reals by columns'
write (MP,9010) (A(K+I),I=1, MIN(10, (NRF)* (NC-KEEP(MBLOCK+3*JB))))
            end if
   80       J1 = J2 + 1
   90    continue
         if (ICNTL(3).gt.3) then
            write (MP,'(A)') ' Contents of KEEP array'
            write (MP,9020) (KEEP(I),I=1,M)
            write (MP,9030) (KEEP(M+I),I=1,N)
9020       format (' Positions of original rows in the permuted matrix' / (10I6))
9030       format (' Positions of columns of permuted matrix ','in or', 'ig','inal matrix '/ (10I6))
            write (MP,'(A)') ' Pointer information from KEEP'
            write (MP,9040) (KEEP(IPTRD+J),J=1,N+1)
 9040       format (' IPTRD ='/ (8X,10I6))
            write (MP,9050) (KEEP(IPTRO+J),J=1,N+1)
 9050       format (' IPTRO ='/ (8X,10I6))
            write (MP,9060) (KEEP(IPTRL+J),J=1,N)
 9060       format (' IPTRL ='/ (8X,10I6))
            write (MP,9070) (KEEP(IPTRU+J),J=1,N)
 9070       format (' IPTRU ='/ (8X,10I6))
            write (MP,'(A)') ' Block structure information from KEEP'
            write (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,NB)
            write (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,NB)
            write (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,NB)
 9080       format (' NBLOCK (order blocks) ='/ (8X,10I6))
9090       format (' MBLOCK (triangular flag and number packed rows) =' / (8X,10I6))
9100       format (' KBLOCK (position of beginning of block) ='/ (8X,10I6))
            if (TRANS) then
               write (MP,9110) (RHS(I),I=1,N)
 9110          format (' RHS =  ',1P,5D12.4/ (8X,5D12.4))
            else
               write (MP,9110) (RHS(I),I=1,M)
            end if
         else
            write (MP,'(A)') ' Contents of KEEP array'
            write (MP,9020) (KEEP(I),I=1,MIN(10,M))
            write (MP,9030) (KEEP(M+I),I=1,MIN(10,N))
            write (MP,'(A)') ' Pointer information from KEEP'
            write (MP,9040) (KEEP(IPTRD+K),K=1,MIN(10,N+1))
            write (MP,9050) (KEEP(IPTRO+K),K=1,MIN(10,N+1))
            write (MP,9060) (KEEP(IPTRL+J),J=1,MIN(10,N))
            write (MP,9070) (KEEP(IPTRU+J),J=1,MIN(10,N))
            write (MP,'(A)') ' Block structure information from KEEP'
            write (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,KB)
            write (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,KB)
            write (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,KB)
            if (TRANS) then
               write (MP,9110) (RHS(I),I=1,MIN(10,N))
            else
               write (MP,9110) (RHS(I),I=1,MIN(10,M))
            end if
         end if
      end if

      if (TRANS) then
         NEQ = N
         NVAR = M
         do 100 I = 1,NEQ
            W(I) = RHS(KEEP(M+I))
  100    continue
      else
         NEQ = M
         NVAR = N
         do 110 I = 1,NEQ
            W(KEEP(I)) = RHS(I)
  110    continue
      end if
      if (JOB.eq.1) then
         if (NB.eq.1 .and. KEEP(MBLOCK+3).ge.0) then
call MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),TRANS,LA-NE, A(NE+1),IRN(NE+1),KEEP(IPTRL+1),KEEP(IPTRU+1),W, X,W(NEQ+1),INFO)
         else
call MA48DD (N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN,KEEP(IPTRD+1), KEEP(IPTRO+1),NB,KEEP(NBLOCK+3),KEEP(IPTRL+1), &
KEEP(IPTRU+1),W,X,TRANS,ICNTL5,W(NEQ+1))
         end if
         GO TO 340
      end if

      do 120 I = 1,NVAR
         X(I) = ZERO
  120 continue
      do 130 I = 1,NEQ
         RHS(I) = W(I)
  130 continue

      OM1 = ZERO
      do 260 K = 1,ICNTL(9)
         if (NB.eq.1 .and. KEEP(MBLOCK+3).ge.0) then
call MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),TRANS,LA-NE, A(NE+1),IRN(NE+1),KEEP(IPTRL+1),KEEP(IPTRU+1), &
W,W(NEQ+1),W(M+N+1),INFO)
         else
call MA48DD (N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN,KEEP(IPTRD+1), KEEP(IPTRO+1),NB,KEEP(NBLOCK+3),KEEP(IPTRL+1), &
KEEP(IPTRU+1),W,W(NEQ+1),TRANS,ICNTL5,W(M+N+1))
         end if
         do 140 I = 1,NVAR
            X(I) = X(I) + W(NEQ+I)
  140    continue
         do 150 I = 1,NEQ
            W(I) = RHS(I)
            W(NEQ+I) = ZERO
            W(2*NEQ+I) = ZERO
  150    continue
         if (TRANS) then
            do 180 J = 1,N
!DIR$ IVDEP
               do 160 JJ = KEEP(IPTRD+J),KEEP(IPTRD+J+1) - 1
                  I = IRN(JJ)
                  W(J) = W(J) - A(JJ)*X(I)
                  W(NEQ+J) = W(NEQ+J) + ABS(A(JJ)*X(I))
                  W(2*NEQ+J) = W(2*NEQ+J) + ABS(A(JJ))
  160          continue
!DIR$ IVDEP
               do 170 JJ = KEEP(IPTRO+J),KEEP(IPTRO+J+1) - 1
                  I = IRN(JJ)
                  W(J) = W(J) - A(JJ)*X(I)
                  W(NEQ+J) = W(NEQ+J) + ABS(A(JJ)*X(I))
                  W(2*NEQ+J) = W(2*NEQ+J) + ABS(A(JJ))
  170          continue
  180       continue
         else
            do 210 J = 1,N
!DIR$ IVDEP
               do 190 JJ = KEEP(IPTRD+J),KEEP(IPTRD+J+1) - 1
                  I = IRN(JJ)
                  W(I) = W(I) - A(JJ)*X(J)
                  W(NEQ+I) = W(NEQ+I) + ABS(A(JJ)*X(J))
                  W(2*NEQ+I) = W(2*NEQ+I) + ABS(A(JJ))
  190          continue
!DIR$ IVDEP
               do 200 JJ = KEEP(IPTRO+J),KEEP(IPTRO+J+1) - 1
                  I = IRN(JJ)
                  W(I) = W(I) - A(JJ)*X(J)
                  W(NEQ+I) = W(NEQ+I) + ABS(A(JJ)*X(J))
                  W(2*NEQ+I) = W(2*NEQ+I) + ABS(A(JJ))
  200          continue
  210       continue
         end if
         DXMAX = ZERO
         do 220 I = 1,NVAR
            DXMAX = MAX(DXMAX,ABS(X(I)))
  220    continue
         OMEGA(1) = ZERO
         OMEGA(2) = ZERO
         do 230 I = 1,NEQ
            TAU = (W(2*NEQ+I)*DXMAX+ABS(RHS(I)))*NVAR*CTAU
            if ((W(NEQ+I)+ABS(RHS(I))).gt.TAU) then
OMEGA(1) = MAX(OMEGA(1),ABS(W(I))/ (W(NEQ+I)+ABS(RHS(I))))
               IW(I) = 1
            else
               if (TAU.gt.ZERO) then
OMEGA(2) = MAX(OMEGA(2),ABS(W(I))/ (W(NEQ+I)+W(2*NEQ+I)*DXMAX))
               end if
               IW(I) = 2
            end if
  230    continue
         if (JOB.eq.2) GO TO 340
         OM2 = OMEGA(1) + OMEGA(2)
         if (OM2.le.EPS) GO TO 270
         if (K.gt.1 .and. OM2.gt.OM1*CNTL(5)) then
            if (OM2.gt.OM1) then
               OMEGA(1) = OLDOMG(1)
               OMEGA(2) = OLDOMG(2)
               do 240 I = 1,NVAR
                  X(I) = W(3*NEQ+I)
  240          continue
            end if
            GO TO 270
         end if
         do 250 I = 1,NVAR
            W(3*NEQ+I) = X(I)
  250    continue
         OLDOMG(1) = OMEGA(1)
         OLDOMG(2) = OMEGA(2)
         OM1 = OM2
  260 continue
      INFO(1) = -8
      if (LP.gt.0 .and. ICNTL(3).ge.1) write (LP,9170) INFO(1),ICNTL(9)
      GO TO 340

  270 if (JOB.le.3) GO TO 340
      if (M.ne.N) GO TO 340
      LCOND(1) = .false.
      LCOND(2) = .false.
      do 280 I = 1,NEQ
         if (IW(I).eq.1) then
            W(I) = W(NEQ+I) + ABS(RHS(I))
            W(NEQ+I) = ZERO
            LCOND(1) = .true.
         else

            W(NEQ+I) = W(NEQ+I) + W(2*NEQ+I)*DXMAX
            W(I) = ZERO
            LCOND(2) = .true.
         end if
  280 continue
      KASE = 0
      do 330 K = 1,2
         if (LCOND(K)) then
            do 310 KK = 1,40
               call MC71AD(N,KASE,W(3*NEQ+1),COND(K),RHS,IW,KEEP71)
               if (KASE.eq.0) GO TO 320
               if (KASE.eq.1) then
                  if (NB.eq.1 .and. KEEP(MBLOCK+3).ge.0) then
call MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3), .NOT.TRANS,LA-NE,A(NE+1),IRN(NE+1), KEEP(IPTRL+1),KEEP(IPTRU+1),W(3*NEQ+1), &
W(2*NEQ+1),RHS,INFO)
                  else
call MA48DD(N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN, KEEP(IPTRD+1),KEEP(IPTRO+1),NB, KEEP(NBLOCK+3),KEEP(IPTRL+1), &
KEEP(IPTRU+1),W(3*NEQ+1),W(2*NEQ+1), .NOT.TRANS,ICNTL5,RHS)
                  end if

                  do 290 I = 1,M
                     W(3*NEQ+I) = W((K-1)*NEQ+I)*W(2*NEQ+I)
  290             continue
               end if
               if (KASE.eq.2) then
                  do 300 I = 1,N
                     W(2*NEQ+I) = W((K-1)*NEQ+I)*W(3*NEQ+I)
  300             continue
                  if (NB.eq.1 .and. KEEP(MBLOCK+3).ge.0) then
call MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),TRANS, LA-NE,A(NE+1),IRN(NE+1),KEEP(IPTRL+1), &
KEEP(IPTRU+1),W(2*NEQ+1),W(3*NEQ+1), RHS,INFO)
                  else
call MA48DD(N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN, KEEP(IPTRD+1),KEEP(IPTRO+1),NB, KEEP(NBLOCK+3),KEEP(IPTRL+1), &
KEEP(IPTRU+1),W(2*NEQ+1),W(3*NEQ+1), TRANS,ICNTL5,RHS)
                  end if
               end if
  310       continue
            INFO(1) = -9
            if (LP.ne.0 .and. ICNTL(3).ge.1) write (LP,9160)
            GO TO 340
  320       if (DXMAX.gt.ZERO) COND(K) = COND(K)/DXMAX
            ERROR(3) = ERROR(3) + OMEGA(K)*COND(K)
         end if
  330 continue

  340 do 350 I = 1,NVAR
         W(I) = X(I)
  350 continue
      if (.not.TRANS) then
         do 360 I = 1,NVAR
            X(KEEP(M+I)) = W(I)
  360    continue
      else
         do 370 I = 1,NVAR
            X(I) = W(KEEP(I))
  370    continue
      end if
      if (JOB.ge.2) then
         ERROR(1) = OMEGA(1)
         ERROR(2) = OMEGA(2)
      end if

      if (MP.gt.0 .and. ICNTL(3).gt.2) then
         write (MP,'(/A,I6)') ' Leaving MA48C/CD with INFO(1) =',INFO(1)
         if (JOB.gt.1) then
            K = 2
            if (JOB.eq.4 .and. INFO(1).ne.-9) K = 3
            write (MP,9120) (ERROR(I),I=1,K)
 9120       format (' ERROR =',1P,3D12.4)
         end if
         if (ICNTL(3).gt.3) then
            write (MP,9130) (X(I),I=1,NVAR)
 9130       format (' X =    ',1P,5D12.4:/ (8X,5D12.4))
         else
            write (MP,9130) (X(I),I=1,MIN(10,NVAR))
         end if
      end if
  380 return

9140 format (' Error return from MA48C/CD because M =',I10,' and N =', I10)
 9150 format (' Error return from MA48C/CD because ','JOB = ',I10)
9160 format (' Error return from MA48C/CD because of ','error in MC71', 'A/AD'/' ERROR(3) not calculated')
9170 format (' Error return from MA48C/CD because of ','nonconvergenc', &
'e of iterative refinement'/' Error INFO(1) = ',I2,'  wit', 'h ICNTL','(9) = ',I10)
      END


subroutine MA48DD(N,NE,LA,A,AA,IRN,IRNA,IPTRD,IPTRO,NB,IBLOCK, IPTRL,IPTRU,RHS,X,TRANS,ICNTL5,W)
      integer N,NE,LA
      double precision A(LA),AA(NE)
      integer IRN(LA),IRNA(NE),IPTRD(N+1),IPTRO(N+1),NB
      integer IBLOCK(3,NB),IPTRL(N),IPTRU(N)
      double precision RHS(N),X(N)
      logical :: TRANS
      integer :: ICNTL5(20)
      double precision W(N)
      integer I,IFLAG(15),IQB(1),J,JB,JJ,J1,K1,K2,NC,NUMB
      external MA50CD
      IQB(1) = 0
      NUMB = IBLOCK(3,1)
      if (.not.TRANS) then
         K1 = N + 1
         do 50 JB = NUMB,1,-1
            NC = IBLOCK(1,JB)
            K2 = K1 - 1
            K1 = K1 - NC
            if (IBLOCK(2,JB).lt.0) then
               do 20 J = K2,K1,-1
                  X(J) = RHS(J)/AA(IPTRD(J))
!DIR$ IVDEP
                  do 10 JJ = IPTRD(J) + 1,IPTRD(J+1) - 1
                     I = IRNA(JJ)
                     RHS(I) = RHS(I) - AA(JJ)*X(J)
   10             continue
   20          continue
            else
               J1 = 1
               if (JB.gt.1) J1 = IBLOCK(3,JB)
call MA50CD(NC,NC,ICNTL5,IQB,IBLOCK(2,JB),TRANS,LA+1-J1, A(J1),IRN(J1),IPTRL(K1),IPTRU(K1),RHS(K1), X(K1),W,IFLAG)
            end if
            if (JB.eq.1) GO TO 50
            do 40 J = K1,K2
!DIR$ IVDEP
               do 30 JJ = IPTRO(J),IPTRO(J+1) - 1
                  I = IRNA(JJ)
                  RHS(I) = RHS(I) - AA(JJ)*X(J)
   30          continue
   40       continue
   50    continue
      else
         K2 = 0
         do 100 JB = 1,NUMB
            NC = IBLOCK(1,JB)
            K1 = K2 + 1
            K2 = K2 + NC
            if (JB.gt.1) then
               do 70 J = K1,K2
                  do 60 JJ = IPTRO(J),IPTRO(J+1) - 1
                     I = IRNA(JJ)
                     RHS(J) = RHS(J) - AA(JJ)*X(I)
   60             continue
   70          continue
            end if
            if (IBLOCK(2,JB).lt.0) then
               do 90 J = K1,K2
                  do 80 JJ = IPTRD(J) + 1,IPTRD(J+1) - 1
                     I = IRNA(JJ)
                     RHS(J) = RHS(J) - AA(JJ)*X(I)
   80             continue
                  X(J) = RHS(J)/AA(IPTRD(J))
   90          continue
            else
               J1 = 1
               if (JB.gt.1) J1 = IBLOCK(3,JB)
call MA50CD(NC,NC,ICNTL5,IQB,IBLOCK(2,JB),TRANS,LA+1-J1, A(J1),IRN(J1),IPTRL(K1),IPTRU(K1),RHS(K1), X(K1),W,IFLAG)
            end if
  100    continue
      end if

      return
      END

      subroutine MA48ID(CNTL,ICNTL)

      double precision CNTL(10)
      integer :: ICNTL(20)


      integer :: I

      do 10 I = 3, 10
         CNTL(I) = 0.0D0
   10 continue
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
      do 20 I = 10, 20
         ICNTL(I) = 0
   20 continue

      END

