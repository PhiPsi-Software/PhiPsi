!-----------------------------------------------------------
! Brief: Solve a sparse linear system using the ITPACK JCG algorithm.
!
! Parameters:
!   Input:  NN     - matrix dimension (N)
!           IA,JA  - CSR row pointer and column index arrays
!           A      - CSR nonzero values
!           RHS    - right-hand side vector
!           IWKSP  - integer workspace (length 3*N)
!           NW     - length of available WKSP workspace
!           IPARM  - integer parameter array (length 12)
!           RPARM  - real parameter array (length 12)
!   In/Out: U     - initial guess on input, solution on output
!   In/Out: WKSP  - real workspace vector
!   Output: IERR  - error flag (0 = success)
!
! Notes:   ITPACK 2C Jacobi Conjugate Gradient driver. Uses IPARM/RPARM
!          to control maximum iterations, stopping tolerance, and storage.
!-----------------------------------------------------------

      subroutine JCG (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,IERR)
      !http://www.netlib.org/itpack/
      !Modified by Fang Shi. 2024-09-08. IA(1) to IA(N+1). BUGFIX2024090801.
      
!       
!     ITPACK 2C MAIN subroutine  JCG  (JACOBI CONJUGATE GRADIENT)   
!     EACH OF THE MAIN SUBROUTINES:   
!           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
!     CAN BE USED INDEPENDENTLY OF THE OTHERS   
!       
! ... function:   
!       
!          THIS subroutine, JCG, DRIVES THE JACOBI CONJUGATE
!          GRADIENT ALGORITHM.
!       
! ... parameter LIST:       
!       
!          N      INPUT integer.  dimension OF THE MATRIX. (= NN)   
!          IA,JA  INPUT integer VECTORS.  THE TWO integer ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  contains THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U contains THE 
!                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT contains 
!                 THE LATEST ESTIMATE TO THE SOLUTION.    
!          IWKSP  integer VECTOR WORKSPACE OF LENGTH 3*N  
!          NW     INPUT integer.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
!                 IPARM(8) IS AMOUNT USED.      
!          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  JACOBI CONJUGATE 
!                 GRADIENT NEEDS THIS TO BE IN LENGTH AT LEAST      
!                 4*N + 2*ITMAX,  IF ISYM = 0  (SYMMETRIC STORAGE)  
!                 4*N + 4*ITMAX,  IF ISYM = 1  (NONSYMMETRIC STORAGE) 
!                 HERE ITMAX = IPARM(1) AND ISYM = IPARM(5) 
!                 (ITMAX IS THE MAXIMUM ALLOWABLE NUMBER OF ITERATIONS) 
!          IPARM  integer VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
!                 SOME integer PARAMETERS WHICH AFFECT THE METHOD.  
!          RPARM  D.P. VECTOR OF LENGTH 12. ALLOWS USER TO SPECIFY SOME 
!                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
!          IER    OUTPUT integer.  ERROR FLAG. (= IERR)   
!       
! ... JCG SUBPROGRAM REFERENCES:      
!       
!          FROM ITPACK    BISRCH, CHGCON, DETERM, DFAULT, ECHALL,   
!                         ECHOUT, EIGVNS, EIGVSS, EQRT1S, ITERM, TIMER, 
!                         ITJCG, IVFILL, PARCON, PERMAT,  
!                         c_PERROR, PERVEC, PJAC2, PMULT, PRBNDX,      
!                         PSTOP, QSORT, DAXPY2, SBELM, SCAL, DCOPY2,  
!                         DDOT2, SUM3, UNSCAL, VEVMW, VFILL, VOUT,   
!                         WEVMW, ZBRENT 
!          SYSTEM         DABS, DLOG10, DBLE(AMAX0), DMAX1, MOD, DSQRT
!       
!     VERSION:  ITPACK 2C (MARCH 1982)
!       
!     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
!                       CENTER FOR NUMERICAL ANALYSIS     
!                       UNIVERSITY OF TEXAS     
!                       AUSTIN, TX  78712       
!                       (512) 471-1242
!       
!     FOR ADDITIONAL DETAILS ON THE   
!          (A) subroutine SEE TOMS ARTICLE 1982 
!          (B) ALGORITHM  SEE CNA REPORT 150    
!       
!     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
!       
!     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
!                          L. HAGEMAN, D. YOUNG 
!                          ACADEMIC PRESS, 1981 
!       
!     **************************************************  
!     *               IMPORTANT NOTE                   *  
!     *                                                *  
!     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
!     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
!     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
!     *                                                *  
!     *   DRELPR      MACHINE RELATIVE PRECISION       *  
!     *   RPARM(1)    STOPPING CRITERION               *  
!     *                                                *  
!     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
!     *   SECOND USED IN TIMER                         *  
!     *                                                *  
!     **************************************************  
!       
!     SPECIFICATIONS FOR ARGUMENTS    
!       
      integer IA(1),JA(1),IWKSP(1),IPARM(12),NN,NW,IERR   
      double precision A(1),RHS(NN),U(NN),WKSP(NW),RPARM(12)
!       
!     SPECIFICATIONS FOR LOCAL VARIABLES
!       
      integer IB1,IB2,IB3,IB4,IB5,IDGTS,IER,IERPER,ITMAX1,LOOP,N,NB,N3
      double precision DIGIT1,DIGIT2,TEMP,TIME1,TIME2,TOL 
!       
! **** BEGIN: ITPACK common 
!       
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
!       
! **** END  : ITPACK common 
!       
! ... VARIABLES IN common BLOCK - ITCOM1
!       
!     IN     - ITERATION NUMBER       
!     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
!     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE format SWITCH 
!     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
!     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
!     NOUT   - OUTPUT element number     
!       
! ... VARIABLES IN common BLOCK - ITCOM2
!       
!     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
!     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
!     CASEII - ADAPTIVE PROCEDURE case SWITCH   
!     HALT   - STOPPING TEST SWITCH   
!     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
!       
! ... VARIABLES IN common BLOCK - ITCOM3
!       
!     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
!     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
!     CME    - ESTIMATE OF LARGEST EIGENVALUE   
!     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
!     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
!     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
!     GAMMA  - ACCELERATION parameter 
!     OMEGA  - OVERRELAXATION parameter FOR SOR AND SSOR  
!     QA     - PSEUDO-RESIDUAL RATIO  
!     QT     - VIRTUAL SPECTRAL RADIUS
!     RHO    - ACCELERATION parameter 
!     RRR    - ADAPTIVE parameter     
!     SIGE   - parameter SIGMA-SUB-E  
!     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
!     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
!     DRELPR - MACHINE RELATIVE PRECISION       
!     STPTST - STOPPING parameter     
!     UDNM   - TWO NORM OF U
!     ZETA   - STOPPING CRITERION     
!       
! ... INITIALIZE common BLOCKS
!     
      
      LEVEL = IPARM(2)    
      
      NOUT = IPARM(4)       
      if (LEVEL.ge.1) write (NOUT,10) 
   10 format ('0'///1X,'BEGINNING OF ITPACK SOLUTION module  JCG')  
      IER = 0     
      if (IPARM(1).le.0) return       
      N = NN      
      if (IPARM(11).eq.0) TIMJ1 = TIMER(DUMMY)  
      if (LEVEL.ge.3) GO TO 20
      call ECHOUT (IPARM,RPARM,1)     
      GO TO 30    
   20 call ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      if (ZETA.ge.TEMP) GO TO 50      
      if (LEVEL.ge.1) write (NOUT,40) ZETA,DRELPR,TEMP    
40 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE JCG'/' ','    RPARM(1) =',D10.3, &
' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/ ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ', &
'    ZETA RESET TO ',D10.3)
      ZETA = TEMP 
   50 continue    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
!       
! ... VERIFY N    
!       
      if (N.gt.0) GO TO 70  
      IER = 11    
      if (LEVEL.ge.0) write (NOUT,60) N 
60 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE JCG '/' ', &
'    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 370   
   70 continue    
!       
! ... REMOVE ROWS AND COLUMNS IF REQUESTED      
!       
      if (IPARM(10).eq.0) GO TO 90    
      TOL = RPARM(8)
      call IVFILL (N,IWKSP,0) 
      call VFILL (N,WKSP,0.0D0)       
      call SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      if (IER.eq.0) GO TO 90
      if (LEVEL.ge.0) write (NOUT,80) IER,TOL   
80 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE JCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  SBELM '/' ', '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ', &
'    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X, ' RPARM(8) = ',D10.3,' (TOL)')
      GO TO 370   
!       
! ... INITIALIZE WKSP BASE ADDRESSES. 
!       
   90 IB1 = 1     
      IB2 = IB1+N 
      IB3 = IB2+N 
      IB4 = IB3+N 
      IB5 = IB4+N 
      IPARM(8) = 4*N+2*ITMAX
      if (ISYM.ne.0) IPARM(8) = IPARM(8)+2*ITMAX
      if (NW.ge.IPARM(8)) GO TO 110   
      IER = 12    
      if (LEVEL.ge.0) write (NOUT,100) NW,IPARM(8)
100 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE JCG '/' ', &
'    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10 ,' (NW)')
      GO TO 370   
!       
! ... PERMUTE TO  RED-BLACK SYSTEM IF REQUESTED 
!       
  110 NB = IPARM(9) 
      if (NB.lt.0) GO TO 170
      N3 = 3*N    
      call IVFILL (N3,IWKSP,0)
      call PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      if (IER.eq.0) GO TO 130 
      if (LEVEL.ge.0) write (NOUT,120) IER,NB   
120 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE JCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ', '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5 &
,' IPARM(9) = ',I5,' (NB)')
      GO TO 370   
!       
! ... PERMUTE MATRIX AND RHS
!       
  130 if (LEVEL.ge.2) write (NOUT,140) NB       
  140 format (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      call PERMAT (N,IA,JA,A,IWKSP,IWKSP(IB3),ISYM,LEVEL,NOUT,IER)  
      if (IER.eq.0) GO TO 160 
      if (LEVEL.ge.0) write (NOUT,150) IER      
150 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE JCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ', '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 370   
  160 call PERVEC (N,RHS,IWKSP)       
      call PERVEC (N,U,IWKSP) 
!       
! ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
! ... DIAGONAL ELEMENTS.    
!       
  170 continue    
      call VFILL (IPARM(8),WKSP,0.0D0)
      call SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)     
      if (IER.eq.0) GO TO 190 
      if (LEVEL.ge.0) write (NOUT,180) IER      
180 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE JCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ', '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)
      GO TO 370   
  190 if (LEVEL.le.2) GO TO 220       
      write (NOUT,200)      
200 format (///1X,'IN THE FOLLOWING, RHO AND GAMMA ARE', ' ACCELERATION PARAMETERS')
      if (ADAPT) write (NOUT,210)     
210 format (1X,'CME IS THE ESTIMATE OF THE LARGEST EIGENVALUE OF', ' THE JACOBI MATRIX')
  220 if (IPARM(11).ne.0) GO TO 230   
      TIMI1 = TIMER(DUMMY)  
!       
! ... COMPUTE INITIAL PSEUDO-RESIDUAL 
!       
  230 continue    
      call DCOPY2 (N,RHS,1,WKSP(IB2),1)
      call PJAC2 (N,IA,JA,A,U,WKSP(IB2)) 
      call VEVMW (N,WKSP(IB2),U)      
!       
! ... ITERATION SEQUENCE    
!       
      ITMAX1 = ITMAX+1      
      do 250 LOOP = 1,ITMAX1
         IN = LOOP-1
         if (MOD(IN,2).eq.1) GO TO 240
!       
! ... CODE FOR THE EVEN ITERATIONS.   
!       
!     U           = U(IN)             WKSP(IB2) = DEL(IN) 
!     WKSP(IB1)   = U(IN-1)           WKSP(IB3) = DEL(IN-1) 
!       
call ITJCG (N,IA,JA,A,U,WKSP(IB1),WKSP(IB2),WKSP(IB3),WKSP(IB4) ,WKSP(IB5))
!       
         if (HALT) GO TO 280
         GO TO 250
!       
! ... CODE FOR THE ODD ITERATIONS.    
!       
!     U           = U(IN-1)           WKSP(IB2) = DEL(IN-1) 
!     WKSP(IB1)   = U(IN)             WKSP(IB3) = DEL(IN) 
!       
240    call ITJCG (N,IA,JA,A,WKSP(IB1),U,WKSP(IB3),WKSP(IB2),WKSP(IB4) ,WKSP(IB5))
!       
         if (HALT) GO TO 280
  250 continue    
!       
! ... ITMAX HAS BEEN REACHED
!       
      if (IPARM(11).ne.0) GO TO 260   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  260 IER = 13    
      if (LEVEL.ge.1) write (NOUT,270) ITMAX    
270 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE JCG'/' ','    FAILURE TO CONVERGE IN ', &
I5 ,' ITERATIONS')
      if (IPARM(3).eq.0) RPARM(1) = STPTST      
      GO TO 310   
!       
! ... METHOD HAS CONVERGED  
!       
  280 if (IPARM(11).ne.0) GO TO 290   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  290 if (LEVEL.ge.1) write (NOUT,300) IN       
  300 format (/1X,'JCG  HAS CONVERGED IN ',I5,' ITERATIONS')
!       
! ... PUT SOLUTION INTO U IF NOT ALREADY THERE. 
!       
  310 continue    
      if (MOD(IN,2).eq.1) call DCOPY2 (N,WKSP(IB1),1,U,1)  
!       
! ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
!       
      call UNSCAL (N,IA,JA,A,RHS,U,WKSP)
!       
! ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
!       
      if (IPARM(9).lt.0) GO TO 340    
call PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(IB3),ISYM,LEVEL,NOUT, IERPER)
      if (IERPER.eq.0) GO TO 330      
      if (LEVEL.ge.0) write (NOUT,320) IERPER   
320 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE JCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ', '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ', '    IER = ',I5)
      if (IER.eq.0) IER = IERPER      
      GO TO 370   
  330 call PERVEC (N,RHS,IWKSP(IB2))  
      call PERVEC (N,U,IWKSP(IB2))    
!       
! ... optional ERROR ANALYSIS 
!       
  340 IDGTS = IPARM(12)     
      if (IDGTS.lt.0) GO TO 350       
      if (IPARM(2).le.0) IDGTS = 0    
      call c_PERROR (N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)
!       
! ... SET return PARAMETERS IN IPARM AND RPARM  
!       
  350 IPARM(8) = IPARM(8)-2*(ITMAX-IN)
      if (IPARM(11).ne.0) GO TO 360   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  360 if (ISYM.ne.0) IPARM(8) = IPARM(8)-2*(ITMAX-IN)     
      if (IPARM(3).ne.0) GO TO 370    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
!       
  370 continue    
      IERR = IER  
      if (LEVEL.ge.3) call ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
!       
      return      
      END 
      subroutine JSI (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,IERR)
!       
!     ITPACK 2C MAIN subroutine  JSI  (JACOBI SEMI-ITERATIVE)       
!     EACH OF THE MAIN SUBROUTINES:   
!           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
!     CAN BE USED INDEPENDENTLY OF THE OTHERS   
!       
! ... function:   
!       
!          THIS subroutine, JSI, DRIVES THE JACOBI SEMI-  
!          ITERATION ALGORITHM.       
!       
! ... parameter LIST:       
!       
!          N      INPUT integer.  dimension OF THE MATRIX. (= NN)   
!          IA,JA  INPUT integer VECTORS.  THE TWO integer ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  contains THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U contains THE 
!                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT contains 
!                 THE LATEST ESTIMATE TO THE SOLUTION.    
!          IWKSP  integer VECTOR WORKSPACE OF LENGTH 3*N  
!          NW     INPUT integer.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
!                 IPARM(8) IS AMOUNT USED.      
!          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  JACOBI SI    
!                 NEEDS THIS TO BE IN LENGTH AT LEAST     
!                 2*N       
!          IPARM  integer VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
!                 SOME integer PARAMETERS WHICH AFFECT THE METHOD.  
!          RPARM  D.P. VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY SOME
!                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
!          IER    OUTPUT integer.  ERROR FLAG. (= IERR)   
!       
! ... JSI SUBPROGRAM REFERENCES:      
!       
!          FROM ITPACK   BISRCH, CHEBY, CHGSI, CHGSME, DFAULT, ECHALL,
!                        ECHOUT, ITERM, TIMER, ITJSI, IVFILL, PAR   
!                        PERMAT, c_PERROR, PERVEC, PJAC2, PMULT, PRBNDX, 
!                        PSTOP, PVTBV, QSORT, DAXPY2, SBELM, SCAL,   
!                        DCOPY2, DDOT2, SUM3, TSTCHG, UNSCAL, VEVMW,  
!                        VFILL, VOUT, WEVMW     
!          SYSTEM        DABS, DLOG10, DBLE(AMAX0), DMAX1, DBLE(FLOAT), 
!                        MOD,DSQRT    
!       
!     VERSION:  ITPACK 2C (MARCH 1982)
!       
!     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
!                       CENTER FOR NUMERICAL ANALYSIS     
!                       UNIVERSITY OF TEXAS     
!                       AUSTIN, TX  78712       
!                       (512) 471-1242
!       
!     FOR ADDITIONAL DETAILS ON THE   
!          (A) subroutine SEE TOMS ARTICLE 1982 
!          (B) ALGORITHM  SEE CNA REPORT 150    
!       
!     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
!       
!     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
!                          L. HAGEMAN, D. YOUNG 
!                          ACADEMIC PRESS, 1981 
!       
!     **************************************************  
!     *               IMPORTANT NOTE                   *  
!     *                                                *  
!     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
!     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
!     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
!     *                                                *  
!     *   DRELPR      MACHINE RELATIVE PRECISION       *  
!     *   RPARM(1)    STOPPING CRITERION               *  
!     *                                                *  
!     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
!     *   SECOND USED IN TIMER                         *  
!     *                                                *  
!     **************************************************  
!       
!     SPECIFICATIONS FOR ARGUMENTS    
!       
      integer IA(1),JA(1),IWKSP(1),IPARM(12),NN,NW,IERR   
      double precision A(1),RHS(NN),U(NN),WKSP(NW),RPARM(12)
!       
!     SPECIFICATIONS FOR LOCAL VARIABLES
!       
      integer IB1,IB2,IB3,ICNT,IDGTS,IER,IERPER,ITMAX1,LOOP,N,NB,N3 
      double precision DIGIT1,DIGIT2,TEMP,TIME1,TIME2,TOL 
!       
! *** BEGIN: ITPACK common  
!       
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
!       
! *** END  : ITPACK common  
!       
! ... VARIABLES IN common BLOCK - ITCOM1
!       
!     IN     - ITERATION NUMBER       
!     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
!     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE format SWITCH 
!     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
!     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
!     NOUT   - OUTPUT element number     
!       
! ... VARIABLES IN common BLOCK - ITCOM2
!       
!     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
!     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
!     CASEII - ADAPTIVE PROCEDURE case SWITCH   
!     HALT   - STOPPING TEST SWITCH   
!     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
!       
! ... VARIABLES IN common BLOCK - ITCOM3
!       
!     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
!     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
!     CME    - ESTIMATE OF LARGEST EIGENVALUE   
!     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
!     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
!     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
!     GAMMA  - ACCELERATION parameter 
!     OMEGA  - OVERRELAXATION parameter FOR SOR AND SSOR  
!     QA     - PSEUDO-RESIDUAL RATIO  
!     QT     - VIRTUAL SPECTRAL RADIUS
!     RHO    - ACCELERATION parameter 
!     RRR    - ADAPTIVE parameter     
!     SIGE   - parameter SIGMA-SUB-E  
!     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
!     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
!     DRELPR - MACHINE RELATIVE PRECISION       
!     STPTST - STOPPING parameter     
!     UDNM   - TWO NORM OF U
!     ZETA   - STOPPING CRITERION     
!       
! ... INITIALIZE common BLOCKS
!       
      LEVEL = IPARM(2)      
      NOUT = IPARM(4)       
!     if (LEVEL.ge.1) write (NOUT,10) 
!  10 format (8X,'BEGINNING OF ITPACK SOLUTION module  JSI')  
      IER = 0     
      if (IPARM(1).le.0) return       
      N = NN      
      if (IPARM(11).eq.0) TIMJ1 = TIMER(DUMMY)  
      if (LEVEL.ge.3) GO TO 20
      call ECHOUT (IPARM,RPARM,2)     
      GO TO 30    
   20 call ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      if (ZETA.ge.TEMP) GO TO 50      
      if (LEVEL.ge.1) write (NOUT,40) ZETA,DRELPR,TEMP    
40 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE JSI'/' ','    RPARM(1) =',D10.3, &
' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/ ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ', &
'    ZETA RESET TO ',D10.3)
      ZETA = TEMP 
   50 continue    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
!       
! ... VERIFY N    
!       
      if (N.gt.0) GO TO 70  
      IER = 21    
      if (LEVEL.ge.0) write (NOUT,60) N 
60 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE JSI '/' ', &
'    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 360   
   70 continue    
!       
! ... REMOVE ROWS AND COLUMNS IF REQUESTED      
!       
      if (IPARM(10).eq.0) GO TO 90    
      TOL = RPARM(8)
      call IVFILL (N,IWKSP,0) 
      call VFILL (N,WKSP,0.0D0)       
      call SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      if (IER.eq.0) GO TO 90
      if (LEVEL.ge.0) write (NOUT,80) IER,TOL   
80 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE JSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  SBELM '/' ', '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ', &
'    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X, ' RPARM(8) = ',D10.3,' (TOL)')
      GO TO 360   
!       
! ... INITIALIZE WKSP BASE ADDRESSES. 
!       
   90 IB1 = 1     
      IB2 = IB1+N 
      IB3 = IB2+N 
      IPARM(8) = 2*N
      if (NW.ge.IPARM(8)) GO TO 110   
      IER = 22    
      if (LEVEL.ge.0) write (NOUT,100) NW,IPARM(8)
100 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE JSI '/' ', &
'    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10 ,' (NW)')
      GO TO 360   
!       
! ... PERMUTE TO  RED-BLACK SYSTEM IF REQUESTED 
!       
  110 NB = IPARM(9) 
      if (NB.lt.0) GO TO 170
      N3 = 3*N    
      call IVFILL (N3,IWKSP,0)
      call PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      if (IER.eq.0) GO TO 130 
      if (LEVEL.ge.0) write (NOUT,120) IER,NB   
120 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE JSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ', '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5 &
,' IPARM(9) = ',I5,' (NB)')
      GO TO 360   
!       
! ... PERMUTE MATRIX AND RHS
!       
  130 if (LEVEL.ge.2) write (NOUT,140) NB       
  140 format (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      call PERMAT (N,IA,JA,A,IWKSP,IWKSP(IB3),ISYM,LEVEL,NOUT,IER)  
      if (IER.eq.0) GO TO 160 
      if (LEVEL.ge.0) write (NOUT,150) IER      
150 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE JSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ', '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 360   
  160 call PERVEC (N,RHS,IWKSP)       
      call PERVEC (N,U,IWKSP) 
!       
! ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
! ... DIAGONAL ELEMENTS.    
!       
  170 continue    
      call VFILL (IPARM(8),WKSP,0.0D0)
      call SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)     
      !SCAL (NN,IA,JA,A,RHS,U,D,LEVEL,NOUT,IER) 
      
      if (IER.eq.0) GO TO 190 
      if (LEVEL.ge.0) write (NOUT,180) IER      
180 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE JSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ', '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)
      GO TO 360   
  190 if (LEVEL.le.2) GO TO 210       
      write (NOUT,200)      
200 format (///1X,'IN THE FOLLOWING, RHO AND GAMMA ARE', ' ACCELERATION PARAMETERS')
  210 if (IPARM(11).ne.0) GO TO 220   
      TIMI1 = TIMER(DUMMY)  
!       
! ... ITERATION SEQUENCE    
!       
  220 ITMAX1 = ITMAX+1      
      do 240 LOOP = 1,ITMAX1
         IN = LOOP-1
         if (MOD(IN,2).eq.1) GO TO 230
!       
! ... CODE FOR THE EVEN ITERATIONS.   
!       
!     U           = U(IN)   
!     WKSP(IB1)   = U(IN-1) 
!       
         call ITJSI (N,IA,JA,A,RHS,U,WKSP(IB1),WKSP(IB2),ICNT)      
!       
         if (HALT) GO TO 270
         GO TO 240
!       
! ... CODE FOR THE ODD ITERATIONS.    
!       
!     U           = U(IN-1) 
!     WKSP(IB1)   = U(IN)   
!       
  230    call ITJSI (N,IA,JA,A,RHS,WKSP(IB1),U,WKSP(IB2),ICNT)      
!       
         if (HALT) GO TO 270
  240 continue    
!       
! ... ITMAX HAS BEEN REACHED
!       
      if (IPARM(11).ne.0) GO TO 250   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  250 IER = 23    
      if (LEVEL.ge.1) write (NOUT,260) ITMAX    
260 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE JSI'/' ','    FAILURE TO CONVERGE IN ',I5 &
,' ITERATIONS')
      if (IPARM(3).eq.0) RPARM(1) = STPTST      
      GO TO 300   
!       
! ... METHOD HAS CONVERGED  
!       
  270 if (IPARM(11).ne.0) GO TO 280   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  280 if (LEVEL.ge.1) write (NOUT,290) IN       
  290 format (9X,'  JSI has converged after ',I5,' interations')
!       
! ... PUT SOLUTION INTO U IF NOT ALREADY THERE. 
!       
  300 continue    
      if (MOD(IN,2).eq.1) call DCOPY2 (N,WKSP(IB1),1,U,1)  
!       
! ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
!       
      call UNSCAL (N,IA,JA,A,RHS,U,WKSP)
!       
! ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
!       
      if (IPARM(9).lt.0) GO TO 330    
call PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(IB3),ISYM,LEVEL,NOUT, IERPER)
      if (IERPER.eq.0) GO TO 320      
      if (LEVEL.ge.0) write (NOUT,310) IERPER   
310 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE JSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ', '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ', '    IER = ',I5)
      if (IER.eq.0) IER = IERPER      
      GO TO 360   
  320 call PERVEC (N,RHS,IWKSP(IB2))  
      call PERVEC (N,U,IWKSP(IB2))    
!       
! ... optional ERROR ANALYSIS 
!       
  330 IDGTS = IPARM(12)     
      if (IDGTS.lt.0) GO TO 340       
      if (IPARM(2).le.0) IDGTS = 0    
      call c_PERROR(N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)

!       
! ... SET return PARAMETERS IN IPARM AND RPARM  
!       
  340 if (IPARM(11).ne.0) GO TO 350   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  350 if (IPARM(3).ne.0) GO TO 360    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
!       
  360 continue    
      IERR = IER  
      if (LEVEL.ge.3) call ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
!       
      return      
      END 
      subroutine SOR (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,IERR)
!       
!     ITPACK 2C MAIN subroutine  SOR  (SUCCESSIVE OVERRELATION)     
!     EACH OF THE MAIN SUBROUTINES:   
!           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
!     CAN BE USED INDEPENDENTLY OF THE OTHERS   
!       
! ... function:   
!       
!          THIS subroutine, SOR, DRIVES THE  SUCCESSIVE   
!          OVERRELAXATION ALGORITHM.  
!       
! ... parameter LIST:       
!       
!          N      INPUT integer.  dimension OF THE MATRIX. (= NN)   
!          IA,JA  INPUT integer VECTORS.  THE TWO integer ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION 
!          RHS    INPUT D.P. VECTOR.  contains THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U contains THE 
!                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT contains 
!                 THE LATEST ESTIMATE TO THE SOLUTION.    
!          IWKSP  integer VECTOR WORKSPACE OF LENGTH 3*N  
!          NW     INPUT integer.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
!                 IPARM(8) IS AMOUNT USED.      
!          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  SOR NEEDS THIS 
!                 TO BE IN LENGTH AT LEAST  N   
!          IPARM  integer VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
!                 SOME integer PARAMETERS WHICH AFFECT THE METHOD.  
!          RPARM  D.P. VECTOR OF LENGTH 12. ALLOWS USER TO SPECIFY SOME 
!                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
!          IER    OUTPUT integer.  ERROR FLAG. (= IERR)   
!       
! ... SOR SUBPROGRAM REFERENCES:      
!       
!          FROM ITPACK   BISRCH, DFAULT, ECHALL, ECHOUT, IPSTR, ITERM,
!                        TIMER, ITSOR, IVFILL, PERMAT, c_PERROR,      
!                        PERVEC, PFSOR1, PMULT, PRBNDX, PSTOP, QSORT, 
!                        SBELM, SCAL, DCOPY2, DDOT2, TAU, UNSCAL, VFILL,
!                        VOUT, WEVMW  
!          SYSTEM        DABS, DLOG10, DBLE(AMAX0), DMAX1, DBLE(FLOAT), 
!                        DSQRT
!       
!     VERSION:  ITPACK 2C (MARCH 1982)
!       
!     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
!                       CENTER FOR NUMERICAL ANALYSIS     
!                       UNIVERSITY OF TEXAS     
!                       AUSTIN, TX  78712       
!                       (512) 471-1242
!       
!     FOR ADDITIONAL DETAILS ON THE   
!          (A) subroutine SEE TOMS ARTICLE 1982 
!          (B) ALGORITHM  SEE CNA REPORT 150    
!       
!     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
!       
!     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
!                          L. HAGEMAN, D. YOUNG 
!                          ACADEMIC PRESS, 1981 
!       
!     **************************************************  
!     *               IMPORTANT NOTE                   *  
!     *                                                *  
!     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
!     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
!     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
!     *                                                *  
!     *   DRELPR      MACHINE RELATIVE PRECISION       *  
!     *   RPARM(1)    STOPPING CRITERION               *  
!     *                                                *  
!     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
!     *   SECOND USED IN TIMER                         *  
!     *                                                *  
!     **************************************************  
!       
!     SPECIFICATIONS FOR ARGUMENTS    
!       
      integer IA(1),JA(1),IWKSP(1),IPARM(12),NN,NW,IERR   
      double precision A(1),RHS(NN),U(NN),WKSP(NW),RPARM(12)
!       
!     SPECIFICATIONS FOR LOCAL VARIABLES
!       
      integer IB1,IB2,IB3,IDGTS,IER,IERPER,ITMAX1,LOOP,N,NB,N3      
      double precision DIGIT1,DIGIT2,TEMP,TIME1,TIME2,TOL 
!       
! *** BEGIN: ITPACK common  
!       
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
!       
! *** END  : ITPACK common  
!       
! ... VARIABLES IN common BLOCK - ITCOM1
!       
!     IN     - ITERATION NUMBER       
!     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
!     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE format SWITCH 
!     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
!     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
!     NOUT   - OUTPUT element number     
!       
! ... VARIABLES IN common BLOCK - ITCOM2
!       
!     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
!     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
!     CASEII - ADAPTIVE PROCEDURE case SWITCH   
!     HALT   - STOPPING TEST SWITCH   
!     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
!       
! ... VARIABLES IN common BLOCK - ITCOM3
!       
!     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
!     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
!     CME    - ESTIMATE OF LARGEST EIGENVALUE   
!     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
!     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
!     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
!     GAMMA  - ACCELERATION parameter 
!     OMEGA  - OVERRELAXATION parameter FOR SOR AND SSOR  
!     QA     - PSEUDO-RESIDUAL RATIO  
!     QT     - VIRTUAL SPECTRAL RADIUS
!     RHO    - ACCELERATION parameter 
!     RRR    - ADAPTIVE parameter     
!     SIGE   - parameter SIGMA-SUB-E  
!     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
!     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
!     DRELPR - MACHINE RELATIVE PRECISION       
!     STPTST - STOPPING parameter     
!     UDNM   - TWO NORM OF U
!     ZETA   - STOPPING CRITERION     
!       
! ... INITIALIZE common BLOCKS
!       
      LEVEL = IPARM(2)      
      NOUT = IPARM(4)       
      if (LEVEL.ge.1) write (NOUT,10) 
   10 format ('0'///1X,'BEGINNING OF ITPACK SOLUTION module  SOR')  
      IER = 0     
      if (IPARM(1).le.0) return       
      N = NN      
      if (IPARM(11).eq.0) TIMJ1 = TIMER(DUMMY)  
      if (LEVEL.ge.3) GO TO 20
      call ECHOUT (IPARM,RPARM,3)     
      GO TO 30    
   20 call ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      if (ZETA.ge.TEMP) GO TO 50      
      if (LEVEL.ge.1) write (NOUT,40) ZETA,DRELPR,TEMP    
40 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE SOR'/' ','    RPARM(1) =',D10.3, &
' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/ ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ', &
'    ZETA RESET TO ',D10.3)
      ZETA = TEMP 
   50 continue    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
!       
! ... VERIFY N    
!       
      if (N.gt.0) GO TO 70  
      IER = 31    
      if (LEVEL.ge.0) write (NOUT,60) N 
60 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SOR '/' ', &
'    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 360   
   70 continue    
!       
! ... REMOVE ROWS AND COLUMNS IF REQUESTED      
!       
      if (IPARM(10).eq.0) GO TO 90    
      TOL = RPARM(8)
      call IVFILL (N,IWKSP,0) 
      call VFILL (N,WKSP,0.0D0)       
      call SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      if (IER.eq.0) GO TO 90
      if (LEVEL.ge.0) write (NOUT,80) IER,TOL   
80 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SOR '/' ', &
'    ERROR DETECTED IN SUBROUTINE  SBELM '/' ', '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ', &
'    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X, ' RPARM(8) = ',D10.3,' (TOL)')
      GO TO 360   
!       
! ... INITIALIZE WKSP BASE ADDRESSES. 
!       
   90 IB1 = 1     
      IB2 = IB1+N 
      IB3 = IB2+N 
      IPARM(8) = N
      if (NW.ge.IPARM(8)) GO TO 110   
      IER = 32    
      if (LEVEL.ge.0) write (NOUT,100) NW,IPARM(8)
100 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SOR '/' ', &
'    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10 ,' (NW)')
      GO TO 360   
!       
! ... PERMUTE TO  RED-BLACK SYSTEM IF REQUESTED 
!       
  110 NB = IPARM(9) 
      if (NB.lt.0) GO TO 170
      N3 = 3*N    
      call IVFILL (N3,IWKSP,0)
      call PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      if (IER.eq.0) GO TO 130 
      if (LEVEL.ge.0) write (NOUT,120) IER,NB   
120 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SOR '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ', '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5 &
,' IPARM(9) = ',I5,' (NB)')
      GO TO 360   
!       
! ... PERMUTE MATRIX AND RHS
!       
  130 if (LEVEL.ge.2) write (NOUT,140) NB       
  140 format (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      call PERMAT (N,IA,JA,A,IWKSP,IWKSP(IB3),ISYM,LEVEL,NOUT,IER)  
      if (IER.eq.0) GO TO 160 
      if (LEVEL.ge.0) write (NOUT,150) IER      
150 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SOR '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ', '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 360   
  160 call PERVEC (N,RHS,IWKSP)       
      call PERVEC (N,U,IWKSP) 
!       
! ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
! ... DIAGONAL ELEMENTS.    
!       
  170 continue    
      call VFILL (IPARM(8),WKSP,0.0D0)
      call SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)     
      if (IER.eq.0) GO TO 190 
      if (LEVEL.ge.0) write (NOUT,180) IER      
180 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SOR '/' ', &
'    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ', '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)
      GO TO 360   
  190 if (LEVEL.le.2) GO TO 220       
      if (ADAPT) write (NOUT,200)     
200 format (///1X,'CME IS THE ESTIMATE OF THE LARGEST EIGENVALUE OF', ' THE JACOBI MATRIX')
      write (NOUT,210)      
  210 format (1X,'OMEGA IS THE RELAXATION FACTOR')
  220 if (IPARM(11).ne.0) GO TO 230   
      TIMI1 = TIMER(DUMMY)  
!       
! ... ITERATION SEQUENCE    
!       
  230 ITMAX1 = ITMAX+1      
      do 240 LOOP = 1,ITMAX1
         IN = LOOP-1
!       
! ... CODE FOR ONE ITERATION. 
!       
!     U           = U(IN)   
!       
         call ITSOR (N,IA,JA,A,RHS,U,WKSP(IB1)) 
!       
         if (HALT) GO TO 270
  240 continue    
!       
! ... ITMAX HAS BEEN REACHED
!       
      if (IPARM(11).ne.0) GO TO 250   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  250 if (LEVEL.ge.1) write (NOUT,260) ITMAX    
260 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE SOR'/' ','    FAILURE TO CONVERGE IN ',I5 &
,' ITERATIONS')
      IER = 33    
      if (IPARM(3).eq.0) RPARM(1) = STPTST      
      GO TO 300   
!       
! ... METHOD HAS CONVERGED  
!       
  270 if (IPARM(11).ne.0) GO TO 280   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  280 if (LEVEL.ge.1) write (NOUT,290) IN       
  290 format (/1X,'SOR  HAS CONVERGED IN ',I5,' ITERATIONS')
!       
! ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
!       
  300 call UNSCAL (N,IA,JA,A,RHS,U,WKSP)
!       
! ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
!       
      if (IPARM(9).lt.0) GO TO 330    
call PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(IB3),ISYM,LEVEL,NOUT, IERPER)
      if (IERPER.eq.0) GO TO 320      
      if (LEVEL.ge.0) write (NOUT,310) IERPER   
310 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SOR '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ', '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ', '    IER = ',I5)
      if (IER.eq.0) IER = IERPER      
      GO TO 360   
  320 call PERVEC (N,RHS,IWKSP(IB2))  
      call PERVEC (N,U,IWKSP(IB2))    
!       
! ... optional ERROR ANALYSIS 
!       
  330 IDGTS = IPARM(12)     
      if (IDGTS.lt.0) GO TO 340       
      if (IPARM(2).le.0) IDGTS = 0    
      call c_PERROR (N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)
!       
! ... SET return PARAMETERS IN IPARM AND RPARM  
!       
  340 if (IPARM(11).ne.0) GO TO 350   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  350 if (IPARM(3).ne.0) GO TO 360    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(5) = OMEGA      
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
!       
  360 continue    
      IERR = IER  
      if (LEVEL.ge.3) call ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
!       
      return      
      END 
subroutine SSORCG (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM, IERR)
!       
!     ITPACK 2C MAIN subroutine  SSORCG  (SYMMETRIC SUCCESSIVE OVER-
!                                        RELAXATION CONJUGATE GRADIENT) 
!     EACH OF THE MAIN SUBROUTINES:   
!           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
!     CAN BE USED INDEPENDENTLY OF THE OTHERS   
!       
! ... function:   
!       
!          THIS subroutine, SSORCG, DRIVES THE  SYMMETRIC SOR-CG    
!          ALGORITHM.       
!       
! ... parameter LIST:       
!       
!          N      INPUT integer.  dimension OF THE MATRIX. (= NN)   
!          IA,JA  INPUT integer VECTORS.  THE TWO integer ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  contains THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U contains THE 
!                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT contains 
!                 THE LATEST ESTIMATE TO THE SOLUTION.    
!          IWKSP  integer VECTOR WORKSPACE OF LENGTH 3*N  
!          NW     INPUT integer.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
!                 IPARM(8) IS AMOUNT USED.      
!          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  SSOR-CG      
!                 NEEDS TO BE IN LENGTH AT LEAST
!                 6*N + 2*ITMAX,  IF IPARM(5)=0  (SYMMETRIC STORAGE)
!                 6*N + 4*ITMAX,  IF IPARM(5)=1  (NONSYMMETRIC STORAGE) 
!          IPARM  integer VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
!                 SOME integer PARAMETERS WHICH AFFECT THE METHOD.  IF
!          RPARM  D.P. VECTOR OF LENGTH 12. ALLOWS USER TO SPECIFY SOME 
!                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
!          IER    OUTPUT integer.  ERROR FLAG. (= IERR)   
!       
! ... SSORCG SUBPROGRAM REFERENCES:   
!       
!          FROM ITPACK    BISRCH, CHGCON, DETERM, DFAULT, ECHALL,   
!                         ECHOUT, EIGVNS, EIGVSS, EQRT1S, ITERM, TIMER, 
!                         ITSRCG, IVFILL, OMEG, OMGCHG, OMGSTR,     
!                         PARCON, PBETA, PBSOR, PERMAT, c_PERROR,     
!                         PERVEC, PFSOR, PJAC2, PMULT, PRBNDX, PSTOP, PVT
!                         QSORT, SBELM, SCAL, DCOPY2, DDOT2, SUM3,    
!                         UNSCAL, VEVMW, VEVPW, VFILL, VOUT, WEVMW, 
!                         ZBRENT      
!          SYSTEM         DABS, DLOG, DLOG10, DBLE(AMAX0), DMAX1, AMIN1,
!                         MOD, DSQRT  
!       
!     VERSION:  ITPACK 2C (MARCH 1982)
!       
!     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
!                       CENTER FOR NUMERICAL ANALYSIS     
!                       UNIVERSITY OF TEXAS     
!                       AUSTIN, TX  78712       
!                       (512) 471-1242
!       
!     FOR ADDITIONAL DETAILS ON THE   
!          (A) subroutine SEE TOMS ARTICLE 1982 
!          (B) ALGORITHM  SEE CNA REPORT 150    
!       
!     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
!       
!     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
!                          L. HAGEMAN, D. YOUNG 
!                          ACADEMIC PRESS, 1981 
!       
!     **************************************************  
!     *               IMPORTANT NOTE                   *  
!     *                                                *  
!     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
!     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
!     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
!     *                                                *  
!     *   DRELPR      MACHINE RELATIVE PRECISION       *  
!     *   RPARM(1)    STOPPING CRITERION               *  
!     *                                                *  
!     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
!     *   SECOND USED IN TIMER                         *  
!     *                                                *  
!     **************************************************  
!       
!     SPECIFICATIONS FOR ARGUMENTS    
!       
      integer IA(1),JA(1),IWKSP(1),IPARM(12),NN,NW,IERR   
      double precision A(1),RHS(NN),U(NN),WKSP(NW),RPARM(12)
!       
!     SPECIFICATIONS FOR LOCAL VARIABLES
!       
integer IB1,IB2,IB3,IB4,IB5,IB6,IB7,IDGTS,IER,IERPER,ITMAX1,LOOP,N ,NB,N3
      double precision BETNEW,DIGIT1,DIGIT2,PBETA,TEMP,TIME1,TIME2,TOL
!       
! *** BEGIN: ITPACK common  
!       
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
!       
! *** END  : ITPACK common  
!       
! ... VARIABLES IN common BLOCK - ITCOM1
!       
!     IN     - ITERATION NUMBER       
!     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
!     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE format SWITCH 
!     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
!     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
!     NOUT   - OUTPUT element number     
!       
! ... VARIABLES IN common BLOCK - ITCOM2
!       
!     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
!     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
!     CASEII - ADAPTIVE PROCEDURE case SWITCH   
!     HALT   - STOPPING TEST SWITCH   
!     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
!       
! ... VARIABLES IN common BLOCK - ITCOM3
!       
!     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
!     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
!     CME    - ESTIMATE OF LARGEST EIGENVALUE   
!     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
!     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
!     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
!     GAMMA  - ACCELERATION parameter 
!     OMEGA  - OVERRELAXATION parameter FOR SOR AND SSOR  
!     QA     - PSEUDO-RESIDUAL RATIO  
!     QT     - VIRTUAL SPECTRAL RADIUS
!     RHO    - ACCELERATION parameter 
!     RRR    - ADAPTIVE parameter     
!     SIGE   - parameter SIGMA-SUB-E  
!     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
!     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
!     DRELPR - MACHINE RELATIVE PRECISION       
!     STPTST - STOPPING parameter     
!     UDNM   - TWO NORM OF U
!     ZETA   - STOPPING CRITERION     
!       
! ... INITIALIZE common BLOCKS
!       
      LEVEL = IPARM(2)      
      NOUT = IPARM(4)       
      if (IPARM(9).ge.0) IPARM(6) = 2 
      if (LEVEL.ge.1) write (NOUT,10) 
   10 format ('0'///1X,'BEGINNING OF ITPACK SOLUTION module  SSORCG') 
      IER = 0     
      if (IPARM(1).le.0) return       
      N = NN      
      if (IPARM(11).eq.0) TIMJ1 = TIMER(DUMMY)  
      if (LEVEL.ge.3) GO TO 20
      call ECHOUT (IPARM,RPARM,4)     
      GO TO 30    
   20 call ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      if (ZETA.ge.TEMP) GO TO 50      
      if (LEVEL.ge.1) write (NOUT,40) ZETA,DRELPR,TEMP    
40 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE SSORCG'/' ','    RPARM(1) =',D10.3, &
' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/ ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ', &
'    ZETA RESET TO ',D10.3)
      ZETA = TEMP 
   50 continue    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
!       
! ... VERIFY N    
!       
      if (N.gt.0) GO TO 70  
      IER = 41    
      if (LEVEL.ge.0) write (NOUT,60) N 
60 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SSORCG '/' ', &
'    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 390   
   70 continue    
!       
! ... REMOVE ROWS AND COLUMNS IF REQUESTED      
!       
      if (IPARM(10).eq.0) GO TO 90    
      TOL = RPARM(8)
      call IVFILL (N,IWKSP,0) 
      call VFILL (N,WKSP,0.0D0)       
      call SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      if (IER.eq.0) GO TO 90
      if (LEVEL.ge.0) write (NOUT,80) IER,TOL   
80 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SSORCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  SBELM '/' ', '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ', &
'    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X, ' RPARM(8) = ',D10.3,' (TOL)')
      GO TO 390   
!       
! ... INITIALIZE WKSP BASE ADDRESSES. 
!       
   90 IB1 = 1     
      IB2 = IB1+N 
      IB3 = IB2+N 
      IB4 = IB3+N 
      IB5 = IB4+N 
      IB6 = IB5+N 
      IB7 = IB6+N 
      IPARM(8) = 6*N+2*ITMAX
      if (ISYM.ne.0) IPARM(8) = IPARM(8)+2*ITMAX
      if (NW.ge.IPARM(8)) GO TO 110   
      IER = 42    
      if (LEVEL.ge.0) write (NOUT,100) NW,IPARM(8)
100 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SSORCG '/' ', &
'    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10 ,' (NW)')
      GO TO 390   
!       
! ... PERMUTE TO  RED-BLACK SYSTEM IF REQUESTED 
!       
  110 NB = IPARM(9) 
      if (NB.lt.0) GO TO 170
      N3 = 3*N    
      call IVFILL (N3,IWKSP,0)
      call PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      if (IER.eq.0) GO TO 130 
      if (LEVEL.ge.0) write (NOUT,120) IER,NB   
120 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SSORCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ', '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5 &
,' IPARM(9) = ',I5,' (NB)')
      GO TO 390   
!       
! ... PERMUTE MATRIX AND RHS
!       
  130 if (LEVEL.ge.2) write (NOUT,140) NB       
  140 format (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      call PERMAT (N,IA,JA,A,IWKSP,IWKSP(IB3),ISYM,LEVEL,NOUT,IER)  
      if (IER.eq.0) GO TO 160 
      if (LEVEL.ge.0) write (NOUT,150) IER      
150 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SSORCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ', '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 390   
  160 call PERVEC (N,RHS,IWKSP)       
      call PERVEC (N,U,IWKSP) 
!       
! ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
! ... DIAGONAL ELEMENTS.    
!       
  170 continue    
      call VFILL (IPARM(8),WKSP,0.0D0)
      call SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)     
      if (IER.eq.0) GO TO 190 
      if (LEVEL.ge.0) write (NOUT,180) IER      
180 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SSORCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ', '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)
      GO TO 390   
  190 if (LEVEL.le.2) GO TO 220       
      write (NOUT,200)      
200 format (///1X,'IN THE FOLLOWING, RHO AND GAMMA ARE', ' ACCELERATION PARAMETERS')
      write (NOUT,210)      
  210 format (1X,'S-PRIME IS AN INITIAL ESTIMATE FOR NEW CME')      
  220 if (IPARM(11).ne.0) GO TO 230   
      TIMI1 = TIMER(DUMMY)  
!       
! ... SPECIAL PROCEDURE FOR FULLY ADAPTIVE case.
!       
  230 continue    
      if (.not.ADAPT) GO TO 250       
      if (.not.BETADT) GO TO 240      
      call VFILL (N,WKSP(IB1),1.D0)   
BETNEW = PBETA(N,IA,JA,A,WKSP(IB1),WKSP(IB2),WKSP(IB3))/ DBLE(FLOAT(N))
      BETAB = DMAX1(BETAB,.25D0,BETNEW) 
  240 call OMEG (0.D0,1)    
      IS = 0      
!       
! ... INITIALIZE FORWARD PSEUDO-RESIDUAL
!       
  250 call DCOPY2 (N,RHS,1,WKSP(IB1),1)
      call DCOPY2 (N,U,1,WKSP(IB2),1)  
      call PFSOR (N,IA,JA,A,WKSP(IB2),WKSP(IB1))
      call VEVMW (N,WKSP(IB2),U)      
!       
! ... ITERATION SEQUENCE    
!       
      ITMAX1 = ITMAX+1      
      do 270 LOOP = 1,ITMAX1
         IN = LOOP-1
         if (MOD(IN,2).eq.1) GO TO 260
!       
! ... CODE FOR THE EVEN ITERATIONS.   
!       
!     U           = U(IN)       WKSP(IB2) = C(IN) 
!     WKSP(IB1)   = U(IN-1)     WKSP(IB3) = C(IN-1)       
!       
call ITSRCG (N,IA,JA,A,RHS,U,WKSP(IB1),WKSP(IB2),WKSP(IB3), WKSP(IB4),WKSP(IB5),WKSP(IB6),WKSP(IB7))
!       
         if (HALT) GO TO 300
         GO TO 270
!       
! ... CODE FOR THE ODD ITERATIONS.    
!       
!     U           = U(IN-1)     WKSP(IB2) = C(IN-1)       
!     WKSP(IB1)   = U(IN)       WKSP(IB3) =C(IN)
!       
260    call ITSRCG (N,IA,JA,A,RHS,WKSP(IB1),U,WKSP(IB3),WKSP(IB2), WKSP(IB4),WKSP(IB5),WKSP(IB6),WKSP(IB7))
!       
         if (HALT) GO TO 300
  270 continue    
!       
! ... ITMAX HAS BEEN REACHED
!       
      if (IPARM(11).ne.0) GO TO 280   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  280 if (LEVEL.ge.1) write (NOUT,290) ITMAX    
290 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE SSORCG'/' ','    FAILURE TO CONVERGE IN' &
,I5,' ITERATIONS')
      IER = 43    
      if (IPARM(3).eq.0) RPARM(1) = STPTST      
      GO TO 330   
!       
! ... METHOD HAS CONVERGED  
!       
  300 if (IPARM(11).ne.0) GO TO 310   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  310 if (LEVEL.ge.1) write (NOUT,320) IN       
  320 format (/1X,'SSORCG  HAS CONVERGED IN ',I5,' ITERATIONS')     
!       
! ... PUT SOLUTION INTO U IF NOT ALREADY THERE. 
!       
  330 continue    
      if (MOD(IN,2).eq.1) call DCOPY2 (N,WKSP(IB1),1,U,1)  
!       
! ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
!       
      call UNSCAL (N,IA,JA,A,RHS,U,WKSP)
!       
! ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
!       
      if (IPARM(9).lt.0) GO TO 360    
call PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(IB3),ISYM,LEVEL,NOUT, IERPER)
      if (IERPER.eq.0) GO TO 350      
      if (LEVEL.ge.0) write (NOUT,340) IERPER   
340 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SSORCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ', '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ', '    IER = ',I5)
      if (IER.eq.0) IER = IERPER      
      GO TO 390   
  350 call PERVEC (N,RHS,IWKSP(IB2))  
      call PERVEC (N,U,IWKSP(IB2))    
!       
! ... optional ERROR ANALYSIS 
!       
  360 IDGTS = IPARM(12)     
      if (IDGTS.lt.0) GO TO 370       
      if (IPARM(2).le.0) IDGTS = 0    
      call c_PERROR (N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)
!       
! ... SET return PARAMETERS IN IPARM AND RPARM  
!       
  370 if (IPARM(11).ne.0) GO TO 380   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  380 IPARM(8) = IPARM(8)-2*(ITMAX-IN)
      if (ISYM.ne.0) IPARM(8) = IPARM(8)-2*(ITMAX-IN)     
      if (IPARM(3).ne.0) GO TO 390    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(5) = OMEGA      
      RPARM(6) = SPECR      
      RPARM(7) = BETAB      
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
!       
  390 continue    
      IERR = IER  
      if (LEVEL.ge.3) call ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
!       
      return      
      END 
subroutine SSORSI (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM, IERR)
!       
!     ITPACK 2C MAIN subroutine  SSORSI  (SYMMETRIC SUCCESSIVE RELAX- 
!                                         ATION SEMI-ITERATION)     
!     EACH OF THE MAIN SUBROUTINES:   
!           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
!     CAN BE USED INDEPENDENTLY OF THE OTHERS   
!       
! ... function:   
!       
!          THIS subroutine, SSORSI, DRIVES THE  SYMMETRIC SOR-SI    
!          ALGORITHM.       
!       
! ... parameter LIST:       
!       
!          N      INPUT integer.  dimension OF THE MATRIX. (= NN)   
!          IA,JA  INPUT integer VECTORS.  THE TWO integer ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  contains THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U contains THE 
!                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT contains 
!                 THE LATEST ESTIMATE TO THE SOLUTION.    
!          IWKSP  integer VECTOR WORKSPACE OF LENGTH 3*N  
!          NW     INPUT integer.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
!                 IPARM(8) IS AMOUNT USED.      
!          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  SSORSI       
!                 NEEDS THIS TO BE IN LENGTH AT LEAST  5*N
!          IPARM  integer VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
!                 SOME integer PARAMETERS WHICH AFFECT THE METHOD.  IF
!          RPARM  D.P. VECTOR OF LENGTH 12. ALLOWS USER TO SPECIFY SOME 
!                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
!          IER    OUTPUT integer.  ERROR FLAG. (= IERR)   
!       
! ... SSORSI SUBPROGRAM REFERENCES:   
!       
!          FROM ITPACK    BISRCH, CHEBY, CHGSI, DFAULT, ECHALL, ECHOUT, 
!                         ITERM, TIMER, ITSRSI, IVFILL, OMEG,       
!                         OMGSTR, PARSI, PBETA, PERMAT, c_PERROR,     
!                         PERVEC, PFSOR, PMULT, PRBNDX, PSSOR1,     
!                         PSTOP, PVTBV, QSORT, SBELM, SCAL, DCOPY2,  
!                         DDOT2, SUM3, TSTCHG, UNSCAL, VEVPW, VFILL, 
!                         VOUT, WEVMW 
!          SYSTEM         DABS, DLOG, DLOG10, DBLE(AMAX0), DMAX1, DBLE(F
!                         MOD, DSQRT  
!       
!     VERSION:  ITPACK 2C (MARCH 1982)
!       
!     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
!                       CENTER FOR NUMERICAL ANALYSIS     
!                       UNIVERSITY OF TEXAS     
!                       AUSTIN, TX  78712       
!                       (512) 471-1242
!       
!     FOR ADDITIONAL DETAILS ON THE   
!          (A) subroutine SEE TOMS ARTICLE 1982 
!          (B) ALGORITHM  SEE CNA REPORT 150    
!       
!     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
!       
!     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
!                          L. HAGEMAN, D. YOUNG 
!                          ACADEMIC PRESS, 1981 
!       
!     **************************************************  
!     *               IMPORTANT NOTE                   *  
!     *                                                *  
!     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
!     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
!     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
!     *                                                *  
!     *   DRELPR      MACHINE RELATIVE PRECISION       *  
!     *   RPARM(1)    STOPPING CRITERION               *  
!     *                                                *  
!     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
!     *   SECOND USED IN TIMER                         *  
!     *                                                *  
!     **************************************************  
!       
!     SPECIFICATIONS FOR ARGUMENTS    
!       
      integer IA(1),JA(1),IWKSP(1),IPARM(12),NN,NW,IERR   
      double precision A(1),RHS(NN),U(NN),WKSP(NW),RPARM(12)
!       
!     SPECIFICATIONS FOR LOCAL VARIABLES
!       
      integer IB1,IB2,IB3,IB4,IB5,IDGTS,IER,IERPER,ITMAX1,LOOP,N,NB,N3
      double precision BETNEW,DIGIT1,DIGIT2,PBETA,TEMP,TIME1,TIME2,TOL
!       
! *** BEGIN: ITPACK common  
!       
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
!       
! *** END  : ITPACK common  
!       
! ... VARIABLES IN common BLOCK - ITCOM1
!       
!     IN     - ITERATION NUMBER       
!     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE format SWITCH 
!     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
!     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
!     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
!     NOUT   - OUTPUT element number     
!       
! ... VARIABLES IN common BLOCK - ITCOM2
!       
!     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
!     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
!     CASEII - ADAPTIVE PROCEDURE case SWITCH   
!     HALT   - STOPPING TEST SWITCH   
!     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
!       
! ... VARIABLES IN common BLOCK - ITCOM3
!       
!     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
!     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
!     CME    - ESTIMATE OF LARGEST EIGENVALUE   
!     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
!     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
!     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
!     GAMMA  - ACCELERATION parameter 
!     OMEGA  - OVERRELAXATION parameter FOR SOR AND SSOR  
!     QA     - PSEUDO-RESIDUAL RATIO  
!     QT     - VIRTUAL SPECTRAL RADIUS
!     RHO    - ACCELERATION parameter 
!     RRR    - ADAPTIVE parameter     
!     SIGE   - parameter SIGMA-SUB-E  
!     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
!     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
!     DRELPR - MACHINE RELATIVE PRECISION       
!     STPTST - STOPPING parameter     
!     UDNM   - TWO NORM OF U
!     ZETA   - STOPPING CRITERION     
!       
! ... INITIALIZE common BLOCKS
!       
      LEVEL = IPARM(2)      
      NOUT = IPARM(4)       
      if (IPARM(9).ge.0) IPARM(6) = 2 
      if (LEVEL.ge.1) write (NOUT,10) 
   10 format ('0'///1X,'BEGINNING OF ITPACK SOLUTION module  SSORSI') 
      IER = 0     
      if (IPARM(1).le.0) return       
      N = NN      
      if (IPARM(11).eq.0) TIMJ1 = TIMER(DUMMY)  
      if (LEVEL.ge.3) GO TO 20
      call ECHOUT (IPARM,RPARM,5)     
      GO TO 30    
   20 call ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      if (ZETA.ge.TEMP) GO TO 50      
      if (LEVEL.ge.1) write (NOUT,40) ZETA,DRELPR,TEMP    
40 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE SSORSI'/' ','    RPARM(1) =',D10.3, &
' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/ ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ', &
'    ZETA RESET TO ',D10.3)
      ZETA = TEMP 
   50 continue    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
!       
! ... VERIFY N    
!       
      if (N.gt.0) GO TO 70  
      IER = 51    
      if (LEVEL.ge.0) write (NOUT,60) N 
60 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SSORSI '/' ', &
'    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 380   
   70 continue    
!       
! ... REMOVE ROWS AND COLUMNS IF REQUESTED      
!       
      if (IPARM(10).eq.0) GO TO 90    
      TOL = RPARM(8)
      call IVFILL (N,IWKSP,0) 
      call VFILL (N,WKSP,0.0D0)       
      call SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      if (IER.eq.0) GO TO 90
      if (LEVEL.ge.0) write (NOUT,80) IER,TOL   
80 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SSORSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  SBELM '/' ', '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ', &
'    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X, ' RPARM(8) = ',D10.3,' (TOL)')
      GO TO 380   
!       
! ... INITIALIZE WKSP BASE ADDRESSES. 
!       
   90 IB1 = 1     
      IB2 = IB1+N 
      IB3 = IB2+N 
      IB4 = IB3+N 
      IB5 = IB4+N 
      IPARM(8) = 5*N
      if (NW.ge.IPARM(8)) GO TO 110   
      IER = 52    
      if (LEVEL.ge.0) write (NOUT,100) NW,IPARM(8)
100 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SSORSI '/' ', &
'    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10 ,' (NW)')
!       
! ... PERMUTE TO  RED-BLACK SYSTEM IF REQUESTED 
!       
  110 NB = IPARM(9) 
      if (NB.lt.0) GO TO 170
      N3 = 3*N    
      call IVFILL (N3,IWKSP,0)
      call PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      if (IER.eq.0) GO TO 130 
      if (LEVEL.ge.0) write (NOUT,120) IER,NB   
120 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SSORSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ', '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5 &
,' IPARM(9) = ',I5,' (NB)')
      GO TO 380   
!       
! ... PERMUTE MATRIX AND RHS
!       
  130 if (LEVEL.ge.2) write (NOUT,140) NB       
  140 format (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      call PERMAT (N,IA,JA,A,IWKSP,IWKSP(IB3),ISYM,LEVEL,NOUT,IER)  
      if (IER.eq.0) GO TO 160 
      if (LEVEL.ge.0) write (NOUT,150) IER      
150 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SSORSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ', '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 380   
  160 call PERVEC (N,RHS,IWKSP)       
      call PERVEC (N,U,IWKSP) 
!       
! ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
! ... DIAGONAL ELEMENTS.    
!       
  170 continue    
      call VFILL (IPARM(8),WKSP,0.0D0)
      call SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)     
      if (IER.eq.0) GO TO 190 
      if (LEVEL.ge.0) write (NOUT,180) IER      
180 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SSORSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ', '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)
      GO TO 380   
  190 if (LEVEL.le.2) GO TO 210       
      write (NOUT,200)      
200 format (///1X,'IN THE FOLLOWING, RHO AND GAMMA ARE', ' ACCELERATION PARAMETERS')
  210 if (IPARM(11).ne.0) GO TO 220   
      TIMI1 = TIMER(DUMMY)  
!       
! ... SPECIAL PROCEDURE FOR FULLY ADAPTIVE case.
!       
  220 continue    
      if (.not.ADAPT) GO TO 240       
      if (.not.BETADT) GO TO 230      
      call VFILL (N,WKSP(IB1),1.D0)   
BETNEW = PBETA(N,IA,JA,A,WKSP(IB1),WKSP(IB2),WKSP(IB3))/ DBLE(FLOAT(N))
      BETAB = DMAX1(BETAB,.25D0,BETNEW) 
  230 call OMEG (0.D0,1)    
      IS = 0      
!       
! ... ITERATION SEQUENCE    
!       
  240 ITMAX1 = ITMAX+1      
      do 260 LOOP = 1,ITMAX1
         IN = LOOP-1
         if (MOD(IN,2).eq.1) GO TO 250
!       
! ... CODE FOR THE EVEN ITERATIONS.   
!       
!     U           = U(IN)   
!     WKSP(IB1)   = U(IN-1) 
!       
call ITSRSI (N,IA,JA,A,RHS,U,WKSP(IB1),WKSP(IB2),WKSP(IB3), WKSP(IB4),WKSP(IB5))
!       
         if (HALT) GO TO 290
         GO TO 260
!       
! ... CODE FOR THE ODD ITERATIONS.    
!       
!     U           = U(IN-1) 
!     WKSP(IB1)   = U(IN)   
!       
250    call ITSRSI (N,IA,JA,A,RHS,WKSP(IB1),U,WKSP(IB2),WKSP(IB3), WKSP(IB4),WKSP(IB5))
!       
         if (HALT) GO TO 290
  260 continue    
!       
! ... ITMAX HAS BEEN REACHED
!       
      if (IPARM(11).ne.0) GO TO 270   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  270 if (LEVEL.ge.1) write (NOUT,280) ITMAX    
280 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE SSORSI'/' ','    FAILURE TO CONVERGE IN ' &
,I5,' ITERATIONS')
      IER = 53    
      if (IPARM(3).eq.0) RPARM(1) = STPTST      
      GO TO 320   
!       
! ... METHOD HAS CONVERGED  
!       
  290 if (IPARM(11).ne.0) GO TO 300   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  300 if (LEVEL.ge.1) write (NOUT,310) IN       
  310 format (/1X,'SSORSI  HAS CONVERGED IN ',I5,' ITERATIONS')     
!       
! ... PUT SOLUTION INTO U IF NOT ALREADY THERE. 
!       
  320 continue    
      if (MOD(IN,2).eq.1) call DCOPY2 (N,WKSP(IB1),1,U,1)  
!       
! ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
!       
      call UNSCAL (N,IA,JA,A,RHS,U,WKSP)
!       
! ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
!       
      if (IPARM(9).lt.0) GO TO 350    
call PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(IB3),ISYM,LEVEL,NOUT, IERPER)
      if (IERPER.eq.0) GO TO 340      
      if (LEVEL.ge.0) write (NOUT,330) IERPER   
330 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE SSORSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ', '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ', '    IER = ',I5)
      if (IER.eq.0) IER = IERPER      
      GO TO 380   
  340 call PERVEC (N,RHS,IWKSP(IB2))  
      call PERVEC (N,U,IWKSP(IB2))    
!       
! ... optional ERROR ANALYSIS 
!       
  350 IDGTS = IPARM(12)     
      if (IDGTS.lt.0) GO TO 360       
      if (IPARM(2).le.0) IDGTS = 0    
      call c_PERROR (N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)
!       
! ... SET return PARAMETERS IN IPARM AND RPARM  
!       
  360 if (IPARM(11).ne.0) GO TO 370   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  370 if (IPARM(3).ne.0) GO TO 380    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(5) = OMEGA      
      RPARM(6) = SPECR      
      RPARM(7) = BETAB      
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
!       
  380 continue    
      IERR = IER  
      if (LEVEL.ge.3) call ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
!       
      return      
      END 
      subroutine RSCG (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,IERR) 
!       
!     ITPACK 2C MAIN subroutine  RSCG  (REDUCED SYSTEM CONJUGATE    
!                                       GRADIENT) 
!     EACH OF THE MAIN SUBROUTINES:   
!           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
!     CAN BE USED INDEPENDENTLY OF THE OTHERS   
!       
! ... function:   
!       
!          THIS subroutine, RSCG, DRIVES THE  REDUCED SYSTEM CG     
!          ALGORITHM.       
!       
! ... parameter LIST:       
!       
!          N     INPUT integer.  dimension OF THE MATRIX. (= NN)    
!                 IN THE RED-BLACK MATRIX.      
!          IA,JA  INPUT integer VECTORS.  THE TWO integer ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  contains THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U contains THE 
!                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT contains 
!                 THE LATEST ESTIMATE TO THE SOLUTION.    
!          IWKSP  integer VECTOR WORKSPACE OF LENGTH 3*N  
!          NW     INPUT integer.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
!                 IPARM(8) IS AMOUNT USED.      
!          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  RSCG NEEDS   
!                 THIS TO BE IN LENGTH AT LEAST 
!                 N+3*NB+2*ITMAX, IF IPARM(5)=0  (SYMMETRIC STORAGE)
!                 N+3*NB+4*ITMAX, IF IPARM(5)=1  (NONSYMMETRIC STORAGE) 
!                 HERE NB IS THE ORDER OF THE BLACK SUBSYSTEM       
!          IPARM  integer VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
!                 SOME integer PARAMETERS WHICH AFFECT THE METHOD.  IF
!          RPARM  D.P. VECTOR OF LENGTH 12. ALLOWS USER TO SPECIFY SOME 
!                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
!          IER    OUTPUT integer. ERROR FLAG. (= IERR)    
!       
! ... RSCG SUBPROGRAM REFERENCES:     
!       
!          FROM ITPACK    BISRCH, CHGCON, DETERM, DFAULT, ECHALL,   
!                         ECHOUT, EIGVNS, EIGVSS, EQRT1S, ITERM, TIMER
!                         ITRSCG, IVFILL, PARCON, PERMAT, 
!                         c_PERROR, PERVEC, PMULT, PRBNDX, PRSBLK,    
!                         PRSRED, PSTOP, QSORT, SBELM, SCAL, DCOPY2, 
!                         DDOT2, SUM3, UNSCAL, VFILL, VOUT, WEVMW,   
!                         ZBRENT      
!          SYSTEM         DABS, DLOG10, DBLE(AMAX0), DMAX1, MOD, DSQRT
!       
!     VERSION:  ITPACK 2C (MARCH 1982)
!       
!     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
!                       CENTER FOR NUMERICAL ANALYSIS     
!                       UNIVERSITY OF TEXAS     
!                       AUSTIN, TX  78712       
!                       (512) 471-1242
!       
!     FOR ADDITIONAL DETAILS ON THE   
!          (A) subroutine SEE TOMS ARTICLE 1982 
!          (B) ALGORITHM  SEE CNA REPORT 150    
!       
!     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
!       
!     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
!                          L. HAGEMAN, D. YOUNG 
!                          ACADEMIC PRESS, 1981 
!       
!     **************************************************  
!     *               IMPORTANT NOTE                   *  
!     *                                                *  
!     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
!     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
!     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
!     *                                                *  
!     *   DRELPR      MACHINE RELATIVE PRECISION       *  
!     *   RPARM(1)    STOPPING CRITERION               *  
!     *                                                *  
!     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
!     *   SECOND USED IN TIMER                         *  
!     *                                                *  
!     **************************************************  
!       
!     SPECIFICATIONS FOR ARGUMENTS    
!       
      integer IA(1),JA(1),IWKSP(1),IPARM(12),NN,NW,IERR   
      double precision A(1),RHS(NN),U(NN),WKSP(NW),RPARM(12)
!       
!     SPECIFICATIONS FOR LOCAL VARIABLES
!       
integer IB1,IB2,IB3,IB4,IB5,IDGTS,IER,IERPER,ITMAX1,JB3,LOOP,N,NB, NR,NRP1,N3
      double precision DIGIT1,DIGIT2,TEMP,TIME1,TIME2,TOL 
!       
! *** BEGIN: ITPACK common  
!       
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
!       
! *** END  : ITPACK common  
!       
! ... VARIABLES IN common BLOCK - ITCOM1
!       
!     IN     - ITERATION NUMBER       
!     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
!     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE format SWITCH 
!     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
!     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
!     NOUT   - OUTPUT element number     
!       
! ... VARIABLES IN common BLOCK - ITCOM2
!       
!     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
!     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
!     CASEII - ADAPTIVE PROCEDURE case SWITCH   
!     HALT   - STOPPING TEST SWITCH   
!     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
!       
! ... VARIABLES IN common BLOCK - ITCOM3
!       
!     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
!     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
!     CME    - ESTIMATE OF LARGEST EIGENVALUE   
!     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
!     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
!     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
!     GAMMA  - ACCELERATION parameter 
!     OMEGA  - OVERRELAXATION parameter FOR SOR AND SSOR  
!     QA     - PSEUDO-RESIDUAL RATIO  
!     QT     - VIRTUAL SPECTRAL RADIUS
!     RHO    - ACCELERATION parameter 
!     RRR    - ADAPTIVE parameter     
!     SIGE   - parameter SIGMA-SUB-E  
!     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
!     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
!     DRELPR - MACHINE RELATIVE PRECISION       
!     STPTST - STOPPING parameter     
!     UDNM   - TWO NORM OF U
!     ZETA   - STOPPING CRITERION     
!       
! ... INITIALIZE common BLOCKS
!       
      LEVEL = IPARM(2)      
      NOUT = IPARM(4)       
      if (LEVEL.ge.1) write (NOUT,10) 
   10 format ('0'///1X,'BEGINNING OF ITPACK SOLUTION module  RSCG') 
      IER = 0     
      if (IPARM(1).le.0) return       
      N = NN      
      if (IPARM(11).eq.0) TIMJ1 = TIMER(DUMMY)  
      if (LEVEL.ge.3) GO TO 20
      call ECHOUT (IPARM,RPARM,6)     
      GO TO 30    
   20 call ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      if (ZETA.ge.TEMP) GO TO 50      
      if (LEVEL.ge.1) write (NOUT,40) ZETA,DRELPR,TEMP    
40 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE RSCG'/' ','    RPARM(1) =',D10.3, &
' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/ ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ', &
'    ZETA RESET TO ',D10.3)
      ZETA = TEMP 
   50 continue    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
!       
! ... VERIFY N    
!       
      if (N.gt.0) GO TO 70  
      IER = 61    
      if (LEVEL.ge.0) write (NOUT,60) N 
60 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE RSCG '/' ', &
'    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 430   
   70 continue    
!       
! ... REMOVE ROWS AND COLUMNS IF REQUESTED      
!       
      if (IPARM(10).eq.0) GO TO 90    
      TOL = RPARM(8)
      call IVFILL (N,IWKSP,0) 
      call VFILL (N,WKSP,0.0D0)       
      call SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      if (IER.eq.0) GO TO 90
      if (LEVEL.ge.0) write (NOUT,80) IER,TOL   
80 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE RSCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  SBELM '/' ', '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ', &
'    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X, ' RPARM(8) = ',D10.3,' (TOL)')
      GO TO 430   
!       
! ... INITIALIZE WKSP BASE ADDRESSES. 
!       
   90 IB1 = 1     
      IB2 = IB1+N 
      JB3 = IB2+N 
!       
! ... PERMUTE TO  RED-BLACK SYSTEM IF POSSIBLE  
!       
      NB = IPARM(9) 
      if (NB.ge.0) GO TO 110
      N3 = 3*N    
      call IVFILL (N3,IWKSP,0)
      call PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      if (IER.eq.0) GO TO 110 
      if (LEVEL.ge.0) write (NOUT,100) IER,NB   
100 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE RSCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ', '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5 &
,' IPARM(9) = ',I5,' (NB)')
      GO TO 430   
  110 if (NB.ge.0.and.NB.le.N) GO TO 130
      IER = 64    
      if (LEVEL.ge.1) write (NOUT,120) IER,NB   
120 format (/10X,'ERROR DETECTED IN RED-BLACK SUBSYSTEM INDEX'/10X, 'IER =',I5,' IPARM(9) =',I5,' (NB)')
      GO TO 430   
  130 if (NB.ne.0.and.NB.ne.N) GO TO 150
      NB = N/2    
      if (LEVEL.ge.2.and.IPARM(9).ge.0) write (NOUT,140) IPARM(9),NB
140 format (/10X,' IPARM(9) = ',I5,' IMPLIES MATRIX IS DIAGONAL'/10X, ' NB RESET TO ',I5)
!       
! ... PERMUTE MATRIX AND RHS
!       
  150 if (IPARM(9).ge.0) GO TO 190    
      if (LEVEL.ge.2) write (NOUT,160) NB       
  160 format (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      call PERMAT (N,IA,JA,A,IWKSP,IWKSP(JB3),ISYM,LEVEL,NOUT,IER)  
      if (IER.eq.0) GO TO 180 
      if (LEVEL.ge.0) write (NOUT,170) IER      
170 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE RSCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ', '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 430   
  180 call PERVEC (N,RHS,IWKSP)       
      call PERVEC (N,U,IWKSP) 
!       
! ... FINISH WKSP BASE ADDRESSES      
!       
  190 IB3 = IB2+NB
      IB4 = IB3+NB
      IB5 = IB4+NB
      NR = N-NB   
      NRP1 = NR+1 
      IPARM(8) = N+3*NB+2*ITMAX       
      if (ISYM.ne.0) IPARM(8) = IPARM(8)+2*ITMAX
      if (NW.ge.IPARM(8)) GO TO 210   
      IER = 62    
      if (LEVEL.ge.0) write (NOUT,200) NW,IPARM(8)
200 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE RSCG '/' ', &
'    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10 ,' (NW)')
      GO TO 430   
!       
! ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
! ... DIAGONAL ELEMENTS.    
!       
  210 continue    
      call VFILL (IPARM(8),WKSP,0.0D0)
      call SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)     
      if (IER.eq.0) GO TO 230 
      if (LEVEL.ge.0) write (NOUT,220) IER      
220 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE RSCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ', '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)
      GO TO 430   
  230 if (LEVEL.le.2) GO TO 260       
      write (NOUT,240)      
240 format (///1X,'IN THE FOLLOWING, RHO AND GAMMA ARE', ' ACCELERATION PARAMETERS')
      if (ADAPT) write (NOUT,250)     
250 format (1X,'CME IS THE ESTIMATE OF THE LARGEST EIGENVALUE OF', ' THE JACOBI MATRIX')
  260 if (IPARM(11).ne.0) GO TO 270   
      TIMI1 = TIMER(DUMMY)  
!       
! ... INITIALIZE FORWARD PSEUDO-RESIDUAL
!       
  270 continue    
      if (N.gt.1) GO TO 280 
      U(1) = RHS(1) 
      GO TO 330   
  280 call DCOPY2 (NR,RHS,1,WKSP(IB1),1) 
      call PRSRED (NB,NR,IA,JA,A,U(NRP1),WKSP(IB1))       
      call DCOPY2 (NB,RHS(NRP1),1,WKSP(IB2),1)   
      call PRSBLK (NB,NR,IA,JA,A,WKSP(IB1),WKSP(IB2))     
      call VEVMW (NB,WKSP(IB2),U(NRP1)) 
!       
! ... ITERATION SEQUENCE    
!       
      ITMAX1 = ITMAX+1      
      do 300 LOOP = 1,ITMAX1
         IN = LOOP-1
         if (MOD(IN,2).eq.1) GO TO 290
!       
! ... CODE FOR THE EVEN ITERATIONS.   
!       
!     U           = U(IN)       WKSP(IB2) = D(IN) 
!     WKSP(IB1)   = U(IN-1)     WKSP(IB3) = D(IN-1)       
!       
call ITRSCG (N,NB,IA,JA,A,U,WKSP(IB1),WKSP(IB2),WKSP(IB3), WKSP(IB4),WKSP(IB5))
!       
         if (HALT) GO TO 330
         GO TO 300
!       
! ... CODE FOR THE ODD ITERATIONS.    
!       
!     U           = U(IN-1)     WKSP(IB2) = D(IN-1)       
!     WKSP(IB1)   = U(IN)       WKSP(IB3) = D(IN) 
!       
290    call ITRSCG (N,NB,IA,JA,A,WKSP(IB1),U,WKSP(IB3),WKSP(IB2), WKSP(IB4),WKSP(IB5))
!       
         if (HALT) GO TO 330
  300 continue    
!       
! ... ITMAX HAS BEEN REACHED
!       
      if (IPARM(11).ne.0) GO TO 310   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  310 if (LEVEL.ge.1) write (NOUT,320) ITMAX    
320 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE RSCG'/' ','    FAILURE TO CONVERGE IN ', &
I5,' ITERATIONS')
      IER = 63    
      if (IPARM(3).eq.0) RPARM(1) = STPTST      
      GO TO 360   
!       
! ... METHOD HAS CONVERGED  
!       
  330 if (IPARM(11).ne.0) GO TO 340   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  340 if (LEVEL.ge.1) write (NOUT,350) IN       
  350 format (/1X,'RSCG  HAS CONVERGED IN ',I5,' ITERATIONS')       
!       
! ... PUT SOLUTION INTO U IF NOT ALREADY THERE. 
!       
  360 continue    
      if (N.eq.1) GO TO 370 
      if (MOD(IN,2).eq.1) call DCOPY2 (N,WKSP(IB1),1,U,1)  
      call DCOPY2 (NR,RHS,1,U,1)       
      call PRSRED (NB,NR,IA,JA,A,U(NRP1),U)     
!       
! ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
!       
  370 call UNSCAL (N,IA,JA,A,RHS,U,WKSP)
!       
! ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
!       
      if (IPARM(9).ge.0) GO TO 400    
call PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(JB3),ISYM,LEVEL,NOUT, IERPER)
      if (IERPER.eq.0) GO TO 390      
      if (LEVEL.ge.0) write (NOUT,380) IERPER   
380 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE RSCG '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ', '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ', '    IER = ',I5)
      if (IER.eq.0) IER = IERPER      
      GO TO 430   
  390 call PERVEC (N,RHS,IWKSP(IB2))  
      call PERVEC (N,U,IWKSP(IB2))    
!       
! ... optional ERROR ANALYSIS 
!       
  400 IDGTS = IPARM(12)     
      if (IDGTS.lt.0) GO TO 410       
      if (IPARM(2).le.0) IDGTS = 0    
      call c_PERROR (N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)
!       
! ... SET return PARAMETERS IN IPARM AND RPARM  
!       
  410 if (IPARM(11).ne.0) GO TO 420   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  420 IPARM(8) = IPARM(8)-2*(ITMAX-IN)
      if (ISYM.ne.0) IPARM(8) = IPARM(8)-2*(ITMAX-IN)     
      if (IPARM(3).ne.0) GO TO 430    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
!       
  430 continue    
      IERR = IER  
      if (LEVEL.ge.3) call ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
!       
      return      
      END 
      subroutine RSSI (NN,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,IERR) 
!       
!     ITPACK 2C MAIN subroutine  RSSI  (REDUCED SYSTEM SEMI-ITERATIVE)
!     EACH OF THE MAIN SUBROUTINES:   
!           JCG, JSI, SOR, SSORCG, SSORSI, RSCG, RSSI     
!     CAN BE USED INDEPENDENTLY OF THE OTHERS   
!       
! ... function:   
!       
!          THIS subroutine, RSSI, DRIVES THE  REDUCED SYSTEM SI     
!          ALGORITHM.       
!       
! ... parameter LIST:       
!       
!          N     INPUT integer.  dimension OF THE MATRIX. (= NN)    
!          IA,JA  INPUT integer VECTORS.  THE TWO integer ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  contains THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U contains THE 
!                 INITIAL GUESS TO THE SOLUTION. ON OUTPUT, IT contains 
!                 THE LATEST ESTIMATE TO THE SOLUTION.    
!          IWKSP  integer VECTOR WORKSPACE OF LENGTH 3*N  
!          NW     INPUT integer.  LENGTH OF AVAILABLE WKSP.  ON OUTPUT, 
!                 IPARM(8) IS AMOUNT USED.      
!          WKSP   D.P. VECTOR USED FOR WORKING SPACE.  RSSI 
!                 NEEDS THIS TO BE IN LENGTH AT LEAST  N + NB       
!                 HERE NB IS THE ORDER OF THE BLACK SUBSYSTEM       
!          IPARM  integer VECTOR OF LENGTH 12.  ALLOWS USER TO SPECIFY
!                 SOME integer PARAMETERS WHICH AFFECT THE METHOD.  IF
!          RPARM  D.P. VECTOR OF LENGTH 12. ALLOWS USER TO SPECIFY SOME 
!                 D.P. PARAMETERS WHICH AFFECT THE METHOD.
!          IER     OUTPUT integer.  ERROR FLAG. (= IERR)  
!       
! ... RSSI SUBPROGRAM REFERENCES:     
!       
!          FROM ITPACK    BISRCH, CHEBY, CHGSI, DFAULT, ECHALL,     
!                         ECHOUT, ITERM, TIMER, ITRSSI, IVFILL,     
!                         PARSI, PERMAT, c_PERROR, PERVEC, PMULT,     
!                         PRBNDX, PRSBLK, PRSRED, PSTOP, QSORT,     
!                         DAXPY2, SBELM, SCAL, DCOPY2, DDOT2, SUM3,    
!                         TSTCHG, UNSCAL, VEVMW, VFILL, VOUT,       
!                         WEVMW       
!          SYSTEM         DABS, DLOG10, DBLE(AMAX0), DMAX1, DBLE(FLOAT),
!                         DSQRT       
!       
!     VERSION:  ITPACK 2C (MARCH 1982)
!       
!     CODE WRITTEN BY:  DAVID KINCAID, ROGER GRIMES, JOHN RESPESS   
!                       CENTER FOR NUMERICAL ANALYSIS     
!                       UNIVERSITY OF TEXAS     
!                       AUSTIN, TX  78712       
!                       (512) 471-1242
!       
!     FOR ADDITIONAL DETAILS ON THE   
!          (A) subroutine SEE TOMS ARTICLE 1982 
!          (B) ALGORITHM  SEE CNA REPORT 150    
!       
!     BASED ON THEORY BY:  DAVID YOUNG, DAVID KINCAID, LOU HAGEMAN  
!       
!     REFERENCE THE BOOK:  APPLIED ITERATIVE METHODS      
!                          L. HAGEMAN, D. YOUNG 
!                          ACADEMIC PRESS, 1981 
!       
!     **************************************************  
!     *               IMPORTANT NOTE                   *  
!     *                                                *  
!     *      WHEN INSTALLING ITPACK ROUTINES ON A      *  
!     *  DIFFERENT COMPUTER, RESET SOME OF THE VALUES  *  
!     *  IN  SUBROUTNE DFAULT.   MOST IMPORTANT ARE    *  
!     *                                                *  
!     *   DRELPR      MACHINE RELATIVE PRECISION       *  
!     *   RPARM(1)    STOPPING CRITERION               *  
!     *                                                *  
!     *   ALSO CHANGE SYSTEM-DEPENDENT ROUTINE         *  
!     *   SECOND USED IN TIMER                         *  
!     *                                                *  
!     **************************************************  
!       
!     SPECIFICATIONS FOR ARGUMENTS    
!       
      integer IA(1),JA(1),IWKSP(1),IPARM(12),NN,NW,IERR   
      double precision A(1),RHS(NN),U(NN),WKSP(NW),RPARM(12)
!       
!     SPECIFICATIONS FOR LOCAL VARIABLES
!       
      integer IB1,IB2,IDGTS,IER,IERPER,ITMAX1,JB3,LOOP,N,NB,NR,NRP1,N3
      double precision DIGIT1,DIGIT2,TEMP,TIME1,TIME2,TOL 
!       
! *** BEGIN: ITPACK common  
!       
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
!       
! *** END  : ITPACK common  
!       
! ... VARIABLES IN common BLOCK - ITCOM1
!       
!     IN     - ITERATION NUMBER       
!     IS     - ITERATION NUMBER WHEN PARAMETERS LAST CHANGED
!     ISYM   - SYMMETRIC/NONSYMMETRIC STORAGE format SWITCH 
!     ITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED       
!     LEVEL  - LEVEL OF OUTPUT CONTROL SWITCH   
!     NOUT   - OUTPUT element number     
!       
! ... VARIABLES IN common BLOCK - ITCOM2
!       
!     ADAPT  - FULLY ADAPTIVE PROCEDURE SWITCH  
!     BETADT - SWITCH FOR ADAPTIVE DETERMINATION OF BETA  
!     CASEII - ADAPTIVE PROCEDURE case SWITCH   
!     HALT   - STOPPING TEST SWITCH   
!     PARTAD - PARTIALLY ADAPTIVE PROCEDURE SWITCH
!       
! ... VARIABLES IN common BLOCK - ITCOM3
!       
!     BDELNM - TWO NORM OF B TIMES DELTA-SUPER-N
!     BETAB  - ESTIMATE FOR THE SPECTRAL RADIUS OF LU MATRIX
!     CME    - ESTIMATE OF LARGEST EIGENVALUE   
!     DELNNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION N      
!     DELSNM - INNER PRODUCT OF PSEUDO-RESIDUAL AT ITERATION S      
!     FF     - ADAPTIVE PROCEDURE DAMPING FACTOR
!     GAMMA  - ACCELERATION parameter 
!     OMEGA  - OVERRELAXATION parameter FOR SOR AND SSOR  
!     QA     - PSEUDO-RESIDUAL RATIO  
!     QT     - VIRTUAL SPECTRAL RADIUS
!     RHO    - ACCELERATION parameter 
!     RRR    - ADAPTIVE parameter     
!     SIGE   - parameter SIGMA-SUB-E  
!     SME    - ESTIMATE OF SMALLEST EIGENVALUE  
!     SPECR  - SPECTRAL RADIUS ESTIMATE FOR SSOR
!     DRELPR - MACHINE RELATIVE PRECISION       
!     STPTST - STOPPING parameter     
!     UDNM   - TWO NORM OF U
!     ZETA   - STOPPING CRITERION     
!       
! ... INITIALIZE common BLOCKS
!       
      LEVEL = IPARM(2)      
      NOUT = IPARM(4)       
      if (LEVEL.ge.1) write (NOUT,10) 
   10 format ('0'///1X,'BEGINNING OF ITPACK SOLUTION module  RSSI') 
      IER = 0     
      if (IPARM(1).le.0) return       
      N = NN      
      if (IPARM(11).eq.0) TIMJ1 = TIMER(DUMMY)  
      if (LEVEL.ge.3) GO TO 20
      call ECHOUT (IPARM,RPARM,7)     
      GO TO 30    
   20 call ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,1) 
   30 TEMP = 5.0D2*DRELPR   
      if (ZETA.ge.TEMP) GO TO 50      
      if (LEVEL.ge.1) write (NOUT,40) ZETA,DRELPR,TEMP    
40 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE RSSI'/' ','    RPARM(1) =',D10.3, &
' (ZETA)'/' ','    A VALUE THIS SMALL MAY HINDER CONVERGENCE '/ ' ','    SINCE MACHINE PRECISION DRELPR =',D10.3/' ', &
'    ZETA RESET TO ',D10.3)
      ZETA = TEMP 
   50 continue    
      TIME1 = RPARM(9)      
      TIME2 = RPARM(10)     
      DIGIT1 = RPARM(11)    
      DIGIT2 = RPARM(12)    
!       
! ... VERIFY N    
!       
      if (N.gt.0) GO TO 70  
      IER = 71    
      if (LEVEL.ge.0) write (NOUT,60) N 
60 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE RSSI '/' ', &
'    INVALID MATRIX DIMENSION, N =',I8)
      GO TO 420   
   70 continue    
!       
! ... REMOVE ROWS AND COLUMNS IF REQUESTED      
!       
      if (IPARM(10).eq.0) GO TO 90    
      TOL = RPARM(8)
      call IVFILL (N,IWKSP,0) 
      call VFILL (N,WKSP,0.0D0)       
      call SBELM (N,IA,JA,A,RHS,IWKSP,WKSP,TOL,ISYM,LEVEL,NOUT,IER) 
      if (IER.eq.0) GO TO 90
      if (LEVEL.ge.0) write (NOUT,80) IER,TOL   
80 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE RSSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  SBELM '/' ', '    WHICH REMOVES ROWS AND COLUMNS OF SYSTEM '/' ', &
'    WHEN DIAGONAL ENTRY TOO LARGE  '/' ','    IER = ',I5,5X, ' RPARM(8) = ',D10.3,' (TOL)')
!       
! ... INITIALIZE WKSP BASE ADDRESSES. 
!       
   90 IB1 = 1     
      IB2 = IB1+N 
      JB3 = IB2+N 
!       
! ... PERMUTE TO  RED-BLACK SYSTEM IF POSSIBLE  
!       
      NB = IPARM(9) 
      if (NB.ge.0) GO TO 110
      N3 = 3*N    
      call IVFILL (N3,IWKSP,0)
      call PRBNDX (N,NB,IA,JA,IWKSP,IWKSP(IB2),LEVEL,NOUT,IER)      
      if (IER.eq.0) GO TO 110 
      if (LEVEL.ge.0) write (NOUT,100) IER,NB   
100 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE RSSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PRBNDX'/' ', '    WHICH COMPUTES THE RED-BLACK INDEXING'/' ','    IER = ',I5 &
,' IPARM(9) = ',I5,' (NB)')
      GO TO 420   
  110 if (NB.ge.0.and.NB.le.N) GO TO 130
      IER = 74    
      if (LEVEL.ge.1) write (NOUT,120) IER,NB   
120 format (/10X,'ERROR DETECTED IN RED-BLACK SUBSYSTEM INDEX'/10X, 'IER =',I5,' IPARM(9) =',I5,' (NB)')
      GO TO 420   
  130 if (NB.ne.0.and.NB.ne.N) GO TO 150
      NB = N/2    
      if (LEVEL.ge.2.and.IPARM(9).ge.0) write (NOUT,140) IPARM(9),NB
140 format (/10X,' IPARM(9) = ',I5,' IMPLIES MATRIX IS DIAGONAL'/10X, ' NB RESET TO ',I5)
!       
! ... PERMUTE MATRIX AND RHS
!       
  150 if (IPARM(9).ge.0) GO TO 190    
      if (LEVEL.ge.2) write (NOUT,160) NB       
  160 format (/10X,'ORDER OF BLACK SUBSYSTEM = ',I5,' (NB)')
      call PERMAT (N,IA,JA,A,IWKSP,IWKSP(JB3),ISYM,LEVEL,NOUT,IER)  
      if (IER.eq.0) GO TO 180 
      if (LEVEL.ge.0) write (NOUT,170) IER      
170 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE RSSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ', '    WHICH DOES THE RED-BLACK PERMUTATION'/' ','    IER = ',I5)
      GO TO 420   
  180 call PERVEC (N,RHS,IWKSP)       
      call PERVEC (N,U,IWKSP) 
!       
! ... INITIALIZE WKSP BASE ADDRESSES  
!       
  190 NR = N-NB   
!       
      NRP1 = NR+1 
      IPARM(8) = N+NB       
      if (NW.ge.IPARM(8)) GO TO 210   
      IER = 72    
      if (LEVEL.ge.0) write (NOUT,200) NW,IPARM(8)
200 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE RSSI '/' ', &
'    NOT ENOUGH WORKSPACE AT ',I10/' ','    SET IPARM(8) =',I10 ,' (NW)')
      GO TO 420   
!       
! ... SCALE LINEAR SYSTEM, U, AND RHS BY THE SQUARE ROOT OF THE     
! ... DIAGONAL ELEMENTS.    
!       
  210 continue    
      call VFILL (IPARM(8),WKSP,0.0D0)
      call SCAL (N,IA,JA,A,RHS,U,WKSP,LEVEL,NOUT,IER)     
      if (IER.eq.0) GO TO 230 
      if (LEVEL.ge.0) write (NOUT,220) IER      
220 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE RSSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  SCAL  '/' ', '    WHICH SCALES THE SYSTEM   '/' ','    IER = ',I5)
      GO TO 420   
  230 if (LEVEL.le.2) GO TO 250       
      write (NOUT,240)      
240 format (///1X,'IN THE FOLLOWING, RHO AND GAMMA ARE', ' ACCELERATION PARAMETERS')
  250 if (IPARM(11).ne.0) GO TO 260   
      TIMI1 = TIMER(DUMMY)  
!       
! ... ITERATION SEQUENCE    
!       
  260 if (N.gt.1) GO TO 270 
      U(1) = RHS(1) 
      GO TO 320   
  270 ITMAX1 = ITMAX+1      
      do 290 LOOP = 1,ITMAX1
         IN = LOOP-1
         if (MOD(IN,2).eq.1) GO TO 280
!       
! ... CODE FOR THE EVEN ITERATIONS.   
!       
!     U           = U(IN)   
!     WKSP(IB1)   = U(IN-1) 
!       
         call ITRSSI (N,NB,IA,JA,A,RHS,U,WKSP(IB1),WKSP(IB2))       
!       
         if (HALT) GO TO 320
         GO TO 290
!       
! ... CODE FOR THE ODD ITERATIONS.    
!       
!     U           = U(IN-1) 
!     WKSP(IB1)   = U(IN)   
!       
  280    call ITRSSI (N,NB,IA,JA,A,RHS,WKSP(IB1),U,WKSP(IB2))       
!       
         if (HALT) GO TO 320
  290 continue    
!       
! ... ITMAX HAS BEEN REACHED
!       
      if (IPARM(11).ne.0) GO TO 300   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  300 if (LEVEL.ge.1) write (NOUT,310) ITMAX    
310 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE RSSI'/' ','    FAILURE TO CONVERGE IN', &
I5,' ITERATIONS')
      IER = 73    
      if (IPARM(3).eq.0) RPARM(1) = STPTST      
      GO TO 350   
!       
! ... METHOD HAS CONVERGED  
!       
  320 if (IPARM(11).ne.0) GO TO 330   
      TIMI2 = TIMER(DUMMY)  
      TIME1 = DBLE(TIMI2-TIMI1)       
  330 if (LEVEL.ge.1) write (NOUT,340) IN       
  340 format (/1X,'RSSI  HAS CONVERGED IN ',I5,' ITERATIONS')       
!       
! ... PUT SOLUTION INTO U IF NOT ALREADY THERE. 
!       
  350 continue    
      if (N.eq.1) GO TO 360 
      if (MOD(IN,2).eq.1) call DCOPY2 (N,WKSP(IB1),1,U,1)  
      call DCOPY2 (NR,RHS,1,U,1)       
      call PRSRED (NB,NR,IA,JA,A,U(NRP1),U)     
!       
! ... UNSCALE THE MATRIX, SOLUTION, AND RHS VECTORS.      
!       
  360 call UNSCAL (N,IA,JA,A,RHS,U,WKSP)
!       
! ... UN-PERMUTE MATRIX,RHS, AND SOLUTION       
!       
      if (IPARM(9).ge.0) GO TO 390    
call PERMAT (N,IA,JA,A,IWKSP(IB2),IWKSP(JB3),ISYM,LEVEL,NOUT, IERPER)
      if (IERPER.eq.0) GO TO 380      
      if (LEVEL.ge.0) write (NOUT,370) IERPER   
370 format ('0','*** F A T A L     E R R O R ************'/'0', '    CALLED FROM ITPACK ROUTINE RSSI '/' ', &
'    ERROR DETECTED IN SUBROUTINE  PERMAT'/' ', '    WHICH UNDOES THE RED-BLACK PERMUTATION   '/' ', '    IER = ',I5)
      if (IER.eq.0) IER = IERPER      
      GO TO 420   
  380 call PERVEC (N,RHS,IWKSP(IB2))  
      call PERVEC (N,U,IWKSP(IB2))    
!       
! ... optional ERROR ANALYSIS 
!       
  390 IDGTS = IPARM(12)     
      if (IDGTS.lt.0) GO TO 400       
      if (IPARM(2).le.0) IDGTS = 0    
      call c_PERROR (N,IA,JA,A,RHS,U,WKSP,DIGIT1,DIGIT2,IDGTS)
!       
! ... SET return PARAMETERS IN IPARM AND RPARM  
!       
  400 if (IPARM(11).ne.0) GO TO 410   
      TIMJ2 = TIMER(DUMMY)  
      TIME2 = DBLE(TIMJ2-TIMJ1)       
  410 if (IPARM(3).ne.0) GO TO 420    
      IPARM(1) = IN 
      IPARM(9) = NB 
      RPARM(2) = CME
      RPARM(3) = SME
      RPARM(9) = TIME1      
      RPARM(10) = TIME2     
      RPARM(11) = DIGIT1    
      RPARM(12) = DIGIT2    
!       
  420 continue    
      IERR = IER  
      if (LEVEL.ge.3) call ECHALL (N,IA,JA,A,RHS,IPARM,RPARM,2)     
!       
      return      
      END 
      subroutine ITJCG (NN,IA,JA,A,U,U1,D,D1,DTWD,TRI)    
!       
! ... function:   
!       
!          THIS subroutine, ITJCG, PERFORMS ONE ITERATION OF THE    
!          JACOBI CONJUGATE GRADIENT ALGORITHM.  IT IS CALLED BY JCG. 
!       
! ... parameter LIST:       
!       
!          N      INPUT integer.  dimension OF THE MATRIX. (= NN)   
!          IA,JA  INPUT integer VECTORS.  contains INFORMATION DEFINING 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR. contains THE NONZERO VALUES OF THE 
!                 LINEAR SYSTEM.      
!          U      INPUT D.P. VECTOR.  contains THE VALUE OF THE     
!                 SOLUTION VECTOR AT THE END OF IN ITERATIONS.      
!          U1     INPUT/OUTPUT D.P. VECTOR.  ON INPUT, IT contains  
!                 THE VALUE OF THE SOLUTION AT THE END OF THE IN-1  
!                 ITERATION.  ON OUTPUT, IT WILL CONTAIN THE NEWEST 
!                 ESTIMATE FOR THE SOLUTION VECTOR.       
!          D      INPUT D.P. VECTOR.  contains THE PSEUDO-RESIDUAL  
!                 VECTOR AFTER IN ITERATIONS.   
!          D1     INPUT/OUTPUT D.P. VECTOR.  ON INPUT, D1 contains  
!                 THE PSEUDO-RESIDUAL VECTOR AFTER IN-1 ITERATIONS.  ON 
!                 OUTPUT, IT WILL CONTAIN THE NEWEST PSEUDO-RESIDUAL
!                 VECTOR.   
!          DTWD   D.P. ARRAY.  USED IN THE COMPUTATIONS OF THE      
!                 ACCELERATION parameter GAMMA AND THE NEW PSEUDO-  
!                 RESIDUAL. 
!          TRI    D.P. ARRAY.  STORES THE TRIDIAGONAL MATRIX associated 
!                 WITH THE EIGENVALUES OF THE CONJUGATE GRADIENT    
!                 POLYNOMIAL. 
!       
! ... SPECIFICATIONS FOR ARGUMENTS    
!       
      integer IA(1),JA(1),NN
      double precision A(1),U(NN),U1(NN),D(NN),D1(NN),DTWD(NN),TRI(2,1) 
!       
! ... SPECIFICATIONS FOR LOCAL VARIABLES
!       
      integer :: N
      double precision CON,C1,C2,C3,C4,DNRM,DTNRM,GAMOLD,RHOOLD,RHOTMP
      logical :: Q1
!       
! ... SPECIFICATIONS FOR function SUBPROGRAMS   
!       
      double precision DDOT2 
!       
! *** BEGIN: ITPACK common  
!       
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
!       
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
!       
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
!       
! *** END  : ITPACK common  
!       
!     DESCRIPTION OF VARIABLES IN common BLOCKS IN subroutine JCG   
!       
! ... COMPUTE NEW ESTIMATE FOR CME IF ADAPT = .true.      
!       
      if (ADAPT) call CHGCON (TRI,GAMOLD,RHOOLD,1)
!       
! ... TEST FOR STOPPING     
!       
      N = NN      
      DELNNM = DDOT2(N,D,1,D,1)
      DNRM = DELNNM 
      CON = CME   
      call PSTOP (N,U,DNRM,CON,1,Q1)  
      if (HALT) GO TO 30    
!       
! ... COMPUTE RHO AND GAMMA - ACCELERATION PARAMETERS     
!       
      call VFILL (N,DTWD,0.D0)
      call PJAC2 (N,IA,JA,A,D,DTWD)    
      DTNRM = DDOT2(N,D,1,DTWD,1)      
      if (ISYM.eq.0) GO TO 10 
      RHOTMP = DDOT2(N,DTWD,1,D1,1)    
      call PARCON (DTNRM,C1,C2,C3,C4,GAMOLD,RHOTMP,1)     
      RHOOLD = RHOTMP       
      GO TO 20    
   10 call PARCON (DTNRM,C1,C2,C3,C4,GAMOLD,RHOOLD,1)     
!       
! ... COMPUTE U(IN+1) AND D(IN+1)     
!       
   20 call SUM3 (N,C1,D,C2,U,C3,U1)   
      call SUM3 (N,C1,DTWD,C4,D,C3,D1)
!       
! ... OUTPUT INTERMEDIATE INFORMATION 
!       
   30 call ITERM (N,A,U,DTWD,1)       
!       
      return      
      END 
      subroutine ITJSI (NN,IA,JA,A,RHS,U,U1,D,ICNT)       
!       
! ... function:   
!       
!          THIS subroutine, ITJSI, PERFORMS ONE ITERATION OF THE    
!          JACOBI SEMI-ITERATIVE ALGORITHM.  IT IS CALLED BY JSI.   
!       
! ... parameter LIST:       
!       
!          N      INPUT integer.  dimension OF THE MATRIX. (= NN)   
!          IA,JA  INPUT integer VECTORS.  THE TWO integer ARRAYS OF 
!                 THE SPARSE MATRIX REPRESENTATION.       
!          A      INPUT D.P. VECTOR.  THE D.P. ARRAY OF THE SPARSE  
!                 MATRIX REPRESENTATION.
!          RHS    INPUT D.P. VECTOR.  contains THE RIGHT HAND SIDE  
!                 OF THE MATRIX PROBLEM.
!          U      INPUT D.P. VECTOR.  contains THE ESTIMATE FOR THE 
!                 SOLUTION VECTOR AFTER IN ITERATIONS.    
!          U1     INPUT/OUTPUT D.P. VECTOR.  ON INPUT, U1 contains THE
!                 SOLUTION VECTOR AFTER IN-1 ITERATIONS.  ON OUTPUT,
!                 IT WILL CONTAIN THE NEWEST ESTIMATE FOR THE SOLUTION
!                 VECTOR.   
!          D      D.P. ARRAY.  D IS USED FOR THE COMPUTATION OF THE 
!                 PSEUDO-RESIDUAL ARRAY FOR THE CURRENT ITERATION.  
!          ICNT   NUMBER OF ITERATIONS SINCE LAST CHANGE OF SME     
!       
! ... SPECIFICATIONS OF ARGUMENTS     
!       
      !double precision A(1),RHS(NN),U(NN),U1(NN),D(NN)    
      
      !BUGFIX2024091302.
      use Global_ITPACK
      integer IA(NN+1),JA(ITPACK_NNZ),NN,ICNT     
      double precision A(ITPACK_NNZ),RHS(NN),U(NN),U1(NN),D(NN)    
      
      integer :: N
      double precision CON,C1,C2,C3,DNRM,DTNRM,OLDNRM     
      logical :: Q1
      double precision DDOT2,PVTBV     
      logical TSTCHG,CHGSME 
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      N = NN      
      if (IN.eq.0) ICNT = 0 
      call DCOPY2 (N,RHS,1,D,1)
      call PJAC2 (N,IA,JA,A,U,D)       
      call VEVMW (N,D,U)    
      OLDNRM = DELNNM       
      DELNNM = DDOT2(N,D,1,D,1)
      DNRM = DELNNM 
      CON = CME   
      call PSTOP (N,U,DNRM,CON,1,Q1)  
      if (HALT) GO TO 40    
      if (.not.ADAPT) GO TO 30
      if (.not.TSTCHG(1)) GO TO 10    
      DTNRM = PVTBV(N,IA,JA,A,D)      
      call CHGSI (DTNRM,1)  
      if (.not.ADAPT) GO TO 30
      GO TO 20    
   10 continue    
      if (CASEII) GO TO 30  
      if (.not.CHGSME(OLDNRM,ICNT)) GO TO 30    
      ICNT = 0    
   20 call DCOPY2 (N,U,1,U1,1) 
      call DAXPY2 (N,GAMMA,D,1,U1,1)   
      GO TO 40    
   30 call PARSI (C1,C2,C3,1) 
      call SUM3 (N,C1,D,C2,U,C3,U1)   
   40 call ITERM (N,A,U,D,2)
      return      
      END 
      subroutine ITSOR (NN,IA,JA,A,RHS,U,WK)    
      integer IA(1),JA(1),NN
      double precision A(1),RHS(NN),U(NN),WK(NN)
      integer IP,IPHAT,IPSTAR,ISS,N   
      double precision DNRM,H,OMEGAP,SPCRM1     
      logical CHANGE,Q1     
      double precision TAU  
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      N = NN      
      if (IN.ne.0) GO TO 20 
      call PSTOP (N,U,0.D0,0.D0,0,Q1) 
      if (ADAPT) GO TO 10   
      CHANGE = .false.      
      IP = 0      
      IPHAT = 2   
      ISS = 0     
      GO TO 30    
   10 CHANGE = .true.       
      IP = 0      
      OMEGAP = OMEGA
      OMEGA = 1.D0
      ISS = 0     
      IPHAT = 2   
      IPSTAR = 4  
      if (OMEGAP.le.1.D0) CHANGE = .false.      
   20 if (.not.CHANGE) GO TO 30       
      CHANGE = .false.      
      IS = IS+1   
      IP = 0      
      ISS = 0     
      OMEGA = DMIN1(OMEGAP,TAU(IS))   
      IPHAT = MAX0(3,IFIX(SNGL((OMEGA-1.D0)/(2.D0-OMEGA)))) 
      IPSTAR = IPSTR(OMEGA) 
   30 continue    
      DELSNM = DELNNM       
      SPCRM1 = SPECR
      call DCOPY2 (N,RHS,1,WK,1)       
      call PFSOR1 (N,IA,JA,A,U,WK)    
      if (DELNNM.eq.0.D0) GO TO 40    
      if (IN.ne.0) SPECR = DELNNM/DELSNM
      if (IP.lt.IPHAT) GO TO 70       
      if (SPECR.ge.1.D0) GO TO 70     
      if (.not.(SPECR.gt.(OMEGA-1.D0))) GO TO 40
      H = SPECR   
      GO TO 50    
   40 ISS = ISS+1 
      H = OMEGA-1.D0
   50 continue    
      DNRM = DELNNM**2      
      call PSTOP (N,U,DNRM,H,1,Q1)    
      if (HALT) GO TO 70    
      if (.not.ADAPT) GO TO 70
      if (IP.lt.IPSTAR) GO TO 70      
      if (OMEGA.gt.1.D0) GO TO 60     
      CME = DSQRT(DABS(SPECR))
      OMEGAP = 2.D0/(1.D0+DSQRT(DABS(1.D0-SPECR)))
      CHANGE = .true.       
      GO TO 70    
   60 if (ISS.ne.0) GO TO 70
      if (SPECR.le.(OMEGA-1.D0)**FF) GO TO 70   
      if ((SPECR+5.D-5).le.SPCRM1) GO TO 70     
      CME = (SPECR+OMEGA-1.D0)/(DSQRT(DABS(SPECR))*OMEGA) 
      OMEGAP = 2.D0/(1.D0+DSQRT(DABS(1.D0-CME*CME)))      
      CHANGE = .true.       
   70 call ITERM (N,A,U,WK,3) 
      IP = IP+1   
      return      
      END 
      subroutine ITSRCG (NN,IA,JA,A,RHS,U,U1,C,C1,D,DL,WK,TRI)      
      integer IA(1),JA(1),NN
double precision A(1),RHS(NN),U(NN),U1(NN),C(NN),C1(NN),D(NN), DL(NN),WK(NN),TRI(2,1)
      integer :: N
      double precision BETNEW,CON,DNRM,GAMOLD,RHOOLD,RHOTMP,T1,T2,T3,T4 
      logical :: Q1
      double precision DDOT2,PBETA,PVTBV 
      logical OMGCHG,OMGSTR 
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      N = NN      
      if (ADAPT.or.PARTAD) call CHGCON (TRI,GAMOLD,RHOOLD,3)
      call DCOPY2 (N,RHS,1,WK,1)       
      call DCOPY2 (N,C,1,D,1)
      call VEVPW (N,D,U)    
      call PBSOR (N,IA,JA,A,D,WK)     
      call VEVMW (N,D,U)    
      call DCOPY2 (N,D,1,DL,1) 
      call VFILL (N,WK,0.D0)
      call PFSOR (N,IA,JA,A,DL,WK)    
      call WEVMW (N,D,DL)   
      DELNNM = DDOT2(N,C,1,C,1)
      if (DELNNM.eq.0.D0) GO TO 30    
      DNRM = DDOT2(N,C,1,DL,1) 
      if (DNRM.eq.0.D0) GO TO 30      
      if (ISYM.eq.0) GO TO 10 
      RHOTMP = DDOT2(N,C,1,C1,1)-DDOT2(N,DL,1,C1,1) 
      call PARCON (DNRM,T1,T2,T3,T4,GAMOLD,RHOTMP,3)      
      RHOOLD = RHOTMP       
      GO TO 20    
   10 call PARCON (DNRM,T1,T2,T3,T4,GAMOLD,RHOOLD,3)      
   20 call SUM3 (N,T1,D,T2,U,T3,U1)   
   30 BDELNM = DDOT2(N,D,1,D,1)
      DNRM = BDELNM 
      CON = SPECR 
      call PSTOP (N,U,DNRM,CON,1,Q1)  
      if (HALT) GO TO 100   
      if (ADAPT) GO TO 40   
      call SUM3 (N,-T1,DL,T2,C,T3,C1) 
      GO TO 100   
   40 continue    
      if (OMGSTR(1)) GO TO 90 
      if (OMGCHG(1)) GO TO 50 
      call SUM3 (N,-T1,DL,T2,C,T3,C1) 
      GO TO 100   
   50 continue    
      if (.not.BETADT) GO TO 60       
      BETNEW = PBETA(N,IA,JA,A,D,WK,C1)/BDELNM  
      BETAB = DMAX1(BETAB,.25D0,BETNEW) 
   60 continue    
      if (CASEII) GO TO 70  
      DNRM = PVTBV(N,IA,JA,A,D)       
      GO TO 80    
   70 call VFILL (N,WK,0.D0)
      call PJAC2 (N,IA,JA,A,D,WK)      
      DNRM = DDOT2(N,WK,1,WK,1)
   80 call OMEG (DNRM,3)    
   90 continue    
      call DCOPY2 (N,RHS,1,WK,1)       
      call DCOPY2 (N,U1,1,C1,1)
      call PFSOR (N,IA,JA,A,C1,WK)    
      call VEVMW (N,C1,U1)  
  100 call ITERM (N,A,U,WK,4) 
      return      
      END 
      subroutine ITSRSI (NN,IA,JA,A,RHS,U,U1,C,D,CTWD,WK) 
      integer IA(1),JA(1),NN
double precision A(1),RHS(NN),U(NN),U1(NN),C(NN),D(NN),WK(NN), CTWD(NN)
      integer :: N
      double precision BETNEW,CON,C1,C2,C3,DNRM 
      logical :: Q1
      double precision DDOT2,PBETA,PVTBV 
      logical OMGSTR,TSTCHG 
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      N = NN      
      call DCOPY2 (N,RHS,1,WK,1)       
      call DCOPY2 (N,U,1,CTWD,1)       
      call PSSOR1 (N,IA,JA,A,CTWD,WK,C,D)       
      call PARSI (C1,C2,C3,3) 
      call SUM3 (N,C1,D,C2,U,C3,U1)   
      BDELNM = DDOT2(N,D,1,D,1)
      DNRM = BDELNM 
      CON = SPECR 
      call PSTOP (N,U,DNRM,CON,1,Q1)  
      if (HALT.or..not.(ADAPT.or.PARTAD)) GO TO 40
      if (OMGSTR(1)) GO TO 40 
      DELNNM = DDOT2(N,C,1,C,1)
      if (IN.eq.IS) DELSNM = DELNNM   
      if (IN.eq.0.or..not.TSTCHG(1)) GO TO 40   
      call DCOPY2 (N,D,1,CTWD,1)       
      call VFILL (N,WK,0.D0)
      call PFSOR (N,IA,JA,A,CTWD,WK)  
      call VEVPW (N,CTWD,C) 
      call VEVMW (N,CTWD,D) 
      DNRM = DDOT2(N,C,1,CTWD,1)       
      call CHGSI (DNRM,3)   
      if (.not.ADAPT) GO TO 40
      if (.not.BETADT) GO TO 10       
      BETNEW = PBETA(N,IA,JA,A,D,WK,CTWD)/BDELNM
      BETAB = DMAX1(BETAB,.25D0,BETNEW) 
   10 continue    
      if (CASEII) GO TO 20  
      DNRM = PVTBV(N,IA,JA,A,D)       
      GO TO 30    
   20 call VFILL (N,WK,0.D0)
      call PJAC2 (N,IA,JA,A,D,WK)      
      DNRM = DDOT2(N,WK,1,WK,1)
   30 call OMEG (DNRM,3)    
   40 call ITERM (N,A,U,WK,5) 
      return      
      END 
      subroutine ITRSCG (N,NNB,IA,JA,A,UB,UB1,DB,DB1,WB,TRI)
      integer IA(1),JA(1),N,NNB       
      double precision A(1),UB(N),UB1(N),DB(NNB),DB1(N),WB(NNB),TRI(2,1)
      integer NB,NR,NRP1    
      double precision CON,C1,C2,C3,C4,DNRM,GAMOLD,RHOOLD,RHOTMP    
      logical :: Q1
      double precision DDOT2 
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      NB = NNB    
      NR = N-NB   
      NRP1 = NR+1 
      if (ADAPT) call CHGCON (TRI,GAMOLD,RHOOLD,2)
      DELNNM = DDOT2(NB,DB,1,DB,1)     
      DNRM = DELNNM 
      CON = CME   
      call PSTOP (NB,UB(NRP1),DNRM,CON,2,Q1)    
      if (HALT) GO TO 30    
      call VFILL (NR,UB1,0.D0)
      call PRSRED (NB,NR,IA,JA,A,DB,UB1)
      call VFILL (NB,WB,0.D0) 
      call PRSBLK (NB,NR,IA,JA,A,UB1,WB)
      DNRM = DDOT2(NB,DB,1,WB,1)       
      if (ISYM.eq.0) GO TO 10 
      RHOTMP = DDOT2(NB,WB,1,DB1,1)    
      call PARCON (DNRM,C1,C2,C3,C4,GAMOLD,RHOTMP,2)      
      RHOOLD = RHOTMP       
      GO TO 20    
   10 call PARCON (DNRM,C1,C2,C3,C4,GAMOLD,RHOOLD,2)      
   20 call SUM3 (NB,C1,DB,C2,UB(NRP1),C3,UB1(NRP1))       
      call SUM3 (NB,C1,WB,C4,DB,C3,DB1) 
   30 call ITERM (NB,A(NRP1),UB(NRP1),WB,6)     
      return      
      END 
      subroutine ITRSSI (N,NNB,IA,JA,A,RHS,UB,UB1,DB)     
      integer IA(1),JA(1),N,NNB       
      double precision A(1),RHS(N),UB(N),UB1(N),DB(NNB)   
      integer NB,NR,NRP1    
      double precision CONST,C1,C2,C3,DNRM      
      logical :: Q1
      double precision DDOT2 
      logical :: TSTCHG
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      NB = NNB    
      NR = N-NB   
      NRP1 = NR+1 
      call DCOPY2 (NR,RHS,1,UB,1)      
      call PRSRED (NB,NR,IA,JA,A,UB(NRP1),UB)   
      call DCOPY2 (NB,RHS(NRP1),1,DB,1)
      call PRSBLK (NB,NR,IA,JA,A,UB,DB) 
      call VEVMW (NB,DB,UB(NRP1))     
      DELNNM = DDOT2(NB,DB,1,DB,1)     
      DNRM = DELNNM 
      CONST = CME 
      call PSTOP (NB,UB(NRP1),DNRM,CONST,2,Q1)  
      if (HALT) GO TO 20    
      if (.not.ADAPT) GO TO 10
      if (.not.TSTCHG(2)) GO TO 10    
      call VFILL (NR,UB1,0.D0)
      call PRSRED (NB,NR,IA,JA,A,DB,UB1)
      DNRM = DDOT2(NR,UB1,1,UB1,1)     
      call CHGSI (DNRM,2)   
      if (.not.ADAPT) GO TO 10
      call DCOPY2 (NB,UB(NRP1),1,UB1(NRP1),1)    
      call DAXPY2 (NB,GAMMA,DB,1,UB1(NRP1),1)    
      GO TO 20    
   10 call PARSI (C1,C2,C3,2) 
      call SUM3 (NB,C1,DB,C2,UB(NRP1),C3,UB1(NRP1))       
   20 call ITERM (NB,A(NRP1),UB(NRP1),DB,7)     
      return      
      END 
      integer function BISRCH (N,K,L) 
      integer N,L,K(N)      
      integer JLEFT,JMID,JRIGHT       
      JLEFT = 1   
      JRIGHT = N  
      if (N.eq.2) GO TO 40  
      JMID = (N+1)/2
   10 if (L.ge.K(JMID)) GO TO 20      
      JRIGHT = JMID 
      GO TO 30    
   20 JLEFT = JMID
   30 if (JRIGHT-JLEFT.eq.1) GO TO 40 
      JMID = JLEFT+(JRIGHT-JLEFT+1)/2 
      GO TO 10    
   40 BISRCH = JLEFT
      return      
      END 
      double precision function CHEBY (QA,QT,RRR,IP,CME,SME)
      integer :: IP
      double precision CME,QA,QT,RRR,SME
      double precision X,Y,Z
      Z = .5D0*(QA+DSQRT(DABS(QA**2-QT**2)))*(1.D0+RRR**IP) 
      X = Z**(1.D0/DBLE(FLOAT(IP)))   
      Y = (X+RRR/X)/(1.D0+RRR)
      CHEBY = .5D0*(CME+SME+Y*(2.D0-CME-SME))   
      return      
      END 
      subroutine CHGCON (TRI,GAMOLD,RHOOLD,IBMTH) 
      integer :: IBMTH
      double precision TRI(2,1),GAMOLD,RHOOLD   
      integer IB2,IB3,IER,IP
      double precision CMOLD,END,START,EIGVSS,EIGVNS      
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      GO TO (10,20,30), IBMTH 
   10 START = CME 
      IP = IN     
      GO TO 40    
   20 START = CME**2
      IP = IN     
      GO TO 40    
   30 if (ADAPT) START = SPR
      if (.not.ADAPT) START = SPECR   
      IP = IN-IS  
   40 if (IP.ge.2) GO TO 60 
      if (IP.eq.1) GO TO 50 
      END = 0.D0  
      CMOLD = 0.D0
      GO TO 110   
   50 END = 1.D0-1.D0/GAMMA 
      TRI(1,1) = END
      TRI(2,1) = 0.D0       
      GO TO 110   
   60 if ((IP.gt.2).and.(DABS(START-CMOLD).le.ZETA*START)) GO TO 120
      CMOLD = START 
      TRI(1,IP) = 1.D0-1.D0/GAMMA     
      TRI(2,IP) = (RHO-1.D0)/(RHO*RHOOLD*GAMMA*GAMOLD)    
      if (ISYM.ne.0) GO TO 80 
      END = EIGVSS(IP,TRI,START,ZETA,ITMAX,IER) 
      if (IER.eq.0) GO TO 100 
      if (LEVEL.ge.2) write (NOUT,70) IER       
70 format (/10X,'DIFFICULTY IN COMPUTATION OF MAXIMUM EIGENVALUE'/15X &
,'OF ITERATION MATRIX'/10X,'SUBROUTINE ZBRENT RETURNED IER =', I5)
      GO TO 100   
   80 IB2 = 1+IP  
      IB3 = IB2+IP/2+1      
      END = EIGVNS(IP,TRI,TRI(1,IB2),TRI(1,IB3),IER)      
      if (IER.eq.0) GO TO 100 
      if (LEVEL.ge.2) write (NOUT,90) IER       
90 format (/10X,'DIFFICULTY IN COMPUTATION OF MAXIMUM EIGENVALUE'/15X &
,'OF ITERATION MATRIX'/10X,'SUBROUTINE EQRT1S RETURNED IER =', I5)
  100 continue    
      if (IER.ne.0) GO TO 130 
  110 if (IBMTH.eq.1) CME = END       
      if (IBMTH.eq.2) CME = DSQRT(DABS(END))    
      if (IBMTH.eq.3.and.ADAPT) SPR = END       
      if (IBMTH.eq.3.and..not.ADAPT) SPECR = END
      return      
  120 ADAPT = .false.       
      PARTAD = .false.      
      return      
  130 ADAPT = .false.       
      PARTAD = .false.      
      if (LEVEL.ge.2) write (NOUT,140) IN,START 
140 format (/10X,'ESTIMATE OF MAXIMUM EIGENVALUE OF JACOBI   '/15X, 'MATRIX (CME) NOT ACCURATE'/10X, &
'ADAPTIVE PROCEDURE TURNED OFF AT ITERATION ',I5/10X, 'FINAL ESTIMATE OF MAXIMUM EIGENVALUE =',D15.7/)
      return      
      END 
      subroutine CHGSI (DTNRM,IBMTH)  
      integer :: IBMTH
      double precision DTNRM
      double precision CMOLD,ZM1,ZM2  
      double precision CHEBY
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      GO TO (10,30,50), IBMTH 
   10 continue    
      if (IN.eq.0) ZM1 = CME
      if (IN.ne.0) ZM1 = CHEBY(QA,QT,RRR,IN-IS,CME,SME)   
      ZM2 = DTNRM/DELNNM    
      CMOLD = CME 
      CME = DMAX1(ZM1,ZM2,CMOLD)      
      if (CME.ge.1.D0) GO TO 20       
      if (CASEII) SME = -CME
      SIGE = (CME-SME)/(2.D0-CME-SME) 
      GAMMA = 2.D0/(2.D0-CME-SME)     
RRR = (1.D0-DSQRT(DABS(1.D0-SIGE*SIGE)))/(1.D0+DSQRT(DABS(1.D0- SIGE*SIGE)))
      IS = IN     
      DELSNM = DELNNM       
      RHO = 1.D0  
      if (LEVEL.ge.2) write (NOUT,90) IN,ZM1,ZM2,CME,GAMMA,CME      
      return      
   20 CME = CMOLD 
      ADAPT = .false.       
      if (LEVEL.ge.2) write (NOUT,110) IN,CME   
      return      
   30 continue    
      if (IN.eq.0) ZM1 = CME
      if (IN.ne.0) ZM1 = CHEBY(QA,QT,RRR,2*(IN-IS),0.D0,0.D0)       
      ZM2 = DSQRT(DABS(DTNRM/DELNNM)) 
      CMOLD = CME 
      CME = DMAX1(ZM1,ZM2,CMOLD)      
      if (CME.ge.1.D0) GO TO 40       
      SIGE = CME*CME/(2.D0-CME*CME)   
      GAMMA = 2.D0/(2.D0-CME*CME)     
RRR = (1.D0-DSQRT(DABS(1.D0-CME*CME)))/(1.D0+DSQRT(DABS(1.D0-CME* CME)))
      IS = IN     
      DELSNM = DELNNM       
      RHO = 1.D0  
      if (LEVEL.ge.2) write (NOUT,90) IN,ZM1,ZM2,CME,GAMMA,CME      
      return      
   40 CME = CMOLD 
      ADAPT = .false.       
      if (LEVEL.ge.2) write (NOUT,110) IN,CME   
      return      
   50 continue    
      if (SPECR.eq.0.D0) SPECR = .171572875D0   
      if (IN.eq.0) GO TO 60 
      ZM1 = CHEBY(QA,QT,RRR,IN-IS,SPECR,0.D0)   
      GO TO 70    
   60 ZM1 = SPECR 
      SPR = SPECR 
   70 ZM2 = DTNRM/DELNNM    
      if (ADAPT) GO TO 80   
      SPECR = DMAX1(ZM1,ZM2,SPECR)    
      IS = IN+1   
      DELSNM = DELNNM       
      if (LEVEL.ge.2) write (NOUT,100) IN,ZM1,ZM2,CME,SPECR 
      return      
   80 SPR = DMAX1(ZM1,ZM2,SPR)
      return      
90 format (/30X,'PARAMETERS WERE CHANGED AT ITERATION NO.',I5/35X, 'SOLUTION TO CHEBYSHEV EQN.       =',D15.7/35X, &
'SOLUTION TO RAYLEIGH QUOTIENT    =',D15.7/35X, 'NEW ESTIMATE FOR CME             =',D15.7/35X, &
'NEW ESTIMATE FOR GAMMA           =',D15.7/35X, 'NEW ESTIMATE FOR SPECTRAL RADIUS =',D15.7/)
100 format (/30X,'PARAMETERS WERE CHANGED AT ITERATION NO.',I5/35X, 'SOLUTION TO CHEBYSHEV EQN.       =',D15.7/35X, &
'SOLUTION TO RAYLEIGH QUOTIENT    =',D15.7/35X, 'NEW ESTIMATE FOR CME             =',D15.7/35X, &
'NEW ESTIMATE FOR SPECTRAL RADIUS =',D15.7/)
110 format (/10X,'ESTIMATE OF MAXIMUM EIGENVALUE OF JACOBI   '/15X, 'MATRIX (CME) TOO LARGE'/10X, &
'ADAPTIVE PROCEDURE TURNED OFF AT ITERATION ',I5/10X, 'FINAL ESTIMATE OF MAXIMUM EIGENVALUE =',D15.7/)
      END 
      logical function CHGSME (OLDNRM,ICNT)     
      integer :: ICNT
      double precision OLDNRM 
      integer :: IP
      double precision Q,RN,SM1,SM2,WP,Z
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      CHGSME = .false.      
      RN = DSQRT(DELNNM/OLDNRM)       
      if (.not.(QA.gt.1.D0.and.RN.gt.1.D0)) return
      if (IN.le.IS+2) return
      ICNT = ICNT+1 
      if (ICNT.lt.3) return 
      CHGSME = .true.       
      SM1 = 0.D0  
      SM2 = 0.D0  
      if (SME.ge.CME) GO TO 10
      IP = IN-IS  
      Q = QA*(1.D0+RRR**IP)/(2.D0*DSQRT(RRR**IP)) 
      Z = (Q+DSQRT(Q**2-1.D0))**(1.D0/DBLE(FLOAT(IP)))    
      WP = (Z**2+1.D0)/(2.D0*Z)       
      SM1 = .5D0*(CME+SME-WP*(CME-SME)) 
      Q = RN*(1.D0+RRR**IP)/((1.D0+RRR**(IP-1))*DSQRT(RRR)) 
      WP = (Q**2+1.D0)/(2.D0*Q)       
      SM2 = .5D0*(CME+SME-WP*(CME-SME)) 
   10 SME = DMIN1(1.25D0*SM1,1.25D0*SM2,SME,-1.D0)
      SIGE = (CME-SME)/(2.D0-CME-SME) 
      GAMMA = 2.D0/(2.D0-CME-SME)     
      RRR = (1.D0-DSQRT(1.D0-SIGE**2))/(1.D0+DSQRT(1.D0-SIGE**2))   
      IS = IN     
      DELSNM = DELNNM       
      RHO = 1.D0  
      if (LEVEL.ge.2) write (NOUT,20) IN,SM1,SM2,SME      
20 format (/30X,'ESTIMATE OF SMALLEST EIGENVALUE OF JACOBI'/37X, 'MATRIX (SME) CHANGED AT ITERATION ',I5/35X, &
'FIRST ESTIMATE OF SME            =',D15.7/35X, 'SECOND ESTIMATE OF SME           =',D15.7/35X, &
'NEW ESTIMATE OF SME              =',D15.7/)
      return      
      END 
      subroutine DAXPY2 (N,DA,DX,INCX,DY,INCY)   
      
      double precision DX(N),DY(N),DA 
      
      if (N.le.0.or.DA.eq.0.D0) return
      if (INCX.eq.INCY) if (INCX-1) 10 , 30 , 70
   10 continue    
      IX = 1      
      IY = 1      
      if (INCX.lt.0) IX = (-N+1)*INCX+1 
      if (INCY.lt.0) IY = (-N+1)*INCY+1 
      do 20 I = 1,N 
         DY(IY) = DY(IY)+DA*DX(IX)    
         IX = IX+INCX       
         IY = IY+INCY       
   20 continue    
      return      
   30 M = N-(N/4)*4 
      if (M.eq.0) GO TO 50  
      do 40 I = 1,M 
         DY(I) = DY(I)+DA*DX(I)       
   40 continue    
      if (N.lt.4) return    
   50 MP1 = M+1   
      do 60 I = MP1,N,4     
         DY(I) = DY(I)+DA*DX(I)       
         DY(I+1) = DY(I+1)+DA*DX(I+1) 
         DY(I+2) = DY(I+2)+DA*DX(I+2) 
         DY(I+3) = DY(I+3)+DA*DX(I+3) 
   60 continue    
      return      
   70 continue    
      NS = N*INCX 
      do 80 I = 1,NS,INCX   
         DY(I) = DA*DX(I)+DY(I)       
   80 continue    
      return      
      END 
      subroutine DCOPY2 (N,DX,INCX,DY,INCY)      
      double precision DX(N),DY(N)
      if (N.le.0) return    
      if (INCX.eq.INCY) if (INCX-1) 10 , 30 , 70
   10 continue    
      IX = 1      
      IY = 1      
      if (INCX.lt.0) IX = (-N+1)*INCX+1 
      if (INCY.lt.0) IY = (-N+1)*INCY+1 
      do 20 I = 1,N 
         DY(IY) = DX(IX)    
         IX = IX+INCX       
         IY = IY+INCY       
   20 continue    
      return      
   30 M = N-(N/7)*7 
      if (M.eq.0) GO TO 50  
      do 40 I = 1,M 
         DY(I) = DX(I)      
   40 continue    
      if (N.lt.7) return    
   50 MP1 = M+1   
      do 60 I = MP1,N,7     
         DY(I) = DX(I)      
         DY(I+1) = DX(I+1)  
         DY(I+2) = DX(I+2)  
         DY(I+3) = DX(I+3)  
         DY(I+4) = DX(I+4)  
         DY(I+5) = DX(I+5)  
         DY(I+6) = DX(I+6)  
   60 continue    
      return      
   70 continue    
      NS = N*INCX 
      do 80 I = 1,NS,INCX   
         DY(I) = DX(I)      
   80 continue    
      return      
      END 
      double precision function DDOT2 (N,DX,INCX,DY,INCY)  
      double precision DX(N),DY(N)
      
      DDOT2 = 0.D0 
      if (N.le.0) return    
      if (INCX.eq.INCY) if (INCX-1) 10 , 30 , 70
   10 continue    
      IX = 1      
      IY = 1      
      if (INCX.lt.0) IX = (-N+1)*INCX+1 
      if (INCY.lt.0) IY = (-N+1)*INCY+1 
      do 20 I = 1,N 
         DDOT2 = DDOT2+DX(IX)*DY(IY)    
         IX = IX+INCX       
         IY = IY+INCY       
   20 continue    
      return      
   30 M = N-(N/5)*5 
      if (M.eq.0) GO TO 50  
      do 40 I = 1,M 
         DDOT2 = DDOT2+DX(I)*DY(I)      
   40 continue    
      if (N.lt.5) return    
   50 MP1 = M+1   
      do 60 I = MP1,N,5     
DDOT2 = DDOT2+DX(I)*DY(I)+DX(I+1)*DY(I+1)+DX(I+2)*DY(I+2)+DX(I+3) *DY(I+3)+DX(I+4)*DY(I+4)
   60 continue    
      return      
   70 continue    
      NS = N*INCX 
      do 80 I = 1,NS,INCX   
         DDOT2 = DDOT2+DX(I)*DY(I)      
   80 continue    
      return      
      END 
      double precision function DETERM (N,TRI,XLMDA)      
      integer :: N
      double precision TRI(2,1),XLMDA 
      integer ICNT,L,NM1    
      double precision D1,D2,D3       
      NM1 = N-1   
      D2 = TRI(1,N)-XLMDA   
      D1 = D2*(TRI(1,NM1)-XLMDA)-TRI(2,N)       
      if (N.eq.2) GO TO 20  
      do 10 ICNT = 2,NM1    
         L = NM1-ICNT+2     
         D3 = D2  
         D2 = D1  
         D1 = (TRI(1,L-1)-XLMDA)*D2-D3*TRI(2,L) 
   10 continue    
   20 DETERM = D1 
      return      
      END 
      subroutine DFAULT (IPARM,RPARM) 
      integer :: IPARM(12)
      double precision RPARM(12)      
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      DRELPR = 7.11D-15
      IPARM(1) = 100
      IPARM(2) = 0
      IPARM(3) = 0
      IPARM(4) = 6
      IPARM(5) = 0
      IPARM(6) = 1
      IPARM(7) = 1
      IPARM(8) = 0
      IPARM(9) = -1 
      IPARM(10) = 0 
      IPARM(11) = 0 
      IPARM(12) = 0 
      RPARM(1) = 0.5D-5     
      RPARM(2) = 0.D0       
      RPARM(3) = 0.D0       
      RPARM(4) = .75D0      
      RPARM(5) = 1.D0       
      RPARM(6) = 0.D0       
      RPARM(7) = .25D0      
      RPARM(8) = 1.D2*DRELPR
      RPARM(9) = 0.D0       
      RPARM(10) = 0.D0      
      RPARM(11) = 0.D0      
      RPARM(12) = 0.D0      
      return      
      END 
      subroutine ECHALL (NN,IA,JA,A,RHS,IPARM,RPARM,ICALL)
      integer IA(1),JA(1),IPARM(12),NN,ICALL    
      double precision A(1),RHS(NN),RPARM(12)   
      integer I,N,NP1,NZRO  
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      if (ICALL.ne.1) GO TO 100       
      N = NN      
      NP1 = N+1   
      NZRO = IA(NP1)-1      
      ZETA = RPARM(1)       
      CME = RPARM(2)
      SME = RPARM(3)
      FF = RPARM(4) 
      OMEGA = RPARM(5)      
      SPECR = RPARM(6)      
      BETAB = RPARM(7)      
      ITMAX = IPARM(1)      
      LEVEL = IPARM(2)      
      ISYM = IPARM(5)       
      ADAPT = .false.       
      PARTAD = .false.      
      BETADT = .false.      
      if (IPARM(6).eq.1.or.IPARM(6).eq.3) ADAPT = .true.  
      if (IPARM(6).eq.1) BETADT = .true.
      if (IPARM(6).eq.2) PARTAD = .true.
      CASEII = .false.      
      if (IPARM(7).eq.2) CASEII = .true.
      if (CASEII) SME = -CME
      if (.not.CASEII.and.SME.eq.0.D0) SME = -1.D0
      SPR = SME   
      IN = 0      
      IS = 0      
      HALT = .false.
      BDELNM = 0.D0 
      DELNNM = 0.D0 
      DELSNM = 0.D0 
      GAMMA = 0.D0
      QA = 0.D0   
      QT = 0.D0   
      RHO = 0.D0  
      RRR = 0.D0  
      SIGE = 0.D0 
      STPTST = 0.D0 
      UDNM = 0.D0 
      if (LEVEL.le.4) GO TO 80
      write (NOUT,10)       
   10 format (///30X,'THE LINEAR SYSTEM IS AS FOLLOWS')   
      write (NOUT,20)       
   20 format (/2X,'IA ARRAY') 
      write (NOUT,30) (IA(I),I=1,NP1) 
   30 format (2X,10(2X,I8)) 
      write (NOUT,40)       
   40 format (/2X,'JA ARRAY') 
      write (NOUT,30) (JA(I),I=1,NZRO)
      write (NOUT,50)       
   50 format (/2X,' A ARRAY') 
      write (NOUT,60) (A(I),I=1,NZRO) 
   60 format (2X,5(2X,D20.13))
      write (NOUT,70)       
   70 format (/2X,'RHS ARRAY')
      write (NOUT,60) (RHS(I),I=1,N)  
   80 write (NOUT,90)       
   90 format (///30X,'INITIAL ITERATIVE PARAMETERS')      
      GO TO 120   
  100 write (NOUT,110)      
  110 format (///30X,'FINAL ITERATIVE PARAMETERS')
  120 write (NOUT,130) IPARM(1),LEVEL,IPARM(3),NOUT,ISYM,IPARM(6)   
130 format (35X,'IPARM(1)  =',I15,4X,'(ITMAX)'/35X,'IPARM(2)  =',I15, &
4X,'(LEVEL) '/35X,'IPARM(3)  =',I15,4X,'(IRESET)'/35X, 'IPARM(4)  =',I15,4X,'(NOUT)  '/35X,'IPARM(5)  =',I15,4X, &
'(ISYM)  '/35X,'IPARM(6)  =',I15,4X,'(IADAPT)')
write (NOUT,140) IPARM(7),IPARM(8),IPARM(9),IPARM(10),IPARM(11), IPARM(12)
140 format (35X,'IPARM(7)  =',I15,4X,'(ICASE)'/35X,'IPARM(8)  =',I15, &
4X,'(NWKSP)'/35X,'IPARM(9)  =',I15,4X,'(NB)    '/35X, 'IPARM(10) =',I15,4X,'(IREMOVE)'/35X,'IPARM(11) =',I15,4X, &
'(ITIME)'/35X,'IPARM(12) =',I15,4X,'(IDGTS)')
      write (NOUT,150) ZETA,CME,SME,FF,OMEGA,SPECR
150 format (35X,'RPARM(1)  =',D15.8,4X,'(ZETA)  '/35X,'RPARM(2)  =', &
D15.8,4X,'(CME)   '/35X,'RPARM(3)  =',D15.8,4X,'(SME)   '/35X, &
'RPARM(4)  =',D15.8,4X,'(FF)    '/35X,'RPARM(5)  =',D15.8,4X, '(OMEGA) '/35X,'RPARM(6)  =',D15.8,4X,'(SPECR) ')
write (NOUT,160) BETAB,RPARM(8),RPARM(9),RPARM(10),RPARM(11), RPARM(12)
160 format (35X,'RPARM(7)  =',D15.8,4X,'(BETAB) '/35X,'RPARM(8)  =', &
D15.8,4X,'(TOL)'/35X,'RPARM(9)  =',D15.8,4X,'(TIME1)'/35X, 'RPARM(10) =',D15.8,4X,'(TIME2)'/35X,'RPARM(11) =',D15.8,4X, &
'(DIGIT1)'/35X,'RPARM(12) =',D15.8,4X,'(DIGIT2)')
      return      
      END 
      subroutine ECHOUT (IPARM,RPARM,IMTHD)     
      integer IPARM(12),IMTHD 
      double precision RPARM(12)      
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      ZETA = RPARM(1)       
      CME = RPARM(2)
      SME = RPARM(3)
      FF = RPARM(4) 
      OMEGA = RPARM(5)      
      SPECR = RPARM(6)      
      BETAB = RPARM(7)      
      ITMAX = IPARM(1)      
      LEVEL = IPARM(2)      
      ISYM = IPARM(5)       
      ADAPT = .false.       
      PARTAD = .false.      
      BETADT = .false.      
      if (IPARM(6).eq.1.or.IPARM(6).eq.3) ADAPT = .true.  
      if (IPARM(6).eq.1) BETADT = .true.
      if (IPARM(6).eq.2) PARTAD = .true.
      CASEII = .false.      
      if (IPARM(7).eq.2) CASEII = .true.
      if (CASEII) SME = -CME
      if (.not.CASEII.and.SME.eq.0.D0) SME = -1.D0
      SPR = SME   
      IN = 0      
      IS = 0      
      HALT = .false.
      BDELNM = 0.D0 
      DELNNM = 0.D0 
      DELSNM = 0.D0 
      GAMMA = 0.D0
      QA = 0.D0   
      QT = 0.D0   
      RHO = 0.D0  
      RRR = 0.D0  
      SIGE = 0.D0 
      STPTST = 0.D0 
      UDNM = 0.D0 
      if (LEVEL.le.2) return
      write (NOUT,10) ISYM,ITMAX,ZETA,ADAPT,CASEII
10 format (///30X,'INITIAL ITERATIVE PARAMETERS',3X, 'RELEVANT SWITCHES'/35X,'ISYM   =',I15,8X,'IPARM(5)'/35X, &
'ITMAX  =',I15,8X,'IPARM(1)'/35X,'ZETA   =',D15.8,8X,'RPARM(1)' /35X,'ADAPT  =',L15,8X,'IPARM(6)'/35X,'CASEII =',L15,8X, &
'IPARM(7)')
      GO TO (80,20,100,60,40,80,20), IMTHD      
   20 write (NOUT,30) FF,CME,SME      
30 format (35X,'FF     =',D15.8,8X,'RPARM(4)'/35X,'CME    =',D15.8,8X ,'RPARM(2)'/35X,'SME    =',D15.8,8X,'RPARM(3)'///)
      return      
   40 write (NOUT,50) PARTAD,FF,CME,OMEGA,SPECR,BETAB,BETADT
50 format (35X,'PARTAD =',L15,8X,'IPARM(6)'/35X,'FF     =',D15.8,8X, &
'RPARM(4)'/35X,'CME    =',D15.8,8X,'RPARM(2)'/35X,'OMEGA  =', &
D15.8,8X,'RPARM(5)'/35X,'SPECR  =',D15.8,8X,'RPARM(6)'/35X, &
'BETAB  =',D15.8,8X,'RPARM(7)'/35X,'BETADT =',L15,8X,'IPARM(6)' ///)
      return      
   60 write (NOUT,70) PARTAD,CME,OMEGA,SPECR,BETAB,BETADT 
70 format (35X,'PARTAD =',L15,8X,'IPARM(6)'/35X,'CME    =',D15.8,8X, &
'RPARM(2)'/35X,'OMEGA  =',D15.8,8X,'RPARM(5)'/35X,'SPECR  =', &
D15.8,8X,'RPARM(6)'/35X,'BETAB  =',D15.8,8X,'RPARM(7)'/35X, 'BETADT =',L15,8X,'IPARM(6)'///)
      return      
   80 if (ADAPT) return     
      write (NOUT,90) CME   
   90 format (35X,'CME    =',D15.8,8X,'RPARM(2)'///)      
  100 continue    
      return      
      END 
      double precision function EIGVNS (N,TRI,D,E2,IER)   
      integer N,IER 
      double precision TRI(2,N),D(N),E2(N)      
      integer :: I
      EIGVNS = 0.D0 
      D(1) = -TRI(1,1)      
      do 10 I = 2,N 
         D(I) = -TRI(1,I)   
         E2(I) = DABS(TRI(2,I))       
   10 continue    
      call EQRT1S (D,E2,N,1,0,IER)    
      EIGVNS = -D(1)
      return      
      END 
      double precision function EIGVSS (N,TRI,START,ZETA,ITMAX,IER) 
      integer N,ITMAX,IER   
      double precision TRI(2,1),START,ZETA,A,B,EPS
      integer MAXFN,NSIG,ITMP 
      EIGVSS = 0.D0 
      ITMP = IFIX(SNGL(-DLOG10(DABS(ZETA))))    
      NSIG = MAX0(ITMP,4)   
      MAXFN = MAX0(ITMAX,50)
      EPS = 0.0D0 
      A = START   
      B = 1.0D0   
      call ZBRENT (N,TRI,EPS,NSIG,A,B,MAXFN,IER)
      EIGVSS = B  
      return      
      END 
      subroutine EQRT1S (D,E2,NN,M,ISW,IERR)    
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      integer NN,M,ISW,IERR 
      double precision D(NN),E2(NN)   
      integer II,I,JJ,J,K1,K,N,IER    
      double precision DELTA,DLAM,EP,ERR,F,P,QP,Q,R,S,TOT 
      N = NN      
      IER = 0     
      DLAM = 0.0D0
      ERR = 0.0D0 
      S = 0.0D0   
      TOT = D(1)  
      Q = 0.0D0   
      J = 0       
      do 30 I = 1,N 
         P = Q    
         if (I.eq.1) GO TO 10 
         if (P.gt.DRELPR*(DABS(D(I))+DABS(D(I-1)))) GO TO 20
   10    E2(I) = 0.0D0      
   20    if (E2(I).eq.0.D0) J = J+1   
         Q = 0.0D0
         if (I.ne.N) Q = DSQRT(DABS(E2(I+1)))   
         TOT = DMIN1(D(I)-P-Q,TOT)    
   30 continue    
      if (ISW.eq.1.and.TOT.lt.0.0D0) GO TO 50   
      do 40 I = 1,N 
         D(I) = D(I)-TOT    
   40 continue    
      GO TO 60    
   50 TOT = 0.0D0 
   60 do 200 K = 1,M
   70    TOT = TOT+S
         DELTA = D(N)-S     
         I = N    
         F = DABS(DRELPR*TOT) 
         if (DLAM.lt.F) DLAM = F      
         if (DELTA.gt.DLAM) GO TO 90  
         if (DELTA.ge.(-DLAM)) GO TO 170
         IER = 602
         if (LEVEL.ge.1) write (NOUT,80)
80    format ('0','*** W A R N I N G ************'/' ', '    IN ITPACK ROUTINE EQRT1S  '/' ', &
'    PARAMETER ISW = 1 BUT MATRIX   '/' ', '    NOT POSITIVE DEFINITE')
         GO TO 210
   90    if (K.eq.N) GO TO 110
         K1 = K+1 
         do 100 J = K1,N    
            if (E2(J).le.(DRELPR*(D(J)+D(J-1)))**2) E2(J) = 0.0D0   
  100    continue 
  110    F = E2(N)/DELTA    
         QP = DELTA+F       
         P = 1.0D0
         if (K.eq.N) GO TO 140
         K1 = N-K 
         do 130 II = 1,K1   
            I = N-II
            Q = D(I)-S-F    
            R = Q/QP
            P = P*R+1.0D0   
            EP = F*R
            D(I+1) = QP+EP  
            DELTA = Q-EP    
            if (DELTA.gt.DLAM) GO TO 120
            if (DELTA.ge.(-DLAM)) GO TO 170     
            IER = 602       
            if (LEVEL.ge.0) write (NOUT,80)     
            GO TO 210       
  120       F = E2(I)/Q     
            QP = DELTA+F    
            E2(I+1) = QP*EP 
  130    continue 
  140    D(K) = QP
         S = QP/P 
         if (TOT+S.gt.TOT) GO TO 70   
         IER = 601
         E2(1) = K
         if (LEVEL.ge.1) write (NOUT,150) K     
150    format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE EQRT1S  '/' ', &
'    SUCCESSIVE ITERATES TO THE',I10/' ', '    EIGENVALUE WERE NOT MONOTONE INCREASING ')
         S = 0.0D0
         DELTA = QP 
         do 160 J = K,N     
            if (D(J).gt.DELTA) GO TO 160
            I = J 
            DELTA = D(J)    
  160    continue 
  170    if (I.lt.N) E2(I+1) = E2(I)*F/QP       
         if (I.eq.K) GO TO 190
         K1 = I-K 
         do 180 JJ = 1,K1   
            J = I-JJ
            D(J+1) = D(J)-S 
            E2(J+1) = E2(J) 
  180    continue 
  190    D(K) = TOT 
         ERR = ERR+DABS(DELTA)
         E2(K) = ERR
  200 continue    
      if (IER.eq.0) GO TO 220 
  210 continue    
  220 IERR = IER  
      return      
      END 
      integer function IPSTR (OMEGA)  
      double precision OMEGA
      integer :: IP
      double precision WM1  
      WM1 = OMEGA-1.D0      
      do 10 IP = 6,940      
         if (DBLE(FLOAT(IP))*(WM1**(IP-1)).gt.0.50D0) GO TO 10      
         IPSTR = IP 
         return   
   10 continue    
      IPSTR = 940 
      return      
      END 
      subroutine ITERM (NN,A,U,WK,IMTHDD)       
      integer NN,IMTHD      
      double precision A(1),U(NN),WK(NN)
      integer I,IMTHDD,IP,N 
      double precision QTFF 
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      N = NN      
      IMTHD = IMTHDD
      if (LEVEL.lt.2) return
      GO TO (10,110,170,210,50,10,110), IMTHD   
   10 if (IN.gt.0) GO TO 30 
      write (NOUT,20)       
20 format (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'// &
' NUMBER OF',5X,'CONVERGENCE',7X,'CME ',11X,'RHO',12X,'GAMMA'/ ' ITERATIONS',4X,'TEST '//)
   30 write (NOUT,40) IN,STPTST,CME,RHO,GAMMA   
   40 format (4X,I5,3X,4D15.7)
      if (LEVEL.ge.4) GO TO 250       
      return      
   50 if (IN.gt.0) GO TO 70 
      write (NOUT,60)       
60 format (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'// &
' NUMBER OF',4X,'CONVERGENCE',7X,'PARAMETER CHANGE TEST',10X, &
'RHO',12X,'GAMMA'/' ITERATIONS',3X,'TEST ',11X,'LHS(QA)',7X, 'RHS(QT**FF)'//)
   70 IP = IN-IS  
      if (IMTHD.eq.7) IP = 2*IP       
      if (IP.lt.3) GO TO 90 
      QTFF = QT**FF 
      write (NOUT,80) IN,STPTST,QA,QTFF,RHO,GAMMA 
   80 format (4X,I5,3X,5D15.7)
      if (LEVEL.ge.4) GO TO 250       
      return      
   90 write (NOUT,100) IN,STPTST,RHO,GAMMA      
  100 format (4X,I5,3X,D15.7,30X,2D15.7)
      if (LEVEL.ge.4) GO TO 250       
      return      
  110 if (IN.gt.0) GO TO 130
      write (NOUT,120)      
120 format (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'// &
' NUMBER OF',4X,'CONVERGENCE',7X,'PARAMETER CHANGE TEST',10X, &
'RHO'/' ITERATIONS',3X,'TEST ',11X,'LHS(QA)',7X,'RHS(QT**FF)'// )
  130 IP = IN-IS  
      if (IMTHD.eq.7) IP = 2*IP       
      if (IP.lt.3) GO TO 150
      QTFF = QT**FF 
      write (NOUT,140) IN,STPTST,QA,QTFF,RHO    
  140 format (4X,I5,3X,5D15.7)
      if (LEVEL.ge.4) GO TO 250       
      return      
  150 write (NOUT,160) IN,STPTST,RHO  
  160 format (4X,I5,3X,D15.7,30X,D15.7) 
      if (LEVEL.ge.4) GO TO 250       
      return      
  170 if (IN.gt.0) GO TO 190
      write (NOUT,180)      
180 format (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'// ' NUMBER OF',4X,'CONVERGENCE',6X,'CME ',9X,'OMEGA',7X, &
'SPECTRAL'/' ITERATIONS',3X,'TEST',38X,'RADIUS'//)
  190 continue    
      write (NOUT,200) IN,STPTST,CME,OMEGA,SPECR
  200 format (4X,I5,3X,4D14.7)
      if (LEVEL.ge.4) GO TO 250       
      return      
  210 if (IN.gt.0) GO TO 230
      write (NOUT,220)      
220 format (////15X,'INTERMEDIATE OUTPUT AFTER EACH ITERATION'// &
' NUMBER OF',4X,'CONVERGENCE',3X,' SPECTRAL',6X,'S-PRIME',9X, 'RHO',10X,'GAMMA'/' ITERATIONS',3X,'TEST ',10X,'RADIUS'//)
  230 continue    
      write (NOUT,240) IN,STPTST,SPECR,SPR,RHO,GAMMA      
  240 format (4X,I5,3X,5D14.7)
      if (LEVEL.ge.4) GO TO 250       
      return      
  250 if (IMTHD.gt.5) GO TO 270       
      write (NOUT,260) IN   
  260 format ('0',2X,'ESTIMATE OF SOLUTION AT ITERATION ',I5)       
      GO TO 290   
  270 write (NOUT,280) IN   
280 format ('0',2X,'ESTIMATE OF SOLUTION AT BLACK POINTS ', 'AT ITERATION ',I5)
  290 do 300 I = 1,N
         WK(I) = U(I)/A(I)  
  300 continue    
      write (NOUT,310) (WK(I),I=1,N)  
  310 format (2X,5(2X,D20.13))
      write (NOUT,320)      
  320 format (//) 
      return      
      END 
      subroutine IVFILL (N,IV,IVAL)   
      integer N,IVAL,IV(N)  
      integer I,M,MP1       
      if (N.le.0) return    
      M = MOD(N,10) 
      if (M.eq.0) GO TO 20  
      do 10 I = 1,M 
         IV(I) = IVAL       
   10 continue    
      if (N.lt.10) return   
   20 MP1 = M+1   
      do 30 I = MP1,N,10    
         IV(I) = IVAL       
         IV(I+1) = IVAL     
         IV(I+2) = IVAL     
         IV(I+3) = IVAL     
         IV(I+4) = IVAL     
         IV(I+5) = IVAL     
         IV(I+6) = IVAL     
         IV(I+7) = IVAL     
         IV(I+8) = IVAL     
         IV(I+9) = IVAL     
   30 continue    
      return      
      END 
      subroutine OMEG (DNRM,IFLAG)    
      integer :: IFLAG
      double precision DNRM 
      double precision TEMP,ZM1,ZM2   
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      ZM1 = 0.D0  
      ZM2 = 0.D0  
      if (IFLAG.eq.1) GO TO 10
ZM1 = ((1.D0-SPR)*(1.D0+BETAB*OMEGA**2)-OMEGA*(2.D0-OMEGA))/(OMEGA *(OMEGA-1.D0-SPR))
      if (.not.CASEII) ZM2 = DNRM/BDELNM
      if (CASEII) ZM2 = DSQRT(DABS(DNRM/BDELNM))
      CME = DMAX1(CME,ZM1,ZM2)
   10 IS = IN+1   
      DELSNM = DELNNM       
      if (CME.ge.(4.D0*BETAB)) GO TO 30 
      TEMP = DSQRT(DABS(1.D0-2.D0*CME+4.D0*BETAB))
      OMEGA = DMAX1((2.D0/(1.D0+TEMP)),1.D0)    
      TEMP = (1.D0-CME)/TEMP
      SPECR = (1.D0-TEMP)/(1.D0+TEMP) 
      if (DABS(OMEGA-1.D0).lt.DRELPR) SPECR = 0.D0
      if (LEVEL.ge.2) write (NOUT,20) IN,BETAB,ZM1,ZM2,CME,OMEGA,SPECR
20 format (/30X,'PARAMETERS WERE CHANGED AT ITERATION NO.',I5/35X, 'NEW ESTIMATE OF BETAB            =',D15.7/35X, &
'SOLUTION TO CHEBYSHEV EQN.       =',D15.7/35X, 'SOLUTION TO RAYLEIGH QUOTIENT    =',D15.7/35X, &
'NEW ESTIMATE FOR CME             =',D15.7/35X, 'NEW ESTIMATE FOR OMEGA           =',D15.7/35X, &
'NEW ESTIMATE FOR SPECTRAL RADIUS =',D15.7/)
      return      
   30 CME = 2.D0*DSQRT(DABS(BETAB))   
      OMEGA = 2.D0/(1.D0+DSQRT(DABS(1.D0-4.D0*BETAB)))    
      SPECR = OMEGA-1.D0    
      ADAPT = .false.       
      PARTAD = .false.      
      if (LEVEL.ge.2) write (NOUT,20) IN,BETAB,ZM1,ZM2,CME,OMEGA,SPECR
      return      
      END 
      logical function OMGCHG (NDUMMY)
      integer :: NDUMMY
      double precision DEL1,DEL2,X    
      double precision PHI  
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      PHI(X) = (1.D0-DSQRT(DABS(1.D0-X)))/(1.D0+DSQRT(DABS(1.D0-X)))
      OMGCHG = .false.      
      if (IN-IS.lt.3) return
      if (SPECR.eq.0.D0) GO TO 10     
      if (SPECR.ge.SPR) return
      DEL1 = -DLOG(DABS(PHI(SPECR)/PHI(SPECR/SPR)))       
      DEL2 = -DLOG(DABS(PHI(SPR)))    
      if ((DEL1/DEL2).ge.FF) return   
   10 OMGCHG = .true.       
      return      
      END 
      logical function OMGSTR (NDUMMY)
      integer :: NDUMMY
      double precision OMSTAR,TEMP,TEMP1,X      
      double precision PHI  
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      PHI(X) = (1.D0-DSQRT(DABS(1.D0-X)))/(1.D0+DSQRT(DABS(1.D0-X)))
      OMGSTR = .false.      
      if (BETAB.ge..25D0.or..not.ADAPT) return  
      OMSTAR = 2.D0/(1.D0+DSQRT(DABS(1.D0-4.D0*BETAB)))   
      if ((OMSTAR.le.1.D0).or.(SPECR.le.0.D0)) GO TO 10   
      TEMP = DLOG(DABS(PHI(OMSTAR-1.D0)))       
      TEMP1 = DLOG(DABS(PHI(SPECR)))  
      if ((TEMP/TEMP1).lt.FF) return  
   10 OMEGA = OMSTAR
      SPECR = OMEGA-1.D0    
      OMGSTR = .true.       
      ADAPT = .false.       
      PARTAD = .false.      
      CME = 2.D0*DSQRT(DABS(BETAB))   
      RRR = PHI(1.D0-SPECR)**2
      GAMMA = 2.D0/(2.D0-SPECR)       
      SIGE = SPECR/(2.D0-SPECR)       
      RHO = 1.D0  
      IS = IN+1   
      DELSNM = DELNNM       
      if (LEVEL.ge.2) write (NOUT,20) IN,CME,OMEGA,SPECR  
20 format (/30X,'OMEGA-STAR, AN ALTERNATE ESTIMATE OF', ' OMEGA, WAS CHOSEN AT ITERATION',I5/35X, &
'NEW ESTIMATE FOR CME             =',D15.7/35X, 'NEW ESTIMATE FOR OMEGA           =',D15.7/35X, &
'NEW ESTIMATE FOR SPECTRAL RADIUS =',D15.7/)
      return      
      END 
      subroutine PARCON (DTNRM,C1,C2,C3,C4,GAMOLD,RHOTMP,IBMTH)     
      integer :: IBMTH
      double precision DTNRM,C1,C2,C3,C4,GAMOLD,RHOTMP    
      integer :: IP
      double precision RHOOLD 
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      IP = IN-IS  
      RHOOLD = RHO
      GAMOLD = GAMMA
      if (IBMTH.le.2) GAMMA = 1.D0/(1.D0-DTNRM/DELNNM)    
      if (IBMTH.eq.3) GAMMA = DELNNM/DTNRM      
      RHO = 1.D0  
      if (IP.eq.0) GO TO 20 
      if (ISYM.eq.0) GO TO 10 
      RHO = 1.D0/(1.D0-GAMMA*RHOTMP/DELSNM)     
      GO TO 20    
   10 RHO = 1.D0/(1.D0-GAMMA*DELNNM/(GAMOLD*DELSNM*RHOOLD)) 
   20 DELSNM = DELNNM       
      RHOTMP = RHOOLD       
      C1 = RHO*GAMMA
      C2 = RHO    
      C3 = 1.D0-RHO 
      C4 = RHO*(1.D0-GAMMA) 
      return      
      END 
      subroutine PARSI (C1,C2,C3,IBMTH) 
      integer :: IBMTH
      double precision C1,C2,C3       
      integer :: IP
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      IP = IN-IS  
      if (IP.eq.0) GO TO 30 
      if (IP.eq.1) GO TO 10 
      RHO = 1.D0/(1.D0-SIGE*SIGE*RHO*.25D0)     
      GO TO 20    
   10 RHO = 1.D0/(1.D0-SIGE*SIGE*.5D0)
   20 C1 = RHO*GAMMA
      C2 = RHO    
      C3 = 1.D0-RHO 
      return      
   30 continue    
      GO TO (40,50,60), IBMTH 
   40 if (CASEII) SME = -CME
      GAMMA = 2.D0/(2.D0-CME-SME)     
      SIGE = (CME-SME)/(2.D0-CME-SME) 
      GO TO 70    
   50 GAMMA = 2.D0/(2.D0-CME*CME)     
      SIGE = CME*CME/(2.D0-CME*CME)   
RRR = (1.D0-DSQRT(DABS(1.D0-CME*CME)))/(1.D0+DSQRT(DABS(1.D0-CME* CME)))
      GO TO 70    
   60 GAMMA = 2.D0/(2.D0-SPECR)       
      SIGE = SPECR/(2.D0-SPECR)       
RRR = (1.D0-DSQRT(DABS(1.D0-SIGE*SIGE)))/(1.D0+DSQRT(DABS(1.D0- SIGE*SIGE)))
   70 RHO = 1.D0  
      C1 = GAMMA  
      C2 = 1.D0   
      C3 = 0.D0   
      return      
      END 
      double precision function PBETA (NN,IA,JA,A,V,W1,W2)
      
      use Global_ITPACK
      integer IA(NN+1),JA(ITPACK_NNZ),NN
      double precision A(ITPACK_NNZ),V(NN),W1(NN),W2(NN) 
      
      integer I,IBGN,IEND,II,ITMP,JAI,JAJJ,JJ,K,N,NM1     
      double precision SUM,TEMP1,TEMP2
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      N = NN      
      PBETA = 0.D0
      if (ISYM.eq.0) GO TO 110
      do 10 I = 1,N 
         W1(I) = V(I)       
   10 continue    
      TEMP1 = 0.D0
      TEMP2 = 0.D0
      ITMP = 2    
      IBGN = IA(1)
      IEND = IA(ITMP)-1     
      if (IEND.lt.IBGN) GO TO 30      
      do 20 I = IBGN,IEND   
         JAI = JA(I)
         TEMP1 = TEMP1-A(I)*W1(JAI)   
   20 continue    
   30 W1(1) = TEMP1 
      W2(1) = 0.D0
      NM1 = N-1   
      do 70 K = 2,NM1       
         TEMP1 = 0.D0       
         TEMP2 = 0.D0       
         IBGN = IA(K)       
         IEND = IA(K+1)-1   
         if (IEND.lt.IBGN) GO TO 60   
         do 50 I = IBGN,IEND
            JAI = JA(I)     
            if (JAI.gt.K) GO TO 40    
            TEMP2 = TEMP2-A(I)*W1(JAI)
            GO TO 50
   40       TEMP1 = TEMP1-A(I)*W1(JAI)
   50    continue 
   60    W1(K) = TEMP1      
         W2(K) = TEMP2      
   70 continue    
      TEMP2 = 0.D0
      IBGN = IA(N)
      IEND = IA(N+1)-1      
      if (IEND.lt.IBGN) GO TO 90      
      do 80 I = IBGN,IEND   
         JAI = JA(I)
         TEMP2 = TEMP2-A(I)*W1(JAI)   
   80 continue    
   90 W2(N) = TEMP2 
      do 100 I = 1,N
         PBETA = PBETA+V(I)*W2(I)     
  100 continue    
      return      
  110 do 130 II = 1,N       
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = 0.D0 
         if (IBGN.gt.IEND) GO TO 130  
         do 120 JJ = IBGN,IEND
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*V(JAJJ)   
  120    continue 
         PBETA = PBETA+SUM*SUM
  130 continue    
      return      
      END 
      subroutine PBSOR (NN,IA,JA,A,U,RHS)       
      
      
      use Global_ITPACK
      integer IA(NN+1),JA(ITPACK_NNZ),NN
      double precision A(ITPACK_NNZ),U(NN),RHS(NN)  
      
      integer I,IBGN,IEND,II,JAJJ,JJ,N,NPL1     
      double precision OMM1,SUM,UI    
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      N = NN      
      NPL1 = N+1  
      OMM1 = OMEGA-1.D0     
      if (ISYM.eq.0) GO TO 40 
      do 30 I = 1,N 
         II = NPL1-I
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         if (IBGN.gt.IEND) GO TO 20   
         do 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   10    continue 
   20    U(II) = OMEGA*SUM-OMM1*U(II) 
   30 continue    
      return      
   40 do 60 II = 1,N
         UI = U(II) 
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         if (IBGN.gt.IEND) GO TO 60   
         do 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UI      
   50    continue 
   60 continue    
      do 90 I = 1,N 
         II = NPL1-I
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         if (IBGN.gt.IEND) GO TO 80   
         do 70 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   70    continue 
   80    U(II) = OMEGA*SUM-OMM1*U(II) 
   90 continue    
      return      
      END 
      subroutine PERMAT (NN,IA,JA,A,P,NEWIA,ISYM,LEVEL,NOUT,IERR)   
      
      
      
      use Global_ITPACK
      
      integer NN,IA(NN+1),JA(ITPACK_NNZ),P(NN),NEWIA(NN),ISYM,IERR    
      double precision A(ITPACK_NNZ) 
      
      integer BISRCH,I,IBGN,IEND,IP,IPP,J,JAJ,JP,IER,K,N,NELS,NEXT,NPL1 
      double precision save,TEMP      
      N = NN      
      IER = 0     
      NPL1 = N+1  
      NELS = IA(NPL1)-1     
      do 10 I = 1,N 
         NEWIA(I) = 0       
   10 continue    
      do 30 I = 1,N 
         IP = P(I)
         IBGN = IA(I)       
         IEND = IA(I+1)-1   
         if (IBGN.gt.IEND) GO TO 90   
         do 20 J = IBGN,IEND
            IPP = IP
            JAJ = JA(J)     
            JP = P(JAJ)     
            if (ISYM.eq.0.and.IP.gt.JP) IPP = JP
            NEWIA(IPP) = NEWIA(IPP)+1 
            JA(J) = -JAJ    
   20    continue 
   30 continue    
      IBGN = 1    
      do 40 I = 1,N 
         K = IBGN+NEWIA(I)  
         NEWIA(I) = IBGN    
         IBGN = K 
   40 continue    
      do 70 J = 1,NELS      
         if (JA(J).gt.0) GO TO 70     
         JAJ = -JA(J)       
         save = A(J)
         NEXT = J 
         JA(J) = JAJ
   50    JP = P(JAJ)
         I = BISRCH(NPL1,IA,NEXT)     
         IP = P(I)
         IPP = IP 
         if (ISYM.ne.0.or.IP.le.JP) GO TO 60    
         IPP = JP 
         JP = IP  
   60    NEXT = NEWIA(IPP)  
         TEMP = save
         save = A(NEXT)     
         A(NEXT) = TEMP     
         JAJ = -JA(NEXT)    
         JA(NEXT) = JP      
         NEWIA(IPP) = NEWIA(IPP)+1    
         if (JAJ.gt.0) GO TO 50       
   70 continue    
      IA(1) = 1   
      do 80 I = 1,N 
         IA(I+1) = NEWIA(I) 
         K = IA(I+1)-IA(I)  
         if (K.eq.1) GO TO 80 
         if (K.lt.1) GO TO 110
         IBGN = IA(I)       
         call QSORT (K,JA(IBGN),A(IBGN),IER)    
         if (IER.ne.0) GO TO 130      
   80 continue    
      GO TO 150   
   90 IER = 301   
      if (LEVEL.ge.0) write (NOUT,100) I
100 format ('0','*** F A T A L     E R R O R ************'/'0', &
'    IN ITPACK ROUTINE PERMAT  '/' ','    NO ENTRY IN ROW ',I10 ,' OF ORIGINAL MATRIX ')
      GO TO 150   
  110 IER = 302   
      if (LEVEL.ge.0) write (NOUT,120) I
120 format ('0','*** F A T A L     E R R O R ************'/'0', &
'    IN ITPACK ROUTINE PRBNDX  '/' ','    NO ENTRY IN ROW ',I10 ,' OF PERMUTED MATRIX ')
      GO TO 150   
  130 IER = 303   
      if (LEVEL.ge.0) write (NOUT,140) I
140 format ('0','*** F A T A L     E R R O R ************'/'0', '    IN ITPACK ROUTINE QSORT   '/' ', &
'    ERROR IN SORTING PERMUTED ROW ',I12/' ', '    CALLED FROM ITPACK ROUTINE PRBNDX   ')
  150 continue    
      IERR = IER  
      return      
      END 
      subroutine c_PERROR (NN,IA,JA,A,RHS,U,W,DIGTT1,DIGTT2,IDGTTS)   
      integer IA(1),JA(1),NN,IDGTTS   
      double precision A(1),RHS(NN),U(NN),W(NN),DIGTT1,DIGTT2       
      integer IDGTS,N       
      double precision BNRM,DIGIT1,DIGIT2,RNRM,TEMP       
      double precision DDOT2 
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      N = NN      
      IDGTS = IDGTTS
      DIGIT1 = 0.D0 
      DIGIT2 = 0.D0 
      if (N.le.0) GO TO 40  
      DIGIT1 = -DLOG10(DABS(DRELPR))  
      if (STPTST.gt.0.D0) DIGIT1 = -DLOG10(DABS(STPTST))  
      BNRM = DDOT2(N,RHS,1,RHS,1)      
      if (BNRM.eq.0.D0) GO TO 10      
      call PMULT (N,IA,JA,A,U,W)      
      call WEVMW (N,RHS,W)  
      RNRM = DDOT2(N,W,1,W,1)
      TEMP = RNRM/BNRM      
      if (TEMP.eq.0.D0) GO TO 10      
      DIGIT2 = -DLOG10(DABS(TEMP))/2.D0 
      GO TO 20    
   10 DIGIT2 = -DLOG10(DABS(DRELPR))  
   20 if ((IDGTS.lt.1).or.(LEVEL.le.0)) GO TO 40
      write (NOUT,30) DIGIT1,DIGIT2   
30 format (/6X,'APPROX. NO. OF DIGITS (EST. REL. ERROR) =',F5.1,2X, &
'(DIGIT1)'/3X,'APPROX. NO. OF DIGITS (EST. REL. RESIDUAL) =', F5.1,2X,'(DIGIT2)')
      if (IDGTS.le.1.or.IDGTS.gt.4) GO TO 40    
      if (IDGTS.ne.3) call VOUT (N,U,2,NOUT)    
      if (IDGTS.ge.3) call VOUT (N,W,1,NOUT)    
   40 continue    
      DIGTT1 = DIGIT1       
      DIGTT2 = DIGIT2       
      return      
      END 
      subroutine PERVEC (N,V,P)       
      integer N,P(N)
      double precision V(N) 
      integer II,NEXT,NOW   
      double precision save,TEMP      
      if (N.le.0) return    
      do 20 II = 1,N
         if (P(II).lt.0) GO TO 20     
         NEXT = P(II)       
         save = V(II)       
   10    continue 
         if (P(NEXT).lt.0) GO TO 20   
         TEMP = save
         save = V(NEXT)     
         V(NEXT) = TEMP     
         NOW = NEXT 
         NEXT = P(NOW)      
         P(NOW) = -NEXT     
         GO TO 10 
   20 continue    
      do 30 II = 1,N
         P(II) = -P(II)     
   30 continue    
      return      
      END 
      subroutine PFSOR (NN,IA,JA,A,U,RHS)       
      
      
      use Global_ITPACK
      
      integer IA(NN+1),JA(ITPACK_NNZ),NN
      double precision A(ITPACK_NNZ),U(NN),RHS(NN)   
      
      integer IBGN,IEND,II,JAJJ,JJ,N  
      double precision OMM1,SUM,UI    
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      N = NN      
      OMM1 = OMEGA-1.D0     
      if (ISYM.eq.0) GO TO 40 
      do 30 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         if (IBGN.gt.IEND) GO TO 20   
         do 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   10    continue 
   20    UI = OMEGA*SUM-OMM1*U(II)    
         U(II) = UI 
   30 continue    
      return      
   40 do 80 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         if (IBGN.gt.IEND) GO TO 60   
         do 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   50    continue 
   60    UI = OMEGA*SUM-OMM1*U(II)    
         U(II) = UI 
         if (IBGN.gt.IEND) GO TO 80   
         do 70 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UI      
   70    continue 
   80 continue    
      return      
      END 
      subroutine PFSOR1 (NN,IA,JA,A,U,RHS)      
      
      use Global_ITPACK
      
      integer IA(NN+1),JA(ITPACK_NNZ),NN
      double precision A(ITPACK_NNZ),U(NN),RHS(NN)   
      
      integer IBGN,IEND,II,JAJJ,JJ,N  
      double precision OMM1,SUM,SUMD,UI 
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      N = NN      
      OMM1 = OMEGA-1.D0     
      SUMD = 0.D0 
      if (ISYM.eq.0) GO TO 40 
      do 30 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         if (IBGN.gt.IEND) GO TO 20   
         do 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   10    continue 
   20    continue 
         UI = OMEGA*SUM-OMM1*U(II)    
         SUMD = SUMD+(UI-U(II))**2    
         U(II) = UI 
   30 continue    
      GO TO 90    
   40 do 80 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         if (IBGN.gt.IEND) GO TO 60   
         do 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   50    continue 
   60    continue 
         UI = OMEGA*SUM-OMM1*U(II)    
         SUMD = SUMD+(UI-U(II))**2    
         U(II) = UI 
         if (IBGN.gt.IEND) GO TO 80   
         do 70 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UI      
   70    continue 
   80 continue    
   90 DELNNM = DSQRT(SUMD)  
      return      
      END 
      subroutine PJAC2 (NN,IA,JA,A,U,RHS)
      
      use Global_ITPACK
      integer IA(NN+1),JA(ITPACK_NNZ),NN
      double precision A(ITPACK_NNZ),U(NN),RHS(NN)     
      integer IBGN,IEND,II,JAJJ,JJ,N  
      double precision RHSII,UII      
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      N = NN      
      if (ISYM.eq.0) GO TO 30 
      do 20 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         if (IBGN.gt.IEND) GO TO 20   
         RHSII = RHS(II)    
         do 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHSII = RHSII-A(JJ)*U(JAJJ) 
   10    continue 
         RHS(II) = RHSII    
   20 continue    
      return      
   30 do 50 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         if (IBGN.gt.IEND) GO TO 50   
         RHSII = RHS(II)    
         UII = U(II)
         do 40 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHSII = RHSII-A(JJ)*U(JAJJ) 
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UII     
   40    continue 
         RHS(II) = RHSII    
   50 continue    
      return      
      END 
      subroutine PMULT (NN,IA,JA,A,U,W) 
      
      use Global_ITPACK

      integer IA(NN+1),JA(ITPACK_NNZ),NN
      double precision A(ITPACK_NNZ),U(NN),W(NN) 
      
      integer IBGN,IEND,II,JJ,N       
      double precision SUM,UII,WII    
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      N = NN      
      if (N.le.0) return    
      if (ISYM.eq.0) GO TO 40 
      do 30 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = 0.0D0
         if (IBGN.gt.IEND) GO TO 20   
         do 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM+A(JJ)*U(JAJJ)   
   10    continue 
   20    W(II) = SUM
   30 continue    
      return      
   40 call VFILL (N,W,0.D0) 
      do 70 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         UII = U(II)
         WII = W(II)
         if (IBGN.gt.IEND) GO TO 60   
         do 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            WII = WII+A(JJ)*U(JAJJ)   
            W(JAJJ) = W(JAJJ)+A(JJ)*UII 
   50    continue 
   60    W(II) = WII
   70 continue    
      return      
      END 
      subroutine PRBNDX (NN,NBLACK,IA,JA,P,IP,LEVEL,NOUT,IER)       
      integer NN,NBLACK,IA(1),JA(1),P(NN),IP(NN),IER      
integer FIRST,NEXT,LAST,I,OLD,YOUNG,IBGN,IEND,J,K,CURTYP,NXTTYP, TYPE,NRED,N
      N = NN      
      IER = 0     
      do 10 I = 1,N 
         P(I) = 0 
         IP(I) = 0
   10 continue    
      FIRST = 1   
   20 P(FIRST) = FIRST      
      if (IA(FIRST+1)-IA(FIRST).gt.1) GO TO 40  
      if (FIRST.eq.N) GO TO 130       
      IBGN = FIRST+1
      do 30 I = IBGN,N      
         if (P(I).ne.0) GO TO 30      
         FIRST = I
         GO TO 20 
   30 continue    
      GO TO 130   
   40 NEXT = 1    
      LAST = 1    
      IP(1) = FIRST 
   50 K = IP(NEXT)
      CURTYP = P(K) 
      NXTTYP = -CURTYP      
      IBGN = IA(K)
      IEND = IA(K+1)-1      
      if (IBGN.gt.IEND) GO TO 110     
      do 100 I = IBGN,IEND  
         J = JA(I)
         type = P(J)
         if (J.eq.K) GO TO 100
         if (type.eq.NXTTYP) GO TO 100
         if (type.ne.0) GO TO 60      
         LAST = LAST+1      
         IP(LAST) = J       
         P(J) = NXTTYP      
         GO TO 100
   60    if (type.eq.CURTYP) GO TO 160
         if (type*NXTTYP.lt.1) GO TO 80 
         OLD = MIN0(IABS(type),IABS(NXTTYP))    
         YOUNG = MAX0(IABS(type),IABS(NXTTYP))  
         do 70 J = YOUNG,N  
            if (IABS(P(J)).eq.YOUNG) P(J) = ISIGN(OLD,P(J)) 
   70    continue 
         CURTYP = P(K)      
         NXTTYP = -CURTYP   
         GO TO 100
   80    OLD = MIN0(IABS(type),IABS(NXTTYP))    
         YOUNG = MAX0(IABS(type),IABS(NXTTYP))  
         do 90 J = YOUNG,N  
            if (IABS(P(J)).eq.YOUNG) P(J) = ISIGN(OLD,-P(J))
   90    continue 
         CURTYP = P(K)      
         NXTTYP = -CURTYP   
  100 continue    
  110 NEXT = NEXT+1 
      if (NEXT.le.LAST) GO TO 50      
      IBGN = FIRST+1
      do 120 I = IBGN,N     
         if (P(I).ne.0) GO TO 120     
         FIRST = I
         GO TO 20 
  120 continue    
  130 NRED = 0    
      NBLACK = 0  
      do 150 I = 1,N
         if (P(I).lt.0) GO TO 140     
         NRED = NRED+1      
         IP(NRED) = I       
         P(I) = NRED
         GO TO 150
  140    NBLACK = NBLACK+1  
         J = N-NBLACK+1     
         IP(J) = I
         P(I) = J 
  150 continue    
      GO TO 180   
  160 IER = 201   
      if (LEVEL.ge.0) write (NOUT,170)
170 format ('0','*** F A T A L     E R R O R ************'/'0', '    IN ITPACK ROUTINE PRBNDX  '/' ', &
'    RED-BLACK INDEXING NOT POSSIBLE')
  180 continue    
      return      
      END 
      subroutine PRSBLK (NNB,NNR,IA,JA,A,UR,VB) 
      integer IA(NNB+NNR+1),JA(1),NNB,NNR     
      double precision A(1),UR(NNR),VB(NNB)     
      integer I,IBGN,IEND,INR,J,JAJ,NB,NR       
      double precision SUM,URI
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      NB = NNB    
      NR = NNR    
      if (ISYM.eq.0) GO TO 30 
      do 20 I = 1,NB
         INR = I+NR 
         IBGN = IA(INR)     
         IEND = IA(INR+1)-1 
         SUM = VB(I)
         if (IBGN.gt.IEND) GO TO 20   
         do 10 J = IBGN,IEND
            JAJ = JA(J)     
            SUM = SUM-A(J)*UR(JAJ)    
   10    continue 
         VB(I) = SUM
   20 continue    
      return      
   30 do 50 I = 1,NR
         IBGN = IA(I)       
         IEND = IA(I+1)-1   
         if (IBGN.gt.IEND) GO TO 50   
         URI = UR(I)
         do 40 J = IBGN,IEND
            JAJ = JA(J)-NR  
            VB(JAJ) = VB(JAJ)-A(J)*URI
   40    continue 
   50 continue    
      return      
      END 
      subroutine PRSRED (NNB,NNR,IA,JA,A,UB,VR) 
      integer IA(NNR+1),JA(1),NNB,NNR     
      double precision A(1),UB(NNB),VR(NNR)     
      integer IBGN,IEND,II,JAJJ,JJ,NB,NR
      double precision SUM  
      NB = NNB    
      NR = NNR    
      do 20 II = 1,NR       
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         if (IBGN.gt.IEND) GO TO 20   
         SUM = VR(II)       
         do 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)-NR
            SUM = SUM-A(JJ)*UB(JAJJ)  
   10    continue 
         VR(II) = SUM       
   20 continue    
      return      
      END 
      subroutine PSSOR1 (NN,IA,JA,A,U,RHS,FR,BR)
      
      use Global_ITPACK
      
      integer IA(NN+1),JA(ITPACK_NNZ),NN
      double precision A(ITPACK_NNZ),U(NN),RHS(NN),FR(NN),BR(NN)   
      
      integer I,IBGN,IEND,II,JAJJ,JJ,N,NPL1     
      double precision OMM1,SUM,UII   
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      N = NN      
      NPL1 = N+1  
      OMM1 = OMEGA-1.D0     
      if (ISYM.eq.0) GO TO 40 
      do 30 II = 1,N
         BR(II) = U(II)     
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         if (IBGN.gt.IEND) GO TO 20   
         do 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   10    continue 
   20    UII = OMEGA*SUM-OMM1*U(II)   
         FR(II) = UII-U(II) 
         U(II) = UII
   30 continue    
      GO TO 90    
   40 do 80 II = 1,N
         BR(II) = U(II)     
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         SUM = RHS(II)      
         if (IBGN.gt.IEND) GO TO 60   
         do 50 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUM = SUM-A(JJ)*U(JAJJ)   
   50    continue 
   60    UII = OMEGA*SUM-OMM1*U(II)   
         FR(II) = UII-U(II) 
         U(II) = UII
         if (IBGN.gt.IEND) GO TO 80   
         do 70 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            RHS(JAJJ) = RHS(JAJJ)-A(JJ)*UII     
   70    continue 
   80 continue    
   90 do 120 I = 1,N
         II = NPL1-I
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         UII = RHS(II)      
         if (IBGN.gt.IEND) GO TO 110  
         do 100 JJ = IBGN,IEND
            JAJJ = JA(JJ)   
            UII = UII-A(JJ)*U(JAJJ)   
  100    continue 
  110    U(II) = OMEGA*UII-OMM1*U(II) 
         BR(II) = U(II)-BR(II)
  120 continue    
      return      
      END 
      subroutine PSTOP (N,U,DNRM,CCON,IFLAG,Q1) 
      integer N,IFLAG       
      double precision U(N),DNRM,CCON 
      logical :: Q1
      double precision CON,TL,TR,UOLD 
      double precision DDOT2 
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      CON = CCON  
      HALT = .false.
      
      Q1 = .false.
      if (IN.ge.1) GO TO 10 
      Q1 = .false.
      UDNM = 1.D0 
      STPTST = 1.D3 
      if (IFLAG.le.0) return
   10 continue    
      if (Q1) GO TO 20      
      if ((IN.gt.5).and.(MOD(IN,5).ne.0)) GO TO 20
      UOLD = UDNM 
      UDNM = DDOT2(N,U,1,U,1)
      if (UDNM.eq.0.D0) UDNM = 1.D0   
      if ((IN.gt.5).and.(DABS(UDNM-UOLD).le.UDNM*ZETA)) Q1 = .true. 
   20 TR = DSQRT(UDNM)      
      TL = 1.D0   
      if (CON.eq.1.D0) GO TO 40       
      if (IFLAG.eq.2) GO TO 30
      TL = DSQRT(DNRM)      
      TR = TR*(1.D0-CON)    
      GO TO 40    
   30 TL = DSQRT(2.D0*DNRM) 
      TR = TR*(1.D0-CON*CON)
   40 STPTST = TL/TR
      if (TL.ge.TR*ZETA) return       
      HALT = .true. 
      return      
      END 
      double precision function PVTBV (N,IA,JA,A,V)       
      
      use Global_ITPACK
      integer IA(N+1),JA(ITPACK_NNZ),N 
      double precision A(ITPACK_NNZ),V(N)  
      
      integer IBGN,IEND,II,JAJJ,JJ    
      double precision SUM,SUMR       
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      PVTBV = 0.D0
      SUM = 0.D0  
      do 20 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         if (IBGN.gt.IEND) GO TO 20   
         SUMR = 0.D0
         do 10 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            SUMR = SUMR-A(JJ)*V(JAJJ) 
   10    continue 
         SUM = SUM+V(II)*SUMR 
   20 continue    
      if (ISYM.eq.0) SUM = 2.D0*SUM   
      PVTBV = SUM 
      return      
      END 
      subroutine QSORT (NN,KEY,data,ERROR)      
      integer NN,ERROR,KEY(NN)
      double precision data(NN)       
      integer TOP,LEFT,RIGHT,I,J,TINY,V,K,IP1,JM1,LLEN,RLEN,N       
      logical :: DONE
      double precision D    
      integer STKLEN,STACK(30)
      data TINY,STKLEN / 9,30 /       
      N = NN      
      if (N.eq.1) return    
      if (N.le.0) GO TO 240 
      ERROR = 0   
      TOP = 1     
      LEFT = 1    
      RIGHT = N   
      DONE = (N.le.TINY)    
      if (DONE) GO TO 150   
      call IVFILL (STKLEN,STACK,0)    
   10 if (DONE) GO TO 150   
      LFRH2 = (LEFT+RIGHT)/2
      K = KEY(LFRH2)
      D = data(LFRH2)       
      KEY(LFRH2) = KEY(LEFT)
      data(LFRH2) = data(LEFT)
      KEY(LEFT) = K 
      data(LEFT) = D
      if (KEY(LEFT+1).le.KEY(RIGHT)) GO TO 20   
      K = KEY(LEFT+1)       
      D = data(LEFT+1)      
      KEY(LEFT+1) = KEY(RIGHT)
      data(LEFT+1) = data(RIGHT)      
      KEY(RIGHT) = K
      data(RIGHT) = D       
   20 if (KEY(LEFT).le.KEY(RIGHT)) GO TO 30     
      K = KEY(LEFT) 
      D = data(LEFT)
      KEY(LEFT) = KEY(RIGHT)
      data(LEFT) = data(RIGHT)
      KEY(RIGHT) = K
      data(RIGHT) = D       
   30 if (KEY(LEFT+1).le.KEY(LEFT)) GO TO 40    
      K = KEY(LEFT+1)       
      D = data(LEFT+1)      
      KEY(LEFT+1) = KEY(LEFT) 
      data(LEFT+1) = data(LEFT)       
      KEY(LEFT) = K 
      data(LEFT) = D
   40 V = KEY(LEFT) 
      I = LEFT+1  
      J = RIGHT   
   50 continue    
   60 I = I+1     
      if (KEY(I).lt.V) GO TO 60       
   70 J = J-1     
      if (KEY(J).gt.V) GO TO 70       
      if (J.lt.I) GO TO 80  
      K = KEY(I)  
      D = data(I) 
      KEY(I) = KEY(J)       
      data(I) = data(J)     
      KEY(J) = K  
      data(J) = D 
      GO TO 50    
   80 K = KEY(LEFT) 
      D = data(LEFT)
      KEY(LEFT) = KEY(J)    
      data(LEFT) = data(J)  
      KEY(J) = K  
      data(J) = D 
      LLEN = J-LEFT 
      RLEN = RIGHT-I+1      
      if (MAX0(LLEN,RLEN).gt.TINY) GO TO 100    
      if (TOP.eq.1) GO TO 90
      TOP = TOP-2 
      LEFT = STACK(TOP)     
      RIGHT = STACK(TOP+1)  
      GO TO 10    
   90 DONE = .true. 
      GO TO 10    
  100 if (MIN0(LLEN,RLEN).gt.TINY) GO TO 120    
      if (LLEN.gt.RLEN) GO TO 110     
      LEFT = I    
      GO TO 10    
  110 RIGHT = J-1 
      GO TO 10    
  120 if (TOP.ge.STKLEN) GO TO 240    
      if (LLEN.gt.RLEN) GO TO 130     
      STACK(TOP) = I
      STACK(TOP+1) = RIGHT  
      RIGHT = J-1 
      GO TO 140   
  130 STACK(TOP) = LEFT     
      STACK(TOP+1) = J-1    
      LEFT = I    
  140 TOP = TOP+2 
      GO TO 10    
  150 I = N-1     
      LEFT = MAX0(0,N-TINY) 
      K = KEY(N)  
      J = N       
  160 if (I.le.LEFT) GO TO 180
      if (KEY(I).le.K) GO TO 170      
      K = KEY(I)  
      J = I       
  170 I = I-1     
      GO TO 160   
  180 if (J.eq.N) GO TO 190 
      KEY(J) = KEY(N)       
      KEY(N) = K  
      D = data(N) 
      data(N) = data(J)     
      data(J) = D 
  190 I = N-1     
      IP1 = N     
  200 if (KEY(I).le.KEY(IP1)) GO TO 220 
      K = KEY(I)  
      D = data(I) 
      J = IP1     
      JM1 = I     
  210 KEY(JM1) = KEY(J)     
      data(JM1) = data(J)   
      JM1 = J     
      J = J+1     
      if (KEY(J).lt.K) GO TO 210      
      KEY(JM1) = K
      data(JM1) = D 
  220 IP1 = I     
      I = I-1     
      if (I.gt.0) GO TO 200 
  230 return      
  240 ERROR = 1   
      GO TO 230   
      END 
      subroutine SBAGN (N,NZ,IA,JA,A,IWORK,LEVELL,NOUTT,IERR)       
      
      integer NZ,IA(N+1),JA(NZ),IWORK(NZ),N,LEVELL,NOUTT,IERR
      double precision A(NZ)
      
      integer I,IER,J,LEVEL,NOUT,NADD,NADDP1,NOW,NP1,NTO,NTN
      NOW = IA(N+1)-1       
      NADD = NZ-NOW 
      IER = 0     
      LEVEL = LEVELL
      NOUT = NOUTT
      if (NADD.le.0) IER = 703
      if (IER.eq.0) GO TO 20
      if (LEVEL.ge.0) write (NOUT,10) IER       
10 format ('0','*** F A T A L     E R R O R ************'/'0', &
'    IN ITPACK ROUTINE SBAGN   '/' ','    IER = ',I10/' ', '    NZ TOO SMALL - NO ROOM FOR NEW ENTRY')
      GO TO 90    
   20 NTO = NOW   
      NTN = NZ    
      do 30 I = 1,NOW       
         JA(NTN) = JA(NTO)  
         A(NTN) = A(NTO)    
         NTO = NTO-1
         NTN = NTN-1
   30 continue    
      do 40 I = 1,NADD      
         JA(I) = 0
         A(I) = 0.D0
   40 continue    
      NP1 = N+1   
      do 50 I = 1,NP1       
         IA(I) = IA(I)+NADD 
   50 continue    
      NADDP1 = NADD+1       
      do 60 I = NADDP1,NZ   
         IWORK(I) = I+1     
   60 continue    
      do 70 I = 1,NADD      
         IWORK(I) = 0       
   70 continue    
      do 80 I = 1,N 
         J = IA(I+1)-1      
         IWORK(J) = -I      
   80 continue    
      IA(N+1) = NADD
      return      
   90 IERR = IER  
      return      
      END 
      subroutine SBELM (NN,IA,JA,A,RHS,IW,RW,TOL,ISYM,LEVEL,NOUT,IER) 
      
      
      use Global_ITPACK
      
      integer NN,IA(NN+1),JA(ITPACK_NNZ),IW(NN),ISYM,LEVEL,NOUT,IER   
      double precision A(ITPACK_NNZ),RHS(NN),RW(NN),TOL  
      
      integer IBGN,ICNT,IEND,JJ,JJDI,KK,N       
      double precision DI   
      N = NN      
      IER = 0     
      ICNT = 0    
      do 10 II = 1,N
         RW(II) = 0.0D0     
         IW(II) = 0 
   10 continue    
      do 20 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         if (IBGN.gt.IEND) GO TO 140  
         do 20 JJ = IBGN,IEND 
            KK = JA(JJ)     
            if (KK.eq.II) GO TO 20    
            RW(II) = DMAX1(RW(II),DABS(A(JJ)))  
            if (ISYM.ne.0) GO TO 20   
            RW(KK) = DMAX1(RW(KK),DABS(A(JJ)))  
   20 continue    
      do 80 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         do 40 JJ = IBGN,IEND 
            if (JA(JJ).ne.II) GO TO 40
            DI = A(JJ)      
            JJDI = JJ       
            if (DI.gt.0.D0) GO TO 50  
            IER = 101       
            if (LEVEL.ge.0) write (NOUT,30) II,DI 
30       format ('0','*** F A T A L     E R R O R ************'/'0', '    IN ITPACK ROUTINE SBELM   '/' ', &
'    DIAGONAL ELEMENT',I10,' NOT POSITIVE  '/' ', '    CURRENT VALUE = ',D15.8)
            return
   40    continue 
         GO TO 140
   50    continue 
         if (RW(II).ne.0.0D0) GO TO 60
         if (1.0D0/DI.le.TOL) GO TO 70
         GO TO 80 
   60    if (RW(II)/DI.gt.TOL) GO TO 80 
   70    ICNT = ICNT+1      
         IW(II) = II
         RW(II) = DI
         A(JJDI) = 1.0D0    
         RHS(II) = RHS(II)/DI 
   80 continue    
      if (ICNT.eq.0) GO TO 130
      do 120 II = 1,N       
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         if (IW(II).eq.0) GO TO 100   
         do 90 JJ = IBGN,IEND 
            KK = JA(JJ)     
            if (KK.eq.II) GO TO 90    
if ((IW(KK).eq.0).and.(ISYM.eq.0)) RHS(KK) = RHS(KK)-A(JJ)* RHS(II)
            A(JJ) = 0.0D0   
   90    continue 
         GO TO 120
  100    do 110 JJ = IBGN,IEND
            KK = JA(JJ)     
            if (KK.eq.II.or.IW(KK).eq.0) GO TO 110
            RHS(II) = RHS(II)-A(JJ)*RHS(KK)     
            A(JJ) = 0.0D0   
  110    continue 
  120 continue    
  130 return      
  140 continue    
      IER = 102   
      if (LEVEL.ge.0) write (NOUT,150) II       
150 format ('0','*** F A T A L     E R R O R ************'/'0', '    IN ITPACK ROUTINE SBELM   '/' ', &
'    NO DIAGONAL ENTRY IN ROW  ',I10)
      return      
      END 
      subroutine SBEND (N,NZ,IA,JA,A,IWORK)     
      integer N,NZ,IA(N+1),JA(NZ),IWORK(NZ)       
      double precision A(NZ)
integer MAXTOP,NEXT,TOP,IDEG,NULINK,JAJ,HLINK,OHLINK,L,I,LINK, MHLINK
      double precision VAL  
      NEXT = 1    
      TOP = IA(N+1)+1       
      MAXTOP = NZ-IA(N+1)+1 
      do 90 I = 1,N 
         IDEG = 0 
         NULINK = IA(I)     
   10    LINK = NULINK      
         if (LINK.le.0) GO TO 80      
         NULINK = IWORK(LINK) 
         JAJ = JA(LINK)     
         VAL = A(LINK)      
         if (NEXT.ge.TOP.and.LINK.ne.TOP) GO TO 20
         JA(LINK) = 0       
         A(LINK) = 0.0D0    
         IWORK(LINK) = 0    
         if (LINK.eq.TOP) GO TO 60    
         GO TO 70 
   20    IA(I) = LINK       
         HLINK = TOP
   30    HLINK = IWORK(HLINK) 
         if (HLINK.gt.0) GO TO 30     
         MHLINK = -HLINK    
         if (IA(MHLINK).ne.TOP) GO TO 40
         IWORK(LINK) = IWORK(TOP)     
         JA(LINK) = JA(TOP) 
         A(LINK) = A(TOP)   
         IA(MHLINK) = LINK  
         if (NULINK.eq.TOP) NULINK = LINK       
         GO TO 60 
   40    HLINK = IA(MHLINK) 
   50    OHLINK = HLINK     
         HLINK = IWORK(OHLINK)
         if (HLINK.ne.TOP) GO TO 50   
         IWORK(LINK) = IWORK(TOP)     
         JA(LINK) = JA(TOP) 
         A(LINK) = A(TOP)   
         if (OHLINK.ne.LINK) IWORK(OHLINK) = LINK 
         if (NULINK.eq.TOP) NULINK = LINK       
   60    TOP = TOP+1
         if (TOP.ge.MAXTOP) GO TO 70  
         if (IWORK(TOP).ne.0) GO TO 70
         GO TO 60 
   70    JA(NEXT) = JAJ     
         A(NEXT) = VAL      
         NEXT = NEXT+1      
         IDEG = IDEG+1      
         GO TO 10 
   80    IA(I) = IDEG       
   90 continue    
      L = IA(1)+1 
      IA(1) = 1   
      do 100 I = 1,N
         IDEG = IA(I+1)     
         IA(I+1) = L
         L = L+IDEG 
  100 continue    
      return      
      END 
      subroutine SBINI (N,NZ,IA,JA,A,IWORK)     
      integer N,NZ,IA(1),JA(NZ),IWORK(NZ),I     
      double precision A(NZ)
      do 10 I = 1,N 
         IA(I) = -I 
   10 continue    
      IA(N+1) = NZ
      call IVFILL (NZ,JA,0) 
      call IVFILL (NZ,IWORK,0)
      call VFILL (NZ,A,0.D0)
      return      
      END 
subroutine SBSIJ (N,NZ,IA,JA,A,IWORK,II,JJ,VALL,MODE,LEVELL,NOUTT, IERR)
      integer N,NZ,IA(1),JA(NZ),IWORK(NZ),II,JJ,MODE,LEVELL,NOUTT,IERR
      double precision A(NZ),VALL     
      integer LINK,NEXT,NPL1,I,J,LEVEL,NOUT,IER 
      double precision VAL,TEMP       
      I = II      
      J = JJ      
      VAL = VALL  
      LEVEL = LEVELL
      NOUT = NOUTT
      IER = 0     
      if (I.le.0.or.I.gt.N) IER = 701 
      if (J.le.0.or.J.gt.N) IER = 701 
      if (IER.eq.0) GO TO 20
      if (LEVEL.ge.0) write (NOUT,10) IER,I,J   
10 format ('0','*** F A T A L     E R R O R ************'/'0', &
'    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', '    ( ',I10,' , ',I10,' )'/' ', &
'    IMPROPER VALUE FOR I OR J ')
      GO TO 130   
   20 NPL1 = N+1  
      LINK = IA(I)
      if (LINK.gt.0) GO TO 30 
      NEXT = IA(NPL1)       
      if (NEXT.lt.1) GO TO 110
      IA(I) = NEXT
      JA(NEXT) = J
      A(NEXT) = VAL 
      IWORK(NEXT) = -I      
      IA(NPL1) = NEXT-1     
      GO TO 130   
   30 if (JA(LINK).eq.J) GO TO 40     
      if (IWORK(LINK).le.0) GO TO 100 
      LINK = IWORK(LINK)    
      GO TO 30    
   40 IER = 700   
      if (MODE.ge.0) GO TO 60 
      if (LEVEL.ge.1) write (NOUT,50) IER,I,J,A(LINK)     
50 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', &
'    ( ',I10,' , ',I10,' )'/' ', '    ENTRY ALREADY SET AND IS LEFT AS ',D15.8)
      GO TO 130   
   60 if (MODE.ge.1) GO TO 80 
      if (LEVEL.ge.1) write (NOUT,70) IER,I,J,A(LINK),VAL 
70 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', &
'    ( ',I10,' , ',I10,' )'/' ', '    ENTRY ALREADY SET - CURRENT VALUE OF',D15.8/' ', &
'                                RESET TO',D15.8)
      A(LINK) = VAL 
      GO TO 130   
   80 TEMP = A(LINK)+VAL    
      if (LEVEL.ge.1) write (NOUT,90) IER,I,J,A(LINK),TEMP
90 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', &
'    ( ',I10,' , ',I10,' )'/' ', '    ENTRY ALREADY SET - CURRENT VALUE OF',D15.8/' ', &
'                                RESET TO',D15.8)
      A(LINK) = TEMP
      GO TO 130   
  100 NEXT = IA(NPL1)       
      if (NEXT.lt.1) GO TO 110
      IWORK(LINK) = NEXT    
      JA(NEXT) = J
      A(NEXT) = VAL 
      IWORK(NEXT) = -I      
      IA(NPL1) = NEXT-1     
      GO TO 130   
  110 IER = 702   
      if (LEVEL.ge.0) write (NOUT,120) IER      
120 format ('0','*** F A T A L     E R R O R ************'/'0', &
'    IN ITPACK ROUTINE SBSIJ   '/' ','    IER = ',I10/' ', '    NZ TOO SMALL - NO ROOM FOR NEW ENTRY')
  130 continue    
      IERR = IER  
      return      
      END 
      subroutine SCAL (NN,IA,JA,A,RHS,U,D,LEVEL,NOUT,IER) 
      use Global_ITPACK
      
      integer IA(NN+1),JA(ITPACK_NNZ),NN,LEVEL,NOUT,IER     
      double precision A(ITPACK_NNZ),RHS(NN),U(NN),D(NN) 
      integer I,IBGN,IEND,II,IM1,J,JADD,JAJJ,JJ,JJPI,N,NP1
      double precision DI   
      N = NN      
      IER = 0   
      do 80 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         if (IBGN.gt.IEND) GO TO 50   
         do 40 JJ = IBGN,IEND 
            if (JA(JJ).ne.II) GO TO 40
            DI = A(JJ)      
            if (DI.gt.0.D0) GO TO 70  
            if (DI.eq.0.D0) GO TO 20  
            IER = 401       
            if (LEVEL.ge.0) write (NOUT,10) II  
10       format ('0','*** F A T A L     E R R O R ************'/'0', '    IN ITPACK ROUTINE SCAL    '/' ', &
'    DIAGONAL ENTRY IN ROW ',I10,' NEGATIVE')
            return
   20       IER = 401       
            if (LEVEL.ge.0) write (NOUT,30)     
30       format ('0','*** F A T A L     E R R O R ************'/'0', '    IN ITPACK ROUTINE SCAL    '/' ', &
'    DIAGONAL ENTRY IN ROW ',I10,' IS ZERO')
            return
   40    continue 
   50    IER = 402
         if (LEVEL.ge.0) write (NOUT,60) II     
60    format ('0','*** F A T A L     E R R O R ************'/'0', '    IN ITPACK ROUTINE SCAL    '/' ', &
'    NO DIAGONAL ENTRY IN ROW',I10)
         return   
   70    continue 
         DI = DSQRT(DABS(DI)) 
         RHS(II) = RHS(II)/DI 
         U(II) = U(II)*DI   
         D(II) = DI 
   80 continue    
      if (N.eq.1) GO TO 110 
      NP1 = N+1   
      do 100 I = 1,N
         IM1 = I-1
         II = NP1-I 
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         JADD = IBGN+IEND   
         do 90 J = IBGN,IEND
            JJ = JADD-J     
            JJPI = JJ+IM1   
            if (JA(JJ).eq.II) IM1 = I 
            A(JJPI) = A(JJ) 
            JA(JJPI) = JA(JJ) 
   90    continue 
         IA(II+1) = IA(II+1)+I-1      
  100 continue    
  110 IA(1) = IA(1)+N       
      do 140 II = 1,N       
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         DI = D(II) 
         if (IBGN.gt.IEND) GO TO 130  
         do 120 JJ = IBGN,IEND
            JAJJ = JA(JJ)   
            A(JJ) = A(JJ)/(DI*D(JAJJ))
  120    continue 
  130    continue 
         A(II) = DI 
  140 continue    
      return      
      END 
      subroutine SUM3 (N,C1,X1,C2,X2,C3,X3)     
      integer :: I
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      double precision X1(N),X2(N),X3(N),C1,C2,C3 
      if (N.le.0) return    
      if (DABS(C3).eq.0.D0) GO TO 20  
      do 10 I = 1,N 
         X3(I) = C1*X1(I)+C2*X2(I)+C3*X3(I)     
   10 continue    
      return      
   20 do 30 I = 1,N 
         X3(I) = C1*X1(I)+C2*X2(I)    
   30 continue    
      return      
      END 
      double precision function TAU (II)
      integer :: II
      double precision T(8) 
data T(1),T(2),T(3),T(4),T(5),T(6),T(7),T(8) / 1.5D0,1.8D0,1.85D0, 1.9D0,1.94D0,1.96D0,1.975D0,1.985D0 /
      TAU = 1.992D0 
      if (II.le.8) TAU = T(II)
      return      
      END 
      
      function TIMER (TIMDMY) 
      real :: TIMDMY
      external ETIME
      real ETIME, TIMER
      TOTAL = 0
      TIMER = TOTAL
      return
      END 
      
      logical function TSTCHG (IBMTH) 
      integer :: IBMTH
      integer :: IP
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      IP = IN-IS  
      if (IBMTH.eq.2) IP = 2*IP       
      if (IN.eq.0) GO TO 10 
      if (IP.lt.3) GO TO 20 
      QA = DSQRT(DABS(DELNNM/DELSNM)) 
      QT = 2.D0*DSQRT(DABS(RRR**IP))/(1.D0+RRR**IP)       
      if ((QA.ge.1.D0).or.(QA.lt.QT**FF)) GO TO 20
   10 TSTCHG = .true.       
      return      
   20 TSTCHG = .false.      
      return      
      END 
      subroutine UNSCAL (N,IA,JA,A,RHS,U,D)     
      
      use Global_ITPACK
      integer IA(N+1),JA(ITPACK_NNZ),N 
      double precision A(ITPACK_NNZ),RHS(N),U(N),D(N)    
      integer IBGN,IEND,II,INEW,IS,JAJJ,JJ,JJPI 
      double precision DI   
      do 10 II = 1,N
         DI = A(II) 
         U(II) = U(II)/DI   
         RHS(II) = RHS(II)*DI 
         D(II) = DI 
   10 continue    
      do 30 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         if (IBGN.gt.IEND) GO TO 30   
         DI = D(II) 
         do 20 JJ = IBGN,IEND 
            JAJJ = JA(JJ)   
            A(JJ) = A(JJ)*DI*D(JAJJ)  
   20    continue 
   30 continue    
      do 60 II = 1,N
         IBGN = IA(II)      
         IEND = IA(II+1)-1  
         IS = N-II
         INEW = IBGN-IS-1   
         A(INEW) = D(II)**2 
         JA(INEW) = II      
         if (IS.eq.0.or.IBGN.gt.IEND) GO TO 50  
         do 40 JJ = IBGN,IEND 
            JJPI = JJ-IS    
            A(JJPI) = A(JJ) 
            JA(JJPI) = JA(JJ) 
   40    continue 
   50    continue 
         IA(II) = INEW      
   60 continue    
      return      
      END 
      subroutine VEVMW (N,V,W)
      integer :: N
      double precision V(N),W(N)      
      integer I,M,MP1       
      if (N.le.0) return    
      M = MOD(N,4)
      if (M.eq.0) GO TO 20  
      do 10 I = 1,M 
         V(I) = V(I)-W(I)   
   10 continue    
      if (N.lt.4) return    
   20 MP1 = M+1   
      do 30 I = MP1,N,4     
         V(I) = V(I)-W(I)   
         V(I+1) = V(I+1)-W(I+1)       
         V(I+2) = V(I+2)-W(I+2)       
         V(I+3) = V(I+3)-W(I+3)       
   30 continue    
      return      
      END 
      subroutine VEVPW (N,V,W)
      integer :: N
      double precision V(N),W(N)      
      integer I,M,MP1       
      if (N.le.0) return    
      M = MOD(N,4)
      if (M.eq.0) GO TO 20  
      do 10 I = 1,M 
         V(I) = V(I)+W(I)   
   10 continue    
      if (N.lt.4) return    
   20 MP1 = M+1   
      do 30 I = MP1,N,4     
         V(I) = V(I)+W(I)   
         V(I+1) = V(I+1)+W(I+1)       
         V(I+2) = V(I+2)+W(I+2)       
         V(I+3) = V(I+3)+W(I+3)       
   30 continue    
      return      
      END 
      subroutine VFILL (N,V,VAL)      
      integer :: N
      double precision V(N),VAL       
      integer I,M,MP1       
      if (N.le.0) return    
      M = MOD(N,10) 
      if (M.eq.0) GO TO 20  
      do 10 I = 1,M 
         V(I) = VAL 
   10 continue    
      if (N.lt.10) return   
   20 MP1 = M+1   
      do 30 I = MP1,N,10    
         V(I) = VAL 
         V(I+1) = VAL       
         V(I+2) = VAL       
         V(I+3) = VAL       
         V(I+4) = VAL       
         V(I+5) = VAL       
         V(I+6) = VAL       
         V(I+7) = VAL       
         V(I+8) = VAL       
         V(I+9) = VAL       
   30 continue    
      return      
      END 
      subroutine VOUT (N,V,ISWT,NOUTT)
      integer N,ISWT,NOUTT  
      double precision V(N) 
      integer I,J,JM1,K,KUPPER,NOUT   
      NOUT = NOUTT
      KUPPER = MIN0(N,8)    
      if (ISWT.eq.1) write (NOUT,10)  
   10 format (//5X,'RESIDUAL VECTOR') 
      if (ISWT.eq.2) write (NOUT,20)  
   20 format (//5X,'SOLUTION VECTOR') 
      write (NOUT,30) (I,I=1,KUPPER)  
   30 format (10X,8I15)     
      write (NOUT,40)       
   40 format (10X,120('-')/)
      do 60 J = 1,N,8       
         KUPPER = MIN0(J+7,N) 
         JM1 = J-1
         write (NOUT,50) JM1,(V(K),K=J,KUPPER)  
   50    format (4X,I5,'+  ',8D15.5)  
   60 continue    
      return      
      END 
      subroutine WEVMW (N,V,W)
      integer :: N
      double precision V(N),W(N)      
      integer I,M,MP1       
      if (N.le.0) return    
      M = MOD(N,4)
      if (M.eq.0) GO TO 20  
      do 10 I = 1,M 
         W(I) = V(I)-W(I)   
   10 continue    
      if (N.lt.4) return    
   20 MP1 = M+1   
      do 30 I = MP1,N,4     
         W(I) = V(I)-W(I)   
         W(I+1) = V(I+1)-W(I+1)       
         W(I+2) = V(I+2)-W(I+2)       
         W(I+3) = V(I+3)-W(I+3)       
   30 continue    
      return      
      END 
      subroutine ZBRENT (N,TRI,EPS,NSIG,AA,BB,MAXFNN,IER) 
      integer IN,IS,ISYM,ITMAX,LEVEL,NOUT       
      common /ITCOM1/ IN,IS,ISYM,ITMAX,LEVEL,NOUT 
      logical ADAPT,BETADT,CASEII,HALT,PARTAD   
      common /ITCOM2/ ADAPT,BETADT,CASEII,HALT,PARTAD     
double precision BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
common /ITCOM3/ BDELNM,BETAB,CME,DELNNM,DELSNM,FF,GAMMA,OMEGA,QA, QT,RHO,RRR,SIGE,SME,SPECR,SPR,DRELPR,STPTST,UDNM,ZETA
      integer NSIG,MAXFNN,IER 
      double precision TRI(2,1),EPS,AA,BB       
      integer IC,MAXFN      
double precision ZERO,HALF,ONE,THREE,TEN,A,B,T,FA,FB,C,FC,D,E,TOL, RM,S,P,Q,R,RONE,TEMP,DETERM
data ZERO / 0.D0 / ,HALF / 5.D-1 / ,ONE / 1.D0 / ,THREE / 3.D0 / , TEN / 10.D0 /
      A = AA      
      B = BB      
      MAXFN = MAXFNN
      IER = 0     
      T = TEN**(-NSIG)      
      IC = 2      
      FA = DETERM(N,TRI,A)  
      FB = DETERM(N,TRI,B)  
      S = B       
      if (FA*FB.gt.ZERO) GO TO 110    
   10 C = A       
      FC = FA     
      D = B-C     
      E = D       
   20 if (DABS(FC).ge.DABS(FB)) GO TO 30
      A = B       
      B = C       
      C = A       
      FA = FB     
      FB = FC     
      FC = FA     
   30 continue    
      TOL = T*DMAX1(DABS(B),0.1D0)    
      RM = (C-B)*HALF       
      if (DABS(FB).le.EPS) GO TO 80   
      if (DABS(C-B).le.TOL) GO TO 80  
      if (IC.ge.MAXFN) GO TO 90       
      if (DABS(E).lt.TOL) GO TO 60    
      if (DABS(FA).le.DABS(FB)) GO TO 60
      S = FB/FA   
      if (A.ne.C) GO TO 40  
      P = (C-B)*S 
      Q = ONE-S   
      GO TO 50    
   40 Q = FA/FC   
      R = FB/FC   
      RONE = R-ONE
      P = S*((C-B)*Q*(Q-R)-(B-A)*RONE)
      Q = (Q-ONE)*RONE*(S-ONE)
   50 if (P.gt.ZERO) Q = -Q 
      if (P.lt.ZERO) P = -P 
      S = E       
      E = D       
      if (P+P.ge.THREE*RM*Q) GO TO 60 
      if (P+P.ge.DABS(S*Q)) GO TO 60  
      D = P/Q     
      GO TO 70    
   60 E = RM      
      D = E       
   70 A = B       
      FA = FB     
      TEMP = D    
      if (DABS(TEMP).le.HALF*TOL) TEMP = DSIGN(HALF*TOL,RM) 
      B = B+TEMP  
      S = B       
      FB = DETERM(N,TRI,S)  
      IC = IC+1   
      if (FB*FC.le.ZERO) GO TO 20     
      GO TO 10    
   80 A = C       
      MAXFN = IC  
      GO TO 130   
   90 IER = 501   
      A = C       
      MAXFN = IC  
      if (LEVEL.ge.1) write (NOUT,100) MAXFN    
100 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE ZBRENT  '/' ', &
'    ALGORITHM FAILED TO CONVERGE   '/' ','    IN',I6, ' ITERATIONS ')
      GO TO 130   
  110 IER = 502   
      MAXFN = IC  
      if (LEVEL.ge.1) write (NOUT,120)
120 format ('0','*** W A R N I N G ************'/'0', '    IN ITPACK ROUTINE ZBRENT  '/' ', &
'    F(A) AND F(B) HAVE SAME SIGN   ')
  130 continue    
      AA = A      
      BB = B      
      MAXFNN = MAXFN
      return      
      END 
