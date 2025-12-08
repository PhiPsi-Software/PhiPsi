 
SUBROUTINE Matrix_Solve_LSOE_Sparse(Key_Indent,Key_LSOE_Sys,c_Key_SLOE,K_CSR_NNZ,&
                 K_CSR_aa,K_CSR_ja,K_CSR_ia,F,D,n)    


use Global_Float_Type
use Global_Common,only:Keywords_Blank,Key_Cond_Number,Space_4,Space_8
use, intrinsic :: ISO_C_BINDING
#ifdef sffortran
#ifndef macos
#ifndef github
use strumpack
#endif
#endif
#endif

implicit none


integer,intent(in)::n,c_Key_SLOE,Key_LSOE_Sys,K_CSR_NNZ
integer,intent(in)::Key_Indent
real(kind=FT),intent(in)::K_CSR_aa(K_CSR_NNZ)
integer,intent(in)::K_CSR_ja(K_CSR_NNZ)
integer,intent(in)::K_CSR_ia(n+1)


real(kind=FT),intent(in)::F(n)
real(kind=FT),intent(out)::D(n)


#ifdef ifort  
include 'mkl_pardiso.fi'
TYPE(MKL_PARDISO_HANDLE) pt(64)
INTEGER maxfct, mnum, mtype, phase, nrhs, error, msglvl
INTEGER iparm(64)
INTEGER idum(1)
REAL*8 ddum(1),pardiso_b(n)
real(kind=8) D_Estimate(n)
integer IWKSP(3*n)
DATA nrhs /1/, maxfct /1/, mnum /1/
integer c_count
integer i
#endif
#ifdef ifx  
include 'mkl_pardiso.fi'
TYPE(MKL_PARDISO_HANDLE) pt(64)
INTEGER maxfct, mnum, mtype, phase, nrhs, error, msglvl
INTEGER iparm(64)
INTEGER idum(1)
REAL*8 ddum(1),pardiso_b(n)
real(kind=8) D_Estimate(n)
integer IWKSP(3*n)
DATA nrhs /1/, maxfct /1/, mnum /1/
integer c_count
integer i
#endif



#ifdef gfortran
#if defined(notlinux) || defined(github)
INCLUDE 'mpif.h'
INCLUDE 'dmumps_struc.h'
INCLUDE 'smumps_struc.h'
TYPE (DMUMPS_STRUC) mumps_par
TYPE (SMUMPS_STRUC) mumps_par_s
INTEGER IERR_6,IERR
#endif
#endif

#ifdef gfortran
#ifndef macos
#ifndef github
integer Ap(n+1),numeric, symbolic, status_7, sys, filenum
real(kind=FT) control_7(20)
real(kind=FT) info_7(90)
real(kind=FT) Ax(K_CSR_NNZ)
integer Ai(K_CSR_NNZ)
#endif
#endif
#endif


#ifdef gfortran
#ifndef github
real(kind=FT),allocatable::K_Value(:)
integer,allocatable::K_Index(:)
integer,allocatable::K_Ptr(:)
real(kind=FT),allocatable::Cond_A(:)
integer,allocatable::Cond_IRN(:),Cond_JCN(:)
integer IERR_8
real(kind=FT) Condition_Num_Norm2
#endif
#endif

#ifdef gfortran
integer SuperLU_nnz, SuperLU_nrhs, SuperLU_ldb, SuperLU_info, SuperLU_iopt
integer,allocatable::SuperLU_rowind(:)
integer,allocatable::SuperLU_colptr(:)
real(kind=FT),allocatable::SuperLU_value(:)
real(kind=FT),allocatable::SuperLU_b(:)
integer*8 SuperLU_factors
integer cc_count
#endif

#ifdef sffortran
#ifndef macos
#ifndef github
type(STRUMPACK_SparseSolver) :: S
integer, target :: STRUMPACK_n
real(kind=FT),allocatable, target::K_Values(:)
integer,allocatable, target::K_Column_indices(:)
integer,allocatable, target::K_Row_Pointers(:)
integer(c_int) :: STRUMPACK_ierr
real(kind=FT),allocatable, target::x_Values(:),F_Values(:)
#endif 
#endif
#endif

1001 FORMAT(5X,'Condition number is',E12.5)
if(Key_Indent==-1)then
elseif(Key_Indent==0)then
  print *, '        Dimension of Linear System:', n
  if (c_Key_SLOE==1) then
      print *, '        Solver: direct method'
  elseif (c_Key_SLOE==2) then
      print *, '        Solver: Gauss elimination'
  elseif (c_Key_SLOE==3) then
      print *, '        Solver: Pardiso'
  elseif (c_Key_SLOE==4) then
      print *, '        Solver: ITPACK'
  elseif (c_Key_SLOE==5) then
      print *, '        Solver: LAPACK'
  elseif (c_Key_SLOE==6) then
      print *, '        Solver: MUMPS'
  elseif (c_Key_SLOE==7) then
      print *, '        Solver: UMFPACK'
  elseif (c_Key_SLOE==8) then
      print *, '        Solver: Lis'
  elseif (c_Key_SLOE==9) then
      print *, '        Solver: SuperLU'
  elseif (c_Key_SLOE==10) then
      print *, '        Solver: y12m'
  elseif (c_Key_SLOE==12) then
      print *, '        Solver: STRUMPACK'
  endif
elseif(Key_Indent==1)then
  print *, '       Dimension of Linear System:', n
  if (c_Key_SLOE==1) then
      print *, '       Solver: direct method'
  elseif (c_Key_SLOE==2) then
      print *, '       Solver: Gauss elimination'
  elseif (c_Key_SLOE==3) then
      print *, '       Solver: Pardiso'
  elseif (c_Key_SLOE==4) then
      print *, '       Solver: ITPACK'
  elseif (c_Key_SLOE==5) then
      print *, '       Solver: LAPACK'
  elseif (c_Key_SLOE==6) then
      print *, '       Solver: MUMPS'
  elseif (c_Key_SLOE==7) then
      print *, '       Solver: UMFPACK'
  elseif (c_Key_SLOE==8) then
      print *, '       Solver: Lis'
  elseif (c_Key_SLOE==9) then
      print *, '       Solver: SuperLU'
  elseif (c_Key_SLOE==10) then
      print *, '       Solver: y12m'
  elseif (c_Key_SLOE==12) then
      print *, '       Solver: STRUMPACK'
  endif
elseif(Key_Indent==2)then
  print *, '      Dimension of Linear System:', n
  if (c_Key_SLOE==1) then
      print *, '      Solver: direct method'
  elseif (c_Key_SLOE==2) then
      print *, '      Solver: Gauss elimination'
  elseif (c_Key_SLOE==3) then
      print *, '      Solver: Pardiso'
  elseif (c_Key_SLOE==4) then
      print *, '      Solver: ITPACK'
  elseif (c_Key_SLOE==5) then
      print *, '      Solver: LAPACK'
  elseif (c_Key_SLOE==6) then
      print *, '      Solver: MUMPS'
  elseif (c_Key_SLOE==7) then
      print *, '      Solver: UMFPACK'
  elseif (c_Key_SLOE==8) then
      print *, '      Solver: Lis'
  elseif (c_Key_SLOE==9) then
      print *, '      Solver: SuperLU'
  elseif (c_Key_SLOE==10) then
      print *, '      Solver: y12m'
  elseif (c_Key_SLOE==12) then
      print *, '      Solver: STRUMPACK'
  endif

elseif(Key_Indent==4)then
  print *, '        Dimension of Linear System:', n
  if (c_Key_SLOE==1) then
      print *, '        Solver: direct method'
  elseif (c_Key_SLOE==2) then
      print *, '        Solver: Gauss elimination'
  elseif (c_Key_SLOE==3) then
      print *, '        Solver: Pardiso'
  elseif (c_Key_SLOE==4) then
      print *, '        Solver: ITPACK'
  elseif (c_Key_SLOE==5) then
      print *, '        Solver: LAPACK'
  elseif (c_Key_SLOE==6) then
      print *, '        Solver: MUMPS'
  elseif (c_Key_SLOE==7) then
      print *, '        Solver: UMFPACK'
  elseif (c_Key_SLOE==8) then
      print *, '        Solver: Lis'
  elseif (c_Key_SLOE==9) then
      print *, '        Solver: SuperLU'
  elseif (c_Key_SLOE==10) then
      print *, '        Solver: y12m'
  elseif (c_Key_SLOE==12) then
      print *, '        Solver: STRUMPACK'
  endif
elseif(Key_Indent==7)then
  print *, '               Dimension of Linear System:', n
  if (c_Key_SLOE==1) then
      print *, '               Solver: direct method'
  elseif (c_Key_SLOE==2) then
      print *, '               Solver: Gauss elimination'
  elseif (c_Key_SLOE==3) then
      print *, '               Solver: Pardiso'
  elseif (c_Key_SLOE==4) then
  elseif (c_Key_SLOE==5) then
      print *, '               Solver: LAPACK'
  elseif (c_Key_SLOE==6) then
      print *, '               Solver: MUMPS'
  elseif (c_Key_SLOE==7) then
      print *, '               Solver: UMFPACK'
  elseif (c_Key_SLOE==8) then
      print *, '               Solver: Lis'
  elseif (c_Key_SLOE==9) then
      print *, '               Solver: SuperLU'
  elseif (c_Key_SLOE==10) then
      print *, '               Solver: y12m'
  elseif (c_Key_SLOE==12) then
      print *, '               Solver: STRUMPACK'
  endif
endif
D(1:n) =ZR

select case(c_Key_SLOE)
case(3)
#if defined(ifort)
    if(Key_Indent==0)then
    elseif (Key_Indent==1)then
    endif
    DO i = 1, 64
        iparm(i) = 0
    END DO
    iparm(1) = 1
    iparm(2) = 2
    iparm(3) = 4
    iparm(4) = 0
    iparm(5) = 0
    iparm(6) = 0
    iparm(8) = 2
    iparm(10) = 13
    iparm(11) = 1
    iparm(13) = 0
    iparm(14) = 0
    iparm(18) = -1
    iparm(19) = -1
    iparm(20) = 0
    error = 0
    msglvl = 0


    mtype = -2


    DO i = 1, 64
    pt(i)%DUMMY = 0
    END DO
    phase = 11
    CALL pardiso (pt, maxfct, mnum, mtype, phase, n, & 
                K_CSR_aa,K_CSR_ia,K_CSR_ja, & 
                idum, nrhs, iparm, msglvl, ddum, ddum, error)

    if(Key_Indent==0)then
    print *,Space_4//'Reordering completed... '
    elseif (Key_Indent==1)then
    print *,Space_8//'Reordering completed... '
    endif
    IF (error .NE. 0) THEN
    WRITE(*,*) 'The following ERROR was detected: ', error
    STOP 1
    END IF

    IF (error .NE. 0) THEN
    WRITE(*,*) 'The following ERROR was detected in Pardiso solver: ', error
    if (error== -4)then
          print *,"Zero pivot, numerical factorization or iterative refinement problem."
          print *,"If the error appears during the solution phase, try to change the   "
          print *,"pivoting perturbation (iparm(10)) and also increase the number of   "
          print *,"iterative refinement steps. If it does not help, consider changing  "
          print *,"the scaling, matching and pivoting options (iparm(11), iparm(13),   "
          print *,"iparm(21))                                                          "
    endif    
    STOP
    END IF
    phase = 22
    CALL pardiso (pt, maxfct, mnum, mtype, phase,n,  & 
                K_CSR_aa,K_CSR_ia, & 
                K_CSR_ja,              & 
                idum, nrhs, iparm, msglvl, ddum, ddum, error)

    if(Key_Indent==0)then
    print *,Space_4//'Factorization completed...'
    elseif (Key_Indent==1)then
    print *,Space_8//'Factorization completed...'
    endif
    IF (error .NE. 0) THEN
    WRITE(*,*) 'The following ERROR was detected: ', error
    STOP 1
    END IF

    IF (error .NE. 0) THEN
    WRITE(*,*) 'The following ERROR was detected in Pardiso solver: ', error
    if (error== -4)then
          print *,"Zero pivot, numerical factorization or iterative refinement problem."
          print *,"If the error appears during the solution phase, try to change the   "
          print *,"pivoting perturbation (iparm(10)) and also increase the number of   "
          print *,"iterative refinement steps. If it does not help, consider changing  "
          print *,"the scaling, matching and pivoting options (iparm(11), iparm(13),   "
          print *,"iparm(21))                                                          "
    endif    
    STOP
    END IF
    phase = 33
    pardiso_b = F
    CALL pardiso (pt, maxfct, mnum, mtype, phase, n, & 
                K_CSR_aa,K_CSR_ia, & 
                K_CSR_ja,              & 
                idum, nrhs, iparm, msglvl, pardiso_b, D, error)


    phase = -1
    call pardiso(pt, maxfct, mnum, mtype, phase, n, &
             ddum, idum, idum, idum, nrhs, &
             iparm, msglvl, ddum, ddum, error)
    if(Key_Indent==0)then
    print *,Space_4//'Solving completed...'
    elseif (Key_Indent==1)then
    print *,Space_8//'Solving completed...'
    endif
#endif      

#if defined(ifx)
    if(Key_Indent==0)then
    elseif (Key_Indent==1)then
    endif
    DO i = 1, 64
    iparm(i) = 0
    END DO
    iparm(1) = 1
    iparm(2) = 2
    iparm(3) = 4
    iparm(4) = 0
    iparm(5) = 0
    iparm(6) = 0
    iparm(8) = 2
    iparm(10) = 13
    iparm(11) = 1
    iparm(13) = 0
    iparm(14) = 0
    iparm(18) = -1
    iparm(19) = -1
    iparm(20) = 0
    error = 0
    msglvl = 0


    mtype = -2


    DO i = 1, 64
    pt(i)%DUMMY = 0
    END DO
    phase = 11
    CALL pardiso (pt, maxfct, mnum, mtype, phase, n, & 
                K_CSR_aa,K_CSR_ia,K_CSR_ja, & 
                idum, nrhs, iparm, msglvl, ddum, ddum, error)

    if(Key_Indent==0)then
    print *,Space_4//'Reordering completed... '
    elseif (Key_Indent==1)then
    print *,Space_8//'Reordering completed... '
    endif
    IF (error .NE. 0) THEN
    WRITE(*,*) 'The following ERROR was detected: ', error
    STOP 1
    END IF

    IF (error .NE. 0) THEN
    WRITE(*,*) 'The following ERROR was detected in Pardiso solver: ', error
    if (error== -4)then
          print *,"Zero pivot, numerical factorization or iterative refinement problem."
          print *,"If the error appears during the solution phase, try to change the   "
          print *,"pivoting perturbation (iparm(10)) and also increase the number of   "
          print *,"iterative refinement steps. If it does not help, consider changing  "
          print *,"the scaling, matching and pivoting options (iparm(11), iparm(13),   "
          print *,"iparm(21))                                                          "
    endif    
    STOP
    END IF
    phase = 22
    CALL pardiso (pt, maxfct, mnum, mtype, phase,n,  & 
                K_CSR_aa,K_CSR_ia, & 
                K_CSR_ja,              & 
                idum, nrhs, iparm, msglvl, ddum, ddum, error)

    if(Key_Indent==0)then
    print *,Space_4//'Factorization completed...'
    elseif (Key_Indent==1)then
    print *,Space_8//'Factorization completed...'
    endif
    IF (error .NE. 0) THEN
    WRITE(*,*) 'The following ERROR was detected: ', error
    STOP 1
    END IF

    IF (error .NE. 0) THEN
    WRITE(*,*) 'The following ERROR was detected in Pardiso solver: ', error
    if (error== -4)then
          print *,"Zero pivot, numerical factorization or iterative refinement problem."
          print *,"If the error appears during the solution phase, try to change the   "
          print *,"pivoting perturbation (iparm(10)) and also increase the number of   "
          print *,"iterative refinement steps. If it does not help, consider changing  "
          print *,"the scaling, matching and pivoting options (iparm(11), iparm(13),   "
          print *,"iparm(21))                                                          "
    endif    
    STOP
    END IF
    phase = 33
    pardiso_b = F
    CALL pardiso (pt, maxfct, mnum, mtype, phase, n, & 
                K_CSR_aa,K_CSR_ia, & 
                K_CSR_ja,              & 
                idum, nrhs, iparm, msglvl, pardiso_b, D, error)


    phase = -1
    call pardiso(pt, maxfct, mnum, mtype, phase, n, &
             ddum, idum, idum, idum, nrhs, &
             iparm, msglvl, ddum, ddum, error)
    if(Key_Indent==0)then
    print *,Space_4//'Solving completed...'
    elseif (Key_Indent==1)then
    print *,Space_8//'Solving completed...'
    endif
#endif   
case(6)
#ifdef gfortran
#if defined(notlinux) || defined(github)
#ifndef SP
    CALL MPI_INIT(IERR_6)
    mumps_par%COMM = MPI_COMM_WORLD
    mumps_par%JOB = -1
    mumps_par%SYM = 0
    mumps_par%PAR = 1
    CALL DMUMPS(mumps_par)
    mumps_par%ICNTL(1) = 6
    mumps_par%ICNTL(2) = 0
    mumps_par%ICNTL(3) = 0
    mumps_par%ICNTL(4) = 2



    mumps_par%N = n
    mumps_par%NZ= K_CSR_NNZ
    ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
    ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
    ALLOCATE( mumps_par%A( mumps_par%NZ ) )
    ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
    call csrcoo(n,3,K_CSR_NNZ,K_CSR_aa,K_CSR_ja,K_CSR_ia,&
           K_CSR_NNZ,mumps_par%A,mumps_par%IRN,mumps_par%JCN,IERR_6)

    mumps_par%RHS(1:n) = F(1:n)
    mumps_par%JOB = 6
    CALL DMUMPS(mumps_par)
    IF (mumps_par%MYID .eq. 0 ) THEN
        D(1:n) = mumps_par%RHS(1:n)
    END IF
    IF (mumps_par%MYID .eq. 0 )THEN
        DEALLOCATE( mumps_par%IRN )
        DEALLOCATE( mumps_par%JCN )
        DEALLOCATE( mumps_par%A   )
        DEALLOCATE( mumps_par%RHS )
    END IF
    mumps_par%JOB = -2
    CALL DMUMPS(mumps_par)
    CALL MPI_FINALIZE(IERR_6)
#endif

#ifdef SP
    CALL MPI_INIT(IERR_6)
    mumps_par_s%COMM = MPI_COMM_WORLD
    mumps_par_s%JOB = -1

    mumps_par_s%SYM = 0
    mumps_par_s%PAR = 1
    CALL SMUMPS(mumps_par_s)
    mumps_par_s%ICNTL(1) = 6
    mumps_par_s%ICNTL(2) = 0
    mumps_par_s%ICNTL(3) = 0
    mumps_par_s%ICNTL(4) = 2

    mumps_par_s%N = n
    mumps_par_s%NZ= K_CSR_NNZ
    ALLOCATE( mumps_par_s%IRN ( mumps_par_s%NZ ) )
    ALLOCATE( mumps_par_s%JCN ( mumps_par_s%NZ ) )
    ALLOCATE( mumps_par_s%A( mumps_par_s%NZ ) )
    ALLOCATE( mumps_par_s%RHS ( mumps_par_s%N  ) )
    call csrcoo(n,3,K_CSR_NNZ,K_CSR_aa,K_CSR_ja,K_CSR_ia,&
           K_CSR_NNZ,mumps_par_s%A,mumps_par_s%IRN,mumps_par_s%JCN,IERR_6)
           mumps_par_s%RHS(1:n) = F(1:n)
    mumps_par_s%JOB = 6
    CALL SMUMPS(mumps_par_s)
    IF (mumps_par_s%MYID .eq. 0 ) THEN
    D(1:n) = mumps_par_s%RHS(1:n)
    END IF
    IF (mumps_par_s%MYID .eq. 0 )THEN
    DEALLOCATE( mumps_par_s%IRN )
    DEALLOCATE( mumps_par_s%JCN )
    DEALLOCATE( mumps_par_s%A   )
    DEALLOCATE( mumps_par_s%RHS )
    END IF
    mumps_par_s%JOB = -2
    CALL SMUMPS(mumps_par_s)
    CALL MPI_FINALIZE(IERR_6)
#endif
#endif
#endif

case(7)
#ifdef gfortran
#ifdef notlinux
#ifndef macos
#ifndef github
    control_7(1:20) = ZR
    info_7(1:90)    = ZR


    call csrcsc(n,1,1,K_CSR_aa(1:K_CSR_NNZ),K_CSR_ja(1:K_CSR_NNZ),K_CSR_ia(1:n+1),&
                Ax(1:K_CSR_NNZ),Ai(1:K_CSR_NNZ),Ap(1:n+1))
    Ap(1:n+1)  = Ap(1:n+1)  -1
    Ai(1:K_CSR_NNZ)  = Ai(1:K_CSR_NNZ)  -1





    call umf4def (control_7)
    control_7(1) = 1
    call umf4pcon (control_7)
    call umf4sym (n, n, Ap, Ai, Ax,symbolic,control_7, info_7)
    if (info_7(1) .lt. 0) then
    print *, 'Error occurred in umf4sym: ', info_7(1)
    stop
    endif
    call umf4num (Ap, Ai, Ax, symbolic, numeric, control_7,info_7)

    if (info_7(1) .lt. 0) then
    print *, 'Error occurred in umf4num: ', info_7 (1)
    stop
    endif
    filenum = 0
    call umf4ssym (symbolic, filenum, status_7)
    if (status_7 .lt. 0) then
    print *, 'Error occurred in umf4ssym: ', status_7
    stop
    endif

    call umf4snum (numeric, filenum, status_7)
    if (status_7 .lt. 0) then
    print *, 'Error occurred in umf4snum: ', status_7
    stop
    endif

    call umf4fsym (symbolic)

    call umf4fnum (numeric)



    call umf4lnum (numeric, filenum, status_7)
    if (status_7 .lt. 0) then
    print *, 'Error occurred in umf4lnum: ', status_7
    stop
    endif

    sys = 0
    call umf4sol (sys, D, F, numeric, control_7, info_7)
    if (info_7 (1) .lt. 0) then
    print *, 'Error occurred in umf4sol: ', info_7 (1)
    stop
    endif
    call umf4fnum (numeric)
#endif
#endif
#endif
#endif


case(8)
#ifdef gfortran
#ifdef notlinux
#ifndef github
    ALLOCATE(K_Value(K_CSR_NNZ))
    ALLOCATE(K_Index(K_CSR_NNZ))
    ALLOCATE(K_Ptr(n+1))

    K_Value(1:K_CSR_NNZ) =K_CSR_aa(1:K_CSR_NNZ)
    K_Index(1:K_CSR_NNZ) = K_CSR_ja(1:K_CSR_NNZ)-1

    K_Ptr(1:n)     = K_CSR_ia(1:n)-1
    K_Ptr(1)       = 0
    K_Ptr(n+1)     = K_CSR_NNZ
    call Matrix_Solve_LSOE_Lis(Key_Indent,K_CSR_NNZ,n,K_Ptr,K_Index,K_Value,F,D)
    if(Key_Cond_Number==1)then
        print *,'    Calculating the condition number...... '
        print *,'    Note: maybe time consuming!            '
        print *,'    Calculator: MC75 of HSL Package        '
        print *,'    http://www.hsl.rl.ac.uk/catalogue/     '
        ALLOCATE(Cond_A(K_CSR_NNZ))
        ALLOCATE(Cond_IRN(K_CSR_NNZ))
        ALLOCATE(Cond_JCN(K_CSR_NNZ))
        call csrcoo(n,3,K_CSR_NNZ,K_CSR_aa,K_CSR_ja,K_CSR_ia,   &
                   K_CSR_NNZ,Cond_A,Cond_IRN,Cond_JCN,IERR_8)
        call Matrix_Condition_Number(n,K_CSR_NNZ,&
                         Cond_A(1:K_CSR_NNZ),Cond_IRN(1:K_CSR_NNZ),&
                         Cond_JCN(1:K_CSR_NNZ),Condition_Num_Norm2)
        write (*,1001) Condition_Num_Norm2
        DEALLOCATE(Cond_A)
        DEALLOCATE(Cond_IRN)
        DEALLOCATE(Cond_JCN)
    endif
    DEALLOCATE(K_Value)
    DEALLOCATE(K_Index)
    DEALLOCATE(K_Ptr)
#endif
#endif
#endif

case(9)
#ifdef gfortran
#ifndef macos
    SuperLU_nnz = K_CSR_NNZ
    allocate(SuperLU_rowind(SuperLU_nnz))
    allocate(SuperLU_value(SuperLU_nnz))
    allocate(SuperLU_colptr(n+1))
    allocate(SuperLU_b(n))

    call CSR_Matrix_Transpose_Real(n, n,SuperLU_nnz,&
                K_CSR_aa(1:SuperLU_nnz),K_CSR_ja(1:SuperLU_nnz),K_CSR_ia(1:n+1), &
                SuperLU_value(1:SuperLU_nnz),SuperLU_rowind(1:SuperLU_nnz),SuperLU_colptr(1:n+1))

    SuperLU_b = F
    SuperLU_ldb = n
    SuperLU_nrhs = 1
    SuperLU_iopt = 1
    call c_fortran_dgssv(SuperLU_iopt,n,  &
                       SuperLU_nnz,SuperLU_nrhs,SuperLU_value,&
                       SuperLU_rowind,SuperLU_colptr,   &
                       SuperLU_b,SuperLU_ldb,SuperLU_factors,SuperLU_info)
    if (SuperLU_info .eq. 0) then
    if(Key_Indent==0)then
          print *,"    Factorization succeeded."
    elseif(Key_Indent==2)then
          print *,"      Factorization succeeded."
    elseif(Key_Indent==4)then
          print *,"        Factorization succeeded."
    elseif(Key_Indent==7)then
          print *,"           Factorization succeeded."
    endif
    else
    if(Key_Indent==0)then
          print *,"    Linear solver SuperLU failed in step 1: error code-",SuperLU_info
    elseif(Key_Indent==2)then
          print *,"      Linear solver SuperLU failed in step 1: error code-",SuperLU_info
    elseif(Key_Indent==4)then
          print *,"        Linear solver SuperLU failed in step 1: error code-",SuperLU_info
    elseif(Key_Indent==7)then
          print *,"           Linear solver SuperLU failed in step 1: error code-",SuperLU_info
    endif
    call Warning_Message('S',Keywords_Blank)
    endif
    SuperLU_iopt = 2
    call c_fortran_dgssv(SuperLU_iopt,n,SuperLU_nnz,SuperLU_nrhs,SuperLU_value, &
                       SuperLU_rowind,SuperLU_colptr,   &
                       SuperLU_b,SuperLU_ldb,SuperLU_factors,SuperLU_info)
    D=SuperLU_b
    if (SuperLU_info .eq. 0) then
    if(Key_Indent==0)then
          print *,"    Solve succeeded."
    elseif(Key_Indent==2)then
          print *,"      Solve succeeded."
    elseif(Key_Indent==4)then
          print *,"        Solve succeeded."
    elseif(Key_Indent==7)then
          print *,"           Solve succeeded."
    endif
    else
    if(Key_Indent==0)then
          print *,"    Linear solver SuperLU failed in step 2: error code-",SuperLU_info
    elseif(Key_Indent==2)then
          print *,"      Linear solver SuperLU failed in step 2: error code-",SuperLU_info
    elseif(Key_Indent==4)then
          print *,"        Linear solver SuperLU failed in step 2: error code-",SuperLU_info
    elseif(Key_Indent==7)then
          print *,"           Linear solver SuperLU failed in step 2: error code-",SuperLU_info
    endif
    call Warning_Message('S',Keywords_Blank)
    endif
    SuperLU_iopt = 3
    call c_fortran_dgssv(SuperLU_iopt,n,SuperLU_nnz,SuperLU_nrhs,SuperLU_value, &
                       SuperLU_rowind,SuperLU_colptr,   &
                       SuperLU_b,SuperLU_ldb,SuperLU_factors,SuperLU_info)
#endif  
#endif        
case(10)
    print *, 'Error :: Solver 10 (y12m) not available for sparse LSOE!'
    call Warning_Message('S',Keywords_Blank)     
    
case(12)
#ifdef sffortran
#ifndef macos
#ifndef github
    1000 FORMAT('###########################################################')
    1100 FORMAT('#              < LSOE Iterative Solver STRUMPACK >         ')
    1300 FORMAT('#              < Exit Iterative Solver STRUMPACK >         ')
    1200 FORMAT('-----------------------------------------------------------')
    WRITE(*,1000)
    WRITE(*,1100)
    WRITE(*,1200)

    call STRUMPACK_init_mt(S, STRUMPACK_DOUBLE, STRUMPACK_MT, 0, c_null_ptr, 1)


    call STRUMPACK_set_matching(S, STRUMPACK_MATCHING_MAX_SMALLEST_DIAGONAL)

    call STRUMPACK_set_reordering_method(S, STRUMPACK_METIS);

    call STRUMPACK_set_compression(S, STRUMPACK_BLR);  

    call STRUMPACK_set_compression_leaf_size(S, 64)
    call STRUMPACK_set_compression_rel_tol(S, dble(1.e-2))

    call STRUMPACK_set_compression_min_sep_size(S, 300)



    ALLOCATE(K_Values(K_CSR_NNZ))
    ALLOCATE(K_Column_indices(K_CSR_NNZ))
    ALLOCATE(K_Row_Pointers(n+1))

    K_Values(1:K_CSR_NNZ) =K_CSR_aa(1:K_CSR_NNZ)
    K_Column_indices(1:K_CSR_NNZ) = K_CSR_ja(1:K_CSR_NNZ)-1

    K_Row_Pointers(1:n)     = K_CSR_ia(1:n)-1
    K_Row_Pointers(1)       = 0
    K_Row_Pointers(n+1)     = K_CSR_NNZ



    ALLOCATE(x_Values(n))
    ALLOCATE(F_Values(n))
    x_Values = ZR
    F_Values = F

    STRUMPACK_n = n

    call STRUMPACK_set_csr_matrix(S, c_loc(STRUMPACK_n), c_loc(K_Row_Pointers), &
                                c_loc(K_Column_indices), c_loc(K_Values), 1)
                                
    STRUMPACK_ierr = STRUMPACK_reorder(S)


    STRUMPACK_ierr = STRUMPACK_solve(S, c_loc(F_Values), c_loc(x_Values), 0)

    D = x_Values

    DEALLOCATE(K_Values)
    DEALLOCATE(K_Column_indices)
    DEALLOCATE(K_Row_Pointers)
    DEALLOCATE(x_Values)
    DEALLOCATE(F_Values)

    WRITE(*,1200)
    WRITE(*,1300)
    WRITE(*,1000)
#endif 
#endif
#endif
end select

RETURN
END SUBROUTINE Matrix_Solve_LSOE_Sparse