 
SUBROUTINE Matrix_Solve_LSOE(Key_Indent,Key_LSOE_Sys,c_Key_SLOE,K,F,D,n)

use Global_Float_Type
use Global_Common,only:Keywords_Blank,Key_Cond_Number,Key_Eigenvalue,Key_Determinant,Space_4,Space_8
use Global_ITPACK
use, intrinsic :: ISO_C_BINDING
#ifdef sffortran
#ifndef macos
#ifndef github
use strumpack
#endif
#endif  
#endif
implicit none
integer,intent(in)::n,c_Key_SLOE,Key_LSOE_Sys
integer,intent(in)::Key_Indent
real(kind=FT),intent(in)::K(n,n),F(n)
real(kind=FT),intent(out)::D(n)
integer i,j,NZ_NUM,IERR
real(kind=FT)K_det

real(kind=FT),ALLOCATABLE::inv_K(:,:)

real(kind=FT),ALLOCATABLE::LIN(:,:)
logical Yes_SINGUL

#if (defined(ifort)) 
include 'mkl_pardiso.fi'
TYPE(MKL_PARDISO_HANDLE) pt(64)
INTEGER maxfct, mnum, mtype, phase, nrhs, error, msglvl
INTEGER iparm(64)
INTEGER idum(1)
REAL*8 ddum(1),pardiso_b(n)
integer I_Sparse_K(n*n),J_Sparse_K(n*n)
real(kind=FT) Sparse_K(n*n),D_Estimate(n)
integer I_Map_Sparse_K(n+1),IWKSP(3*n)
DATA nrhs /1/, maxfct /1/, mnum /1/
integer c_count
#endif
#if (defined(ifx)) 
include 'mkl_pardiso.fi'
TYPE(MKL_PARDISO_HANDLE) pt(64)
INTEGER maxfct, mnum, mtype, phase, nrhs, error, msglvl
INTEGER iparm(64)
INTEGER idum(1)
REAL*8 ddum(1),pardiso_b(n)
integer I_Sparse_K(n*n),J_Sparse_K(n*n)
real(kind=FT) Sparse_K(n*n),D_Estimate(n)
integer I_Map_Sparse_K(n+1),IWKSP(3*n)
DATA nrhs /1/, maxfct /1/, mnum /1/
integer c_count
#endif

integer,ALLOCATABLE::ITPACK_J_Sparse_K(:)
real(kind=FT),ALLOCATABLE::ITPACK_Sparse_K(:)
real(kind=FT),ALLOCATABLE::ITPACK_D_Estimate(:)
integer,ALLOCATABLE:: ITPACK_I_Map_Sparse_K(:)
integer ITPACK_c_count
integer,ALLOCATABLE::ITPACK_IWKSP(:)
integer ITPACK_IPARM(12),ITPACK_ITMAX,ITPACK_NW
parameter(ITPACK_ITMAX=90000)
real(kind=FT) ITPACK_RPARM(12)
real(kind=FT),ALLOCATABLE::ITPACK_WKSP(:)

#ifndef Silverfrost
integer info,ipiv(n),lda
real(kind=FT),ALLOCATABLE::K_inv(:,:)
real(kind=FT) Lapack_work(n)
real(kind=FT) Norm_2_K,Norm_2_K_inv
real(kind=FT) Cond_K
#endif

#ifdef gfortran
#if defined(notlinux) || defined(github)
#ifndef sffortranmpi
INCLUDE 'mpif.h'
INCLUDE 'smumps_struc.h'
INCLUDE 'dmumps_struc.h'
TYPE (SMUMPS_STRUC) s_mumps_par
TYPE (DMUMPS_STRUC) d_mumps_par
INTEGER IERR_6
#endif
#endif
#endif

#ifdef gfortran
integer Ap(n+1),numeric, symbolic, status_7, sys, filenum
real(kind=FT) control_7(20), info_7(90)
real(kind=FT),allocatable:: Ax(:)
integer,allocatable:: Ai(:)
#endif

#ifdef gfortran
real(kind=FT),allocatable::K_Value(:)
integer,allocatable::K_Index(:)
integer,allocatable::K_Ptr(:)
real(kind=FT),allocatable::Cond_A(:)
integer,allocatable::Cond_IRN(:),Cond_JCN(:)
real(kind=FT) Condition_Num_Norm2
#endif

#ifdef gfortran
#ifndef macos
integer SuperLU_nnz, SuperLU_nrhs, SuperLU_ldb, SuperLU_info, SuperLU_iopt
integer,allocatable::SuperLU_rowind(:)
integer SuperLU_colptr(n+1)
real(kind=FT),allocatable::SuperLU_value(:)
real(kind=FT) SuperLU_b(n)
integer*8 SuperLU_factors
integer cc_count
#endif
#endif

#ifdef cbfortran
#else
integer y12m_NN,y12m_NN1,y12m_IHA
real(kind=FT),allocatable::K_Value_y12m(:)
integer,allocatable::K_Index_y12m(:),K_Index_Row_y12m(:)
real(kind=FT) PIVOT(n)
integer HA(n,11),IFLAG(10),IFAIL
real(kind=FT) AFLAG(8)
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
  
#ifndef Silverfrost
integer det_INFO
#endif

#ifndef Silverfrost
real(kind=FT) c_WR(n),c_WI(n)
real(kind=FT),allocatable::K_tem(:,:)
#endif

1001 FORMAT(5X,'Condition number is',E12.5)

if (Key_Determinant==1) then
#ifndef Silverfrost
  print *,'    Attention :: Calculating Determinant of K using LAPACK...'
  print *,'                 Could be time consuming!'
  call Matrix_Det_by_LAPACK(n,K,K_det,det_INFO)
#endif
endif

if (Key_Eigenvalue==1) then
#ifndef Silverfrost
  print *,'    Attention :: Calculating Eigenvalue of K using LAPACK...'
  print *,'                 Could be very time consuming!'
  allocate(K_tem(n,n))
  K_tem = K
  call Matrix_Eigenvalues_by_LAPACK(n,K_tem,c_WR,c_WI)
  print *,'                 Minval of real parts:', minval(c_WR)
  print *,'                 Maxval of real parts:', maxval(c_WR)
  print *,'                 Minval of imag parts:', minval(c_WI)
  print *,'                 Maxval of imag parts:', maxval(c_WI)
  deallocate(K_tem)
#endif
endif

if(Key_Indent==-1)then
elseif(Key_Indent==0)then
  print *,'    Dimension of Linear System:', n
  if (c_Key_SLOE==1) then
      print *, '    Solver: direct method (inverse matrix)' 
  elseif (c_Key_SLOE==2) then
      print *, '    Solver: Gauss elimination' 
  elseif (c_Key_SLOE==3) then
      print *, '    Solver: Pardiso' 
  elseif (c_Key_SLOE==4) then
      print *, '    Solver: ITPACK' 
  elseif (c_Key_SLOE==5) then
      print *, '    Solver: LAPACK' 
  elseif (c_Key_SLOE==6) then
      print *, '    Solver: MUMPS' 
  elseif (c_Key_SLOE==7) then
      print *, '    Solver: UMFPACK' 
  elseif (c_Key_SLOE==8) then
      print *, '    Solver: Lis' 
  elseif (c_Key_SLOE==9) then
      print *, '    Solver: SuperLU' 
  elseif (c_Key_SLOE==10) then
      print *, '    Solver: y12m' 
  elseif (c_Key_SLOE==12) then
      print *, '    Solver: STRUMPACK' 
  endif
elseif(Key_Indent==2)then
  print *,'      Dimension of Linear System:', n
  if (c_Key_SLOE==1) then
      print *, '      Solver: direct method (inverse matrix)' 
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
  print *,'        Dimension of Linear System:', n
  if (c_Key_SLOE==1) then
      print *, '        Solver: direct method (inverse matrix)' 
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
  print *,'           Dimension of Linear System:', n
  if (c_Key_SLOE==1) then
      print *, '           Solver: direct method (inverse matrix)' 
  elseif (c_Key_SLOE==2) then
      print *, '           Solver: Gauss elimination' 
  elseif (c_Key_SLOE==3) then
      print *, '           Solver: Pardiso' 
  elseif (c_Key_SLOE==4) then
      print *, '           Solver: ITPACK' 
  elseif (c_Key_SLOE==5) then
      print *, '           Solver: LAPACK' 
  elseif (c_Key_SLOE==6) then
      print *, '           Solver: MUMPS' 
  elseif (c_Key_SLOE==7) then
      print *, '           Solver: UMFPACK' 
  elseif (c_Key_SLOE==8) then
      print *, '           Solver: Lis' 
  elseif (c_Key_SLOE==9) then
      print *, '           Solver: SuperLU' 
  elseif (c_Key_SLOE==10) then
      print *, '           Solver: y12m' 
  elseif (c_Key_SLOE==12) then
      print *, '           Solver: STRUMPACK' 
  endif
endif
D(1:n) =ZR

select case(c_Key_SLOE)
case(1)
    ALLOCATE(inv_K(n,n))
    call Matrix_Inverse(K,inv_K,n)
    D = MATMUL(inv_K,F)
    DEALLOCATE(inv_K)
case(2)
    ALLOCATE(LIN(n,n+1))
    LIN(:,1:n) = K
    LIN(:,n+1) = F
    call Solver_Gauss(LIN, n, n+1, n, D, Yes_SINGUL)
    if(Yes_SINGUL)then
    print *, '      Error :: Singular matrix!'
    return
    end if
    DEALLOCATE(LIN)
case(3)
#if defined(ifort)
    if(Key_Indent==0)then
    elseif (Key_Indent==1)then
    endif
    if (Key_LSOE_Sys==1) then
        NZ_NUM = 0
        do i=1,n
            c_count = 0
            do j=i,n
                if(K(i,j).ne.ZR) then
                    c_count = c_count +1
                    NZ_NUM             = NZ_NUM +1
                    Sparse_K(NZ_NUM)   = K(i,j)
                    J_Sparse_K(NZ_NUM) = j
                    if (c_count.eq.1) then
                        I_Map_Sparse_K(i)  = NZ_NUM
                    end if
                end if
            end do
            I_Map_Sparse_K(i+1) = NZ_NUM+1
        end do
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
            Sparse_K(1:NZ_NUM),I_Map_Sparse_K,J_Sparse_K(1:NZ_NUM), & 
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
    phase = 22
    CALL pardiso (pt, maxfct, mnum, mtype, phase,n,  & 
            Sparse_K(1:NZ_NUM),I_Map_Sparse_K, & 
            J_Sparse_K(1:NZ_NUM),              & 
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
    phase = 33
    pardiso_b = F
    CALL pardiso (pt, maxfct, mnum, mtype, phase, n, & 
            Sparse_K(1:NZ_NUM),I_Map_Sparse_K, & 
            J_Sparse_K(1:NZ_NUM),              & 
            idum, nrhs, iparm, msglvl, pardiso_b, D, error)
    phase = -1
    call pardiso(pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, idum, nrhs, &
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
    if (Key_LSOE_Sys==1) then
        NZ_NUM = 0
        do i=1,n
            c_count = 0
            do j=i,n
                if(K(i,j).ne.ZR) then
                    c_count = c_count +1
                    NZ_NUM             = NZ_NUM +1
                    Sparse_K(NZ_NUM)   = K(i,j)
                    J_Sparse_K(NZ_NUM) = j
                    if (c_count.eq.1) then
                        I_Map_Sparse_K(i)  = NZ_NUM
                    end if
                end if
            end do
            I_Map_Sparse_K(i+1) = NZ_NUM+1
        end do
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
            Sparse_K(1:NZ_NUM),I_Map_Sparse_K,J_Sparse_K(1:NZ_NUM), & 
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
    phase = 22
    CALL pardiso (pt, maxfct, mnum, mtype, phase,n,  & 
            Sparse_K(1:NZ_NUM),I_Map_Sparse_K, & 
            J_Sparse_K(1:NZ_NUM),              & 
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
    phase = 33
    pardiso_b = F
    CALL pardiso (pt, maxfct, mnum, mtype, phase, n, & 
            Sparse_K(1:NZ_NUM),I_Map_Sparse_K, & 
            J_Sparse_K(1:NZ_NUM),              & 
            idum, nrhs, iparm, msglvl, pardiso_b, D, error)
    phase = -1
    call pardiso(pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, idum, nrhs, &
         iparm, msglvl, ddum, ddum, error)

    if(Key_Indent==0)then
        print *,Space_4//'Solving completed...'
    elseif (Key_Indent==1)then
        print *,Space_8//'Solving completed...'
    endif
#endif

case(4)

    
    print *, '          ****************************************'
    print *, '               < LSOE Iterative Solver ITPACK >   '
    print *, '          ----------------------------------------'
    print *, '          Preparing the sparse format matrix......'
    ALLOCATE(ITPACK_J_Sparse_K(n*n))
    ALLOCATE(ITPACK_Sparse_K(n*n))
    ALLOCATE(ITPACK_WKSP(4*n + 4*ITPACK_ITMAX))
    allocate(ITPACK_D_Estimate(n))
    allocate(ITPACK_I_Map_Sparse_K(n+1))
    allocate(ITPACK_IWKSP(3*n))

    if (Key_LSOE_Sys==1) then
        NZ_NUM = 0
        do i=1,n
            ITPACK_c_count = 0
            do j=i,n
                if(K(i,j).ne.ZR) then
                    ITPACK_c_count = ITPACK_c_count +1
                    NZ_NUM             = NZ_NUM +1
                    ITPACK_Sparse_K(NZ_NUM)   = K(i,j)
                    ITPACK_J_Sparse_K(NZ_NUM) = j
                    if (ITPACK_c_count.eq.1) then
                        ITPACK_I_Map_Sparse_K(i)  = NZ_NUM
                    end if
                end if
            end do
            ITPACK_I_Map_Sparse_K(i+1) = NZ_NUM+1
        end do
        ITPACK_IPARM(5)  = 0
        ITPACK_NW = 4*n + 2*ITPACK_ITMAX
    elseif (Key_LSOE_Sys==0) then
        NZ_NUM = 0
        do i=1,n
            ITPACK_c_count = 0
            do j=1,n
                if(K(i,j).ne.ZR) then
                    ITPACK_c_count = ITPACK_c_count +1
                    NZ_NUM             = NZ_NUM +1
                    ITPACK_Sparse_K(NZ_NUM)   = K(i,j)
                    ITPACK_J_Sparse_K(NZ_NUM) = j
                    if (ITPACK_c_count.eq.1) then
                        ITPACK_I_Map_Sparse_K(i)  = NZ_NUM
                    end if
                end if
            end do
            ITPACK_I_Map_Sparse_K(i+1) = NZ_NUM+1
        end do
        ITPACK_IPARM(5)  = 1
        ITPACK_NW = 4*n + 4*ITPACK_ITMAX
    end if
    print *, "          Iterative solving......"
    CALL DFAULT(ITPACK_IPARM,ITPACK_RPARM)
    ITPACK_IPARM(1)  = ITPACK_ITMAX
    ITPACK_IPARM(2)  = 1
    ITPACK_IPARM(9)  = -1
    ITPACK_IPARM(12) = 0
    ITPACK_RPARM(1)  = 1.0D-7
    ITPACK_NNZ = NZ_NUM
    ITPACK_D_Estimate(1:n) = ZR
    call JSI(n,ITPACK_I_Map_Sparse_K,ITPACK_J_Sparse_K(1:NZ_NUM),&
            ITPACK_Sparse_K(1:NZ_NUM),F,&
            ITPACK_D_Estimate,ITPACK_IWKSP,ITPACK_NW,&
            ITPACK_WKSP,ITPACK_IPARM,ITPACK_RPARM,IERR)
    D = ITPACK_D_Estimate

    DEALLOCATE(ITPACK_J_Sparse_K)
    DEALLOCATE(ITPACK_Sparse_K)
    DEALLOCATE(ITPACK_WKSP)
    deallocate(ITPACK_D_Estimate)
    deallocate(ITPACK_I_Map_Sparse_K)
    deallocate(ITPACK_IWKSP)

    print *, '          ----------------------------------------'
    print *, '                  < Exit Solver ITPACK >          '
    print *, '          ****************************************'
      
case(5)
#ifndef Silverfrost
        lda = n
        if(FT==8)then
            call dgetrf( n, n, K,lda, ipiv, info)
            D = F
            call dgetrs( 'n', n, 1, K, lda, ipiv, D, n, info )
        elseif(FT==4)then
            call sgetrf( n, n, K,lda, ipiv, info)
            D = F
            call sgetrs( 'n', n, 1, K, lda, ipiv, D, n, info )
        endif
        if(Key_Cond_Number==1) then
            print *,'    Calculating the condition number...... '
            print *,'    Note: maybe time consuming!            '
            ALLOCATE(K_inv(n,n))
            K_inv = K
            call dgetri(n,K_inv,lda,ipiv,Lapack_work,n,info)
            call Matrix_Norm2(n,K,Norm_2_K)
            call Matrix_Norm2(n,K_inv,Norm_2_K_inv)
            Cond_K =Norm_2_K*Norm_2_K_inv
            DEALLOCATE(K_inv)
            if(Key_Indent==0)then
                print *, '    Conditional number of K:', Cond_K
            elseif(Key_Indent==4)then
                print *, '        Conditional number of K:', Cond_K
            elseif(Key_Indent==7)then
                print *, '           Conditional number of K:', Cond_K
            endif
        endif
#endif

case(6)
#ifdef gfortran
#if defined(notlinux) || defined(github)
#ifndef sffortranmpi
    if(FT==8)then
        CALL MPI_INIT(IERR_6)
        d_mumps_par%COMM = MPI_COMM_WORLD
        d_mumps_par%JOB = -1
        if (Key_LSOE_Sys==1) then
        d_mumps_par%SYM = 1
        else
        d_mumps_par%SYM = 0
        endif
        d_mumps_par%PAR = 1
        CALL DMUMPS(d_mumps_par)
        d_mumps_par%ICNTL(1) = 6
        d_mumps_par%ICNTL(2) = 0
        d_mumps_par%ICNTL(3) = 0
        d_mumps_par%ICNTL(4) = 2
        if (Key_LSOE_Sys==1) then
            NZ_NUM = 0
            do i=1,n
                do j=i,n
                    if(K(i,j).ne.ZR) then
                        NZ_NUM  = NZ_NUM +1
                    end if
                end do
            end do
        elseif (Key_LSOE_Sys==0) then
            NZ_NUM = 0
            do i=1,n
                do j=1,n
                        if(K(i,j).ne.ZR) then
                            NZ_NUM  = NZ_NUM +1
                        end if
                end do
            end do
        endif
        d_mumps_par%N = n
        d_mumps_par%NZ= NZ_NUM
        ALLOCATE(d_mumps_par%IRN (d_mumps_par%NZ))
        ALLOCATE(d_mumps_par%JCN (d_mumps_par%NZ))
        ALLOCATE(d_mumps_par%A (d_mumps_par%NZ))
        ALLOCATE(d_mumps_par%RHS (d_mumps_par%N))
        NZ_NUM = 0
        if (Key_LSOE_Sys==1) then
            do i=1,n
                do j=i,n
                        if(K(i,j).ne.ZR) then
                            NZ_NUM  = NZ_NUM +1
                            d_mumps_par%IRN(NZ_NUM) = i
                            d_mumps_par%JCN(NZ_NUM) = j
                            d_mumps_par%A(NZ_NUM)   = K(i,j)
                        end if
                end do
            end do
        elseif (Key_LSOE_Sys==0) then
            do i=1,n
                do j=1,n
                    if(K(i,j).ne.ZR) then
                        NZ_NUM  = NZ_NUM +1
                        d_mumps_par%IRN(NZ_NUM) = i
                        d_mumps_par%JCN(NZ_NUM) = j
                        d_mumps_par%A(NZ_NUM)   = K(i,j)
                    end if
                end do
            end do
        endif
        d_mumps_par%RHS(1:n) = F(1:n)
        d_mumps_par%JOB = 6
        CALL DMUMPS(d_mumps_par)
        IF (d_mumps_par%MYID .eq. 0) THEN
            D(1:n) = d_mumps_par%RHS(1:n)
        END IF
        IF (d_mumps_par%MYID .eq. 0)THEN
            DEALLOCATE(d_mumps_par%IRN)
            DEALLOCATE(d_mumps_par%JCN)
            DEALLOCATE(d_mumps_par%A  )
            DEALLOCATE(d_mumps_par%RHS)
        END IF
        d_mumps_par%JOB = -2
        CALL DMUMPS(d_mumps_par)
        CALL MPI_FINALIZE(IERR)

    elseif(FT==4)then

        CALL MPI_INIT(IERR_6)
        s_mumps_par%COMM = MPI_COMM_WORLD
        s_mumps_par%JOB = -1
        if (Key_LSOE_Sys==1) then
        s_mumps_par%SYM = 1
        else
        s_mumps_par%SYM = 0
        endif
        s_mumps_par%PAR = 1
        CALL SMUMPS(s_mumps_par)
        s_mumps_par%ICNTL(1) = 6
        s_mumps_par%ICNTL(2) = 0
        s_mumps_par%ICNTL(3) = 0
        s_mumps_par%ICNTL(4) = 2
        if (Key_LSOE_Sys==1) then
            NZ_NUM = 0
            do i=1,n
                do j=i,n
                    if(K(i,j).ne.ZR) then
                        NZ_NUM  = NZ_NUM +1
                    end if
                end do
            end do
        elseif (Key_LSOE_Sys==0) then
            NZ_NUM = 0
            do i=1,n
                do j=1,n
                        if(K(i,j).ne.ZR) then
                            NZ_NUM  = NZ_NUM +1
                        end if
                end do
            end do
        endif

        s_mumps_par%N = n
        s_mumps_par%NZ= NZ_NUM
        ALLOCATE( s_mumps_par%IRN ( s_mumps_par%NZ ) )
        ALLOCATE( s_mumps_par%JCN ( s_mumps_par%NZ ) )
        ALLOCATE( s_mumps_par%A( s_mumps_par%NZ ) )
        ALLOCATE( s_mumps_par%RHS ( s_mumps_par%N  ) )
        NZ_NUM = 0
        if (Key_LSOE_Sys==1) then
            do i=1,n
                do j=i,n
                        if(K(i,j).ne.ZR) then
                            NZ_NUM  = NZ_NUM +1
                            s_mumps_par%IRN(NZ_NUM) = i
                            s_mumps_par%JCN(NZ_NUM) = j
                            s_mumps_par%A(NZ_NUM)   = real(K(i,j))
                        end if
                end do
            end do
        elseif(Key_LSOE_Sys==0)then
            do i=1,n
                do j=1,n
                        if(K(i,j).ne.ZR) then
                            NZ_NUM  = NZ_NUM +1
                            s_mumps_par%IRN(NZ_NUM) = i
                            s_mumps_par%JCN(NZ_NUM) = j
                            s_mumps_par%A(NZ_NUM)   = real(K(i,j))
                        end if
                end do
            end do
        endif
        s_mumps_par%RHS(1:n) = real(F(1:n))
        s_mumps_par%JOB = 6
        CALL SMUMPS(s_mumps_par)
        IF ( s_mumps_par%MYID .eq. 0 ) THEN
            D(1:n) = s_mumps_par%RHS(1:n)
        END IF
        IF (s_mumps_par%MYID .eq. 0)THEN
            DEALLOCATE( s_mumps_par%IRN )
            DEALLOCATE( s_mumps_par%JCN )
            DEALLOCATE( s_mumps_par%A   )
            DEALLOCATE( s_mumps_par%RHS )
        END IF
        s_mumps_par%JOB = -2
        CALL SMUMPS(s_mumps_par)
        CALL MPI_FINALIZE(IERR)
    endif    
#endif
#endif
#endif

case(7)
#ifdef gfortran
#ifdef notlinux
#ifndef macos
#ifndef github
    Ap(1:n+1) = 0

    NZ_NUM = COUNT(K.ne.ZR)

    ALLOCATE(Ax(NZ_NUM))
    ALLOCATE(Ai(NZ_NUM))
    NZ_NUM = 0
    do j=1,n
        do i=1,n
            if(K(i,j).ne.ZR) then
                NZ_NUM  = NZ_NUM +1
                Ai(NZ_NUM) = i-1
                Ax(NZ_NUM) = K(i,j)
            end if
        end do
        Ap(j+1) = NZ_NUM
    end do
    call umf4def (control_7)
    control_7(1) = 1
    call umf4pcon (control_7)
    call umf4sym (n, n, Ap, Ai, Ax, symbolic, control_7, info_7)
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
    DEALLOCATE(Ax)
    DEALLOCATE(Ai)
#endif
#endif
#endif
#endif

case(8)
#ifdef gfortran
#ifdef notlinux
#ifndef github
    NZ_NUM = COUNT(K.ne.ZR)
    ALLOCATE(K_Value(NZ_NUM))
    ALLOCATE(K_Index(NZ_NUM))
    ALLOCATE(K_Ptr(n+1))

    NZ_NUM = 0
    K_Ptr(1) = 0
    do i=1,n
        do j=1, n
        if (K(i,j) .ne. ZR) then
            NZ_NUM = NZ_NUM+1
            K_Value(NZ_NUM) = K(i,j)
            K_Index(NZ_NUM) = j-1
        endif
        enddo
        K_Ptr(i+1) = NZ_NUM
    enddo

    call Matrix_Solve_LSOE_Lis(Key_Indent,NZ_NUM,n,K_Ptr,K_Index,K_Value,F,D)

    if(Key_Cond_Number==1)then
        print *,'    Calculating the condition number...... '
        print *,'    Note: maybe time consuming!            '
        print *,'    Calculator: MC75 of HSL Package        '
        print *,'    http://www.hsl.rl.ac.uk/catalogue/     '
        ALLOCATE(Cond_A(NZ_NUM))
        ALLOCATE(Cond_IRN(NZ_NUM))
        ALLOCATE(Cond_JCN(NZ_NUM))
        NZ_NUM = 0
        do i=1,n
            do j=1, n
                if (K(i,j) .ne. ZR) then
                    NZ_NUM = NZ_NUM+1
                    Cond_A(NZ_NUM) = K(i,j)
                    Cond_IRN(NZ_NUM) = i
                    Cond_JCN(NZ_NUM) = j
                endif
            enddo
            K_Ptr(i+1) = NZ_NUM
        enddo
        call Matrix_Condition_Number(n,NZ_NUM,Cond_A(1:NZ_NUM),Cond_IRN(1:NZ_NUM), &
                     Cond_JCN(1:NZ_NUM),Condition_Num_Norm2)
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
    SuperLU_nnz = COUNT(abs(K).GE.(Tol_30))
    allocate(SuperLU_rowind(SuperLU_nnz))
    allocate(SuperLU_value(SuperLU_nnz))

    SuperLU_nnz = 0
    do j=1,n
        cc_count=0
        do i=1, n
            if (abs(K(i,j)) .GE.(Tol_30)) then
                SuperLU_nnz = SuperLU_nnz+1
                cc_count=cc_count+1
                if(cc_count==1)then
                    SuperLU_colptr(j) = SuperLU_nnz
                endif
                SuperLU_value(SuperLU_nnz) = K(i,j)
                SuperLU_rowind(SuperLU_nnz) = i
            endif
        enddo
    enddo
    SuperLU_colptr(n+1) = SuperLU_nnz +1

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

#ifdef cbfortran
#else
    NZ_NUM = COUNT(K.ne.ZR)
    y12m_NN  = 10*NZ_NUM
    y12m_NN1 = 10*NZ_NUM
    y12m_IHA = n
    ALLOCATE(K_Value_y12m(y12m_NN))
    ALLOCATE(K_Index_y12m(y12m_NN))
    ALLOCATE(K_Index_Row_y12m(y12m_NN1))
    NZ_NUM = 0
    do i=1,n
        do j=1, n
            if (K(i,j) .ne. ZR) then
                NZ_NUM = NZ_NUM+1
                K_Value_y12m(NZ_NUM) = K(i,j)
                K_Index_y12m(NZ_NUM) = j
                K_Index_Row_y12m(NZ_NUM) = i
            endif
        enddo
    enddo


    call Solver_Y12MAF(n, NZ_NUM, K_Value_y12m,  &
    K_Index_y12m, y12m_NN, K_Index_Row_y12m, y12m_NN1, &
    PIVOT, HA, y12m_IHA, AFLAG, IFLAG, F, IFAIL)
    D=F
    if(ifail.gt.0)then
        if(Key_Indent==0)then
            print *,"    Linear solver y12m failed: error code-",ifail
        elseif(Key_Indent==2)then
            print *,"      Linear solver y12m failed: error code-",ifail
        elseif(Key_Indent==4)then
            print *,"        Linear solver y12m failed: error code-",ifail
        elseif(Key_Indent==7)then
            print *,"           Linear solver y12m failed: error code-",ifail
        endif
        call Warning_Message('S',Keywords_Blank)
    endif
    DEALLOCATE(K_Value_y12m)
    DEALLOCATE(K_Index_y12m)
    DEALLOCATE(K_Index_Row_y12m)
#endif 

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

    NZ_NUM = COUNT(K.ne.ZR)
    ALLOCATE(K_Values(NZ_NUM))
    ALLOCATE(K_Column_indices(NZ_NUM))
    ALLOCATE(K_Row_Pointers(n+1))

    NZ_NUM = 0
    K_Row_Pointers(1) = 0
    do i=1,n
        do j=1, n
            if (K(i,j) .ne. ZR) then
                NZ_NUM = NZ_NUM+1
                K_Values(NZ_NUM) = K(i,j)
                K_Column_indices(NZ_NUM) = j-1
            endif
        enddo
        K_Row_Pointers(i+1) = NZ_NUM
    enddo

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
END SUBROUTINE Matrix_Solve_LSOE
