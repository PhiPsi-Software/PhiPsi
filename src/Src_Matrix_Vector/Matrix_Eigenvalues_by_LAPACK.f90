!-----------------------------------------------------------
! Brief: Compute eigenvalues of a general real square matrix via LAPACK.
!
! Parameters:
!   Input:  n      - order of the matrix
!           Matrix - input matrix (n,n)
!   Output: WR    - real parts of the eigenvalues (size n)
!           WI    - imaginary parts of the eigenvalues (size n)
!
! Notes:   Thin wrapper around DGEEV with 'N' for both left and
!          right eigenvector flags, so only eigenvalues are returned.
!-----------------------------------------------------------

subroutine Matrix_Eigenvalues_by_LAPACK(n,Matrix,WR,WI)
!     The eigenvalues of a square matrix can be calculated, and the eigenvectors can also be calculated (can be calculated but not calculated).
!     Call the LAPACK library DGEEV subroutine to calculate eigenvalues and eigenvectors
!
!      (1) Official online documentation for the DGEEV subroutine in the LAPACK library:
!         http://www.netlib.org/lapack/explore-html/d9/d28/dgeev_8f.html
!         DGEEV computes for an N-by-N real nonsymmetric matrix A, the
!         eigenvalues and, optionally, the left and/or right eigenvectors.
!         The right eigenvector v(j) of A satisfies:
!                       A * v(j) = lambda(j) * v(j),where lambda(j) is its eigenvalue.
!     (2) Example of using the DGEEV subroutine from the LAPACK library:
!         https://www.nag.com/lapack-ex/node87.html
!         call DGEEV('No left vectors','Vectors (right)',N,A,LDA,WR,WI,
!                     DUMMY,1,VR,LDVR,WORK,LWORK,INFO)
use Global_Float_Type
implicit none

integer,intent(in)::n
real(kind=FT),intent(in)::Matrix(n,n)
real(kind=FT),intent(out)::WR(n),WI(n)

integer out_INFO,LWORK
real(kind=FT) VR(n,n)
integer :: LDVR
real(kind=FT) DUMMY(1,N)
integer :: lda
real(kind=FT),allocatable::Lapack_work(:)


lda     = n
LDVR    = n
WR(1:n) = ZR
WI(1:n) = ZR
LWORK = 100*n
allocate(Lapack_work(LWORK))

call dgeev ('N','N', n,     Matrix, lda, &
WR, WI, DUMMY, 1, VR, LDVR,Lapack_work, LWORK, out_INFO)

deallocate(Lapack_work)


return
END subroutine Matrix_Eigenvalues_by_LAPACK



