 
      SUBROUTINE Matrix_Eigenvalues_by_LAPACK(n,Matrix,WR,WI)
      use Global_Float_Type
      implicit none

      integer,intent(in)::n
      real(kind=FT),intent(in)::Matrix(n,n)
      real(kind=FT),intent(out)::WR(n),WI(n)

      integer out_INFO,LWORK
      real(kind=FT) VR(n,n)
      integer LDVR
      real(kind=FT) DUMMY(1,N)
      integer lda
      real(kind=FT),ALLOCATABLE::Lapack_work(:)
      
      
      lda     = n
      LDVR    = n
      WR(1:n) = ZR
      WI(1:n) = ZR
      LWORK = 100*n
      ALLOCATE(Lapack_work(LWORK))
      
      call dgeev ('N','N', n,     Matrix, lda,
     &             WR, WI, DUMMY, 1, VR,
     &             LDVR,Lapack_work, LWORK, out_INFO)
     
      DEALLOCATE(Lapack_work)
      
      
      return
      END SUBROUTINE Matrix_Eigenvalues_by_LAPACK



