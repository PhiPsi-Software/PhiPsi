 
      SUBROUTINE Matrix_Inverse_Lapack(a,n)   
      use Global_Float_Type
      implicit none 
      integer,intent(in)::n
      real(kind=FT),intent(inout):: a(n,n)
      
      integer lda,info,ipiv(n),lwork
      real(kind=FT) work(n)
      
      lda = n
      lwork = n

      if(FT==8)then
          call dgetrf(n,n,a,lda,ipiv,info)
          call dgetri(n,a,lda,ipiv,work,lwork,info)
      elseif(FT==4)then
          call sgetrf(n,n,a,lda,ipiv,info)
          call sgetri(n,a,lda,ipiv,work,lwork,info) 
      endif
      return
      END SUBROUTINE Matrix_Inverse_Lapack
    


