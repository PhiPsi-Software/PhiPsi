      subroutine SPARSKIT2_amux_complex(n, x, y, nnz,a,ja,ia) 
      implicit none
      integer(8),intent(in)::nnz
      integer,intent(in)::n
      complex*16 ,intent(in)::x(n),a(nnz) 
      integer,intent(in)::ja(nnz)
      integer(kind=8),intent(in)::ia(n+1)
      complex*16 ,intent(inout)::y(n)
      
      complex*16 t
      integer i
      integer(8) k
      
      
      do i = 1,n
         t = (0.0d0,0.0d0)
         do k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
         enddo
         y(i) = t
      enddo
      return
      end subroutine SPARSKIT2_amux_complex
      
      
      
      
      
      
      
      
      

