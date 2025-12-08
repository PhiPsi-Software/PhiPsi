 
      SUBROUTINE Matrixes_Equal_Is_Int(Matrix_A,Matrix_B,m,n,Yes)   
      use Global_Float_Type
      implicit none
      integer,intent(in):: m,n
      integer,intent(in):: Matrix_A(m,n),Matrix_B(m,n)
      logical,intent(out)::Yes
      
      integer i,j
      
      Yes = .False.
      
      do i=1,m
          do j=1,n
              if(Matrix_A(i,j) .ne. Matrix_B(i,j)) then
                  return
              end if
          end do
      end do
      
      Yes = .True.

      return
      END SUBROUTINE Matrixes_Equal_Is_Int
    


