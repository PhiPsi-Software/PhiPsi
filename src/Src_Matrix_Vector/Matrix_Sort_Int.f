 
      SUBROUTINE Matrix_Sort_Int(m,n,Matrix)   
      use Global_Float_Type    
      implicit none
      integer i,j,m,n
      integer p
      integer Matrix(m,n),a
          
      do j=2,n
          do p=1,m
              a=Matrix(p,j)
              do i=j-1,1,-1
                  if (Matrix(p,i).le.a) goto 10
                  Matrix(p,i+1)=Matrix(p,i)
              end do
              i=0
   10         Matrix(p,i+1)=a
          end do
      end do
      
      return
      END SUBROUTINE Matrix_Sort_Int
    


