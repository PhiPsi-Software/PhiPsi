 
      SUBROUTINE Matrix_Count_Row_Int(m,n,Matrix,Count_Row,Num_Once)   
      use Global_Float_Type      
      implicit none
      integer i,j,m,n
      integer Num_Once
      integer Matrix(m,n),Count_Row(m)
      integer Vector_1(n),Vector_2(n)
      
      do i=1,m
          Count_Row(i)=1
          Vector_1 = Matrix(i,1:n)
          do j=1,m
              if(j.ne.i)then
                  Vector_2 = Matrix(j,1:n)
                  if (ALL(Vector_1.EQ.Vector_2)) then
                      Count_Row(i) = Count_Row(i) + 1
                  end if
                  
              end if
          end do
      end do
      
      Num_Once=count(Count_Row.eq.1)
      
      return
      END SUBROUTINE Matrix_Count_Row_Int
    


