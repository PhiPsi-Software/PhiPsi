 
      SUBROUTINE Matrix_Sort_Row_Int(m,n,Start_m,Finish_m,Matrix)   
      use Global_Float_Type
      implicit none
      integer,intent(in)::m,n,Start_m,Finish_m
      integer,intent(inout) :: Matrix(m,n)
      
      integer i,j,iii
      integer Tem_Vector(n),a
      
      do iii=Start_m,Finish_m
          Tem_Vector = Matrix(iii,:)
          do j=2, n
              a=Tem_Vector(j)
              do i=j-1,1,-1
                  if (Tem_Vector(i).le.a) goto 10
                  Tem_Vector(i+1)=Tem_Vector(i)
              end do
              i=0
   10         Tem_Vector(i+1) = a
          end do
          Matrix(iii,:) = Tem_Vector
      end do
      
      return
      END SUBROUTINE Matrix_Sort_Row_Int
    


