 
      SUBROUTINE Matrix_Unique_Row_Dou(m,n,m_Finish,Matrix,
     &                                 Uniqued_Matrix,Uniqued_m,
     &                                 Uni_Mat_Count)   
      use Global_Float_Type
      implicit none
      integer,intent(in):: m,n,m_Finish
      real(kind=FT),intent(in):: Matrix(m,n)
      integer,intent(out):: Uniqued_m
      real(kind=FT),intent(out):: Uniqued_Matrix(m,n)
      integer,intent(out):: Uni_Mat_Count(m)
      
      real(kind=FT) tem(n)
      integer i,j,k
      
      Uniqued_Matrix(1:m,1:n) = ZR

      Uni_Mat_Count(1:m)=1
      
      k = 1
      Uniqued_Matrix(1,:) = Matrix(1,:)
      
      outer: do i=2,m_Finish
          do j=1,k
              tem = Uniqued_Matrix(j,:) - Matrix(i,:)
              if ((abs(maxval(tem))<=Tol_11) .and.
     &            (abs(minval(tem))<=Tol_11)) then
                  Uni_Mat_Count(j) = Uni_Mat_Count(j)+1
                  cycle outer
              end if
          end do
          k = k + 1
          Uniqued_Matrix(k,:) = Matrix(i,:)
      end do outer
                  
      Uniqued_m = k            

      return
      END SUBROUTINE Matrix_Unique_Row_Dou
    


