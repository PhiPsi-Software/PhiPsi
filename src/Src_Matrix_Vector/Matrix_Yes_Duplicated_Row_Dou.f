 
      SUBROUTINE Matrix_Yes_Duplicated_Row_Dou(m,n,Matrix,Yes_Dup,
     &                                         num_Dup)   

      use Global_Float_Type
      implicit none
      integer,intent(in):: m,n
      real(kind=FT),intent(in):: Matrix(m,n)
      logical,intent(out):: Yes_Dup
      integer,intent(out):: num_Dup
      integer Uniqued_m,Uni_Mat_Count(m)
      real(kind=FT) Uniqued_Matrix(m,n),tem(n)
      integer i,j,k
      
      Uniqued_Matrix(1:m,1:n) = ZR
      Yes_Dup = .False.
      Uni_Mat_Count(1:m)=1
      
      k = 1
      Uniqued_Matrix(1,:) = Matrix(1,:)
      
      outer: do i=2,m
          do j=1,k
              tem = Uniqued_Matrix(j,:) - Matrix(i,:)
              if ((abs(maxval(tem))<=Tol_11) .and.
     &            (abs(minval(tem))<=Tol_11)) then
                  Uni_Mat_Count(j) = Uni_Mat_Count(j)+1
                  Yes_Dup = .True.
                  cycle outer
              end if
          end do
          k = k + 1
          Uniqued_Matrix(k,:) = Matrix(i,:)
      end do outer
                  
      Uniqued_m = k    
      
      num_Dup = m-Uniqued_m

      return
      END SUBROUTINE Matrix_Yes_Duplicated_Row_Dou
    


