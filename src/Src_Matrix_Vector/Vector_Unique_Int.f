 
      SUBROUTINE Vector_Unique_Int(n,n_Op_F,Vector,
     &                             Uniqued_Vec,Uniqued_n)   
      use Global_Float_Type
      implicit none
      integer,intent(in):: n,n_Op_F
      integer,intent(in):: Vector(n)
      integer,intent(out)::Uniqued_Vec(n)
      integer,intent(out)::Uniqued_n
      
      integer i,j,k
      Uniqued_Vec(1:n) = 0
      k = 1
      Uniqued_Vec(1) = Vector(1)
          
      outer: do i=2,n_Op_F
          do j=1,k
              if (Uniqued_Vec(j) .eq. Vector(i)) then
                  cycle outer
              end if
          end do
          k = k + 1
          Uniqued_Vec(k) = Vector(i)
      end do outer
             
      Uniqued_n   = k

      return
      END SUBROUTINE Vector_Unique_Int
    


