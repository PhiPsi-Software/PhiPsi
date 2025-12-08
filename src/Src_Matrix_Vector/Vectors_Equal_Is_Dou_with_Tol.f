 
      SUBROUTINE Vectors_Equal_Is_Dou_with_Tol(Vector_A,Vector_B,n,Yes)   

      use Global_Float_Type
      implicit none
      integer,intent(in):: n
      real(kind=FT),intent(in):: Vector_A(n),Vector_B(n)
      logical,intent(out)::Yes
      
      integer i
      
      Yes = .False.
      
      do i=1,n
          if(abs(Vector_A(i)-Vector_B(i))>Tol_11) then
              return
          end if
      end do
      
      Yes = .True.

      return
      END SUBROUTINE Vectors_Equal_Is_Dou_with_Tol
    


