 
      SUBROUTINE Vector_Location_Int(n,Vector,Variable,location,Yes_In)   
      
      use Global_Float_Type          
      implicit none
      integer,intent(in)::n,Vector(n),Variable
      integer,intent(out)::location
      logical,intent(out):: Yes_In
      integer i
      
      Yes_In =.False.
      location = 0
      do i=1,n
          if (Vector(i) == Variable) then
              location = i
              Yes_In = .True.
              exit
          end if
      end do  
      
      return
      END SUBROUTINE Vector_Location_Int
    


