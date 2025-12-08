 
      SUBROUTINE Vector_Location_Int_v2(n,Vector,Variable,location)   
      
      use Global_Float_Type          
      implicit none
      integer,intent(in)::n,Variable
      integer, intent(in) :: Vector(n)
      integer,intent(out)::location
      integer i
      
      location = 0
      do i=1,n
          if (Vector(i) == Variable) then
              location = i
              return
          end if
      end do  
      
      return
      END SUBROUTINE Vector_Location_Int_v2
      
    


