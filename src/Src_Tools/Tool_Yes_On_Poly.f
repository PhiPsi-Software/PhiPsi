 
      subroutine Tool_Yes_On_Poly(x,y,xpol,ypol,npol,Yes_ONOUT)
      use Global_Float_Type
      integer i,npol
      real(kind=FT) x,y,xpol(npol),ypol(npol)
      logical,intent(out)::Yes_ONOUT

      logical Yes_ON
      

      Yes_ONOUT = .False.
      
      do i=1,npol-1 
          
          call Tool_Yes_On_Line(x,y,
     &                  [xpol(i),ypol(i)],
     &                  [xpol(i+1),ypol(i+1)],Yes_ON)
          if (Yes_ON.eqv..True.) then
              Yes_ONOUT = .True.
              exit
          endif

      end do 
      
      return 
      end SUBROUTINE Tool_Yes_On_Poly                          
