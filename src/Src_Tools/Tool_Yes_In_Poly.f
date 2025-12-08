 
      subroutine Tool_Yes_In_Poly(x,y,xpol,ypol,npol,Yes_INOUT)
      use Global_Float_Type
      implicit none
      integer i,j,npol
      real(kind=FT) x,y,xpol(npol),ypol(npol)
      logical Yes_INOUT

      Yes_INOUT = .False.
      if((x.gt.maxval(xpol)) .or.
     &   (x.lt.minval(xpol)) .or.
     &   (y.gt.maxval(ypol)) .or.
     &   (y.lt.minval(ypol))) then
             Yes_INOUT = .False.
         goto 100
         
      else
          j = npol-1 
          do i=1,npol-1
              if ( ((ypol(i).gt.y).neqv. (ypol(j).gt.y)) .and.
     &               (x .lt. (xpol(j)-xpol(i)) * (y-ypol(i)) 
     &                / (ypol(j)-ypol(i)) + xpol(i)) ) then
                  Yes_INOUT = .NOT. Yes_INOUT
              end if
              j = i
          end do
      end if
      
  100 continue
  
      return 
      end SUBROUTINE Tool_Yes_In_Poly                          
