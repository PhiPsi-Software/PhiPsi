 
      subroutine Tool_Dis_Point_to_Seg(a,b,c,Dis)

      use Global_Float_Type 
      implicit none
      real(kind=FT), intent(in) :: a(2), b(2), c(2)
      real(kind=FT), intent(out) :: Dis
      real(kind=FT):: t(2), n(2), ac(2)
      real(kind=FT):: dd
      t   = b - a
      dd  = sqrt(t(1)**2+t(2)**2)
      t   = t/dd
      n   = (/-t(2), t(1)/)
      ac  = c - a
      Dis = abs(ac(1)*n(1)+ac(2)*n(2))
      
      return 
      end SUBROUTINE Tool_Dis_Point_to_Seg                       
