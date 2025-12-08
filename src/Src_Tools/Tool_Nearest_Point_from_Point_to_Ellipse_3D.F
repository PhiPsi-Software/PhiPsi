 
      subroutine Tool_Nearest_Point_from_Point_to_Ellipse_3D(
     &               In_Point,Center,a,b,vec_U,vec_V,Out_Point)

      use Global_Float_Type
      use fmin_module
      
      implicit none
      real(kind=FT),intent(in)::In_Point(3),Center(3),a,b
      real(kind=FT),intent(in)::vec_U(3),vec_V(3)
      real(kind=FT),intent(out)::Out_Point(3)
      real(kind=FT) :: ax, bx, xmin
      real(kind=FT),parameter :: tol = 1.0D-8
#ifndef Silverfrost
      ax = ZR
      bx = TWO*pi

      xmin = fmin(func,ax,bx,tol)
      
      Out_Point(1:3) = Center(1:3) + 
     &                 a*cos(xmin)*vec_U(1:3)+b*sin(xmin)*vec_V(1:3)
     
 
      contains
      
      function func(x) result(f)
      use Global_Float_Type 
      implicit none
 
      real(kind=FT),intent(in) :: x
      real(kind=FT)            :: f
 
      f = 
     & (In_Point(1)-Center(1)-a*cos(x)*vec_U(1)-b*sin(x)*vec_V(1))**2 +
     & (In_Point(2)-Center(2)-a*cos(x)*vec_U(2)-b*sin(x)*vec_V(2))**2 +
     & (In_Point(3)-Center(3)-a*cos(x)*vec_U(3)-b*sin(x)*vec_V(3))**2
 
      end function func
 
      
#endif
      end SUBROUTINE Tool_Nearest_Point_from_Point_to_Ellipse_3D                         
