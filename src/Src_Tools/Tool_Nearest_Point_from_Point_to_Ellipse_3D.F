!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
      subroutine Tool_Nearest_Point_from_Point_to_Ellipse_3D(
     &               In_Point,Center,a,b,vec_U,vec_V,Out_Point)
C     Calculate the point on the ellipse that is closest to the vertex In_Point.
c     2022-05-07.

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
