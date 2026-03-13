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
 
      subroutine Tool_Yes_Point_in_Oblique_Ellipse(point,
     &                  x0,y0,a,b,theta,
     &                  Return_Statu,Tol)
C     Determine whether a point is inside, outside, or on a tilted ellipse.
c     Return_Status=0 indicates on the ellipse; =1 indicates inside; =2 indicates outside.
c     See the equation of an inclined ellipse in \Theoretical_Documents\011 How to Draw the Equation of Any Inclined Ellipse (Baidu Wenku)
c     2020-08-9.

      use Global_Float_Type
      use Global_Inclusion
      use Global_Common
      
      implicit none
      real(kind=FT),intent(in)::point(2),x0,y0,a,b,theta
      real(kind=FT),intent(in):: Tol
      integer,intent(out):: Return_Statu
      real(kind=FT) x,y,temp1,temp2,c,c_theta
      
      c_theta = theta*pi/180.0D0
      c = sqrt(a**2-b**2)
      x = point(1)
      y = point(2)
      temp1 =(a**2-c**2*cos(c_theta)**2)*(x-x0)**2 +
     &       (a**2-c**2*sin(c_theta)**2)*(y-y0)**2 -
     &       c**2*sin(TWO*c_theta)*(x-x0)*(y-y0)
      temp2 = (a**2)*(b**2)
      
      if((temp1-temp2)>=Tol)then
          Return_Statu =2
      elseif((temp1-temp2)<=-Tol)then
          Return_Statu =1
      else 
          Return_Statu =0
      endif
      
      return 
      end SUBROUTINE Tool_Yes_Point_in_Oblique_Ellipse                          
