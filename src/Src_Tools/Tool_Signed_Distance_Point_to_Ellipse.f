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
 
      subroutine Tool_Signed_Distance_Point_to_Ellipse(
     &                      x0,y0,a,b,theta,Point_C,S_Distance)
C     Calculate the signed distance from a point to a tilted ellipse
c     Negative inside the ellipse, positive outside the ellipse, zero on the ellipse.
c     2020-08-09.

      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(in)::x0,y0,a,b,theta,Point_C(2)
      real(kind=FT),intent(out)::S_Distance
      
      real(kind=FT) c_Dis,Tol
      integer Return_Statu
      
      Tol = 1.0D-10
      
      call Tool_Yes_Point_in_Oblique_Ellipse(Point_C,
     &                  x0,y0,a,b,theta,
     &                  Return_Statu,Tol)
      if(Return_Statu==1)then
          S_Distance =-ONE
      elseif(Return_Statu==2)then
          S_Distance =ONE
      else
          S_Distance =ZR
      endif
      return 
      end SUBROUTINE Tool_Signed_Distance_Point_to_Ellipse                        
