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
 
      subroutine Tool_2_Points_of_Cir_give_Point_L(
     &                      x,y,r,
     &                      Point_I,Length_of_Arc,Point_Status,
     &                      P1_P2)
C     Given a circle and a point I on the circle, calculate the coordinates of the two points on the circle where the arc length from I is Length_of_Arc.
c     Usage: See the diagram at the lower left corner of note V5-58_P117.
c     Algorithm principle: The essence of this problem is to calculate the intersection points of two circles.

      use Global_Float_Type
      use Global_Elem_Area_Vol
      
      implicit none
      real(kind=FT),intent(in)::x,y,r,Point_I(2),Length_of_Arc
      real(kind=FT),intent(out)::P1_P2(2,2)
      integer,intent(out)::Point_Status
      real(kind=FT) Tool_Function_2Point_Dis
      real(kind=FT) c_r
      real(kind=FT) x1,y1,r1
      real(kind=FT) Radian_Angle
      real(kind=FT) c_Inters(2,2)
      integer c_Circles_Status,c_num_Inters
      
      
      P1_P2(1:2,1:2) = ZR      
      c_r = Tool_Function_2Point_Dis(Point_I,[x,y])
      if(abs(c_r - r)<=Tol_11)then
          Point_Status = 1
          Radian_Angle = Length_of_Arc/(TWO*pi*r)*(TWO*pi)
          x1 = Point_I(1)
          y1 = Point_I(2)
          r1 = TWO*(r*sin(Radian_Angle/TWO))
          call Tool_Intersection_Circle_and_Circle(
     &                      x,y,r,
     &                      x1,y1,r1,c_Circles_Status,
     &                      c_num_Inters,c_Inters)  
          P1_P2 = c_Inters
      else
          Point_Status =-1
          P1_P2(1:2,1:2) = ZR
      endif
      
      
      return 
      end SUBROUTINE Tool_2_Points_of_Cir_give_Point_L                         
