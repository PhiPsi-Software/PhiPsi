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
 
      subroutine Tool_Shorten_or_Extend_Line(Line_AB,delta_L,
     &                                      Point_String,
     &                                      new_Line_AB,new_Point)
C     Calculate the new endpoints based on the existing line segments and offset direction
      ! Shorten or extend line_AB at point a or b by the increment of offset_L.
      ! delta_L can be negative.
      ! Point_String ='A' or 'B'.
      !                A                       B
      !
      !                 ---------------------- 
      !
      ! Point_String ='A',delta_L < 0:
      !
      !                |<---delta_L--->|
      !                                 ------
      !              
      !              
      ! Point_String ='B',delta_L > 0:
      !
      !                                        |<---delta_L--->|
      !                -----------------------------------------

      use Global_Float_Type      
      implicit none
      character*1,intent(in) :: Point_String
      real(kind=FT),intent(in) ::Line_AB(2,2),delta_L
      real(kind=FT),intent(out)::new_Line_AB(2,2),new_Point(2)
      
      real(kind=FT) a_x,a_y,b_x,b_y,theta
      
      a_x = Line_AB(1,1)
      a_y = Line_AB(1,2)
      b_x = Line_AB(2,1)
      b_y = Line_AB(2,2)
      

      theta = atan2(b_y-a_y,b_x-a_x)
      
      
      select case(Point_String)
      case('A')
          new_Line_AB(1,1) = Line_AB(1,1)-delta_L*cos(theta)
          new_Line_AB(1,2) = Line_AB(1,2)-delta_L*sin(theta)
          new_Line_AB(2,:) = Line_AB(2,:)
          new_Point = [new_Line_AB(1,1),new_Line_AB(1,2)]
      case('B')  
          new_Line_AB(2,1) = Line_AB(2,1)+delta_L*cos(theta)
          new_Line_AB(2,2) = Line_AB(2,2)+delta_L*sin(theta)
          new_Line_AB(1,:) = Line_AB(1,:)
          new_Point = [new_Line_AB(2,1),new_Line_AB(2,2)]
      end select
      
      return 
      end SUBROUTINE Tool_Shorten_or_Extend_Line                          
