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
 
      subroutine Tool_Point_Given_Point_Normal_and_Distance_3D(
     &                                      Point,
     &                                      Normal_Vector,Distance,
     &                                      New_Point)
      ! Given a spatial point, a normal direction, and a distance, calculate the coordinates of a new
      ! spatial point.
      
      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(in) ::Point(3),Normal_Vector(3),Distance
      real(kind=FT),intent(out)::New_Point(3)
      real(kind=FT) theta1,theta2,theta3,norm_Normal_Vector
      
      norm_Normal_Vector = sqrt(Normal_Vector(1)**2 + 
     &                      Normal_Vector(2)**2 + Normal_Vector(3)**2)
      theta1 = acos(Normal_Vector(1)/norm_Normal_Vector)
      theta2 = acos(Normal_Vector(2)/norm_Normal_Vector) 
      theta3 = acos(Normal_Vector(3)/norm_Normal_Vector) 
      
      New_Point(1) = Point(1) + Distance*cos(theta1)
      New_Point(2) = Point(2) + Distance*cos(theta2)
      New_Point(3) = Point(3) + Distance*cos(theta3)
      
      return 
      end SUBROUTINE Tool_Point_Given_Point_Normal_and_Distance_3D                      
