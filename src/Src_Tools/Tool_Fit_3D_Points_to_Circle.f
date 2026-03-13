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
 
      subroutine Tool_Fit_3D_Points_to_Circle(num_Point,
     &                                    In_Points,Out_Points)
C     Fit 3D spatial points into a circle, and then output the points on the fitted circle.
c     Idea: (1) Calculate the average coordinates to find the center of the circle;
c     (2) Calculate the average distance from the points to the center of the circle, which is the radius.
c     (3) Calculate the outward normal vector of the plane, thereby obtaining the equation of the circle;
c     (4) Loop through each point, find the point on the circle closest to each point, and mark it as the new coordinate.
c     ---------------Reference materials for obtaining out-of-plane normal vectors-----------------
c     Theory Ref:https://www.freesion.com/article/7157945105/
c                 or
c     \theory_documents\026 Optimal Spatial Circle Fitting for 3D Discrete Points and Implementation-2021-10-30.pdf
c                 
c     Firstly written on 2021-10-30.
      !......................
      ! Variable Declaration
      !......................
      use Global_Float_Type     
      implicit none
      integer,intent(in)::num_Point
      real(kind=FT),intent(in)::In_Points(num_Point,3)
      real(kind=FT),intent(out)::Out_Points(num_Point,3)
      real(kind=FT) Center(3),n_Vector(3)
      real(kind=FT) Ave_Radius,Tool_Function_2Point_Dis_3D
      real(kind=FT) Nearest_Point(3)
      integer i_P
       real(kind=FT) sum_1,sum_2,sum_3
      
      
      sum_1 = ZR
      sum_2 = ZR
      sum_3 = ZR
      do i_P = 1,num_Point
          sum_1 = sum_1 + In_Points(i_P,1)
          sum_2 = sum_2 + In_Points(i_P,2)
          sum_3 = sum_3 + In_Points(i_P,3)
      enddo
      Center(1) = sum_1/dble(num_Point)
      Center(2) = sum_2/dble(num_Point)
      Center(3) = sum_3/dble(num_Point)
      
      
      Ave_Radius = ZR
      do i_P = 1,num_Point
          Ave_Radius =  Ave_Radius + 
     &      Tool_Function_2Point_Dis_3D(In_Points(i_P,1:3),Center)
      enddo
      Ave_Radius = Ave_Radius/dble(num_Point)
      

      
      call Tool_Get_Normal_Vector_of_Points_3D(
     &                 In_Points,num_Point,n_Vector)
     
      
      Out_Points = ZR
      do i_P = 1,num_Point
          call Tool_Nearest_Point_from_Point_to_Cricle_3D(
     &                 In_Points(i_P,1:3),Center,Ave_Radius,
     &                 n_Vector,Nearest_Point)     
          Out_Points(i_P,1:3) = Nearest_Point(1:3)
      enddo

      
      
      return 
      end SUBROUTINE Tool_Fit_3D_Points_to_Circle                 
