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
 
      subroutine Tool_3D_Yes_Points_on_Same_Plane(num_points,
     &                                            Points,Yes_on)
c     Check whether the points are in the same plane.
c     Added on 2022-04-23.

      !......................
      ! Variable Declaration
      !......................
      use Global_Float_Type
      use Global_Common
      implicit none
      integer,intent(in)::num_points
      real(kind=FT),intent(in)::Points(num_points,3)
      logical,intent(out)::Yes_on
      real(kind=FT) Tri_P1(3),Tri_P2(3),Tri_P3(3)
      real(kind=FT) c_Distance,c_PER(3)
      integer i_Point
      logical c_Yes,c_Yes_PER_in,c_Yes_PER_on
      
      Yes_on = .True.
      
      if(num_points<=3)then
          print *,'    Error :: Number of points is less than 4!'
          print *,'             in Tool_3D_Yes_Points_on_Same_Plane.f!'
          call Warning_Message('S',Keywords_Blank)   
          return
      endif      
      
      Tri_P1(1:3) = Points(1,1:3)
      Tri_P2(1:3) = Points(2,1:3)
      Tri_P3(1:3) = Points(3,1:3)
      call Tool_Yes_Point_on_Line_Segment_3D(Tri_P1,Tri_P2,Tri_P3,c_Yes)
      if(c_Yes)then
          print *,'    Error :: Points 1-3 are collinear!'
          print *,'             in Tool_3D_Yes_Points_on_Same_Plane.f!'
          call Warning_Message('S',Keywords_Blank)   
          return
      endif
      do i_Point = 1,num_points-3
          call Tool_Dis_Point_to_3D_Tri(Points(i_Point+3,1:3),
     &                            Tri_P1,Tri_P2,Tri_P3,
     &                            c_Distance,c_PER(1:3),
     &                            c_Yes_PER_in,c_Yes_PER_on)
          if(c_Distance > Tol_10) then
              Yes_on = .False.
              return
          endif
      enddo

      
      return 
      end SUBROUTINE Tool_3D_Yes_Points_on_Same_Plane                  
