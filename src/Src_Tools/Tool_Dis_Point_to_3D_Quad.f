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
 
      subroutine Tool_Dis_Point_to_3D_Quad(Point,
     &                            Quad_P1,Quad_P2,Quad_P3,Quad_P4,
     &                            Distance,PER,Yes_PER_in,Yes_PER_on)
C     To calculate the distance from a point in space to a spatial quadrilateral (1,2,3,4), the idea is to divide 
c     it into two triangles (1,2,3) and (1,3,4), and then call Tool_Cal_Dis_Point_to_3D_Tri for each.
c     PER indicates the coordinates of the foot of the perpendicular.
c     Regarding the determination of the distance sign: it is positive if it is consistent with the orientation 
c     from the origin O to the plane (c_Vector_o_Orient(1:3))

      !......................
      ! Variable Declaration
      !......................
      use Global_Float_Type   
      use Global_Common
      
      implicit none
      real(kind=FT),intent(in)::Point(3),
     &                 Quad_P1(3),Quad_P2(3),Quad_P3(3),Quad_P4(3)
      real(kind=FT),intent(out):: Distance,PER(3)
      logical,intent(out):: Yes_PER_in,Yes_PER_on

      
      real(kind=FT) Tri123_Distance,Tri123_PER(3)
      real(kind=FT) Tri134_Distance,Tri134_PER(3)
      logical Tri123_Yes_PER_in,Tri134_Yes_PER_in
      logical Tri123_Yes_PER_on,Tri134_Yes_PER_on

      Yes_PER_in = .False.
      Yes_PER_on = .False.
      Distance   =  ZR
      PER        =  ZR
      
      call Tool_Dis_Point_to_3D_Tri
     &               (Point,Quad_P1,Quad_P2,Quad_P3,
     &                Tri123_Distance,Tri123_PER,
     &                Tri123_Yes_PER_in,Tri123_Yes_PER_on)
      call Tool_Dis_Point_to_3D_Tri
     &               (Point,Quad_P1,Quad_P3,Quad_P4,
     &                Tri134_Distance,Tri134_PER,
     &                Tri134_Yes_PER_in,Tri134_Yes_PER_on)
      
      if(abs(Tri123_Distance-Tri134_Distance) > Tol_10)then
          print *,'    Error:: non-planar space quadrilateral!'
          print *,'            in Tool_Dis_Point_to_3D_Quad.f'
          call Warning_Message('S',Keywords_Blank)
      endif
      
      Distance = Tri123_Distance
      PER      = Tri123_PER
      if(Tri123_Yes_PER_in .or. Tri134_Yes_PER_in)then
          Yes_PER_in =.True.
      endif
      
      return 
      end SUBROUTINE Tool_Dis_Point_to_3D_Quad                  
