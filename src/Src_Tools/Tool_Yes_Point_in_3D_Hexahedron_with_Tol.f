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
 
      subroutine Tool_Yes_Point_in_3D_Hexahedron_with_Tol(Point,
     &                            A,B,C,D,E,F,G,H,Tol,
     &                            Yes_in,Yes_on)
C     Determine whether a point in space is inside a spatial hexahedron (split into 3 spatial tetrahedrons for determination: H-BCGF, H-EABF, H-ABCD)
c     With tolerance.

      !......................
      ! Variable Declaration
      !......................
      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::Point(3),
     &                 A(3),B(3),C(3),D(3),E(3),F(3),G(3),H(3),
     &                 Tol
      logical,intent(out):: Yes_in,Yes_on
      logical c_Yes_in_1,c_Yes_on_1,c_Yes_in_2,c_Yes_on_2
      logical c_Yes_in_3,c_Yes_on_3
      real(kind=FT) coor_x_max,coor_x_min,coor_y_max,coor_y_min
      real(kind=FT) coor_z_max,coor_z_min
      
      Yes_in = .False.
      Yes_on = .False.
      
      coor_x_max = max(A(1),B(1),C(1),D(1),E(1),F(1),G(1),H(1))
      if(Point(1) > (coor_x_max+Tol_8))then
          return
      endif
      coor_x_min = min(A(1),B(1),C(1),D(1),E(1),F(1),G(1),H(1))
      if(Point(1) < (coor_x_min-Tol_8))then
          return
      endif
      
      coor_y_max = max(A(2),B(2),C(2),D(2),E(2),F(2),G(2),H(2))
      if(Point(2) > (coor_y_max+Tol_8))then
          return
      endif
      coor_y_min = min(A(2),B(2),C(2),D(2),E(2),F(2),G(2),H(2))
      if(Point(2) < (coor_y_min-Tol_8))then
          return
      endif
      
      coor_z_max = max(A(3),B(3),C(3),D(3),E(3),F(3),G(3),H(3))
      if(Point(3) > (coor_z_max+Tol_8))then
          return
      endif
      coor_z_min = min(A(3),B(3),C(3),D(3),E(3),F(3),G(3),H(3))
      if(Point(3) < (coor_z_min-Tol_8))then
          return
      endif

      call Tool_Yes_Point_in_3D_Pyramid_with_Tol(Point,
     &                            H,B,C,G,F,Tol,
     &                            c_Yes_in_1,c_Yes_on_1)
      if (c_Yes_in_1) then
          Yes_in = .True.
          return
      endif
      if (c_Yes_on_1) then
          Yes_on = .True.
          return
      endif
      
      call Tool_Yes_Point_in_3D_Pyramid_with_Tol(Point,
     &                            H,E,A,B,F,Tol,
     &                            c_Yes_in_2,c_Yes_on_2)
      if (c_Yes_in_2) then
          Yes_in = .True.
          return
      endif
      if (c_Yes_on_2) then
          Yes_on = .True.
          return
      endif
      
      call Tool_Yes_Point_in_3D_Pyramid_with_Tol(Point,
     &                            H,A,B,C,D,Tol,
     &                            c_Yes_in_3,c_Yes_on_3)
      if (c_Yes_in_3) then
          Yes_in = .True.
          return
      endif
      if (c_Yes_on_3) then
          Yes_on = .True.
          return
      endif
      
      
      return 
      end SUBROUTINE Tool_Yes_Point_in_3D_Hexahedron_with_Tol               
