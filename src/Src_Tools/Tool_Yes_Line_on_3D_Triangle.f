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
 
      subroutine Tool_Yes_Line_on_3D_Triangle(P1,P2,A,B,C,Yes_on,
     &                                        Point)
c     Determine whether a spatial line segment lies on a spatial triangle.
c     It is possible that both points are inside the triangle; it is also possible that one point is 
c     inside the triangle and the other is outside, but both points are in the plane of the triangle.
c     NEWFTU2022053001.
c     2022-05-30.

      !......................
      ! Variable Declaration
      !......................
      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::P1(3),P2(3),A(3),B(3),C(3)
      logical,intent(out):: Yes_on
      real(kind=FT),intent(out)::Point(3)
      logical Yes_P1_on,Yes_P2_on
      real(kind=FT) Tool_Function_Dis_Point_to_3D_Tri_v2,c_Dis
      
      Yes_on = .False.
      Point(1:3)  = -TEN_15
      
      call Tool_Yes_Point_on_3D_Triangle(P1,A,B,C,Yes_P1_on)
      call Tool_Yes_Point_on_3D_Triangle(P2,A,B,C,Yes_P2_on)
      
      if((Yes_P1_on .eqv. .False.) .and. (Yes_P2_on .eqv. .False.))then
          return
      endif
      
      if((Yes_P1_on .eqv. .True.) .and. (Yes_P2_on .eqv. .True.))then
          Yes_on = .True.
          Point  = P1
          return
      endif      
      
      if((Yes_P1_on .eqv. .True.) .and. (Yes_P2_on .eqv. .False.))then
          c_Dis = Tool_Function_Dis_Point_to_3D_Tri_v2(P2,A,B,C)
          if(abs(c_Dis)<=Tol_11) then
              Yes_on = .True.
              Point  = P1
              return
          endif
      endif       
      
      if((Yes_P2_on .eqv. .True.) .and. (Yes_P1_on .eqv. .False.))then
          c_Dis = Tool_Function_Dis_Point_to_3D_Tri_v2(P1,A,B,C)
          if(abs(c_Dis)<=Tol_11) then
              Yes_on = .True.
              Point  = P2
              return
          endif
      endif         
      return 
      end SUBROUTINE Tool_Yes_Line_on_3D_Triangle         
