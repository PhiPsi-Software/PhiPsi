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
 
      subroutine Tool_Intersection_Line_and_Arc(
     &                                   A,B,
     &                                   c_Arc_Crack_Coor,
     &                                   num_Inter,Inter,Yes_Cross) 
C     Calculate the intersection points of line segment AB and the arc segment (Added on 2017-07-16)
c     A: Coordinates of point A of line segment AB;
c     B: Coordinates of point B of line segment AB;
c     c_Arc_Crack_Coor(1:10),
c     ---------------------------
c     The calculation approach is to first find the intersection points between the line segment and the circle on which the arc lies. If there are intersection points, then determine whether each intersection point lies within the range of the arc.

      use Global_Float_Type
      use Global_Elem_Area_Vol
      
      implicit none
      real(kind=FT),intent(in)::A(2),B(2),c_Arc_Crack_Coor(11)
      integer,intent(out)::num_Inter
      real(kind=FT),intent(out)::Inter(2,2)
      logical,intent(out)::Yes_Cross
      
      real(kind=FT) o_x,o_y
      real(kind=FT) r
      integer c_num_Inter,State
      real(kind=FT) c_Inter(2,2)
      integer i_Inter
      logical c_Yes_ON
      
      Yes_Cross = .False.

      o_x               = c_Arc_Crack_Coor(1)
      o_y               = c_Arc_Crack_Coor(2)
      r                 = c_Arc_Crack_Coor(4)
     
      c_num_Inter = 0
      call Tool_Intersection_Line_and_Circle(o_x,o_y,
     &                                       r,A,B,
     &                                       c_num_Inter,State,
     &                                       c_Inter)

      num_Inter = 0
      do i_Inter = 1,c_num_Inter
          call Tool_Yes_On_Arc(c_Inter(i_Inter,1),c_Inter(i_Inter,2),
     &                         c_Arc_Crack_Coor,c_Yes_ON)
          if(c_Yes_ON)then
              Yes_Cross = .True.
              num_Inter = num_Inter +1
              Inter(num_Inter,1:2) = c_Inter(i_Inter,1:2)
          endif
      enddo
      
      return 
      end SUBROUTINE Tool_Intersection_Line_and_Arc                         
