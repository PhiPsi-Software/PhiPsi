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
 
      subroutine Tool_Intersection_of_AB_and_Brick_Ele(A,B,c_El,
     &                 Yes_Inter,InterSection_P) 
C     Calculate the intersection of line segment AB and the element with number 'Element' (2020-01-01)
c     Current functionality: There may be multiple intersections, but this program calculates only one and then stops.

 
      use Global_Float_Type
      use Global_Elem_Area_Vol
      use Global_Model
      
      implicit none
      real(kind=FT),intent(in)::A(3),B(3)
      integer,intent(in)::c_El
      real(kind=FT),intent(out)::InterSection_P(3)
      logical,intent(out)::Yes_Inter
      real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
      integer Triangle(12,3)
      real(kind=FT) Point1(3),Point2(3),Point3(3)
      logical c_Yes_Inter
      real(kind=FT) c_InterSection_P(3)
      integer i_Tri
      
      Yes_Inter = .False.
      InterSection_P(1:3)=ZR 
      
      c_X_NODES = G_X_NODES(1:8,c_El)
      c_Y_NODES = G_Y_NODES(1:8,c_El)  
      c_Z_NODES = G_Z_NODES(1:8,c_El)  
      
      Triangle(1,1:3) = [1,2,5]
      Triangle(2,1:3) = [2,6,5]
      Triangle(3,1:3) = [2,3,6]
      Triangle(4,1:3) = [3,7,6]
      Triangle(5,1:3) = [3,7,8]
      Triangle(6,1:3) = [3,8,4]
      Triangle(7,1:3) = [4,8,5]
      Triangle(8,1:3) = [4,5,1]
      Triangle(9,1:3) = [5,6,8]
      Triangle(10,1:3) = [6,7,8]
      Triangle(11,1:3) = [1,2,4]
      Triangle(12,1:3) = [2,3,4]    
      
      do i_Tri = 1,12
          Point1 =  [c_X_NODES(Triangle(i_Tri,1)),
     &               c_Y_NODES(Triangle(i_Tri,1)),
     &               c_Z_NODES(Triangle(i_Tri,1))]
          Point2 =  [c_X_NODES(Triangle(i_Tri,2)),
     &               c_Y_NODES(Triangle(i_Tri,2)),
     &               c_Z_NODES(Triangle(i_Tri,2))]
          Point3 =  [c_X_NODES(Triangle(i_Tri,3)),
     &               c_Y_NODES(Triangle(i_Tri,3)),
     &               c_Z_NODES(Triangle(i_Tri,3))] 
          c_Yes_Inter = .False.
          call Tool_Intersection_of_AB_and_Triangle_3D(A,B,
     &                 Point1,Point2,Point3,
     &                 c_Yes_Inter,c_InterSection_P) 
          if (c_Yes_Inter .eqv. .True.) then
              Yes_Inter = .True.
              InterSection_P =c_InterSection_P
              return
          endif
      enddo
      
      return 
      end SUBROUTINE Tool_Intersection_of_AB_and_Brick_Ele                         
