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
 
      subroutine Tool_Yes_Crack_Polygon_Intersects(c_Crack,
     &                                    Clo_Supp_Domain,num_Domain_P,
     &                                    Yes_Intersect)
C     Determine whether crack c_Crack intersects with the polygon
      use Global_Float_Type
      use Global_Elem_Area_Vol
      use Global_Crack
      use Global_Crack_Common
      
      implicit none
      integer,intent(in)::c_Crack,num_Domain_P
      real(kind=FT),intent(in)::Clo_Supp_Domain(num_Domain_P,2)
      logical,intent(out)::Yes_Intersect
      integer i_S,i_Edge
      real(kind=FT) crack_p1(2),crack_p2(2),edge_p1(2),edge_p2(2)
      real(kind=FT) c_X,c_Y
      logical c_Yes_Cross
      
      Yes_Intersect = .False.
      do i_S = 1,Each_Cr_Poi_Num(c_Crack)-1
          crack_p1 = Crack_Coor(c_Crack,i_S,  1:2)
          crack_p2 = Crack_Coor(c_Crack,i_S+1,1:2)
          do i_Edge = 1,num_Domain_P-1
              edge_p1 = Clo_Supp_Domain(i_Edge,1:2)
              edge_p2 = Clo_Supp_Domain(i_Edge+1,1:2)
              call Tool_Intersection(crack_p1,crack_p2,
     &                                   edge_p1,edge_p2 ,
     &                                   c_X,c_Y,c_Yes_Cross)  
              if(c_Yes_Cross)then
                  Yes_Intersect = .True.
                  return
              endif 
          enddo
      enddo
      
      return 
      end SUBROUTINE Tool_Yes_Crack_Polygon_Intersects                          
