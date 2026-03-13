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
 
      subroutine Tool_Rank_LineSeg_Points(Line_AB,num_Points)
C     Reorder the unordered points that make up line segment AB from start to finish so that they connect end to end.
c     Sort according to the distance from each point to endpoint 1
      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(inout):: Line_AB(num_Points,2)
      integer,intent(in):: num_Points
      
      real(kind=FT) new_Line_AB(num_Points,2),L(num_Points),
     &                 A(2),B(2)
      integer i_P,Location_Min
      
      A(1) = Line_AB(1,1)
      A(2) = Line_AB(1,2)
      B(1) = Line_AB(num_Points,1)
      B(2) = Line_AB(num_Points,2) 
      
      new_Line_AB(1,1:2) = A
      new_Line_AB(num_Points,1:2) = B
      
      
      L(1:num_Points) =  1.0D8
      do i_P=2,num_Points-1
          L(i_P) = sqrt((Line_AB(i_P,2)-A(2))**2+
     &                  (Line_AB(i_P,1)-A(1))**2)
      end do
      
      do i_P=2,num_Points-1
          Location_Min = minLoc(L(1:num_Points),1) 
          new_Line_AB(i_P,1:2) = Line_AB(Location_Min,1:2)
          L(Location_Min) = 1.0D8 
      end do
      
      
      Line_AB = new_Line_AB
      
      return 
      end SUBROUTINE Tool_Rank_LineSeg_Points              
