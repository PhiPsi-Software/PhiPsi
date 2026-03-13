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
 
      subroutine Tool_3D_Yes_Points_Lines_Cross(num_points,
     &                                          Points,Yes_Cross)
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
      logical,intent(out)::Yes_Cross
      real(kind=FT) Lines(num_points,2,3)
      integer i_Line,j_Line
      LOGICAL C_Yes_Inter
      real(kind=FT) C_InterSection_P(3),A(3),B(3),C(3),D(3)
      
      Yes_Cross = .False.
      
      if(num_points<=3)then
          print *,'     Error :: Number of points is less than 4!'
          print *,'              in Tool_3D_Yes_Points_Lines_Cross.f!'
          call Warning_Message('S',Keywords_Blank)   
          return
      endif      
      
      do i_Line=1,num_points-1
          Lines(i_Line,1,1:3) = Points(i_Line,1:3)
          Lines(i_Line,2,1:3) = Points(i_Line+1,1:3)
      enddo
      Lines(num_points,1,1:3) = Points(num_points,1:3)
      Lines(num_points,2,1:3) = Points(1,1:3)
      
      do i_Line=1,num_points-1
          A=Lines(i_Line,1,1:3)
          B=Lines(i_Line,2,1:3)
          do j_Line=i_Line+1,num_points
              if(abs(i_Line-j_Line)==1)then
                  cycle
              endif
              if((i_Line==1 .and. j_Line == num_points) .or. 
     &           (j_Line==1 .and. i_Line == num_points))then
                  cycle
              endif              
              C=Lines(j_Line,1,1:3)
              D=Lines(j_Line,2,1:3)  
              CALL Tool_Intersection_of_AB_and_CD_3D(A,B,C,D,
     &                                   C_Yes_Inter,C_InterSection_P) 
              if(C_Yes_Inter)then
                  Yes_Cross = .True.
                  return
              endif
          enddo
      enddo
      
      return 
      end SUBROUTINE Tool_3D_Yes_Points_Lines_Cross                 
