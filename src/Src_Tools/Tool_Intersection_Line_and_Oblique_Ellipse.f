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
 
      subroutine Tool_Intersection_Line_and_Oblique_Ellipse(
     &                      x0,y0,a,b,theta,Seg_A,Seg_B,
     &                      num_Inter,Inter)
C     Calculate the intersection points of line segment AB and the ellipse, assuming that one endpoint of the line segment is inside the ellipse and the other is outside.
c     2020-08-09.

      use Global_Float_Type      
      use Global_Common
      implicit none
      real(kind=FT),intent(in):: x0,y0,a,b,theta,Seg_A(2),Seg_B(2)
      integer,intent(out)::num_Inter
      real(kind=FT),intent(out)::Inter(2)
      real(kind=FT) Tol
      integer Return_Statu_A,Return_Statu_B,Return_Statu_P
      
      real(kind=FT) in_P(2),out_P(2),Possi_P(2)
      
      integer i_Try
      
      
      Tol = 1.0D-10
      call Tool_Yes_Point_in_Oblique_Ellipse(Seg_A,
     &                  x0,y0,a,b,theta,
     &                  Return_Statu_A,Tol)
      call Tool_Yes_Point_in_Oblique_Ellipse(Seg_B,
     &                  x0,y0,a,b,theta,
     &                  Return_Statu_B,Tol)
      if((Return_Statu_A+Return_Statu_B)/=3)then
           num_Inter = 0
           return
      endif
      

      
      if(Return_Statu_A ==1)then
          in_P(1:2)  =Seg_A(1:2)
          out_P(1:2) =Seg_B(1:2)
      else
          in_P(1:2)  =Seg_B(1:2)
          out_P(1:2) =Seg_A(1:2)
      endif
      do i_Try = 1,5000
          Possi_P(1:2) = (in_P(1:2) + out_P(1:2))/TWO
          call Tool_Yes_Point_in_Oblique_Ellipse(Possi_P,
     &                  x0,y0,a,b,theta,
     &                  Return_Statu_P,Tol)
          if(Return_Statu_P==1)then   
              in_P = Possi_P
          elseif(Return_Statu_P==2)then
              out_P = Possi_P
          elseif(Return_Statu_P==0)then
              Inter(1:2) =Possi_P
              num_Inter = 1
              exit
          endif
      enddo

      
 
      return 
      end SUBROUTINE Tool_Intersection_Line_and_Oblique_Ellipse  
