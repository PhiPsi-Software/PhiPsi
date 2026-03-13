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
 
      subroutine Tool_Arc_Extend_or_Shorten(
     &                      delta_L_Arc,
     &                      Old_Arc_Cr,
     &                      New_Arc_Cr,Extend_Status)
C     Extend or shorten the arc (along the original direction (clockwise or counterclockwise) delta_L_Arc).
c     delta_L_Arc > 0 means extension; delta_L_Arc < 0 means shortening
      use Global_Float_Type
      use Global_Common
      
      implicit none
      real(kind=FT),intent(in)::delta_L_Arc,Old_Arc_Cr(11)
      real(kind=FT),intent(out)::New_Arc_Cr(11)
      logical,intent(out)::Extend_Status
      
      real(kind=FT) o_x,o_y
      real(kind=FT) r
      real(kind=FT) c_Arc_Center(2)
      real(kind=FT) c_Arc_Direction
      real(kind=FT) Point_I(2)
      real(kind=FT) P1_P2(2,2)
      integer       Point_Status
      real(kind=FT) c_Arc_Start_Point(2)
      logical       c_Yes_feasible
      real(kind=FT) Arc_Radian_P1,Arc_Radian_P2
      real(kind=FT) Old_Radian,Old_Arc_L
      real(kind=FT) c_Arc_Radian_S_P1,c_Arc_Radian_E_P1
      real(kind=FT) c_Arc_Radian_S_P2,c_Arc_Radian_E_P2
      
 1001 FORMAT(5X,'-- Error :: abs(-delta_L_Arc) is > Old_Arc_L!') 
 1002 FORMAT(5X,'            Error in Tool_Arc_Extend_or_Shorten.f!')  
 
      Extend_Status = .False.
      o_x               = Old_Arc_Cr(1)
      o_y               = Old_Arc_Cr(2)
      c_Arc_Center(1:2) = [o_x,o_Y]
      c_Arc_Direction   = Old_Arc_Cr(3)
      r                 = Old_Arc_Cr(4)
      Old_Radian        = Old_Arc_Cr(7)
      Old_Arc_L         = Old_Radian*r
      c_Arc_Start_Point = [Old_Arc_Cr(8),Old_Arc_Cr(9)]
      Point_I(1)        = Old_Arc_Cr(10)
      Point_I(2)        = Old_Arc_Cr(11)
      
      if(delta_L_Arc < Tol_11)then
          if(abs(delta_L_Arc) > Old_Arc_L)then
           write (*,1001)
           write (*,1002)
           call Warning_Message('S',Keywords_Blank)
          endif
      endif
      
      
      call Tool_2_Points_of_Cir_give_Point_L(
     &                      o_x,o_y,r,
     &                      Point_I,delta_L_Arc,Point_Status,
     &                      P1_P2)
      if(Point_Status ==1 .and.(sum(P1_P2(2,1:2))>Tol_11)) then       
          call Tool_Arc_r_and_Radian_Given_Coors(c_Arc_Direction,
     &                           c_Arc_Start_Point,
     &                           P1_P2(1,1:2),c_Arc_Center,
     &                           c_Yes_feasible,
     &                           r,c_Arc_Radian_S_P1,c_Arc_Radian_E_P1,
     &                           Arc_Radian_P1)
          call Tool_Arc_r_and_Radian_Given_Coors(c_Arc_Direction,
     &                           c_Arc_Start_Point,
     &                           P1_P2(2,1:2),c_Arc_Center,
     &                           c_Yes_feasible,
     &                           r,c_Arc_Radian_S_P2,c_Arc_Radian_E_P2,
     &                           Arc_Radian_P2)
          if(delta_L_Arc >= Tol_11)then
              if(Arc_Radian_P1>=Arc_Radian_P2)then
                  New_Arc_Cr     = Old_Arc_Cr
                  New_Arc_Cr(5)  = c_Arc_Radian_S_P1
                  New_Arc_Cr(6)  = c_Arc_Radian_E_P1
                  New_Arc_Cr(7)  = Arc_Radian_P1
                  New_Arc_Cr(10) = P1_P2(1,1)
                  New_Arc_Cr(11) = P1_P2(1,2)
                  Extend_Status  = .True.
              elseif(Arc_Radian_P1 < Arc_Radian_P2)then
                  New_Arc_Cr     = Old_Arc_Cr
                  New_Arc_Cr(5)  = c_Arc_Radian_S_P2
                  New_Arc_Cr(6)  = c_Arc_Radian_E_P2
                  New_Arc_Cr(7)  = Arc_Radian_P2
                  New_Arc_Cr(10) = P1_P2(2,1)
                  New_Arc_Cr(11) = P1_P2(2,2) 
                  Extend_Status  = .True.
              endif
          endif
          if(delta_L_Arc <= (-Tol_11))then
              if(Arc_Radian_P1 <= Arc_Radian_P2)then
                  New_Arc_Cr     = Old_Arc_Cr
                  New_Arc_Cr(5)  = c_Arc_Radian_S_P1
                  New_Arc_Cr(6)  = c_Arc_Radian_E_P1
                  New_Arc_Cr(7)  = Arc_Radian_P1
                  New_Arc_Cr(10) = P1_P2(1,1)
                  New_Arc_Cr(11) = P1_P2(1,2)
                  Extend_Status  = .True.
              elseif(Arc_Radian_P1 > Arc_Radian_P2)then
                  New_Arc_Cr     = Old_Arc_Cr
                  New_Arc_Cr(5)  = c_Arc_Radian_S_P2
                  New_Arc_Cr(6)  = c_Arc_Radian_E_P2
                  New_Arc_Cr(7)  = Arc_Radian_P2
                  New_Arc_Cr(10) = P1_P2(2,1)
                  New_Arc_Cr(11) = P1_P2(2,2) 
                  Extend_Status  = .True.
              endif
          endif
      endif
      return 
      end SUBROUTINE Tool_Arc_Extend_or_Shorten                     
