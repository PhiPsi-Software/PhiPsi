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
 
      subroutine Cal_Angle_of_AB_and_BC(Line_AB,Line_BC,angle_AB_BC)
c     Calculates the included angle of line AB and BC.
c     If the direction is clockwise or anti-clockwise, controlled by Direc_String.
c
c                                         B                       A
c                                         ------------------------
c                                        -
c           clockwise                  -        angle_AB_B
c                \   angle_AB_B      -
c                 \                -    
c                  \             -
c                   \          -    
c                    \       -    
c                    -/     C
      use Global_Float_Type
      use Global_Common
      implicit none
      real(kind=FT),intent(in)::Line_AB(2,2),Line_BC(2,2)
      real(kind=FT),intent(out)::angle_AB_BC
      real(kind=FT) a_x,a_y,b_x,b_y,c_x,c_y,tem1,tem2,angle
      real(kind=FT) tt1
      real(kind=FT) a(3),b(3),Direction(3)
      
      if  ((Line_AB(2,1) .ne. Line_BC(1,1)) .or.
     &     (Line_AB(2,2) .ne. Line_BC(1,2))) then
          print *,'    ********************************'
     &                 //'*********************************'    
          print *,'    Attention should be paid here!'
          print *,'    This situation rarely' 
     &              //' appears in subroutine Cal_Angle_of_AB_and_BC!'
          print *,'    ********************************'
     &                 //'*********************************'  
      end if     
      a_x = Line_AB(1,1)
      a_y = Line_AB(1,2)
      b_x = Line_AB(2,1)
      b_y = Line_AB(2,2)
      c_x = Line_BC(2,1)
      c_y = Line_BC(2,2)

      tem1 = sqrt((a_x-b_x)**2 + (a_y-b_y)**2)
      tem2 = sqrt((c_x-b_x)**2 + (c_y-b_y)**2)
      tt1 = DOT_PRODUCT([a_x-b_x,a_y-b_y],[c_x-b_x,c_y-b_y])/(tem1*tem2)
      if (tt1 < -ONE) then
          tt1 = -ONE
      elseif(tt1 > ONE)then
          tt1 = ONE
      endif
      angle  = acos(tt1)
      
      a = [a_x-b_x,a_y-b_y,ZR]
      b = [c_x-b_x,c_y-b_y,ZR]
      
      call Vector_Cross_Product_3(a,b,Direction)
      if (Direction(3) >= ZR) then
          angle_AB_BC = angle
      elseif (Direction(3) < ZR) then
          angle_AB_BC = TWO*pi-angle
      end if
      return 
      end SUBROUTINE Cal_Angle_of_AB_and_BC              
