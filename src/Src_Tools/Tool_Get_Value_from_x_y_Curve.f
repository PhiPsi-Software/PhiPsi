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
 
      subroutine Tool_Get_Value_from_x_y_Curve(Curve_x,Curve_y,
     &                                         num_Point,x,y)
C     Given an x value, obtain the corresponding y value from the xy curve, linearly.

      use Global_Float_Type
      implicit none
      integer,intent(in)::  num_Point
      real(kind=FT),intent(in):: Curve_x(num_Point),
     &                              Curve_y(num_Point),x
      real(kind=FT),intent(out):: y
      
      integer i
      integer location_x
      real(kind=FT) L_x,R_x,L_y,R_y
      real(kind=FT) temp
      
      do i =1,num_Point-1
          if(x>=Curve_x(i). and. x<=Curve_x(i+1))then
              location_x = i
              exit
          end if
      end do
      
      L_x = Curve_x(location_x)
      R_x = Curve_x(location_x+1)
      L_y = Curve_y(location_x)
      R_y = Curve_y(location_x+1)
      
      temp = R_x - L_x
      if(R_x==L_x)then
          temp = 1.0D30
      endif
      y = L_y + (x-L_x)*(R_y-L_y)/temp
      
      if (x < minval(Curve_x))then
          y = Curve_y(1)
      endif
      if (x > maxval(Curve_x))then
          y = Curve_y(num_Point)
      endif
      
      return 
      end SUBROUTINE Tool_Get_Value_from_x_y_Curve                          
