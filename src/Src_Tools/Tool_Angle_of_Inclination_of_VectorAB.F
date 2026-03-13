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
 
      subroutine Tool_Angle_of_Inclination_of_VectorAB(
     &                                     A,B,Angle)
c     Calculate the inclination of line segment AB, the inclination of the vector pointing from A to B (in degrees), with the range of the inclination being 0-360 degrees.
c     The arctangent function directly gives a tilt angle range of (-90, 90) degrees.
c     https://www.easycalculation.com/analytical/slope-length-equation.php
C     Formula:
C     m = (Y2-Y1) / (X2-X1) 
C     Distance = Square Root ( (X2-X1)2 + (Y2-Y1)2 ) 
C     Angle = arctan ( m ) 
C     Line of Equation is y = mx + b 
C     Where, m = slope
      
      use Global_Float_Type
      
      implicit none
      real(kind=FT),intent(in):: A(2),B(2)
      real(kind=FT),intent(out)::Angle
      real(kind=FT) X1,Y1,X2,Y2,m
      real(kind=FT) c_Angle
      
      X1 = A(1)
      Y1 = A(2)
      X2 = B(1)
      Y2 = B(2)

      if (abs(X2-X1)<=Tol_11)then
          m = (Y2-Y1) / Tol_11
      else
          m = (Y2-Y1) / (X2-X1)
      endif
      c_Angle = atan(m) *Con_180/pi
      
      if (abs(X2-X1)<=Tol_11)then
          if(Y1<=Y2)then
              Angle = Con_90
          else
              Angle = Con_270
          endif
      elseif(abs(X2-X1)>Tol_11)then
          if(X1 <= X2) then
              if(Y1<=Y2)then
                  Angle = c_Angle
              else
                  Angle = Con_360 - abs(c_Angle)
              endif
          elseif(X1 > X2) then
              if(Y1<=Y2)then
                  Angle = Con_180 - abs(c_Angle)
              else
                  Angle = Con_180 + abs(c_Angle)
              endif
          endif
      endif
      
#ifndef Silverfrost
      if (isnan(Angle)) then
          print *, '    Error :: Tool_Angle_of_Inclination'
     &           //'_of_VectorAB.f'
      endif
#endif     
      
      return 
      end SUBROUTINE Tool_Angle_of_Inclination_of_VectorAB                  
