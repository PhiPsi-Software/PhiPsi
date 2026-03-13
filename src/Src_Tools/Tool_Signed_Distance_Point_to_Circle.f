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
 
      subroutine Tool_Signed_Distance_Point_to_Circle(
     &                             x0,y0,r,Point_C,S_Distance)
C     Calculate Symbol Distance
c     This function calculates the signed distance from the Point_C to the circle.
c     Negative inside the circle, positive outside the circle, zero on the circle

      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(in)::x0,y0,r,Point_C(2)
      real(kind=FT),intent(out)::S_Distance
      
      real(kind=FT) c_Dis
      
      c_Dis = sqrt((x0-Point_C(1))**2 + (y0-Point_C(2))**2)
      
      if(c_Dis < r)then
          S_Distance = -abs(c_Dis-r)
      elseif(c_Dis > R)then
          S_Distance = abs(c_Dis-r)
      endif
      
      if(abs(c_Dis-R)<1.0D-8*r)then
          S_Distance = ZR
      endif
      
      return 
      end SUBROUTINE Tool_Signed_Distance_Point_to_Circle                        
