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
 
      subroutine Tool_Volume_Pyramid(A,B,C,D,E,Volume)
       
      ! Calculate the volume of a quadrilateral pyramid
      ! Theory, see my notes V3-P138
      use Global_Float_Type
      IMPLICIT NONE

      real(kind=FT),intent(in)::A(3),B(3),C(3),D(3),E(3)
      real(kind=FT),intent(out)::Volume
      
      real(kind=FT) Area_BCDE,Area_BCD,Area_CDE
      real(kind=FT) Distance
      
      call Tool_Area_Tri_3D(B,C,D,Area_BCD)
      call Tool_Area_Tri_3D(C,D,E,Area_CDE)
      
      Area_BCDE = Area_BCD + Area_CDE
      
      
      
      if(Area_BCD>=Tol_20) then
          call Tool_Dis_Point_to_3D_Tri_only_Dis(A,B,C,D,Distance)
      else
          call Tool_Dis_Point_to_3D_Tri_only_Dis(A,C,D,E,Distance)
      endif
      
     
      Volume = ONE/THR*Area_BCDE*abs(Distance)
      
      RETURN
      END subroutine Tool_Volume_Pyramid
