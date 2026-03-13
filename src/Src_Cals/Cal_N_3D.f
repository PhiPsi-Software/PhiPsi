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
 
      subroutine Cal_N_3D(kesi,yita,zeta,N)
     
c     This function calculates N, dNdkesi, J and the determinant of Jacobian matrix.
      use Global_Float_Type      
      implicit none
      
      real(kind=FT),intent(in)::kesi,yita,zeta
      real(kind=FT),intent(out):: N(3,24)
      
      real(kind=FT) N1,N2,N3,N4,N5,N6,N7,N8
      real(kind=FT) o
      
      o      = ZR
      
      N1 = (ONE-kesi)*(ONE-yita)*(ONE-zeta)/EIG
      N2 = (ONE+kesi)*(ONE-yita)*(ONE-zeta)/EIG
      N3 = (ONE+kesi)*(ONE+yita)*(ONE-zeta)/EIG
      N4 = (ONE-kesi)*(ONE+yita)*(ONE-zeta)/EIG
      N5 = (ONE-kesi)*(ONE-yita)*(ONE+zeta)/EIG
      N6 = (ONE+kesi)*(ONE-yita)*(ONE+zeta)/EIG
      N7 = (ONE+kesi)*(ONE+yita)*(ONE+zeta)/EIG
      N8 = (ONE-kesi)*(ONE+yita)*(ONE+zeta)/EIG
      
      N(1,1:24) =
     &        [N1,o,o,N2,o,o,N3,o,o,N4,o,o,N5,o,o,N6,o,o,N7,o,o,N8,o,o]
      N(2,1:24) =
     &        [o,N1,o,o,N2,o,o,N3,o,o,N4,o,o,N5,o,o,N6,o,o,N7,o,o,N8,o]
      N(3,1:24) =
     &        [o,o,N1,o,o,N2,o,o,N3,o,o,N4,o,o,N5,o,o,N6,o,o,N7,o,o,N8]
      
      return 
      end SUBROUTINE Cal_N_3D               
