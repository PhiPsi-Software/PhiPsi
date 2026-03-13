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
 
      subroutine Cal_Coor_by_KesiYita(kesi,yita,X_NODES,Y_NODES,
     &                                Out_x,Out_y)
C     Calculate the coordinates in the global coordinate system based on the coordinates in the local coordinate system
      use Global_Float_Type
      implicit none
      real(kind=FT),intent(in)::kesi,yita,X_NODES(4),Y_NODES(4)
      real(kind=FT),intent(out)::Out_x,Out_y
      real(kind=FT) N(4)
      
      N  = 0.25D0 * [(1-kesi)*(1-yita),(1+kesi)*(1-yita),
     &               (1+kesi)*(1+yita),(1-kesi)*(1+yita)]
      
      Out_x  = N(1)*X_NODES(1)+N(2)*X_NODES(2)+
     &         N(3)*X_NODES(3)+N(4)*X_NODES(4)
      Out_y  = N(1)*Y_NODES(1)+N(2)*Y_NODES(2)+
     &         N(3)*Y_NODES(3)+N(4)*Y_NODES(4)
     
      return 
      end SUBROUTINE Cal_Coor_by_KesiYita                          
