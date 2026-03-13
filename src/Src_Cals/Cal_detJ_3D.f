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
 
      subroutine Cal_detJ_3D(kesi,yita,zeta,
     &                    X_NODES,Y_NODES,Z_NODES,
     &                    detJ)
     
c     This function calculates determinant of Jacobian matrix.

      use Global_Float_Type
      implicit none
      
      real(kind=FT),intent(in)::kesi,yita,zeta
      real(kind=FT),intent(in)::X_NODES(8),Y_NODES(8),Z_NODES(8)
      real(kind=FT),intent(out)::detJ       
      real(kind=FT) temp(3,8),Coor(3,8),J(3,3),dNdkesi(8,3)
      
      Coor(1,:) = [X_NODES(1),X_NODES(2),X_NODES(3),X_NODES(4),
     &             X_NODES(5),X_NODES(6),X_NODES(7),X_NODES(8)]
      Coor(2,:) = [Y_NODES(1),Y_NODES(2),Y_NODES(3),Y_NODES(4),
     &             Y_NODES(5),Y_NODES(6),Y_NODES(7),Y_NODES(8)]
      Coor(3,:) = [Z_NODES(1),Z_NODES(2),Z_NODES(3),Z_NODES(4),
     &             Z_NODES(5),Z_NODES(6),Z_NODES(7),Z_NODES(8)]
      
      temp(1,1:8)=[ -(ONE - yita)*(ONE - zeta)/EIG,
     &               (ONE - yita)*(ONE - zeta)/EIG,
     &               (ONE + yita)*(ONE - zeta)/EIG,
     &              -(ONE + yita)*(ONE - zeta)/EIG,
     &              -(ONE - yita)*(ONE + zeta)/EIG,
     &               (ONE - yita)*(ONE + zeta)/EIG,
     &               (ONE + yita)*(ONE + zeta)/EIG,
     &              -(ONE + yita)*(ONE + zeta)/EIG]  
     
      temp(2,1:8)=[ -(ONE - kesi)*(ONE - zeta)/EIG,
     &              -(ONE + kesi)*(ONE - zeta)/EIG,
     &               (ONE + kesi)*(ONE - zeta)/EIG,
     &               (ONE - kesi)*(ONE - zeta)/EIG,
     &              -(ONE - kesi)*(ONE + zeta)/EIG,
     &              -(ONE + kesi)*(ONE + zeta)/EIG,
     &               (ONE + kesi)*(ONE + zeta)/EIG,
     &               (ONE - kesi)*(ONE + zeta)/EIG] 
     
      temp(3,1:8)=[ -(ONE - kesi)*(ONE - yita)/EIG,
     &              -(ONE + kesi)*(ONE - yita)/EIG,
     &              -(ONE + kesi)*(ONE + yita)/EIG,
     &              -(ONE - kesi)*(ONE + yita)/EIG,
     &               (ONE - kesi)*(ONE - yita)/EIG,
     &               (ONE + kesi)*(ONE - yita)/EIG,
     &               (ONE + kesi)*(ONE + yita)/EIG,
     &               (ONE - kesi)*(ONE + yita)/EIG]   
              
      dNdkesi = transpose(temp)

      J = MATMUL(Coor,dNdkesi)     

      call Matrix_Det_3x3(J,detJ) 
      
      
      return 
      end SUBROUTINE Cal_detJ_3D                  
