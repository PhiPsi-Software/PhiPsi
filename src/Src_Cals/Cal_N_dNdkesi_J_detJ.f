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
 
      subroutine Cal_N_dNdkesi_J_detJ(kesi,yita,
     &                                X_NODES,Y_NODES,
     &                                detJ,J,N,dNdkesi)
     
c     This function calculates N, dNdkesi, J and the determinant of Jacobian matrix.
      use Global_Float_Type      
      implicit none
      
      real(kind=FT),intent(in)::kesi,yita,X_NODES(4),Y_NODES(4)
      real(kind=FT),intent(out)::detJ,dNdkesi(4,2)
      real(kind=FT),intent(out):: N(2,8)
      real(kind=FT) N1,N2,N3,N4
      real(kind=FT) o
      
      real(kind=FT) temp(2,4),Coor(2,4),J(2,2) 
      
      o   = ZR
      
      Coor(1,:) = [X_NODES(1),X_NODES(2),X_NODES(3),X_NODES(4)]
      Coor(2,:) = [Y_NODES(1),Y_NODES(2),Y_NODES(3),Y_NODES(4)]
      
      N1 = (ONE-kesi)*(ONE-yita)/FOU
      N2 = (ONE+kesi)*(ONE-yita)/FOU
      N3 = (ONE+kesi)*(ONE+yita)/FOU
      N4 = (ONE-kesi)*(ONE+yita)/FOU
      N(1,:) = [N1,o,N2,o,N3,o,N4,o]
      N(2,:) = [o,N1,o,N2,o,N3,o,N4]

      temp(1,1:4) =[(yita-ONE)/FOU,(-yita+ONE)/FOU,
     &                    (yita+ONE)/FOU,(-yita-ONE)/FOU]
      temp(2,1:4) =[(kesi-ONE)/FOU,(-kesi-ONE)/FOU,
     &                    (kesi+ONE)/FOU,(-kesi+ONE)/FOU]
              
      dNdkesi = transpose(temp)

      J = MATMUL(Coor,dNdkesi)
      
      call Matrix_Det_2x2(J,detJ)  
      
      return 
      end SUBROUTINE Cal_N_dNdkesi_J_detJ                     
