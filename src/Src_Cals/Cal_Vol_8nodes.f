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
 
      subroutine Cal_Vol_8nodes(i_Elem,x,y,z,vol)
       
      ! Calculate the volume of a 3D 8-node hexahedral element
      use Global_Float_Type
      IMPLICIT NONE
      integer ,intent(in)::i_Elem
      real(kind=FT),intent(in):: x(8)
      real(kind=FT),intent(in):: y(8)
      real(kind=FT),intent(in):: z(8)
      real(kind=FT) ,intent(out)::vol
      real(kind=FT) A(3),B(3),C(3),D(3),E(3),F(3),G(3),H(3)
     
      A(1) = x(8);A(2) = y(8);A(3) = z(8)
      B(1) = x(5);B(2) = y(5);B(3) = z(5)
      C(1) = x(6);C(2) = y(6);C(3) = z(6)
      D(1) = x(7);D(2) = y(7);D(3) = z(7)
      
      E(1) = x(4);E(2) = y(4);E(3) = z(4)
      F(1) = x(1);F(2) = y(1);F(3) = z(1)      
      G(1) = x(2);G(2) = y(2);G(3) = z(2)
      H(1) = x(3);H(2) = y(3);H(3) = z(3)  
      
      call Tool_Volume_Hexahedron(A,B,C,D,E,F,G,H,vol)

      RETURN
      END subroutine Cal_Vol_8nodes
