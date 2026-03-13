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
 
      subroutine Tool_Point_Giv_2P_and_L(A,B,L,Out_C,Statu,L_BC)
C     Calculate the coordinates of the point, where point C is at a distance L from A, along the direction of AB.

c     case1:Statu = -1
c  
C     *-------------------*----------*
C     A                   C          B

c     case2:Statu = 1
c  
C     *---------*---------*
C     A         B         C

c     case3: Status = 0, point B and point C coincide
c  
C     *-------------------*
C     A                  B(C)

      use Global_Float_Type     
      implicit none
      real(kind=FT),intent(in)::A(2),B(2),L
      real(kind=FT),intent(out)::Out_C(2),L_BC
      integer,intent(out)::Statu
      real(kind=FT) L_AB,theta_AB
     
      L_AB = sqrt((A(1)-B(1))**2+(A(2)-B(2))**2)
      if(L_AB > L)then
          Statu = -1
      elseif(L_AB < L)then
          Statu = 1
      else
          Statu = 0
      endif
      theta_AB = atan2(B(2)-A(2),B(1)-A(1))

      Out_C(1) = A(1)+L*cos(theta_AB)
      Out_C(2) = A(2)+L*sin(theta_AB)
      
      L_BC = sqrt((Out_C(1)-B(1))**2+(Out_C(2)-B(2))**2)
      
      return 
      end SUBROUTINE Tool_Point_Giv_2P_and_L                          
