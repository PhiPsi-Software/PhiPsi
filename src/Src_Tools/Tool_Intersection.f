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
 
      subroutine Tool_Intersection(A,B,C,D,X,Y,Yes_Cross)
C     Calculate the intersection point of a line segment
c     A: Coordinates of point A of line segment AB;
c     B: Coordinates of point B of line segment AB;
c     C: Coordinates of point C of line segment CD;
c     D: Coordinates of point D of line segment CD.
      use Global_Float_Type
      use Global_Elem_Area_Vol
      
      implicit none
      real(kind=FT),intent(in)::A(2),B(2),C(2),D(2)
      real(kind=FT),intent(out)::X,Y
      logical,intent(out)::Yes_Cross
      
      real(kind=FT) L_AB,L_CD,tem,k2,k1,b1,b2
      real(kind=FT) L_AX,L_BX,L_CX,L_DX
      logical Yes_in_AB,Yes_in_CD

      logical Yes_A1_eq_B1,Yes_C1_eq_D1
      real(kind=FT) Tiny_Value
      
      Tiny_Value = Tol_10
      
      Yes_A1_eq_B1 = .False.
      if(abs(A(1)-B(1)) <= Tiny_Value)then
          Yes_A1_eq_B1 = .True.
      endif
      
      Yes_C1_eq_D1 = .False.
      if(abs(C(1)-D(1)) <= Tiny_Value)then
          Yes_C1_eq_D1 = .True.
      endif
      
      L_AB = sqrt((A(1)-B(1))**2+(A(2)-B(2))**2)
      L_CD = sqrt((C(1)-D(1))**2+(C(2)-D(2))**2)
      
      if (Yes_A1_eq_B1) then
          X=A(1)
          
          tem = D(1)-C(1)
          
          if(tem==ZR) then
              tem = Tol_30
          endif
          
          k2=(D(2)-C(2))/tem
          b2=C(2)-k2*C(1) 
          Y=k2*X+b2
      end if
      if (Yes_C1_eq_D1) then
          X=C(1)
          
          tem =B(1)-A(1)
          
          if(tem==ZR) then
              tem = Tol_30
          endif
          
          k1=(B(2)-A(2))/tem
          b1=A(2)-k1*A(1)
          Y=k1*X+b1
      end if

      if ((Yes_A1_eq_B1.eqv..False.) .and.
     &    (Yes_C1_eq_D1.eqv..False.)) then
          k1=(B(2)-A(2))/(B(1)-A(1))
          k2=(D(2)-C(2))/(D(1)-C(1))
          b1=A(2)-k1*A(1)
          b2=C(2)-k2*C(1)
          if (k1.eq.k2) then
              X=1.0D9
              Y=1.0D9
          else
              X=(b2-b1)/(k1-k2)
              Y=k1*X+b1
          end if
      end if

      Yes_Cross = .False.
      Yes_in_AB = .False.
      Yes_in_CD = .False.

      L_AX = sqrt((A(1)-X)**2+(A(2)-Y)**2)
      L_BX = sqrt((B(1)-X)**2+(B(2)-Y)**2)
      if (((L_AX + L_BX)-L_AB) .le. Tiny_Value) then
           Yes_in_AB = .True.
      end if


      L_CX = sqrt((C(1)-X)**2+(C(2)-Y)**2)
      L_DX = sqrt((D(1)-X)**2+(D(2)-Y)**2)
      if (((L_CX + L_DX) - L_CD) .le.Tiny_Value) then
           Yes_in_CD = .True.
      end if

      if ((Yes_in_AB .eqv..True.).and.(Yes_in_CD.eqv..True.)) then 
          Yes_Cross = .True.
      else
          Yes_Cross = .False.
      end if
      
      return 
      end SUBROUTINE Tool_Intersection                          
