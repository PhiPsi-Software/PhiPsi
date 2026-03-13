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
 
      subroutine Cal_Weighted_Ave_Stress_of_Point(c_Key_Ave_Stress,
     &                       Point,Search_R,a_Weight,l_Ordina,
     &                       WA_S_xx,WA_S_yy,WA_S_xy)
C     Calculate the weighted average stress of the point.
C     Stresses of all the Gauss points inside the circle of radius Search_R
c     are calculated.
c     For SHI Fang's Formation, see page 25 of my doctoral dissertation for details.
c     For ordinary Formation, See p498 of book:Extended Finite Element Method: Theory and Applications  
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Stress
      
      implicit none
      integer,intent(in)::c_Key_Ave_Stress
      real(kind=FT),intent(in)::Point(2),Search_R,a_Weight,
     &                             l_Ordina
      real(kind=FT),intent(out)::WA_S_xx,WA_S_yy,WA_S_xy
      integer Tip_Elem,num_Surround_Ele,c_E
      integer i_E,i_G,c_num_Gauss,c_Global_G
      real(kind=FT) x_GP,y_GP,Sxx_GP,Syy_GP,Sxy_GP,c_L
      real(kind=FT)  All_weight,weight,tem1
      integer num_G
      
      call Cal_Ele_Num_by_Coors(Point(1),Point(2),Tip_Elem)
      
      num_Surround_Ele = num_Ele_Eles(Tip_Elem)
      
      All_weight  = ZR
      WA_S_xx     = ZR
      WA_S_yy     = ZR
      WA_S_xy     = ZR
      tem1 = ONE/((TWO*pi)**(1.5D0))/(l_Ordina**3)
      num_G = 0
      do i_E = 1,num_Surround_Ele
          c_E = Ele_Elements(Tip_Elem,i_E)
          do i_G = 1,num_GP_Elem(c_E)
              c_Global_G = Ele_GP_Start_Num(c_E) + i_G -1
              x_GP = Gauss_CoorX(c_Global_G)
              y_GP = Gauss_CoorY(c_Global_G)
              Sxx_GP = Stress_xx_Gauss(c_Global_G)
              Syy_GP = Stress_yy_Gauss(c_Global_G)
              Sxy_GP = Stress_xy_Gauss(c_Global_G)
              c_L = sqrt((x_GP-Point(1))**2 + (y_GP-Point(2))**2)
              if (c_L <= Search_R) then
                  if(c_Key_Ave_Stress==1)then
                      weight = (ONE-((c_L/Search_R)**a_Weight))**3
                  elseif(c_Key_Ave_Stress==2)then
                      weight = tem1*exp(-c_L**2/(TWO*l_Ordina**2))
                  endif
                  All_weight = All_weight + weight
                  num_G = num_G + 1
                  WA_S_xx     = WA_S_xx + weight * Sxx_GP
                  WA_S_yy     = WA_S_yy + weight * Syy_GP
                  WA_S_xy     = WA_S_xy + weight * Sxy_GP
              endif
          enddo
      enddo
      WA_S_xx = WA_S_xx /All_weight
      WA_S_yy = WA_S_yy /All_weight
      WA_S_xy = WA_S_xy /All_weight
      RETURN
      END SUBROUTINE Cal_Weighted_Ave_Stress_of_Point