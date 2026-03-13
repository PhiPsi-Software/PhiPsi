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
 
      subroutine Cal_Point_Aperture_and_Tangent_disp_NoFEM(c_Crack,
     &                 Point,U_e,Cr_Omega,Aperture,Tangent_disp)
      ! Calculate the aperture at a certain point on the crack (that is, the normal displacement
      ! difference) and the tangential displacement difference (Tangent_disp)
      !--------------
      ! The input degrees of freedom displacements only include the displacements of the enhanced degrees
      ! of freedom.
      
      !+++++++++++++++++++++++++++++
      ! Read public variable module
      !+++++++++++++++++++++++++++++
      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Model
      use Global_Elem_Area_Vol
      
      implicit none
      real(kind=FT),intent(in)::Point(2),U_e(Enrich_Freedom),Cr_Omega
      integer,intent(in)::c_Crack
      real(kind=FT),intent(out)::Aperture,Tangent_disp
      real(kind=FT) Up_Point(2),Low_Point(2),Cr_p1(2),Cr_p2(2),
     &                 Up_Disp(2),Low_Disp(2),delta_L,
     &                 Delta_Disp(2),cos_angle,x1,y1,x2,y2
      real(kind=FT) Vector_a_b(2),Vector_ap_bp(2),sign_Aper,
     &                 Magni_Aper,Vector_n_Up(2),Vector_n_Low(2),
     &                 Up_Disp_Pro(2),Low_Disp_Pro(2),Disp_Normal(2)
      integer c_Elem
      real(kind=FT) c_Kesi,c_Yita
      real(kind=FT) Vector_t(2),Disp_Tang(2),sign_Tang_u,Magni_Tang_u
      real(kind=FT) tem_value1,tem_value2
      
      delta_L = Delta_Factor_Aper*Ave_Elem_L_Enrich
      
      Up_Point(1) = Point(1)  -  delta_L*sin(Cr_Omega)
      Up_Point(2) = Point(2)  +  delta_L*cos(Cr_Omega)
      
      Low_Point(1) = Point(1)  +  delta_L*sin(Cr_Omega)
      Low_Point(2) = Point(2)  -  delta_L*cos(Cr_Omega) 
      
      call Cal_Ele_Num_by_Coors(Up_Point(1),Up_Point(2),c_Elem)
      call Cal_KesiYita_by_Coor(Up_Point,c_Elem,c_Kesi,c_Yita)
      call Cal_Any_Point_Disp_KesiYita_NoFEM(c_Elem,c_Kesi,c_Yita,0,U_e,
     &                                 Up_Disp)
      
      call Cal_Ele_Num_by_Coors(Low_Point(1),Low_Point(2),c_Elem)
      call Cal_KesiYita_by_Coor(Low_Point,c_Elem,c_Kesi,c_Yita)
      call Cal_Any_Point_Disp_KesiYita_NoFEM(c_Elem,c_Kesi,c_Yita,0,U_e,
     &                                 Low_Disp)

      Vector_n_Up  = [-sin(Cr_Omega),cos(Cr_Omega)]
      Vector_n_Low = -Vector_n_Up
      Up_Disp_Pro(1:2) = dot_product(Up_Disp,Vector_n_Up)*Vector_n_Up
      Low_Disp_Pro(1:2)= dot_product(Low_Disp,Vector_n_Low)*Vector_n_Low
      Disp_Normal = Up_Disp_Pro - Low_Disp_Pro
      Magni_Aper = sqrt(Disp_Normal(1)**2+Disp_Normal(2)**2)
      x1 = Disp_Normal(1) 
      y1 = Disp_Normal(2) 
      x2 = Vector_n_Up(1) 
      y2 = Vector_n_Up(2) 
      
      
      tem_value1 = x1**2+y1**2
      tem_value2 = x2**2+y2**2
      if (tem_value1==ZR) tem_value1 = Tol_20
      if (tem_value2==ZR) tem_value2 = Tol_20
      cos_angle = (x1*x2+y1*y2)/(sqrt(tem_value1)*sqrt(tem_value2))
      
      if (cos_angle > 0) then
          sign_Aper = ONE
      else
          sign_Aper =-ONE
      end if    
      Aperture = sign_Aper*Magni_Aper
      
      Vector_t  = [cos(Cr_Omega),sin(Cr_Omega)]
      Up_Disp_Pro(1:2) = dot_product(Up_Disp,Vector_t)*Vector_t
      Low_Disp_Pro(1:2)= dot_product(Low_Disp,Vector_t)*Vector_t
      Disp_Tang = Up_Disp_Pro - Low_Disp_Pro
      Magni_Tang_u = sqrt(Disp_Tang(1)**2+Disp_Tang(2)**2)
      x1 = Disp_Tang(1) 
      y1 = Disp_Tang(2) 
      x2 = Vector_t(1) 
      y2 = Vector_t(2) 
      
      
      tem_value1 = x1**2+y1**2
      tem_value2 = x2**2+y2**2
      if (tem_value1==ZR) tem_value1 = Tol_20
      if (tem_value2==ZR) tem_value2 = Tol_20
      cos_angle = (x1*x2+y1*y2)/(sqrt(tem_value1)*sqrt(tem_value2))
      
      if (cos_angle > 0) then
          sign_Tang_u = ONE
      else
          sign_Tang_u =-ONE
      end if    
      Tangent_disp = sign_Tang_u*Magni_Tang_u
      
      return 
      end SUBROUTINE Cal_Point_Aperture_and_Tangent_disp_NoFEM               
