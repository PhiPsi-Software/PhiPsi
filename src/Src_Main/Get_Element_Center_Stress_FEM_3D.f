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
 
      SUBROUTINE Get_Element_Center_Stress_FEM_3D(Yes_Add_Insitu,
     &        isub,c_DISP,
     &        Stress_El_Center)
c     Stress at the center of the computational unit (FEM)
c     2020-03-13
c     Stress_El_Center(num_Ele,1) - xx
c     Stress_El_Center(num_Ele,2) - yy
c     Stress_El_Center(num_Ele,3) - zz
c     Stress_El_Center(num_Ele,4) - xy
c     Stress_El_Center(num_Ele,5) - yz
c     Stress_El_Center(num_Ele,6) - xz
c     Stress_El_Center(num_Ele,7) - vm

c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Dynamic
      use Global_Material
      use Global_Stress
      use Global_Disp
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      integer, intent(in)::isub
      logical, intent(in)::Yes_Add_Insitu
      real(kind=FT), intent(in)::c_DISP(num_Node*3)
      real(kind=FT), intent(out)::Stress_El_Center(num_Elem,7)       
      real(kind=FT) c_T_Alpha,c_TStress(6)     
      integer i_E
      real(kind=FT) c_D(6,6),U(24)
      real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8),
     &              c_kesi,c_yita,c_zeta,c_Stress(6)
      integer c_NN(8)
      real(kind=FT) Str_InSitu_El(6)
      real(kind=FT) kesi(Num_Gauss_P_FEM_3D),yita(Num_Gauss_P_FEM_3D),
     &               zeta(Num_Gauss_P_FEM_3D),weight(Num_Gauss_P_FEM_3D)
     
      
c     -----------------
C     Initialized to 0
c     -----------------
      Stress_El_Center(1:num_Elem,1:7) = ZR
     
c     -----------------
C     Inter-unit cycle
c     -----------------
      do i_E = 1,Num_Elem
          ! Not related to New Ao.
          if(Key_XA/=2) then
              c_D     = D(Elem_Mat(i_E),1:6,1:6)
              c_T_Alpha = T_Alpha(Elem_Mat(i_E))
          ! Xin Ao. Key_XA==2. IMPROV2023031905.
          else
              c_D     = Elem_D_XA(i_E,1:6,1:6)
              c_T_Alpha = Elem_TEC_XA(i_E)
          endif
          !c_NN    = G_NN(i_E,1:8)
          !c_X_NODES = G_X_NODES(i_E,1:8)
          !c_Y_NODES = G_Y_NODES(i_E,1:8)    
          !c_Z_NODES = G_Z_NODES(i_E,1:8)   
          c_NN    = G_NN(1:8,i_E)
          c_X_NODES = G_X_NODES(1:8,i_E)
          c_Y_NODES = G_Y_NODES(1:8,i_E)    
          c_Z_NODES = G_Z_NODES(1:8,i_E)            
          U =
     &      [c_DISP(c_NN(1)*3-2),c_DISP(c_NN(1)*3-1),c_DISP(c_NN(1)*3),
     &       c_DISP(c_NN(2)*3-2),c_DISP(c_NN(2)*3-1),c_DISP(c_NN(2)*3),
     &       c_DISP(c_NN(3)*3-2),c_DISP(c_NN(3)*3-1),c_DISP(c_NN(3)*3),
     &       c_DISP(c_NN(4)*3-2),c_DISP(c_NN(4)*3-1),c_DISP(c_NN(4)*3),
     &       c_DISP(c_NN(5)*3-2),c_DISP(c_NN(5)*3-1),c_DISP(c_NN(5)*3),
     &       c_DISP(c_NN(6)*3-2),c_DISP(c_NN(6)*3-1),c_DISP(c_NN(6)*3),
     &       c_DISP(c_NN(7)*3-2),c_DISP(c_NN(7)*3-1),c_DISP(c_NN(7)*3),
     &       c_DISP(c_NN(8)*3-2),c_DISP(c_NN(8)*3-1),c_DISP(c_NN(8)*3)]    
          c_kesi = ZR                                              
          c_yita = ZR
          c_zeta = ZR
     
          call Cal_Ele_Str_N8_3D(i_E,0,1,1,
     &                           c_X_NODES,c_Y_NODES,c_Z_NODES,
     &                           c_D,c_kesi,c_yita,c_zeta,U,
     &                           c_Stress)   
          Stress_El_Center(i_E,1) =Stress_El_Center(i_E,1)+c_Stress(1)
          Stress_El_Center(i_E,2) =Stress_El_Center(i_E,2)+c_Stress(2)
          Stress_El_Center(i_E,3) =Stress_El_Center(i_E,3)+c_Stress(3)
          Stress_El_Center(i_E,4) =Stress_El_Center(i_E,4)+c_Stress(4)
          Stress_El_Center(i_E,5) =Stress_El_Center(i_E,5)+c_Stress(5)
          Stress_El_Center(i_E,6) =Stress_El_Center(i_E,6)+c_Stress(6)
             
          ! Subtract thermal expansion stress, Theory: Equation 15.1.98 from 'Fundamentals of Finite Element
          ! Method (5th Edition)', 2019-09-24
          if(Key_Thermal_Stress==1)then
              !IMPROV2023031302.
              c_TStress = c_T_Alpha*Elem_T_for_Stress(i_E)*
     &               MATMUL(c_D,[ONE,ONE,ONE,ZR,ZR,ZR]) 
              Stress_El_Center(i_E,1)=Stress_El_Center(i_E,1)-
     &                                c_TStress(1)
              Stress_El_Center(i_E,2)=Stress_El_Center(i_E,2)-
     &                                c_TStress(2)
              Stress_El_Center(i_E,3)=Stress_El_Center(i_E,3)-
     &                                c_TStress(3)
              Stress_El_Center(i_E,4)=Stress_El_Center(i_E,4)-
     &                                c_TStress(4)
              Stress_El_Center(i_E,5)=Stress_El_Center(i_E,5)-
     &                                c_TStress(5)
              Stress_El_Center(i_E,6)=Stress_El_Center(i_E,6)-
     &                                c_TStress(6)
          endif


          ! Add the initial stress field
          if(Key_InSitu_Strategy==2 .and. Yes_Add_Insitu)then
              Str_InSitu_El(1) = sum(Str_xx_InSitu(c_NN))/8
              Str_InSitu_El(2) = sum(Str_yy_InSitu(c_NN))/8
              Str_InSitu_El(3) = sum(Str_zz_InSitu(c_NN))/8
              Str_InSitu_El(4) = sum(Str_xy_InSitu(c_NN))/8
              Str_InSitu_El(5) = sum(Str_yz_InSitu(c_NN))/8
              Str_InSitu_El(6) = sum(Str_xz_InSitu(c_NN))/8
              Stress_El_Center(i_E,1)=Stress_El_Center(i_E,1)+
     &                                Str_InSitu_El(1)
              Stress_El_Center(i_E,2)=Stress_El_Center(i_E,2)+
     &                                Str_InSitu_El(2)
              Stress_El_Center(i_E,3)=Stress_El_Center(i_E,3)+
     &                                Str_InSitu_El(3)
              Stress_El_Center(i_E,4)=Stress_El_Center(i_E,4)+
     &                                Str_InSitu_El(4)
              Stress_El_Center(i_E,5)=Stress_El_Center(i_E,5)+
     &                                Str_InSitu_El(5)
              Stress_El_Center(i_E,6)=Stress_El_Center(i_E,6)+
     &                                Str_InSitu_El(6)  
          endif
          
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Subtract the stress corresponding to the initial stress field, Theory: Equation 15.1.98
          ! from 'Fundamentals of the Finite Element Method (5th Edition)'
          !2022-06-03
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(Key_InSitu_Strategy==4)then
              ! Obtain the Gauss point stress based on the initial strain. Similar to the handling method of
              ! thermal stress.
              c_Stress = MATMUL(c_D,[ InSitu_Strain_Gaus_xx(i_E,1),
     &                                InSitu_Strain_Gaus_yy(i_E,1),
     &                                InSitu_Strain_Gaus_zz(i_E,1),
     &                                InSitu_Strain_Gaus_xy(i_E,1),
     &                                InSitu_Strain_Gaus_yz(i_E,1),
     &                                InSitu_Strain_Gaus_xz(i_E,1)])   
              Stress_El_Center(i_E,1) = Stress_El_Center(i_E,1) - 
     &                                        c_Stress(1)
              Stress_El_Center(i_E,2) = Stress_El_Center(i_E,2) - 
     &                                        c_Stress(2)  
              Stress_El_Center(i_E,3) = Stress_El_Center(i_E,3) - 
     &                                        c_Stress(3)
              Stress_El_Center(i_E,4) = Stress_El_Center(i_E,4) - 
     &                                        c_Stress(4)
              Stress_El_Center(i_E,5) = Stress_El_Center(i_E,5) - 
     &                                        c_Stress(5)
              Stress_El_Center(i_E,6) = Stress_El_Center(i_E,6) - 
     &                                        c_Stress(6)     
          endif
          
          ! Calculate Mises equivalent stress
          call Tool_von_Mises_3D(Stress_El_Center(i_E,1),
     &                           Stress_El_Center(i_E,2),
     &                           Stress_El_Center(i_E,3),
     &                           Stress_El_Center(i_E,4),
     &                           Stress_El_Center(i_E,5),
     &                           Stress_El_Center(i_E,6),     
     &                           Stress_El_Center(i_E,7))
      end do  

     
      RETURN
      END SUBROUTINE Get_Element_Center_Stress_FEM_3D
