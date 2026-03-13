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
 
SUBROUTINE Get_Gauss_Stress_FEM_3D(isub,c_DISP)
!Computational Gauss Stress for 3D FEM problems.
!Store in global variables: Stress_xx_Gauss, Stress_yy_Gauss, Stress_zz_Gauss, 
!                           Stress_xy_Gauss, Stress_yz_Gauss, Stress_xz_Gauss, 
!                           Stress_vm_Gauss
!2026-02-02.

!----------------------------
!Read public variable module
!----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Dynamic
use Global_Material
use Global_Stress
use Global_Disp
use omp_lib

!--------------------------
!Variable Type Declaration
!--------------------------
implicit none
integer, intent(in)::isub
real(kind=FT), intent(in)::c_DISP(num_Node*3)
!include 'omp_lib.h'
integer i_E
real(kind=FT) c_D(6,6),U(24)
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
real(kind=FT) T_Kesi(8),T_Yita(8),T_Zeta(8)
real(kind=FT) c_kesi,c_yita,c_zeta
real(kind=FT) c_Stress(6)
integer c_NN(8),i_N,C_Node,c_Count(num_Node) 
real(kind=FT) kesi_N_Enr(Num_Gauss_P_FEM_3D),yita_N_Enr(Num_Gauss_P_FEM_3D)
real(kind=FT) zeta_N_Enr(Num_Gauss_P_FEM_3D),weight_N_Enr(Num_Gauss_P_FEM_3D)    
integer G_Counter,i_G
real(kind=FT) kesi(Num_Gauss_P_FEM_3D),yita(Num_Gauss_P_FEM_3D),zeta(Num_Gauss_P_FEM_3D)
real(kind=FT) c_S_1,c_S_2,c_S_3,c_Ele_Gauss(6)
real(kind=FT) Vector_S1(3),Vector_S2(3),Vector_S3(3)

print *,'    Calculating stress of all Gauss points...'


call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,kesi_N_Enr,yita_N_Enr,zeta_N_Enr,weight_N_Enr)

G_Counter = 0

!-----------------
!Inter-unit cycle
!-----------------
!.............................
! OpenMP multi-core computing
!.............................
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,i_G,c_Stress,Vector_S1,Vector_S2,Vector_S3, &
!$OMP              c_D,c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES,U,kesi,yita,zeta, &
!$OMP              c_kesi,c_yita,c_zeta,G_Counter,c_Ele_Gauss,c_S_1,c_S_2,c_S_3) 
do i_E = 1,Num_Elem
  c_D(1:6,1:6)     = D(Elem_Mat(i_E),1:6,1:6)    
  !c_v     = Material_Para(Elem_Mat(i_E),2)                 ! Poisson's ratio

  !Traditional index locations
  c_NN    = G_NN(1:8,i_E)
  c_X_NODES = G_X_NODES(1:8,i_E)
  c_Y_NODES = G_Y_NODES(1:8,i_E)    
  c_Z_NODES = G_Z_NODES(1:8,i_E)        
  
  U = &
    [c_DISP(c_NN(1)*3-2),c_DISP(c_NN(1)*3-1),c_DISP(c_NN(1)*3), &
     c_DISP(c_NN(2)*3-2),c_DISP(c_NN(2)*3-1),c_DISP(c_NN(2)*3),&
     c_DISP(c_NN(3)*3-2),c_DISP(c_NN(3)*3-1),c_DISP(c_NN(3)*3),&
     c_DISP(c_NN(4)*3-2),c_DISP(c_NN(4)*3-1),c_DISP(c_NN(4)*3),&
     c_DISP(c_NN(5)*3-2),c_DISP(c_NN(5)*3-1),c_DISP(c_NN(5)*3),&
     c_DISP(c_NN(6)*3-2),c_DISP(c_NN(6)*3-1),c_DISP(c_NN(6)*3),&
     c_DISP(c_NN(7)*3-2),c_DISP(c_NN(7)*3-1),c_DISP(c_NN(7)*3),&
     c_DISP(c_NN(8)*3-2),c_DISP(c_NN(8)*3-1),c_DISP(c_NN(8)*3)]
     
  c_Ele_Gauss(1:6) = ZR
  
  kesi(1:Num_Gauss_P_FEM_3D)   = kesi_N_Enr
  yita(1:Num_Gauss_P_FEM_3D)   = yita_N_Enr
  zeta(1:Num_Gauss_P_FEM_3D)   = zeta_N_Enr
  
  ! Loop over each Gauss point.
  do i_G = 1,Num_Gauss_P_FEM_3D
      G_Counter =  (i_E-1)*Num_Gauss_P_FEM_3D +i_G 
      c_kesi = kesi(i_G)                                          
      c_yita = yita(i_G)
      c_zeta = zeta(i_G)

    !First 1 denotes i_N, second 1 indicates calculating stress, third 1 indicates the Cartesian
    !coordinate system.
    call Cal_Ele_Str_N8_3D(i_E,1,1,1,     &    
                             c_X_NODES,c_Y_NODES,c_Z_NODES, &
                             c_D,c_kesi,c_yita,c_zeta,U, &
                             c_Stress)   
                             
      Stress_xx_Gauss(G_Counter) = c_Stress(1)
      Stress_yy_Gauss(G_Counter) = c_Stress(2)
      Stress_zz_Gauss(G_Counter) = c_Stress(3)
      Stress_xy_Gauss(G_Counter) = c_Stress(4)
      Stress_yz_Gauss(G_Counter) = c_Stress(5)
      Stress_xz_Gauss(G_Counter) = c_Stress(6)
      
      ! Accumulation of stresses at each Gauss point of the current element (used to calculate the average
      ! Gauss stress of each element)
      c_Ele_Gauss(1) = c_Ele_Gauss(1)+c_Stress(1)
      c_Ele_Gauss(2) = c_Ele_Gauss(2)+c_Stress(2)
      c_Ele_Gauss(3) = c_Ele_Gauss(3)+c_Stress(3)
      c_Ele_Gauss(4) = c_Ele_Gauss(4)+c_Stress(4)
      c_Ele_Gauss(5) = c_Ele_Gauss(5)+c_Stress(5)
      c_Ele_Gauss(6) = c_Ele_Gauss(6)+c_Stress(6)
      
      ! Thermal stress, Theory: Equation 15.1.98 from the Finite Element Method Basic Tutorial (5th
      ! Edition), 2019-09-24
      if(Key_Thermal_Stress==1)then
           !To be done!
      endif

      
      ! Add initial stress field.
      if(Key_InSitu_Strategy==2)then
           !To be done.
      endif
      
      ! Subtract the stress corresponding to the initial strain field, theory: Equation 15.1.98 in
      ! 'Fundamentals of the Finite Element Method (5th Edition)'.
      if(Key_InSitu_Strategy==4)then
          !To be done.  
      endif      
      
      ! If the element is killed or destroyed, the Gauss stress is 0.
      if (allocated(Elem_Break)) then
          if(Elem_Break(i_E)) then
              Stress_xx_Gauss(G_Counter) = ZR
              Stress_yy_Gauss(G_Counter) = ZR
              Stress_zz_Gauss(G_Counter) = ZR
              Stress_xy_Gauss(G_Counter) = ZR
              Stress_yz_Gauss(G_Counter) = ZR
              Stress_xz_Gauss(G_Counter) = ZR
              c_Ele_Gauss(1:6) = ZR
          endif  
      endif
      call Tool_von_Mises_3D(Stress_xx_Gauss(G_Counter),Stress_yy_Gauss(G_Counter),Stress_zz_Gauss(G_Counter), &
                         Stress_xy_Gauss(G_Counter),Stress_yz_Gauss(G_Counter),Stress_xz_Gauss(G_Counter), &    
                         Stress_vm_Gauss(G_Counter))

  end do  
  ! Average Gauss stress of each element
  if (allocated(Elem_Ave_Gauss_Stress)) then
    Elem_Ave_Gauss_Stress(i_E,1:6) = c_Ele_Gauss(1:6)/Num_Gauss_P_FEM_3D
    ! Calculate the average principal stress for each element
    call Tool_Principal_Stresses_3D(Elem_Ave_Gauss_Stress(i_E,1),Elem_Ave_Gauss_Stress(i_E,2), &
                                      Elem_Ave_Gauss_Stress(i_E,3),Elem_Ave_Gauss_Stress(i_E,4), &
                                      Elem_Ave_Gauss_Stress(i_E,5),Elem_Ave_Gauss_Stress(i_E,6), &
                                      c_S_1,c_S_2,c_S_3, &
                                      Vector_S1,Vector_S2,Vector_S3)
    Elem_Ave_Gauss_S1(i_E) =  c_S_1          
  endif
end do  
!$omp end parallel do        

RETURN
END SUBROUTINE Get_Gauss_Stress_FEM_3D
