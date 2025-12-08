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
 
SUBROUTINE Get_Node_Stress_FEM_IN_OUT_3D(Yes_Add_Insitu, isub,c_DISP, &
                Stress_xx_N,Stress_yy_N,Stress_zz_N, &
                Stress_xy_N,Stress_yz_N,Stress_xz_N,Stress_vm_N) 
! Computational Node Stress
! Store into global variables: Stress_xx_N, Stress_yy_N, Stress_xy_N, Stress_vm_N
!                              Stress_zz_N,Stress_yz_N,Stress_xz_N

!-----------------------------          
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Dynamic
use Global_Material
use Global_Stress
use Global_Disp
use Global_POST
use omp_lib

!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
integer, intent(in)::isub
logical, intent(in)::Yes_Add_Insitu
real(kind=FT), intent(in)::c_DISP(num_Node*3)
real(kind=FT), intent(out)::Stress_xx_N(num_Node), &
                  Stress_yy_N(num_Node), &
                  Stress_zz_N(num_Node), &
                  Stress_xy_N(num_Node), &
                  Stress_yz_N(num_Node), &
                  Stress_xz_N(num_Node), &
                  Stress_vm_N(num_Node)     
real(kind=FT) c_T_Alpha,c_TStress(6)     
integer i_E
real(kind=FT) c_D(6,6),U(24)
real(kind=FT) c_v   
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8), &
       T_Kesi(8),T_Yita(8),T_Zeta(8),  &
       c_kesi,c_yita,c_zeta,c_Stress(6)
integer c_NN(8),i_N,C_Node,c_Count(num_Node) 
real(kind=FT) c_v_all(Num_Node)
!2024-03-19.
integer i_Thread,max_threads
integer,ALLOCATABLE::c_Count_Threads(:,:) 
real(kind=FT),ALLOCATABLE::Stress_xx_N_Threads(:,:), &
                  Stress_yy_N_Threads(:,:), &
                  Stress_zz_N_Threads(:,:), &
                  Stress_xy_N_Threads(:,:), &
                  Stress_yz_N_Threads(:,:), &
                  Stress_xz_N_Threads(:,:), &
                  Stress_vm_N_Threads(:,:)  


!2023-02-17.
if (Key_Simple_Post==1) return

print *,'    Calculating stress of all nodes...'

!------------------
! Initialized to 0
!------------------
Stress_xx_N(1:num_Node) = ZR
Stress_yy_N(1:num_Node) = ZR
Stress_zz_N(1:num_Node) = ZR
Stress_xy_N(1:num_Node) = ZR
Stress_yz_N(1:num_Node) = ZR
Stress_xz_N(1:num_Node) = ZR
Stress_vm_N(1:num_Node) = ZR      
c_Count(1:num_Node)     = 0

!------------------
! Inter-unit cycle
!------------------
T_Kesi = [-ONE,  ONE,  ONE, -ONE, -ONE,  ONE,  ONE, -ONE]
T_Yita = [-ONE, -ONE,  ONE,  ONE, -ONE, -ONE,  ONE,  ONE]
T_Zeta = [-ONE, -ONE, -ONE, -ONE,  ONE,  ONE,  ONE,  ONE]

!2024-03-19.
!------------------------------------
! OpenMP parallel computing version.
!IMPROV2024031901.   
!------------------------------------
max_threads = omp_get_max_threads()
if (allocated(c_Count_Threads)) deallocate(c_Count_Threads)
if (allocated(Stress_xx_N_Threads)) deallocate(Stress_xx_N_Threads)
if (allocated(Stress_yy_N_Threads)) deallocate(Stress_yy_N_Threads)
if (allocated(Stress_zz_N_Threads)) deallocate(Stress_zz_N_Threads)
if (allocated(Stress_xy_N_Threads)) deallocate(Stress_xy_N_Threads)
if (allocated(Stress_yz_N_Threads)) deallocate(Stress_yz_N_Threads)
if (allocated(Stress_xz_N_Threads)) deallocate(Stress_xz_N_Threads)
if (allocated(Stress_vm_N_Threads)) deallocate(Stress_vm_N_Threads)
ALLOCATE(c_Count_Threads(num_Node,max_threads))   
c_Count_Threads = 0
ALLOCATE(Stress_xx_N_Threads(num_Node,max_threads))   
ALLOCATE(Stress_yy_N_Threads(num_Node,max_threads))   
ALLOCATE(Stress_zz_N_Threads(num_Node,max_threads))   
ALLOCATE(Stress_xy_N_Threads(num_Node,max_threads))   
ALLOCATE(Stress_yz_N_Threads(num_Node,max_threads))   
ALLOCATE(Stress_xz_N_Threads(num_Node,max_threads))   
ALLOCATE(Stress_vm_N_Threads(num_Node,max_threads)) 
Stress_xx_N_Threads = ZR  
Stress_yy_N_Threads = ZR  
Stress_zz_N_Threads = ZR  
Stress_xy_N_Threads = ZR  
Stress_yz_N_Threads = ZR  
Stress_xz_N_Threads = ZR  
Stress_vm_N_Threads = ZR  

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i_Thread,i_E,c_D,c_T_Alpha,c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES, &
!$OMP                                  U,i_N,C_Node,c_kesi,c_yita,c_zeta,c_Stress,c_TStress) 
i_Thread = omp_get_thread_num()+1
!$OMP DO SCHEDULE(static) 
do i_E = 1,Num_Elem
  c_D     = D(Elem_Mat(i_E),1:6,1:6)
  c_T_Alpha = T_Alpha(Elem_Mat(i_E))
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
    ! Node Loop
  do i_N = 1,8
    C_Node = c_NN(i_N)
    c_kesi = T_Kesi(i_N)                                               
    c_yita = T_Yita(i_N)
    c_zeta = T_Zeta(i_N)
    ! For degenerate elements, use the stress at the element center as the stress for the 8 nodes.
    ! Otherwise, a determinant of zero may occur, making it impossible to calculate Inverse_J.
    ! IMPROV2023061403.
    if(Yes_Degenarated_Elem(i_E) .eqv. .True.)then
        c_kesi = ZR
        c_yita = ZR
        c_zeta = ZR
    endif
    call Cal_Ele_Str_N8_3D(i_E,i_N,1,1,     &
                             c_X_NODES,c_Y_NODES,c_Z_NODES, &
                             c_D,c_kesi,c_yita,c_zeta,U, &
                             c_Stress)   
    c_Count_Threads(C_Node,i_Thread) = c_Count_Threads(C_Node,i_Thread) + 1
    Stress_xx_N_Threads(C_Node,i_Thread) = Stress_xx_N_Threads(C_Node,i_Thread) + c_Stress(1)
    Stress_yy_N_Threads(C_Node,i_Thread) = Stress_yy_N_Threads(C_Node,i_Thread) + c_Stress(2)
    Stress_zz_N_Threads(C_Node,i_Thread) = Stress_zz_N_Threads(C_Node,i_Thread) + c_Stress(3)    
    Stress_xy_N_Threads(C_Node,i_Thread) = Stress_xy_N_Threads(C_Node,i_Thread) + c_Stress(4)
    Stress_yz_N_Threads(C_Node,i_Thread) = Stress_yz_N_Threads(C_Node,i_Thread) + c_Stress(5)
    Stress_xz_N_Threads(C_Node,i_Thread) = Stress_xz_N_Threads(C_Node,i_Thread) + c_Stress(6)                
      ! Subtract thermal expansion stress, Theory: Equation 15.1.98 from 'Fundamentals of Finite Element
      ! Method (5th Edition)', 2019-09-24
      if(Key_Thermal_Stress==1)then
         !c_TStress=c_T_Alpha*Thermal_Str_Temper(Elem_Mat(i_E))*MATMUL(c_D,[ONE,ONE,ONE,ZR,ZR,ZR]) 
          !IMPROV2023031302.
          c_TStress=c_T_Alpha*Elem_T_for_Stress(i_E)*MATMUL(c_D,[ONE,ONE,ONE,ZR,ZR,ZR]) 
          Stress_xx_N_Threads(C_Node,i_Thread) = Stress_xx_N_Threads(C_Node,i_Thread)- c_TStress(1)
          Stress_yy_N_Threads(C_Node,i_Thread) = Stress_yy_N_Threads(C_Node,i_Thread)- c_TStress(2)
          Stress_zz_N_Threads(C_Node,i_Thread) = Stress_zz_N_Threads(C_Node,i_Thread)- c_TStress(3)     
          Stress_xy_N_Threads(C_Node,i_Thread) = Stress_xy_N_Threads(C_Node,i_Thread)- c_TStress(4)
          Stress_yz_N_Threads(C_Node,i_Thread) = Stress_yz_N_Threads(C_Node,i_Thread)- c_TStress(5)
          Stress_xz_N_Threads(C_Node,i_Thread) = Stress_xz_N_Threads(C_Node,i_Thread)- c_TStress(6)       
      endif
      ! Add initial stress field
      if(Key_InSitu_Strategy==2 .and. Yes_Add_Insitu)then
          Stress_xx_N_Threads(C_Node,i_Thread) = Stress_xx_N_Threads(C_Node,i_Thread) +Str_xx_InSitu(C_Node)
          Stress_yy_N_Threads(C_Node,i_Thread) = Stress_yy_N_Threads(C_Node,i_Thread) +Str_yy_InSitu(C_Node)
          Stress_zz_N_Threads(C_Node,i_Thread) = Stress_zz_N_Threads(C_Node,i_Thread) +Str_zz_InSitu(C_Node)
          Stress_xy_N_Threads(C_Node,i_Thread) = Stress_xy_N_Threads(C_Node,i_Thread) +Str_xy_InSitu(C_Node)
          Stress_yz_N_Threads(C_Node,i_Thread) = Stress_yz_N_Threads(C_Node,i_Thread) +Str_yz_InSitu(C_Node)
          Stress_xz_N_Threads(C_Node,i_Thread) = Stress_xz_N_Threads(C_Node,i_Thread) +Str_xz_InSitu(C_Node)     
      endif
      ! Subtract the stress corresponding to the initial strain field, theory: Equation 15.1.98 in
      ! 'Fundamentals of the Finite Element Method (5th Edition)', 2022-06-03
      if(Key_InSitu_Strategy==4)then
          ! Obtain the Gauss point stress based on the initial strain. Similar to the handling method of
          ! thermal stress.
          c_Stress = MATMUL(c_D,[ &
                 InSitu_Strain_Gaus_xx(i_E,1),&
                 InSitu_Strain_Gaus_yy(i_E,1),&
                 InSitu_Strain_Gaus_zz(i_E,1),&
                 InSitu_Strain_Gaus_xy(i_E,1),&
                 InSitu_Strain_Gaus_yz(i_E,1),&
                 InSitu_Strain_Gaus_xz(i_E,1)])   
          Stress_xx_N_Threads(C_Node,i_Thread) = Stress_xx_N_Threads(C_Node,i_Thread)- c_Stress(1)
          Stress_yy_N_Threads(C_Node,i_Thread) = Stress_yy_N_Threads(C_Node,i_Thread)- c_Stress(2)
          Stress_zz_N_Threads(C_Node,i_Thread) = Stress_zz_N_Threads(C_Node,i_Thread)- c_Stress(3)     
          Stress_xy_N_Threads(C_Node,i_Thread) = Stress_xy_N_Threads(C_Node,i_Thread)- c_Stress(4)
          Stress_yz_N_Threads(C_Node,i_Thread) = Stress_yz_N_Threads(C_Node,i_Thread)- c_Stress(5)
          Stress_xz_N_Threads(C_Node,i_Thread) = Stress_xz_N_Threads(C_Node,i_Thread)- c_Stress(6)       
      endif              
    end do
end do  
!$omp end do
!$omp end parallel   
  
! Summarize the calculation results of each thread.
DO i_Thread = 1,omp_get_max_threads()
    c_Count = c_Count + c_Count_Threads(:,i_Thread)
    Stress_xx_N = Stress_xx_N + Stress_xx_N_Threads(:,i_Thread)
    Stress_yy_N = Stress_yy_N + Stress_yy_N_Threads(:,i_Thread)
    Stress_zz_N = Stress_zz_N + Stress_zz_N_Threads(:,i_Thread)
    Stress_xy_N = Stress_xy_N + Stress_xy_N_Threads(:,i_Thread)
    Stress_yz_N = Stress_yz_N + Stress_yz_N_Threads(:,i_Thread)
    Stress_xz_N = Stress_xz_N + Stress_xz_N_Threads(:,i_Thread)
    Stress_vm_N = Stress_vm_N + Stress_vm_N_Threads(:,i_Thread)
ENDDO
  
!---------------------
! Node Stress Average
!---------------------
! OpenMP multi-core computing
!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_N)   &    
!$OMP            SCHEDULE(static)   
do i_N=1,num_Node 
    Stress_xx_N(i_N) = Stress_xx_N(i_N)/c_Count(i_N)
    Stress_yy_N(i_N) = Stress_yy_N(i_N)/c_Count(i_N)
    Stress_zz_N(i_N) = Stress_zz_N(i_N)/c_Count(i_N)
    Stress_xy_N(i_N) = Stress_xy_N(i_N)/c_Count(i_N)
    Stress_yz_N(i_N) = Stress_yz_N(i_N)/c_Count(i_N)
    Stress_xz_N(i_N) = Stress_xz_N(i_N)/c_Count(i_N)          
    call Tool_von_Mises_3D(Stress_xx_N(i_N),Stress_yy_N(i_N),Stress_zz_N(i_N), &
                         Stress_xy_N(i_N),Stress_yz_N(i_N),Stress_xz_N(i_N), &    
                         Stress_vm_N(i_N))
end do
!$omp end parallel do 

RETURN
END SUBROUTINE Get_Node_Stress_FEM_IN_OUT_3D
