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
 
SUBROUTINE Get_Node_Str_XFEM_IN_OUT_3D(Yes_Add_Insitu,isub,c_DISP,Key_to_cal,Key_CoorSys,Str_xx_N,Str_yy_N,Str_zz_N,&
                                     Str_xy_N,Str_yz_N,Str_xz_N,Str_vm_N)
! Calculate the stress or strain at the nodes. Key_to_cal = 1 to calculate stress, Key_to_cal = 2 to
! calculate strain.
! Key_CoorSys=1, Cartesian coordinate system, Key_CoorSys==2, Cylindrical coordinate system
! Store in variables: Str_xx_N, Str_yy_N, Str_xy_N, Str_vm_N, Str_zz_N, Str_yz_N, Str_xz_N

!-----------------------------          
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Dynamic
use Global_Material
use Global_Stress
use Global_Strain
use Global_Disp
use Global_Crack_Common
use Global_Crack_3D
use Global_POST
use Global_XFEM_Elements
use omp_lib
      
!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
integer, intent(in)::isub,Key_to_cal,Key_CoorSys
logical, intent(in)::Yes_Add_Insitu
real(kind=FT), intent(in)::c_DISP(num_Node*3)
real(kind=FT), intent(out)::Str_xx_N(num_Node),Str_yy_N(num_Node),Str_zz_N(num_Node), &
                          Str_xy_N(num_Node),Str_yz_N(num_Node),Str_xz_N(num_Node), &
                          Str_vm_N(num_Node)     
real(kind=FT) c_T_Alpha,c_TStr(6)     
integer i_E,mat_num
real(kind=FT) c_D(6,6),U(24)
real(kind=FT) c_v   
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8),T_Kesi(8),T_Yita(8),T_Zeta(8), c_kesi,c_yita,c_zeta,c_Str(6)
integer c_NN(8),i_N,C_Node,c_Count(num_Node) 
real(kind=FT) c_v_all(Num_Node)
logical Yes_Enrich_Ele,Yes_Enrich_Node
real(kind=FT) Rot_c_D_Comp(6,6),c_D_Comp(6,6),Volume_Ratio
real(kind=FT) T_Matrix(6,6),TT_Matrix(6,6)
real(kind=FT) Ele_Str_xx_N(Num_Elem,8),Ele_Str_yy_N(Num_Elem,8),Ele_Str_zz_N(Num_Elem,8), &
            Ele_Str_xy_N(Num_Elem,8),Ele_Str_yz_N(Num_Elem,8),Ele_Str_xz_N(Num_Elem,8), &
            Ele_Str_vm_N(Num_Elem,8)   
integer c_Elem
            
!2023-02-17.
if (Key_Simple_Post==1) return
      
!------------------
! Initialized to 0
!------------------
Ele_Str_xx_N(1:Num_Elem,1:8) = ZR
Ele_Str_yy_N(1:Num_Elem,1:8) = ZR
Ele_Str_zz_N(1:Num_Elem,1:8) = ZR
Ele_Str_xy_N(1:Num_Elem,1:8) = ZR
Ele_Str_yz_N(1:Num_Elem,1:8) = ZR
Ele_Str_xz_N(1:Num_Elem,1:8) = ZR
Ele_Str_vm_N(1:Num_Elem,1:8) = ZR

!------------------
! Inter-unit cycle
!------------------
T_Kesi = [-ONE,  ONE,  ONE, -ONE, -ONE,  ONE,  ONE, -ONE]
T_Yita = [-ONE, -ONE,  ONE,  ONE, -ONE, -ONE,  ONE,  ONE]
T_Zeta = [-ONE, -ONE, -ONE, -ONE,  ONE,  ONE,  ONE,  ONE]
      
      
!.....................................................
! OpenMP Multi-Core Computing: Option 2.
! Change to save the stress of 8 nodes for each element
! to avoid data conflicts.
! At the same time, it is divided into FEM elements and
! XFEM elements to maintain load balance.
!--------------
! FEM element.
!.....................................................
!2023-07-11.
if (Key_to_cal==1) then
    if(Key_CoorSys==1)then
          print *,'    Calculating stress of nodes of FEM elements...'
    elseif(Key_CoorSys==2)then
          print *,'    Calculating stress of nodes of FEM elements in cylindrical CS...'          
    endif
    elseif(Key_to_cal==2)then
    if(Key_CoorSys==1)then
          print *,'    Calculating strain of nodes of FEM elements...'
    elseif(Key_CoorSys==2)then
          print *,'    Calculating strain of nodes of FEM elements in cylindrical CS...'          
    endif          
endif
      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_Elem,mat_num,c_D,Volume_Ratio, &
!$OMP            c_T_Alpha,c_X_NODES,c_D_Comp,T_Matrix,TT_Matrix,Rot_c_D_Comp, &
!$OMP            c_Y_NODES,c_Z_NODES,c_NN,U,C_Node,c_kesi,c_yita,c_zeta,i_N, &
!$OMP            Yes_Enrich_Node,c_Str,c_TStr)  &
!$OMP            SCHEDULE(static)  
do i_E = 1,num_FEM_Elem
  c_Elem = FEM_Elem_List(i_E)
  mat_num = Elem_Mat(c_Elem)
  c_D     = D(Elem_Mat(c_Elem),1:6,1:6)
  c_T_Alpha = T_Alpha(Elem_Mat(c_Elem))
  !If the mat is composite material.
  if (Material_Type(mat_num)==5)then
      Volume_Ratio = Material_Para_Added(mat_num,10)
      c_D_Comp = D_Comp(mat_num,1:6,1:6)
      T_Matrix = Ele_ComMat_RotMatrix(c_Elem,1:6,1:6)
      TT_Matrix= TRANSPOSE(T_Matrix)
      Rot_c_D_Comp = MATMUL(TT_Matrix,c_D_Comp)
      Rot_c_D_Comp = MATMUL(Rot_c_D_Comp,T_Matrix)
      c_D =(ONE-Volume_Ratio)*c_D + Volume_Ratio*Rot_c_D_Comp
  endif                
  c_NN    = G_NN(1:8,c_Elem)
  c_X_NODES = G_X_NODES(1:8,c_Elem)
  c_Y_NODES = G_Y_NODES(1:8,c_Elem)    
  c_Z_NODES = G_Z_NODES(1:8,c_Elem)  
  U =[c_DISP(c_NN(1)*3-2),c_DISP(c_NN(1)*3-1),c_DISP(c_NN(1)*3),&
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
      !Check if the node is an enriched node.
      Yes_Enrich_Node = .False.
      if(sum(Enriched_Node_Type_3D(C_Node,1:num_Crack)) >=1)then
          Yes_Enrich_Node = .True.
      endif
      
      ! For degenerate elements, use the stress at the element center as the stress for the 8 nodes.
      ! Otherwise, a determinant of zero may occur, making it impossible to calculate Inverse_J.
      ! IMPROV2023061403.
      if(Yes_Degenarated_Elem(c_Elem) .eqv. .True.)then
          c_kesi = ZR
          c_yita = ZR
          c_zeta = ZR
      endif
      !If not an enriched node.
      if (Yes_Enrich_Node .eqv. .False.) then
        call Cal_Ele_Str_N8_3D(c_Elem,i_N,Key_to_cal,Key_CoorSys,c_X_NODES,c_Y_NODES,&
            c_Z_NODES,c_D,c_kesi,c_yita,c_zeta,U,c_Str)                  
      !If an enriched node.
      elseif  (Yes_Enrich_Node .eqv. .True.) then
          call Cal_Any_Point_Str_KesiYita_3D(c_Elem,i_N,Key_to_cal,Key_CoorSys, & 
                 c_kesi,c_yita,c_zeta,1,DISP,  &
                   c_Str(1),c_Str(2),c_Str(3), & 
                   c_Str(4),c_Str(5),c_Str(6))
      endif
      Ele_Str_xx_N(c_Elem,i_N) = c_Str(1)
      Ele_Str_yy_N(c_Elem,i_N) = c_Str(2)
      Ele_Str_zz_N(c_Elem,i_N) = c_Str(3)
      Ele_Str_xy_N(c_Elem,i_N) = c_Str(4)
      Ele_Str_yz_N(c_Elem,i_N) = c_Str(5)
      Ele_Str_xz_N(c_Elem,i_N) = c_Str(6)
      ! Subtract thermal expansion stress, Theory: Equation 15.1.98 from 'Fundamentals of Finite Element
      ! Method (5th Edition)', 2019-09-24
      if(Key_Thermal_Stress==1 .and. Key_to_cal==1)then
          c_TStr = c_T_Alpha*Elem_T_for_Stress(c_Elem)*MATMUL(c_D,[ONE,ONE,ONE,ZR,ZR,ZR]) 
          Ele_Str_xx_N(c_Elem,i_N)=Ele_Str_xx_N(c_Elem,i_N)-c_TStr(1)
          Ele_Str_yy_N(c_Elem,i_N)=Ele_Str_yy_N(c_Elem,i_N)-c_TStr(2)
          Ele_Str_zz_N(c_Elem,i_N)=Ele_Str_zz_N(c_Elem,i_N)-c_TStr(3)
          Ele_Str_xy_N(c_Elem,i_N)=Ele_Str_xy_N(c_Elem,i_N)-c_TStr(4)
          Ele_Str_yz_N(c_Elem,i_N)=Ele_Str_yz_N(c_Elem,i_N)-c_TStr(5)
          Ele_Str_xz_N(c_Elem,i_N)=Ele_Str_xz_N(c_Elem,i_N)-c_TStr(6)
      endif
      ! Obtain Gauss point stress based on the initial strain. Similar to the treatment method for thermal
      ! stress. 2022-06-03.
      if(Key_InSitu_Strategy==4 )then
          c_TStr = MATMUL(c_D,[ &
                 InSitu_Strain_Gaus_xx(c_Elem,1), &
                 InSitu_Strain_Gaus_yy(c_Elem,1), &
                 InSitu_Strain_Gaus_zz(c_Elem,1), &
                 InSitu_Strain_Gaus_xy(c_Elem,1), &
                 InSitu_Strain_Gaus_yz(c_Elem,1), &
                 InSitu_Strain_Gaus_xz(c_Elem,1)])          
          Ele_Str_xx_N(c_Elem,i_N)=Ele_Str_xx_N(c_Elem,i_N)-c_TStr(1)
          Ele_Str_yy_N(c_Elem,i_N)=Ele_Str_yy_N(c_Elem,i_N)-c_TStr(2)
          Ele_Str_zz_N(c_Elem,i_N)=Ele_Str_zz_N(c_Elem,i_N)-c_TStr(3)
          Ele_Str_xy_N(c_Elem,i_N)=Ele_Str_xy_N(c_Elem,i_N)-c_TStr(4)
          Ele_Str_yz_N(c_Elem,i_N)=Ele_Str_yz_N(c_Elem,i_N)-c_TStr(5)
          Ele_Str_xz_N(c_Elem,i_N)=Ele_Str_xz_N(c_Elem,i_N)-c_TStr(6)
      endif              
    end do
end do  
!$omp end parallel do          
      
if(Enrich_Freedom==0) goto 100

!........................................................
! OpenMP Multi-Core Computing: Option 2.
! Change to save the stress of 8 nodes for each element 
! to avoid data conflicts.
! At the same time, it is divided into FEM elements and
! XFEM elements to maintain load balance.
!----------------
!
! XFEM element.
!
!........................................................
!2023-07-11.
if (Key_to_cal==1) then
    if(Key_CoorSys==1)then
          print *,'    Calculating stress of nodes of XFEM elements...'
    elseif(Key_CoorSys==2)then
          print *,'    Calculating stress of nodes of XFEM elements in cylindrical CS...'          
    endif
    elseif(Key_to_cal==2)then
    if(Key_CoorSys==1)then
          print *,'    Calculating strain of nodes of XFEM elements...'
    elseif(Key_CoorSys==2)then
          print *,'    Calculating strain of nodes of XFEM elements in cylindrical CS...'          
    endif          
endif
      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_Elem,mat_num,c_D,Volume_Ratio, &
!$OMP            c_T_Alpha,c_X_NODES,c_D_Comp,T_Matrix,TT_Matrix,Rot_c_D_Comp, &
!$OMP            c_Y_NODES,c_Z_NODES,c_NN,U,C_Node,c_kesi,c_yita,c_zeta,i_N, &
!$OMP            Yes_Enrich_Node,c_Str,c_TStr)  &
!$OMP            SCHEDULE(static)  
do i_E = 1,num_XFEM_Elem
  c_Elem = XFEM_Elem_List(i_E)
  mat_num = Elem_Mat(c_Elem)
  c_D     = D(Elem_Mat(c_Elem),1:6,1:6)
  c_T_Alpha = T_Alpha(Elem_Mat(c_Elem))
  !If the mat is composite material.
  if (Material_Type(mat_num)==5)then
      Volume_Ratio = Material_Para_Added(mat_num,10)
      c_D_Comp = D_Comp(mat_num,1:6,1:6)
      T_Matrix = Ele_ComMat_RotMatrix(c_Elem,1:6,1:6)
      TT_Matrix= TRANSPOSE(T_Matrix)
      Rot_c_D_Comp = MATMUL(TT_Matrix,c_D_Comp)
      Rot_c_D_Comp = MATMUL(Rot_c_D_Comp,T_Matrix)
      c_D =(ONE-Volume_Ratio)*c_D + Volume_Ratio*Rot_c_D_Comp
  endif                
  c_NN    = G_NN(1:8,c_Elem)
  c_X_NODES = G_X_NODES(1:8,c_Elem)
  c_Y_NODES = G_Y_NODES(1:8,c_Elem)    
  c_Z_NODES = G_Z_NODES(1:8,c_Elem)  
  U =[c_DISP(c_NN(1)*3-2),c_DISP(c_NN(1)*3-1),c_DISP(c_NN(1)*3),&
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
      !Check if the node is an enriched node.
      Yes_Enrich_Node = .False.
      if(sum(Enriched_Node_Type_3D(C_Node,1:num_Crack)) >=1)then
          Yes_Enrich_Node = .True.
      endif
      
      ! For degraded elements, the stress at the element center is used as the stress for the 8 nodes.
      ! Otherwise, a determinant of zero may occur, making it impossible to calculate Inverse_J.
      ! IMPROV2023061403.
      if(Yes_Degenarated_Elem(c_Elem) .eqv. .True.)then
          c_kesi = ZR
          c_yita = ZR
          c_zeta = ZR
      endif
      !If not an enriched node.
      if (Yes_Enrich_Node .eqv. .False.) then
        call Cal_Ele_Str_N8_3D(c_Elem,i_N,Key_to_cal,Key_CoorSys,c_X_NODES,c_Y_NODES,&
            c_Z_NODES,c_D,c_kesi,c_yita,c_zeta,U,c_Str)                  
      !If an enriched node.
      elseif  (Yes_Enrich_Node .eqv. .True.) then
          call Cal_Any_Point_Str_KesiYita_3D(c_Elem,i_N,Key_to_cal,Key_CoorSys, & 
                   c_kesi,c_yita,c_zeta,1,DISP,  &
                   c_Str(1),c_Str(2),c_Str(3), & 
                   c_Str(4),c_Str(5),c_Str(6))
      endif
      Ele_Str_xx_N(c_Elem,i_N) = c_Str(1)
      Ele_Str_yy_N(c_Elem,i_N) = c_Str(2)
      Ele_Str_zz_N(c_Elem,i_N) = c_Str(3)
      Ele_Str_xy_N(c_Elem,i_N) = c_Str(4)
      Ele_Str_yz_N(c_Elem,i_N) = c_Str(5)
      Ele_Str_xz_N(c_Elem,i_N) = c_Str(6)
      
      ! Subtract thermal expansion stress, Theory: Equation 15.1.98 from 'Fundamentals of Finite Element
      ! Method (5th Edition)', 2019-09-24
      if(Key_Thermal_Stress==1 .and. Key_to_cal==1)then
          c_TStr = c_T_Alpha*Elem_T_for_Stress(c_Elem)*MATMUL(c_D,[ONE,ONE,ONE,ZR,ZR,ZR]) 
          Ele_Str_xx_N(c_Elem,i_N)=Ele_Str_xx_N(c_Elem,i_N)-c_TStr(1)
          Ele_Str_yy_N(c_Elem,i_N)=Ele_Str_yy_N(c_Elem,i_N)-c_TStr(2)
          Ele_Str_zz_N(c_Elem,i_N)=Ele_Str_zz_N(c_Elem,i_N)-c_TStr(3)
          Ele_Str_xy_N(c_Elem,i_N)=Ele_Str_xy_N(c_Elem,i_N)-c_TStr(4)
          Ele_Str_yz_N(c_Elem,i_N)=Ele_Str_yz_N(c_Elem,i_N)-c_TStr(5)
          Ele_Str_xz_N(c_Elem,i_N)=Ele_Str_xz_N(c_Elem,i_N)-c_TStr(6)
      endif
      
      ! Obtain Gauss point stress based on the initial strain. Similar to the treatment method for thermal
      ! stress. 2022-06-03.
      if(Key_InSitu_Strategy==4 )then
          c_TStr = MATMUL(c_D,[ &
                 InSitu_Strain_Gaus_xx(c_Elem,1), &
                 InSitu_Strain_Gaus_yy(c_Elem,1), &
                 InSitu_Strain_Gaus_zz(c_Elem,1), &
                 InSitu_Strain_Gaus_xy(c_Elem,1), &
                 InSitu_Strain_Gaus_yz(c_Elem,1), &
                 InSitu_Strain_Gaus_xz(c_Elem,1)])          
          Ele_Str_xx_N(c_Elem,i_N)=Ele_Str_xx_N(c_Elem,i_N)-c_TStr(1)
          Ele_Str_yy_N(c_Elem,i_N)=Ele_Str_yy_N(c_Elem,i_N)-c_TStr(2)
          Ele_Str_zz_N(c_Elem,i_N)=Ele_Str_zz_N(c_Elem,i_N)-c_TStr(3)
          Ele_Str_xy_N(c_Elem,i_N)=Ele_Str_xy_N(c_Elem,i_N)-c_TStr(4)
          Ele_Str_yz_N(c_Elem,i_N)=Ele_Str_yz_N(c_Elem,i_N)-c_TStr(5)
          Ele_Str_xz_N(c_Elem,i_N)=Ele_Str_xz_N(c_Elem,i_N)-c_TStr(6)
      endif              
    end do
end do  
!$omp end parallel do           

100 continue
      
!------------------------------------------------------------------------------------------------
! Convert to nodal stress. 2023-07-11.
! Not suitable for parallel computing. Besides, it takes very little time, so parallelization is
! unnecessary.
!------------------------------------------------------------------------------------------------
Str_xx_N(1:num_Node) = ZR
Str_yy_N(1:num_Node) = ZR
Str_zz_N(1:num_Node) = ZR
Str_xy_N(1:num_Node) = ZR
Str_yz_N(1:num_Node) = ZR
Str_xz_N(1:num_Node) = ZR
Str_vm_N(1:num_Node) = ZR      
c_Count(1:num_Node)  = 0
do i_E = 1,Num_Elem
    c_NN   = G_NN(1:8,i_E)
    do i_N = 1,8
        C_Node = c_NN(i_N)   
        c_Count(C_Node) = c_Count(C_Node) + 1  
        Str_xx_N(C_Node) = Str_xx_N(C_Node) + Ele_Str_xx_N(i_E,i_N)
        Str_yy_N(C_Node) = Str_yy_N(C_Node) + Ele_Str_yy_N(i_E,i_N)
        Str_zz_N(C_Node) = Str_zz_N(C_Node) + Ele_Str_zz_N(i_E,i_N)
        Str_xy_N(C_Node) = Str_xy_N(C_Node) + Ele_Str_xy_N(i_E,i_N)
        Str_yz_N(C_Node) = Str_yz_N(C_Node) + Ele_Str_yz_N(i_E,i_N)
        Str_xz_N(C_Node) = Str_xz_N(C_Node) + Ele_Str_xz_N(i_E,i_N)
    enddo
enddo

!-----------------------
! Average nodal stress.
!-----------------------
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_N) SCHEDULE(STATIC)     
do i_N=1,num_Node 
    Str_xx_N(i_N) = Str_xx_N(i_N)/c_Count(i_N)
    Str_yy_N(i_N) = Str_yy_N(i_N)/c_Count(i_N)
    Str_zz_N(i_N) = Str_zz_N(i_N)/c_Count(i_N)
    Str_xy_N(i_N) = Str_xy_N(i_N)/c_Count(i_N)
    Str_yz_N(i_N) = Str_yz_N(i_N)/c_Count(i_N)
    Str_xz_N(i_N) = Str_xz_N(i_N)/c_Count(i_N)          
    call Tool_von_Mises_3D(Str_xx_N(i_N),Str_yy_N(i_N),Str_zz_N(i_N),Str_xy_N(i_N), &
                         Str_yz_N(i_N),Str_xz_N(i_N),Str_vm_N(i_N))
end do
!$omp end parallel do      
      
RETURN
END SUBROUTINE Get_Node_Str_XFEM_IN_OUT_3D
