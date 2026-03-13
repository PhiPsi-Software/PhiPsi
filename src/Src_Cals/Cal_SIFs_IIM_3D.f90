
subroutine Cal_SIFs_IIM_3D(iter,c_DISP,Tab_Num)

use Global_Float_Type
use Global_Common
use Global_Model
use Global_Material
use Global_DISP
use Global_Crack_Common
use Global_Crack_3D
use Global_Ragged_Array_Real_Classs
use Global_Cal_Ele_Num_by_Coors_3D
use Global_Elem_Area_Vol
use Global_POST
use Global_Filename
use Global_Crack
use Global_Stress
use module_Cal_3D_SIFs_IIM_Auxiliary_Fields_Part1
use module_Cal_3D_SIFs_IIM_Auxiliary_Fields_Part2
use module_Cal_3D_SIFs_IIM_q_Function
use module_Cal_3D_SIFs_IIM_Integration_Domain
use module_Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_3_Points

implicit none
integer, intent(in) :: iter, Tab_Num
real(kind=FT), intent(in) :: c_DISP(Total_FD)
real(kind=FT) :: SMALL_R
integer :: i_Crack, i_Vertex, i_Elem, i_GP, i_Mode
integer :: i, j, k
integer :: num_Front_Vertices
integer :: Vertex_Node_ID
real(kind=FT) :: Tip_Coords(3)
real(kind=FT) :: e1_growth(3), e2_normal(3), e3_tangent(3)
real(kind=FT) :: Rotation_Global_to_Local(3,3)
real(kind=FT) :: R_outer
real(kind=FT) :: R_inner
real(kind=FT) :: L_front
integer :: Num_Domain_Elems
integer, allocatable :: Domain_Elem_List(:)
real(kind=FT), allocatable :: q_nodal(:,:)
integer :: Elem_at_Tip, Material_ID
real(kind=FT) :: Young_Modulus, Poisson_Ratio, Shear_Modulus
real(kind=FT) :: E_prime, Kappa_3D
integer :: Current_Elem
integer :: Node_IDs(8)
real(kind=FT) :: Node_X(8), Node_Y(8), Node_Z(8)
real(kind=FT) :: Constitutive_Matrix(6,6)
integer :: nGP
real(kind=FT), allocatable :: GP_xi(:), GP_eta(:), GP_zeta(:), GP_weights(:)
real(kind=FT) :: Shape_Functions(8)
real(kind=FT) :: dN_dxi(8), dN_deta(8), dN_dzeta(8)
real(kind=FT) :: Jacobian_Matrix(3,3), Inv_Jacobian(3,3), Det_Jacobian
real(kind=FT) :: dN_dx(8), dN_dy(8), dN_dz(8)
real(kind=FT) :: GP_coords_global(3)
real(kind=FT) :: GP_coords_local(3)
real(kind=FT) :: q_value_at_GP
real(kind=FT) :: grad_q_global(3)
real(kind=FT) :: grad_q_local(3)
real(kind=FT) :: Strain_Voigt(6)
real(kind=FT) :: Stress_Voigt(6)
real(kind=FT) :: Temp_Stress_Voigt(6) 
real(kind=FT) :: Stress_Tensor(3,3)
real(kind=FT) :: Strain_Tensor(3,3)
real(kind=FT) :: Grad_u_Tensor(3,3)
real(kind=FT) :: Stress_Local(3,3)
real(kind=FT) :: Strain_Local(3,3)
real(kind=FT) :: Grad_u_Local(3,3)
real(kind=FT) :: r_polar, theta_polar
real(kind=FT) :: Aux_Stress_Local(3,3)
real(kind=FT) :: Aux_Strain_Local(3,3)
real(kind=FT) :: Aux_Grad_u_Local(3,3)
real(kind=FT) :: M_integral(3)
real(kind=FT), allocatable :: M_integral_Elements(:,:)
real(kind=FT) :: Integrand_value
real(kind=FT) :: Term_W_delta1j
real(kind=FT) :: Term_sigma_gradu_aux
real(kind=FT) :: Term_sigma_aux_gradu
real(kind=FT) :: Element_DOF_values(24)
integer :: DOF_indices(24)
integer :: Node_i, dof_index
character(5) :: iter_string
real(kind=FT) :: p_fluid
real(kind=FT) :: Stress_from_pressure(3,3)
integer :: Element_Search_Cache
logical :: Is_Boundary_Crack
integer :: Sigma_for_Denoise

character(5) temp,temp1,temp2
character(200) c_File_name_test,c_File_name_test_1,c_File_name_test_2,c_File_name_test_3
character(200) c_File_name_test_4,c_File_name_test_5,c_File_name_test_6
integer i_C,num_Cr_Edges,i_V
integer ii,jj,kk,c_elem

real(kind=FT) :: Wint
real(kind=FT) :: A(3)
real(kind=FT) :: dV
real(kind=FT) :: Term_W,L_eff

integer:: Location_ESM(MDOF_3D)    
integer c_POS_3D_c_Ele(8)  
integer::Location_ESM_C_Crack(MDOF_3D),Location_ESM_C_Cr_NoFEM(MDOF_3D)
integer num_Loc_ESM_C_Cr_NoFEM
integer num_Loc_ESM
integer num_Loc_ESM_C_Crack
real(kind=FT) U_elem(MDOF_3D)
real(kind=FT) B(6,MDOF_3D),tem_B(6,MDOF_3D)
integer num_B,num_tem_B
real(kind=FT) q_line_integral
logical :: Has_Fluid_Pressure
real(kind=FT) :: s_trac
real(kind=FT) :: t_global(3), t_local(3)

integer :: nGP_face, i_GP_face
real(kind=FT), allocatable :: GP_face_global(:,:), GP_face_w(:), GP_face_n_global(:,:)
real(kind=FT) :: q_face
real(kind=FT) :: Integrand_face
real(kind=FT) :: nrm_n, n_unit(3)
real(kind=FT) :: fac_two_faces
real(kind=FT) :: q_integral_denom
real(kind=FT) :: q_volume_integral
real(kind=FT) :: elem_vol_temp
real(kind=FT) :: J_mat(3,3), detJ
real(kind=FT) :: X(8), Y(8), Z(8),q_value
integer iii
integer c_OUT_Elem_Old(5000)
real(kind=FT) Final_Point(5000,3) 
integer Final_Point_Prop(5000) 
real(kind=FT) c_KI_3D,c_KII_3D,c_KIII_3D
real(kind=FT),ALLOCATABLE::Check_Theta(:)
real(kind=FT),ALLOCATABLE::KI_eq(:),Theta_All(:)
real(kind=FT),ALLOCATABLE::Spe_Pri_Stress(:)      
real(kind=FT) tem_part1,tem_part2,Stress_Theta,Tao_Theta
real(kind=FT) tem_root_2pir
real(kind=FT) c_Theta 
integer i_Theta_min
real(kind=FT) Schollm_Max_Theta_in_pi
real(kind=FT) c_x_old,c_y_old,c_z_old
real(kind=FT),ALLOCATABLE::Sai_Cr(:),Ld_Cr(:),Ks_Cr(:),D_Cr(:) 
real(kind=FT) r_for_c_Theta_finding
integer num_of_check,i_Check_Theta 
integer Num_CrMesh_Outlines,i_Out_Node,c_Mesh_Node
integer c_OUT_Elem,num_Ver_in_Model
integer n_Sigma,mat_num
logical logical_close
integer in_Ele_Num
integer Ele_Num_Cache 
real(kind=FT) cc_Part1,cc_Part2,cc_Part3
real(kind=FT) c_x,c_y,c_z
integer c_Mesh_Node_Next
real(kind=FT) Tool_Function_2Point_Dis_3D
real(kind=FT) c_KIc 
real(kind=FT) Normal_Stress
real(kind=FT) Aux_dSig11_dX1, Aux_dSig12_dX2, Aux_dSig13_dX3
real(kind=FT) Aux_dSig12_dX1, Aux_dSig22_dX2, Aux_dSig23_dX3
real(kind=FT) Aux_dSig13_dX1, Aux_dSig23_dX2, Aux_dSig33_dX3
real(kind=FT) Aux_dEps11_dX1, Aux_dEps12_dX1, Aux_dEps13_dX1
real(kind=FT) Aux_dEps22_dX1, Aux_dEps23_dX1, Aux_dEps33_dX1
real(kind=FT) Aux_ddu1_dX1dX1, Aux_ddu2_dX1dX1, Aux_ddu3_dX1dX1
real(kind=FT) Aux_ddu1_dX2dX1, Aux_ddu2_dX2dX1, Aux_ddu3_dX2dX1 
real(kind=FT) Aux_ddu1_dX3dX1, Aux_ddu2_dX3dX1, Aux_ddu3_dX3dX1
real(kind=FT) Integrand_value1,Integrand_value2

1001 FORMAT(5X,'-- Calculating Vertex ',I4,' /',I4,' of crack',I4,'...')  
1002 FORMAT(8X,'-- Calculating Vertex ',I4,' /',I4,' of crack',I4,'...')  
1101 FORMAT(5X,'-- Calculating SIFs using IIM for crack',I4,' /',I4,'...')  
1102 FORMAT(8X,'-- Calculating SIFs using IIM for crack',I4,' /',I4,'...')  

2001 FORMAT(5X,'   SIFs of Vertex ',I4,' of crack ',I4, ':',F11.5,' / ', F11.5,' / ', F11.5, ' MPa.m^1/2.')  


R_inner = 0.5D0 * Ave_Elem_L_Enrich  
R_outer = 3.0D0 * Ave_Elem_L_Enrich                  
L_front = Ave_Elem_L_Enrich
SMALL_R = 1.0D-12 



nGP = Key_3D_Sifs_IIM_num_Gauss_Points
allocate(GP_xi(nGP))
allocate(GP_eta(nGP))
allocate(GP_zeta(nGP))
allocate(GP_weights(nGP))

call Cal_Gauss_Points_3D_8nodes(nGP,GP_xi, GP_eta, GP_zeta, GP_weights)

do i_Crack = 1, num_Crack
    if(Tab_Num==5)then
        write(*,1101) i_Crack, num_Crack
    elseif(Tab_Num==8)then
        write(*,1102) i_Crack, num_Crack
    endif
    Element_Search_Cache = 1
    num_Front_Vertices = Crack3D_Meshed_Outline_num(i_Crack)
    
    Has_Fluid_Pressure = .false.
    if (Crack_Type_Status_3D(i_Crack,1) == 1) then
        Has_Fluid_Pressure = .true.
    endif
    if(Crack_Type_Status_3D(i_Crack,1) == 1) then
        p_fluid = Crack_Pressure(i_Crack)  
    else
        p_fluid = ZR
    endif
    s_trac = ONE

    if(allocated(KI_3D(i_Crack)%row))    deallocate(KI_3D(i_Crack)%row)
    if(allocated(KII_3D(i_Crack)%row))   deallocate(KII_3D(i_Crack)%row)
    if(allocated(KIII_3D(i_Crack)%row))  deallocate(KIII_3D(i_Crack)%row)
    if(allocated(KI_eq_3D(i_Crack)%row)) deallocate(KI_eq_3D(i_Crack)%row)

    allocate(KI_3D(i_Crack)%row(num_Front_Vertices))
    allocate(KII_3D(i_Crack)%row(num_Front_Vertices))
    allocate(KIII_3D(i_Crack)%row(num_Front_Vertices))
    allocate(KI_eq_3D(i_Crack)%row(num_Front_Vertices))

    KI_3D(i_Crack)%row(:) = ZR
    KII_3D(i_Crack)%row(:) = ZR
    KIII_3D(i_Crack)%row(:) = ZR
    KI_eq_3D(i_Crack)%row(:) = ZR

    do i_Vertex = 1, num_Front_Vertices
    

        Vertex_Node_ID = Crack3D_Meshed_Outline(i_Crack)%row(i_Vertex, 1)
        Tip_Coords(1:3) = Crack3D_Meshed_Node(i_Crack)%row(Vertex_Node_ID, 1:3)

        e1_growth(1:3)  =  Crack3D_Meshed_Vertex_x_Vector(i_Crack)%row(i_Vertex, 1:3)
        e2_normal(1:3)  =  Crack3D_Meshed_Vertex_y_Vector(i_Crack)%row(i_Vertex, 1:3)
        e3_tangent(1:3) =  Crack3D_Meshed_Vertex_z_Vector(i_Crack)%row(i_Vertex, 1:3)
        
        call Vector_Normalize(3,e1_growth(1:3))   
        call Vector_Normalize(3,e2_normal(1:3))   
        call Vector_Normalize(3,e3_tangent(1:3))   

        Rotation_Global_to_Local(1, 1:3) = e1_growth
        Rotation_Global_to_Local(2, 1:3) = e2_normal
        Rotation_Global_to_Local(3, 1:3) = e3_tangent
        
        if(Key_Save_Nothing /= 1) then 
            if (Key_Save_3D_SIFs_IIM_Integral_cylinders == 1) then 
                write(temp,'(I5)') iter
                write(temp1,'(I5)') i_Crack
                write(temp2,'(I5)') i_Vertex
                c_File_name_test=trim(Full_Pathname)//'.integral_cylinder_'//ADJUSTL(temp)//'_'//ADJUSTL(temp1)//'_'//ADJUSTL(temp2)     
                call Tool_chrpak_s_blank_delete(c_File_name_test)
                open(101,file=c_File_name_test,status='unknown')  
                write(101, '(15E20.12)') Tip_Coords,e1_growth,e2_normal,e3_tangent,L_front,R_inner,R_outer
                close(101)
            endif
        endif    
        

        call Cal_Ele_Num_by_Coors_3D(Tip_Coords(1), Tip_Coords(2), Tip_Coords(3), &
                                     Element_Search_Cache, Elem_at_Tip)

        if(Elem_at_Tip == 0) cycle

        Material_ID = Elem_Mat(Elem_at_Tip)
        if(Key_XA /= 2) then
            Young_Modulus = Material_Para(Material_ID, 1)
            Poisson_Ratio = Material_Para(Material_ID, 2)
        else
            Young_Modulus = Elem_E_XA(Elem_at_Tip)
            Poisson_Ratio = Elem_Mu_XA(Elem_at_Tip)
        endif

        Shear_Modulus = Young_Modulus / (TWO * (ONE + Poisson_Ratio))
        E_prime = Young_Modulus / (ONE - Poisson_Ratio**2)
        Kappa_3D = THR - FOU * Poisson_Ratio

        call Cal_3D_SIFs_IIM_Integration_Domain(Tip_Coords, e1_growth, e2_normal, e3_tangent, &
                                      R_inner, R_outer, L_front, &
                                      Num_Domain_Elems, Domain_Elem_List, q_nodal)
        

        allocate(M_integral_Elements(Num_Domain_Elems,3))
        M_integral_Elements= ZR

!$OMP PARALLEL DO default(none) SHARED(Num_Domain_Elems,Domain_Elem_List,G_NN,gp_zeta,i_Vertex,&
!$OMP      Elem_Mat,D,c_POS_3D,SMALL_R,s_trac,i_Crack,&
!$OMP      c_DISP,Tip_Coords,e1_growth, e2_normal,e3_tangent,R_inner,R_outer,L_front,&
!$OMP      Rotation_Global_to_Local,Shear_Modulus, Poisson_Ratio, Kappa_3D,elegaus_yes_fem_asemd, &
!$OMP      M_integral_Elements,Has_Fluid_Pressure,p_fluid,TWO,gp_weights,num_gau_points_3d,&
!$OMP      q_nodal,ZR,gp_xi,gp_eta,ngp,num_crack,mdof_3d,G_X_NODES,G_Y_NODES,g_z_nodes) &
!$OMP      PRIVATE(i_Elem,Current_Elem,Node_IDs,Node_X,Node_Y,Node_Z,&
!$OMP      Constitutive_Matrix,Node_i,DOF_indices,Element_DOF_values,&
!$OMP      Location_ESM,num_Loc_ESM,c_POS_3D_c_Ele,&
!$OMP      Location_ESM_C_Crack,num_Loc_ESM_C_Crack,Location_ESM_C_Cr_NoFEM,num_Loc_ESM_C_Cr_NoFEM,&
!$OMP      U_elem,Shape_Functions, dN_dxi, dN_deta, dN_dzeta,&
!$OMP      Jacobian_Matrix, Det_Jacobian,Inv_Jacobian,dN_dx, dN_dy, dN_dz,GP_coords_global,&
!$OMP      q_value_at_GP, grad_q_global,GP_coords_local,grad_q_local,& 
!$OMP      Strain_Voigt,&
!$OMP      B, num_B,i_C,tem_B,num_tem_B,Stress_Voigt,Grad_u_Tensor,& 
!$OMP      Stress_Local,Strain_Local,Grad_u_Local,r_polar,theta_polar,i_Mode,&
!$OMP      Aux_Stress_Local,Aux_Strain_Local,Aux_Grad_u_Local,Stress_from_pressure,& 
!$OMP      dV,Wint,jj,i,j,A,Integrand_value,i_GP,Material_ID,Stress_Tensor,Strain_Tensor,& 
!$OMP      fac_two_faces,nGP_face,GP_face_global,GP_face_w,GP_face_n_global,& 
!$OMP      i_GP_face,n_unit,t_global,t_local,q_face,Integrand_face,Normal_Stress,&
!$OMP      Aux_dSig11_dX1, Aux_dSig12_dX2, Aux_dSig13_dX3,    &
!$OMP      Aux_dSig12_dX1, Aux_dSig22_dX2, Aux_dSig23_dX3,    &
!$OMP      Aux_dSig13_dX1, Aux_dSig23_dX2, Aux_dSig33_dX3,    &
!$OMP      Aux_dEps11_dX1, Aux_dEps12_dX1, Aux_dEps13_dX1,    &
!$OMP      Aux_dEps22_dX1, Aux_dEps23_dX1, Aux_dEps33_dX1,    &
!$OMP      Aux_ddu1_dX1dX1, Aux_ddu2_dX1dX1, Aux_ddu3_dX1dX1, &
!$OMP      Aux_ddu1_dX2dX1, Aux_ddu2_dX2dX1, Aux_ddu3_dX2dX1, & 
!$OMP      Aux_ddu1_dX3dX1, Aux_ddu2_dX3dX1, Aux_ddu3_dX3dX1,Integrand_value1,Integrand_value2)
        do i_Elem = 1, Num_Domain_Elems
           
            Current_Elem = Domain_Elem_List(i_Elem)
            if(Current_Elem <= 0) cycle

            Node_IDs(1:8) = G_NN(1:8, Current_Elem)
            Node_X(1:8) = G_X_NODES(1:8, Current_Elem)
            Node_Y(1:8) = G_Y_NODES(1:8, Current_Elem)
            Node_Z(1:8) = G_Z_NODES(1:8, Current_Elem)
            

            Material_ID = Elem_Mat(Current_Elem)
            Constitutive_Matrix = D(Material_ID, 1:6, 1:6)
            
            

            Location_ESM(1:MDOF_3D)   = 0
            num_Loc_ESM               = 0
            if(num_Crack/=0)then
                do i_C =1,num_Crack
                  c_POS_3D_c_Ele(1:8) = c_POS_3D(Node_IDs,i_C) 
                  call Location_Element_Stiff_Matrix_3D(Current_Elem,i_C, &
                                              c_POS_3D_c_Ele(1:8), &
                                              Location_ESM_C_Crack, &
                                              num_Loc_ESM_C_Crack, &
                                              Location_ESM_C_Cr_NoFEM, &
                                              num_Loc_ESM_C_Cr_NoFEM)
                  Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
                  num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
                end do
            endif            

            U_elem(1:num_Loc_ESM) = c_DISP(Location_ESM(1:num_Loc_ESM))
                

            do i_GP = 1, nGP
                call Cal_Shape_Functions_Hex8(GP_xi(i_GP), GP_eta(i_GP), GP_zeta(i_GP), &
                                                 Shape_Functions, dN_dxi, dN_deta, dN_dzeta)
                call Cal_Jacobian_3D_Hex8(Node_X, Node_Y, Node_Z, dN_dxi, dN_deta, dN_dzeta, &
                                        Jacobian_Matrix, Det_Jacobian)
                if(Det_Jacobian <= ZR) cycle
                call Matrix_Inverse_3x3(Jacobian_Matrix, Inv_Jacobian)  
                call Cal_Transform_Derivatives_Hex8(dN_dxi, dN_deta, dN_dzeta, Inv_Jacobian, &
                                           dN_dx, dN_dy, dN_dz)
                GP_coords_global(1) = sum(Shape_Functions * Node_X)
                GP_coords_global(2) = sum(Shape_Functions * Node_Y)
                GP_coords_global(3) = sum(Shape_Functions * Node_Z)
                call Cal_3D_SIFs_IIM_q_Function(Current_Elem, q_nodal, Shape_Functions, &
                                           dN_dx, dN_dy, dN_dz, &
                                           GP_coords_global, Tip_Coords, &
                                           e1_growth, e2_normal, e3_tangent, &
                                           R_inner, R_outer, L_front, &
                                           q_value_at_GP, grad_q_global)
                                           
                

                if(abs(q_value_at_GP) < 1.0D-15) cycle

                GP_coords_local = matmul(Rotation_Global_to_Local, GP_coords_global - Tip_Coords)
                grad_q_local = matmul(Rotation_Global_to_Local, grad_q_global)
                


                
                
                
                
                
                EleGaus_yes_FEM_asemd(Current_Elem,1:Num_Gau_Points_3D)= .False. 
                B(1:6,1:MDOF_3D) = ZR
                num_B = 0 
                if(num_Crack/=0)then
                    do i_C =1,num_Crack 
                        call Cal_B_Matrix_Crack_3D(GP_xi(i_GP), GP_eta(i_GP), GP_zeta(i_GP), &
                                    i_C,Current_Elem,i_GP, &
                                    Node_IDs,Node_X,Node_Y,Node_Z, &
                                    tem_B,num_tem_B)     
                        B(1:6,num_B+1:num_B+num_tem_B)=tem_B(1:6,1:num_tem_B)
                        num_B = num_B + num_tem_B     
                    end do
                endif
                
                
                Strain_Voigt = MATMUL(B(1:6,1:num_Loc_ESM),U_elem(1:num_Loc_ESM))
                
                Stress_Voigt = matmul(Constitutive_Matrix, Strain_Voigt)
                
                
                call Tool_Voigt_to_Tensor_Symmetric(Stress_Voigt, Stress_Tensor)
                call Tool_Voigt_to_Tensor_Symmetric(Strain_Voigt, Strain_Tensor)
                
                                
                Grad_u_Tensor(1,1) = DOT_PRODUCT(B(1, 1:num_B:3), U_elem(1:num_B:3))
                Grad_u_Tensor(1,2) = DOT_PRODUCT(B(2, 2:num_B:3), U_elem(1:num_B:3))
                Grad_u_Tensor(1,3) = DOT_PRODUCT(B(3, 3:num_B:3), U_elem(1:num_B:3))
                Grad_u_Tensor(2,1) = DOT_PRODUCT(B(1, 1:num_B:3), U_elem(2:num_B:3))
                Grad_u_Tensor(2,2) = DOT_PRODUCT(B(2, 2:num_B:3), U_elem(2:num_B:3))
                Grad_u_Tensor(2,3) = DOT_PRODUCT(B(3, 3:num_B:3), U_elem(2:num_B:3))
                Grad_u_Tensor(3,1) = DOT_PRODUCT(B(1, 1:num_B:3), U_elem(3:num_B:3))  
                Grad_u_Tensor(3,2) = DOT_PRODUCT(B(2, 2:num_B:3), U_elem(3:num_B:3))  
                Grad_u_Tensor(3,3) = DOT_PRODUCT(B(3, 3:num_B:3), U_elem(3:num_B:3))
               
               
                Stress_Local = matmul(matmul(Rotation_Global_to_Local, Stress_Tensor),transpose(Rotation_Global_to_Local))
                Strain_Local = matmul(matmul(Rotation_Global_to_Local, Strain_Tensor),transpose(Rotation_Global_to_Local))
                Grad_u_Local = matmul(matmul(Rotation_Global_to_Local, Grad_u_Tensor),transpose(Rotation_Global_to_Local))
                
                
                
                
                
                

                r_polar = sqrt(GP_coords_local(1)**2 + GP_coords_local(2)**2)
                if(r_polar < SMALL_R) r_polar = SMALL_R
                theta_polar = atan2(GP_coords_local(2), GP_coords_local(1))
                

                do i_Mode = 1, 3
                    call Cal_3D_SIFs_IIM_Auxiliary_Fields_Part1(i_Mode, r_polar, theta_polar, &
                                                  Shear_Modulus, Poisson_Ratio, Kappa_3D, &
                                                  Aux_Stress_Local, Aux_Strain_Local, Aux_Grad_u_Local)
                    
                    dV = abs(Det_Jacobian) * GP_weights(i_GP)
                    
                    Wint = ZR
                    do i = 1, 3
                        do j = 1, 3
                            Wint = Wint + Stress_Local(i,j) * Aux_Strain_Local(i,j) 
                        enddo
                    enddo 
                    
                    A(1:3) = ZR
                    do jj = 1, 3
                        do i = 1, 3
                            A(jj) = A(jj) + Stress_Local(i,jj)     * Aux_Grad_u_Local(i,1) &
                                          + Aux_Stress_Local(i,jj) * Grad_u_Local(i,1)
                        enddo
                    enddo
                    
                    
                    Integrand_value = (A(1) - Wint) * grad_q_local(1) &
                                    + (A(2)) * grad_q_local(2)        &
                                    + (A(3)) * grad_q_local(3)
                                    
                                    
                    M_integral_Elements(i_Elem,i_Mode) = M_integral_Elements(i_Elem,i_Mode) + Integrand_value * dV
                enddo
                
                
                do i_Mode = 1, 3
                    call Cal_3D_SIFs_IIM_Auxiliary_Fields_Part2(i_Mode, r_polar, theta_polar,       &
                                         Shear_Modulus, Poisson_Ratio, Kappa_3D,            &
                                         Aux_dSig11_dX1, Aux_dSig12_dX2, Aux_dSig13_dX3,    &
                                         Aux_dSig12_dX1, Aux_dSig22_dX2, Aux_dSig23_dX3,    &
                                         Aux_dSig13_dX1, Aux_dSig23_dX2, Aux_dSig33_dX3,    &
                                         Aux_dEps11_dX1, Aux_dEps12_dX1, Aux_dEps13_dX1,    &
                                         Aux_dEps22_dX1, Aux_dEps23_dX1, Aux_dEps33_dX1,    &
                                         Aux_ddu1_dX1dX1, Aux_ddu2_dX1dX1, Aux_ddu3_dX1dX1, &
                                         Aux_ddu1_dX2dX1, Aux_ddu2_dX2dX1, Aux_ddu3_dX2dX1, & 
                                         Aux_ddu1_dX3dX1, Aux_ddu2_dX3dX1, Aux_ddu3_dX3dX1) 
                    
                    dV = abs(Det_Jacobian) * GP_weights(i_GP)
                    
                    Integrand_value1 =  Stress_Local(1,1)*(Aux_ddu1_dX1dX1 - Aux_dEps11_dX1) + &
                                        Stress_Local(1,2)*(Aux_ddu2_dX1dX1 - Aux_dEps12_dX1) + &
                                        Stress_Local(1,3)*(Aux_ddu3_dX1dX1 - Aux_dEps13_dX1) + &
                                        Stress_Local(1,2)*(Aux_ddu1_dX2dX1 - Aux_dEps12_dX1) + &
                                        Stress_Local(2,2)*(Aux_ddu2_dX2dX1 - Aux_dEps22_dX1) + &
                                        Stress_Local(2,3)*(Aux_ddu3_dX2dX1 - Aux_dEps23_dX1) + &
                                        Stress_Local(1,3)*(Aux_ddu1_dX3dX1 - Aux_dEps13_dX1) + &
                                        Stress_Local(2,3)*(Aux_ddu2_dX3dX1 - Aux_dEps23_dX1) + &
                                        Stress_Local(3,3)*(Aux_ddu3_dX3dX1 - Aux_dEps33_dX1) 

                    Integrand_value2 = (Aux_dSig11_dX1+Aux_dSig12_dX2+Aux_dSig13_dX3)*Grad_u_Local(1,1) + &
                                       (Aux_dSig12_dX1+Aux_dSig22_dX2+Aux_dSig23_dX3)*Grad_u_Local(2,1) + &
                                       (Aux_dSig13_dX1+Aux_dSig23_dX2+Aux_dSig33_dX3)*Grad_u_Local(3,1) 
                    
                                       
                    
                                                 
                    M_integral_Elements(i_Elem,i_Mode) = M_integral_Elements(i_Elem,i_Mode) + &
                               ((Integrand_value1 + Integrand_value2) * q_value_at_GP) * dV
                enddo






            enddo




            if (Has_Fluid_Pressure .and. abs(p_fluid) > 1.0D-16) then

                fac_two_faces = TWO
                
                call Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_3_Points(i_Crack, Current_Elem, nGP_face, &
                                             GP_face_global, GP_face_w, GP_face_n_global)  
                                             
                                             
                                             
                if (nGP_face > 0) then
                    do i_GP_face = 1, nGP_face
                        
                        n_unit(1:3) = GP_face_n_global(i_GP_face,1:3)
                        
                        
                        call Get_Normal_Insitu_Stress(Current_Elem,n_unit(1:3),Normal_Stress)
                        
                        t_global(1:3) = (p_fluid-Normal_Stress) * n_unit(1:3)
                        
                        t_local(1:3)  = matmul(Rotation_Global_to_Local, t_global)
                        
                        call Cal_3D_SIFs_IIM_q_Function(Current_Elem, q_nodal, Shape_Functions, &
                                                   dN_dx, dN_dy, dN_dz, &
                                                   GP_face_global(i_GP_face,1:3), Tip_Coords, &
                                                   e1_growth, e2_normal, e3_tangent, &
                                                   R_inner, R_outer, L_front, &
                                                   q_face, grad_q_global)

                        if (abs(q_face) < 1.0D-15) cycle
                        
                        


                        GP_coords_local = matmul(Rotation_Global_to_Local, GP_face_global(i_GP_face,1:3) - Tip_Coords)
                        
                        
                        r_polar = sqrt(GP_coords_local(1)**2 + GP_coords_local(2)**2)
                        
                        
                        if (r_polar < SMALL_R) r_polar = SMALL_R
                        theta_polar = atan2(GP_coords_local(2), GP_coords_local(1))
                        

                        do i_Mode = 1,3
                            call Cal_3D_SIFs_IIM_Auxiliary_Fields_Part1(i_Mode, r_polar, theta_polar, &
                                                          Shear_Modulus, Poisson_Ratio, Kappa_3D, &
                                                          Aux_Stress_Local, Aux_Strain_Local, Aux_Grad_u_Local)
                            Integrand_face = t_local(1) * Aux_Grad_u_Local(1,1) &
                                           + t_local(2) * Aux_Grad_u_Local(2,1) &
                                           + t_local(3) * Aux_Grad_u_Local(3,1)
                            M_integral_Elements(i_Elem,i_Mode) = M_integral_Elements(i_Elem,i_Mode) - &
                                fac_two_faces * s_trac * Integrand_face * q_face * GP_face_w(i_GP_face)
                        enddo

                    enddo
                endif

                if (allocated(GP_face_global))    deallocate(GP_face_global)
                if (allocated(GP_face_w))         deallocate(GP_face_w)
                if (allocated(GP_face_n_global))  deallocate(GP_face_n_global)

            endif
        enddo
!$OMP END PARALLEL DO

        if (L_front > 1.0d-10) then
            
            call Cal_3D_SIFs_IIM_q_Line_Integral(i_Crack, i_Vertex, L_front, q_line_integral)
            L_eff = q_line_integral
        else
            L_eff = 1.0d0
        endif
        

        
        
        
        M_integral(1) = sum(M_integral_Elements(1:Num_Domain_Elems,1))
        M_integral(2) = sum(M_integral_Elements(1:Num_Domain_Elems,2))
        M_integral(3) = sum(M_integral_Elements(1:Num_Domain_Elems,3))
        
        deallocate(M_integral_Elements)
        
        KI_3D(i_Crack)%row(i_Vertex)   = (E_prime / TWO * M_integral(1))/ L_eff
        KII_3D(i_Crack)%row(i_Vertex)  = (E_prime / TWO * M_integral(2))/ L_eff
        KIII_3D(i_Crack)%row(i_Vertex) = (Shear_Modulus * M_integral(3))/ L_eff

        if(allocated(Domain_Elem_List)) deallocate(Domain_Elem_List)
        if(allocated(q_nodal)) deallocate(q_nodal)
        
     
        if (Key_3D_SIFs_Print==1) then
            write(*,2001) i_Vertex,i_Crack,KI_3D(i_Crack)%row(i_Vertex)/1.0D6, &
                                           KII_3D(i_Crack)%row(i_Vertex)/1.0D6, &
                                           KIII_3D(i_Crack)%row(i_Vertex)/1.0D6
        endif
       

    enddo

enddo

deallocate(GP_xi, GP_eta, GP_zeta, GP_weights)

if (Key_Denoise_Vertex_Value>=1)then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,num_Cr_Edges,n_Sigma,logical_close)  
  do i_C = 1,num_Crack   
      num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
      if (i_C==1) then
          if(Tab_Num==5)then
              print *,'    Denoising vertex values...'   
          elseif(Tab_Num==8)then
              print *,'       Denoising vertex values...'   
          endif                  
      endif
      n_Sigma =2 
      
      
      logical_close = .True.
      if (Boundary_Cracks(i_C) .eqv. .True.) logical_close = .False.
      
      call Tool_Denoise_Data(KI_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,Key_Denoise_Vertex_Value,n_Sigma,logical_close) 
      call Tool_Denoise_Data(KII_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,Key_Denoise_Vertex_Value,n_Sigma,logical_close)       
      call Tool_Denoise_Data(KIII_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,Key_Denoise_Vertex_Value,n_Sigma,logical_close)      
  enddo
!$OMP END PARALLEL DO
endif
      
if (Key_Smooth_Vertex_Value>=1)then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,num_Cr_Edges,n_Sigma,logical_close)        
      do i_C = 1,num_Crack     
          num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
          if (i_C==1) then
              if(Tab_Num==5)then
                  print *,'    Smoothing vertex values...' 
              elseif(Tab_Num==8)then
                  print *,'       Smoothing vertex values...'   
              endif                
          endif
          
          logical_close = .True.
          if (Boundary_Cracks(i_C) .eqv. .True.) logical_close = .False.
          
          call Tool_Smooth_Data(KI_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,&
                  Key_Smooth_Vertex_Value,Smooth_Vertex_n, logical_close)
          call Tool_Smooth_Data(KII_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,&
                  Key_Smooth_Vertex_Value,Smooth_Vertex_n,logical_close)   
          call Tool_Smooth_Data(KIII_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,&
                  Key_Smooth_Vertex_Value,Smooth_Vertex_n,logical_close)      
      enddo
!$OMP END PARALLEL DO 
             
      if (Key_Smooth_Vertex_Value2>=1)then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,num_Cr_Edges,n_Sigma,logical_close)           
          do i_C = 1,num_Crack     
              num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
              if (i_C==1) then
                  if(Tab_Num==5)then
                      print *,'    Smoothing vertex values...'
                  elseif(Tab_Num==8)then
                      print *,'       Smoothing vertex values...'   
                  endif 
              endif
              
             logical_close = .True.
             if (Boundary_Cracks(i_C) .eqv. .True.) logical_close = .False.
          
          
              call Tool_Smooth_Data(KI_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,&
                                Key_Smooth_Vertex_Value2,Smooth_Vertex_n2,logical_close) 
              call Tool_Smooth_Data(KII_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,&
                                Key_Smooth_Vertex_Value2,Smooth_Vertex_n2,logical_close) 
              call Tool_Smooth_Data(KIII_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,&
                                Key_Smooth_Vertex_Value2,Smooth_Vertex_n2,logical_close) 
          enddo
!$OMP END PARALLEL DO               
      endif           
endif  

 
      
write(temp,'(I5)') iter
if(Key_Save_Nothing /= 1) then 
    if(Tab_Num==5)then
          print *,'    Saving cvk1,cvk2,cvk3 files...'
    elseif(Tab_Num==8)then
          print *,'       Saving cvk1,cvk2,cvk3 files...'
    endif 
    c_File_name_test   = trim(Full_Pathname)//'.cvk1_'//ADJUSTL(temp)  
    open(101,file=c_File_name_test,status='unknown')  
    c_File_name_test   = trim(Full_Pathname)//'.cvk2_'//ADJUSTL(temp)  
    open(102,file=c_File_name_test,status='unknown')  
    c_File_name_test   = trim(Full_Pathname)//'.cvk3_'//ADJUSTL(temp)  
    open(103,file=c_File_name_test,status='unknown')  
    do i_C = 1,num_Crack   
          num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
          write(101, '(50000E20.8)') (KI_3D(i_C)%row(i_V),i_v=1,num_Cr_Edges)
          write(102, '(50000E20.8)') (KII_3D(i_C)%row(i_V),i_v=1,num_Cr_Edges)
          write(103, '(50000E20.8)') (KIII_3D(i_C)%row(i_V),i_v=1,num_Cr_Edges)     
    enddo
    close(101)
    close(102)
    close(103)
endif      
      
      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,Num_CrMesh_Outlines,    &
!$OMP      c_OUT_Elem_Old,Final_Point,Final_Point_Prop,&
!$OMP      i_Out_Node,c_Mesh_Node,c_x_old,c_y_old,c_z_old,&
!$OMP      num_of_check,Check_Theta,Spe_Pri_Stress,Sai_Cr,Ld_Cr,&
!$OMP      Ks_Cr,D_Cr,KI_eq,Theta_All,r_for_c_Theta_finding,&
!$OMP      tem_root_2pir,c_Mesh_Node_Next,in_Ele_Num,mat_num,c_KIc,&
!$OMP      c_x,c_y,c_z,c_OUT_Elem,c_KI_3D,c_KII_3D,c_KIII_3D,&
!$OMP      i_Check_Theta,Schollm_Max_Theta_in_pi,tem_part1,tem_part2,&
!$OMP      cc_Part1,cc_Part2,cc_Part3,Stress_Theta,Tao_Theta,&
!$OMP      i_Theta_min,c_Theta,Ele_Num_Cache)    
do i_C =1,num_Crack  
    Ele_Num_Cache = 1
    Num_CrMesh_Outlines = Crack3D_Meshed_Outline_num(i_C)
    c_OUT_Elem_Old(1:5000)    = 1
    Final_Point(1:5000,1:3)   = ZR 
    Final_Point_Prop(1:5000)  = 0
    do i_Out_Node = 1,Num_CrMesh_Outlines
      c_Mesh_Node = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
      c_x_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1) 
      c_y_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,2) 
      c_z_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,3)
      Final_Point(i_Out_Node,1)=c_x_old
      Final_Point(i_Out_Node,2)=c_y_old
      Final_Point(i_Out_Node,3)=c_z_old
    enddo
    num_of_check = 720
    if(allocated(Check_Theta)) deallocate(Check_Theta)
    if(allocated(Spe_Pri_Stress)) deallocate(Spe_Pri_Stress)
    if(allocated(Sai_Cr)) deallocate(Sai_Cr)
    if(allocated(Ld_Cr)) deallocate(Ld_Cr)
    if(allocated(Ks_Cr)) deallocate(Ks_Cr)
    if(allocated(D_Cr)) deallocate(D_Cr)
    if(allocated(KI_eq)) deallocate(KI_eq)
    if(allocated(Theta_All)) deallocate(Theta_All)
    allocate(Check_Theta(num_of_check+1))
    allocate(Spe_Pri_Stress(num_of_check+1))
    allocate(Sai_Cr(Num_CrMesh_Outlines))
    allocate(Ld_Cr(Num_CrMesh_Outlines))
    allocate(Ks_Cr(Num_CrMesh_Outlines))  
    allocate(D_Cr(Num_CrMesh_Outlines))  
    allocate(KI_eq(Num_CrMesh_Outlines))
    allocate(Theta_All(Num_CrMesh_Outlines))
    
    r_for_c_Theta_finding = 1.0*Ave_Elem_L
    tem_root_2pir = sqrt(TWO*pi*r_for_c_Theta_finding)
    
    
    do i_Out_Node = 1,Num_CrMesh_Outlines
      c_Mesh_Node = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
      c_Mesh_Node_Next =Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,2)
      in_Ele_Num  = Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(c_Mesh_Node)
      if(in_Ele_Num ==0 ) then
          cycle
      endif
  
      if(Key_XA/=2) then
          mat_num = Elem_Mat(in_Ele_Num)
          c_KIc  = Material_Para(mat_num,6)
      else
          c_KIc  = Elem_KIc_XA(in_Ele_Num)
      endif
          
      c_x  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1) 
      c_y  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,2) 
      c_z  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,3) 
      
      call Cal_Ele_Num_by_Coors_3D(c_x,c_y,c_z,Ele_Num_Cache,c_OUT_Elem)
      if(c_OUT_Elem ==0) then
          cycle
      endif
  
      Ld_Cr(i_Out_Node) = Tool_Function_2Point_Dis_3D(Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1:3),&
                                                      Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node_Next,1:3))

      Ks_Cr(i_Out_Node) = ONE/Ld_Cr(i_Out_Node)
      D_Cr(i_Out_Node) = Ld_Cr(i_Out_Node)**3
      
      c_KI_3D   = KI_3D(i_C)%row(i_Out_Node)
      c_KII_3D  = KII_3D(i_C)%row(i_Out_Node)
      c_KIII_3D = KIII_3D(i_C)%row(i_Out_Node)
      do i_Check_Theta =1,num_of_check+1
        Schollm_Max_Theta_in_pi = Schollm_Max_Theta*pi/Con_180
        Check_Theta(i_Check_Theta)=-Schollm_Max_Theta_in_pi+TWO*Schollm_Max_Theta_in_pi/num_of_check*(i_Check_Theta-1)   
 
        tem_part1 = THR*cos(Check_Theta(i_Check_Theta)/TWO) +cos(THR*Check_Theta(i_Check_Theta)/TWO)
        tem_part2 = THR*sin(Check_Theta(i_Check_Theta)/TWO) +THR*sin(THR*Check_Theta(i_Check_Theta)/TWO)     
        Stress_Theta = c_KI_3D/FOU/tem_root_2pir*tem_part1 -c_KII_3D/FOU/tem_root_2pir*tem_part2
        Tao_Theta=c_KIII_3D*cos(Check_Theta(i_Check_Theta)/TWO)/tem_root_2pir
        Spe_Pri_Stress(i_Check_Theta) = Stress_Theta/TWO + ZP5*sqrt(Stress_Theta**2 + FOU*Tao_Theta**2)
      enddo 
      i_Theta_min = maxloc(Spe_Pri_Stress(1:num_of_check+1),1)
      c_Theta = Check_Theta(i_Theta_min)
      Theta_All(i_Out_Node) = c_Theta
      
      cc_Part1 = c_KI_3D*(cos(c_Theta/TWO))**2
      cc_Part2 = THR/TWO*c_KII_3D*sin(c_Theta)
      cc_Part3 =sqrt((cc_Part1 - cc_Part2)**2 + FOU*c_KIII_3D**2)
      KI_eq(i_Out_Node)  =ZP5*cos(c_Theta/TWO)*(cc_Part1-cc_Part2+cc_Part3)
      KI_eq_3D(i_C)%row(i_Out_Node) = KI_eq(i_Out_Node)
    end do
enddo      
!$OMP END PARALLEL DO   
      
if(Key_Save_Nothing  == 0) then
    if(Tab_Num==5)then
        print *,'    Saving cvke file for cracks...'
    elseif(Tab_Num==8)then
        print *,'       Saving cvke file for cracks...'
    endif

    write(temp,'(I5)') iter
    c_File_name_test = trim(Full_Pathname)//'.cvke_'//ADJUSTL(temp)  
    open(104,file=c_File_name_test,status='unknown')  
    do i_C = 1,num_Crack   
          num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
          write(104, '(50000E20.12)') (KI_eq_3D(i_C)%row(i_v) ,i_v=1,num_Cr_Edges)   
    enddo
    close(104)
endif    

return

end subroutine Cal_SIFs_IIM_3D



