!-----------------------------------------------------------
! Brief: Declare the interface for Cal_HF_Initial_Pressure which
!        computes the initial crack-water pressure for the first
!        iteration of an HF step.
!
! Notes:   Interface declaration only; passed crack-tip growth flags,
!          previous-step pressure and point counts.
!-----------------------------------------------------------

module Global_Inter_Cal_HF_Initial_Pressure
INTERFACE
    subroutine Cal_HF_Initial_Pressure(iter,ifra,Counter, Yes_Growth, Last_Cracks_CalP_Num_ifra, Last_Cr_CalP_Pres_ifra, &
    Initial_Cr_CalP_Pres)
    use Global_Float_Type
    use Global_Common
    use Global_Crack
    use Global_Crack_Common
    use Global_HF
    implicit none
    integer,intent(in)::iter,ifra,Counter
    real(kind=FT),intent(in)::Last_Cr_CalP_Pres_ifra(Max_Num_Cr,Max_Num_Cr_CalP)
    integer,intent(in)::Last_Cracks_CalP_Num_ifra(Max_Num_Cr)
    real(kind=FT),intent(out)::Initial_Cr_CalP_Pres(Max_Num_Cr,Max_Num_Cr_CalP)
    logical,intent(in)::Yes_Growth(Max_Num_Cr,2)
END subroutine Cal_HF_Initial_Pressure
END INTERFACE
end module Global_Inter_Cal_HF_Initial_Pressure

!-----------------------------------------------------------
! Brief: Declare the interface for Cal_HF_Resid which assembles the
!        dense coupled residual for the HF mechanical-fluid system.
!
! Notes:   Interface declaration only; receives global K, Q, F_U,
!          last-step displacements, and HF state arrays.
!-----------------------------------------------------------

module Global_Inter_Cal_HF_Resid
INTERFACE
    subroutine Cal_HF_Resid(ifra,iter,Counter_Iter,R,Total_FD, num_FreeD,num_free_CalP, globalK,Coupled_Q,F_U, &
    Last_DISP,Last_Last_DISP, freeDOF,freeDOF_HF,Local_freeDOF_HF, Last_CalP_Pres,delta_Time,H,c_S)
    use Global_Float_Type
    use Global_Crack
    use Global_Crack_Common
    use Global_HF
    use Global_Material
    implicit none
    integer, intent(in)::ifra,iter,Counter_Iter,Total_FD,num_FreeD,num_free_CalP
    real(kind=FT),intent(in):: Coupled_Q(Total_FD,num_Tol_CalP_Water), globalK(Total_FD,Total_FD), &
    F_U(Total_FD),Last_DISP(Total_FD), Last_Last_DISP(Total_FD), Last_CalP_Pres(num_Tol_CalP_Water),delta_Time, &
    H(num_Tol_CalP_Water,num_Tol_CalP_Water), c_S(num_Tol_CalP_Water)
    integer, intent(in)::freeDOF(Total_FD), freeDOF_HF(num_Tol_CalP_Water), Local_freeDOF_HF(num_Tol_CalP_Water)
    real(kind=FT),intent(out)::R(Total_FD+num_Tol_CalP_Water)
end subroutine Cal_HF_Resid
END INTERFACE
end module Global_Inter_Cal_HF_Resid

!-----------------------------------------------------------
! Brief: Declare the interface for Cal_HF_Jacobian_NR which builds
!        the dense NR Jacobian of the coupled HF system.
!
! Notes:   Interface declaration only; uses globalK, Coupled_Q, free
!          DOF lists and the fluid leak-off matrix H.
!-----------------------------------------------------------

module Global_Inter_Cal_HF_Jacobian_NR
INTERFACE
    subroutine Cal_HF_Jacobian_NR(Counter,NR_Deri,Total_FD, num_FreeD,num_free_CalP,globalK,Coupled_Q, &
    freeDOF,freeDOF_HF,Local_freeDOF_HF,delta_Time,H)
    use Global_Float_Type
    use Global_Crack
    use Global_Crack_Common
    implicit none
    integer, intent(in)::Counter,Total_FD,num_FreeD,num_free_CalP
    real(kind=FT),intent(in):: Coupled_Q(Total_FD,num_Tol_CalP_Water), globalK(Total_FD,Total_FD),delta_Time, &
    H(num_Tol_CalP_Water,num_Tol_CalP_Water)
    integer, intent(in)::freeDOF(Total_FD), freeDOF_HF(num_Tol_CalP_Water), Local_freeDOF_HF(num_Tol_CalP_Water)
    real(kind=FT),intent(out):: NR_Deri(Total_FD+num_Tol_CalP_Water, Total_FD+num_Tol_CalP_Water)
end subroutine Cal_HF_Jacobian_NR
END INTERFACE
end module Global_Inter_Cal_HF_Jacobian_NR

!-----------------------------------------------------------
! Brief: Declare the interface for the dense linear-system-of-
!        equations solver Matrix_Solve_LSOE.
!
! Notes:   Interface declaration only; solves D = K^-1 F for a dense
!          K matrix.
!-----------------------------------------------------------

module Global_Inter_Matrix_Solve_LSOE
INTERFACE
    SUBROUTINE Matrix_Solve_LSOE(Key_Indent,Key_LSOE_Sys,c_Key_SLOE,K,F,D,n)
    use Global_Float_Type
    implicit none
    integer,intent(in)::n,c_Key_SLOE,Key_LSOE_Sys
    integer,intent(in)::Key_Indent
    real(kind=FT),intent(in)::K(n,n),F(n)
    real(kind=FT),intent(out)::D(n)
END SUBROUTINE Matrix_Solve_LSOE
END INTERFACE
end module Global_Inter_Matrix_Solve_LSOE

!-----------------------------------------------------------
! Brief: Declare the interface for Get_Contact_State_Eles_IN which
!        classifies elements as contacting or free based on the
!        crack-opening distribution.
!
! Notes:   Interface declaration only; receives per-crack opening
!          and the previous-step contact state.
!-----------------------------------------------------------

module Global_Inter_Get_Contact_State_Eles_IN
INTERFACE
    SUBROUTINE Get_Contact_State_Eles_IN(i_Contact, c_Cracks_CalP_Aper, tem_Elem_Conta_Sta, Yes_Contact, &
    Yes_Contact_Con,num_Contact_Ele)
    use Global_Float_Type
    use Global_Crack
    use Global_Crack_Common
    use Global_Model
    use Global_Common
    use Global_Contact
    use Global_HF
    use Global_Elem_Area_Vol
    implicit none
    integer,intent(in)::i_Contact
    integer,intent(inout)::tem_Elem_Conta_Sta(Num_Elem,Max_Num_Cr)
    real(kind=FT),intent(in)::c_Cracks_CalP_Aper(Max_Num_Cr,Max_Num_Cr_CalP)
    logical,intent(out)::Yes_Contact_Con,Yes_Contact
    integer,intent(out)::num_Contact_Ele
END SUBROUTINE Get_Contact_State_Eles_IN
END INTERFACE
end module Global_Inter_Get_Contact_State_Eles_IN

!-----------------------------------------------------------
! Brief: Declare the interface for Cal_HF_Line_Searching which
!        performs a backtracking line search to enforce the
!        coupling-convergence condition for the HF NR step.
!
! Notes:   Interface declaration only; uses the NR direction delta_x
!          and the current residual.
!-----------------------------------------------------------

module Global_Inter_Cal_HF_Line_Searching
INTERFACE
    subroutine Cal_HF_Line_Searching(ifra,iter,Counter_Iter, num_ALlDOF,c_Total_FD,num_FreeD,num_free_CalP, &
    Total_freeDOF,c_num_Tol_CalP_Water,num_Total_FD,freeDOF_HF, Local_freeDOF_HF,NR_Deri, delta_x, c_R,Initial_DISP, &
    c_DISP,c_freeDOF,c_CalP_Pres,c_globalK,Coupled_Q,F_U, c_Temp_Cr_CalP_Aper,c_Temp_total_Time, &
    Num_Line_Search,total_Time,x)
    use Global_Float_Type
    use Global_Common
    use Global_Crack
    use Global_Crack_Common
    use Global_HF
    use Global_Material
    implicit none
    integer,intent(in)::ifra,iter,Counter_Iter,c_Total_FD, num_FreeD,num_free_CalP,num_ALlDOF,c_num_Tol_CalP_Water, &
    num_Total_FD
    integer,intent(in)::Local_freeDOF_HF(c_num_Tol_CalP_Water)
    integer,intent(in)::Total_freeDOF(num_ALlDOF)
    real(kind=FT),intent(in)::c_R(num_ALlDOF)
    real(kind=FT),intent(in)::Coupled_Q(c_Total_FD,c_num_Tol_CalP_Water), c_globalK(c_Total_FD,c_Total_FD),F_U(c_Total_FD)
    real(kind=FT),intent(in)::delta_x(num_FreeD+num_free_CalP)
    real(kind=FT),intent(in)::c_DISP(c_Total_FD)
    real(kind=FT),intent(in)::Initial_DISP(c_Total_FD)
    integer,intent(in)::c_freeDOF(c_Total_FD)
    integer,intent(in)::freeDOF_HF(c_num_Tol_CalP_Water)
    real(kind=FT),intent(in)::c_CalP_Pres(c_num_Tol_CalP_Water)
    real(kind=FT),intent(in)::c_Temp_Cr_CalP_Aper(Max_Num_Cr,Max_Num_Cr_CalP)
    real(kind=FT),intent(in)::c_Temp_total_Time
    real(kind=FT),intent(in)::NR_Deri(c_Total_FD+c_num_Tol_CalP_Water,c_Total_FD+c_num_Tol_CalP_Water)
    integer,intent(inout)::Num_Line_Search
    real(kind=FT),intent(inout)::total_Time
    real(kind=FT),intent(out)::x(num_FreeD+num_free_CalP)
end subroutine Cal_HF_Line_Searching
END INTERFACE
end module Global_Inter_Cal_HF_Line_Searching

!-----------------------------------------------------------
! Brief: Declare the interface for the iterative 2D contact-state
!        driver Determine_Contact_State_by_Iteration.
!
! Notes:   Interface declaration only; updates the stiffness matrix
!          and displacement vector after the contact pass.
!-----------------------------------------------------------

module Global_Inter_Determine_Contact_State_by_Iteration
INTERFACE
    SUBROUTINE Determine_Contact_State_by_Iteration(iter,ifra,Counter_Iter, Contact_DISP,c_Total_Freedom,usual_FD,enrich_FD, &
    c_freeDOF,c_num_freeDOF,c_F,ori_globalK,c_globalK)
    use Global_Float_Type
    use Global_Common
    use Global_Filename
    use Global_Model
    use Global_Elem_Area_Vol
    use Global_Crack
    use Global_Crack_Common
    use Global_HF
    use Global_Contact
    implicit none
    integer,intent(in)::iter,ifra,Counter_Iter
    integer,intent(in)::c_Total_Freedom,c_num_freeDOF
    integer,intent(in)::usual_FD,enrich_FD
    real(kind=FT),intent(inout)::Contact_DISP(c_Total_Freedom)
    integer,intent(in)::c_freeDOF(c_Total_Freedom)
    real(kind=FT),intent(in)::c_F(c_Total_Freedom)
    real(kind=FT),intent(in)::ori_globalK(c_Total_Freedom,c_Total_Freedom)
    real(kind=FT),intent(out)::c_globalK(c_Total_Freedom,c_Total_Freedom)
END SUBROUTINE Determine_Contact_State_by_Iteration
END INTERFACE
end module Global_Inter_Determine_Contact_State_by_Iteration

!-----------------------------------------------------------
! Brief: Declare the interface for the CSR-format sparse FEM
!        stiffness assembly Assemble_Stiffness_Matrix_SPARS_FEM.
!
! Notes:   Interface declaration only; outputs the CSR row/column
!          pointers and value arrays plus the total Gauss count.
!-----------------------------------------------------------

module Global_Inter_Assemble_Stiffness_Matrix_SPARS_FEM
INTERFACE
    SUBROUTINE Assemble_Stiffness_Matrix_SPARS_FEM(isub,freeDOF,num_FreeD, &
    K_CSR_aa,K_CSR_ja,K_CSR_ia,K_CSR_NNZ_Max,K_CSR_NNZ, T_Freedom,Total_Num_G_P)
    use Global_Float_Type
    use Global_Model
    use Global_Filename
    use Global_Common
    use Global_Material
    implicit none
    integer,intent(in)::isub,num_FreeD,T_Freedom,K_CSR_NNZ_Max
    integer,intent(in)::freeDOF(1:num_FreeD)
    real(kind=FT),intent(out)::K_CSR_aa(K_CSR_NNZ_Max)
    integer,intent(out)::K_CSR_ja(K_CSR_NNZ_Max)
    integer,intent(out)::K_CSR_ia(num_FreeD+1)
    integer,intent(out)::K_CSR_NNZ
    integer,intent(out)::Total_Num_G_P
END SUBROUTINE Assemble_Stiffness_Matrix_SPARS_FEM
END INTERFACE
end module Global_Inter_Assemble_Stiffness_Matrix_SPARS_FEM

!-----------------------------------------------------------
! Brief: Declare the interface for the CSR-format sparse XFEM
!        stiffness assembly Assemble_Stiffness_Matrix_SPARS_XFEM.
!
! Notes:   Interface declaration only; uses OpenMP and a progress
!          bar; supports both free and fixed DOF lists.
!-----------------------------------------------------------

module Global_Inter_Assemble_Stiffness_Matrix_SPARS_XFEM
INTERFACE
    SUBROUTINE Assemble_Stiffness_Matrix_SPARS_XFEM(isub, freeDOF,num_FreeD, fixedDOF,num_FixedD, &
    K_CSR_aa,K_CSR_ja,K_CSR_ia, K_CSR_NNZ_Max,K_CSR_NNZ, T_Freedom,Total_Num_G_P)
    use Global_Float_Type
    use Global_Crack
    use Global_Crack_Common
    use Global_Model
    use Global_Filename
    use Global_Common
    use Global_Material
    use Global_Contact
    use omp_lib
    use CMD_Progress
    use ISO_FORTRAN_ENV

    implicit none
    integer,intent(in)::isub,num_FreeD,T_Freedom,K_CSR_NNZ_Max
    integer,intent(in)::freeDOF(1:num_FreeD)
    integer,intent(out)::Total_Num_G_P
    real(kind=FT),intent(out)::K_CSR_aa(K_CSR_NNZ_Max)
    integer,intent(out)::K_CSR_ja(K_CSR_NNZ_Max)
    integer,intent(out)::K_CSR_ia(num_FreeD+1)
    integer,intent(out)::K_CSR_NNZ
    integer,intent(in)::num_FixedD
    integer,intent(in)::fixedDOF(1:num_FixedD)

END SUBROUTINE Assemble_Stiffness_Matrix_SPARS_XFEM
END INTERFACE
end module Global_Inter_Assemble_Stiffness_Matrix_SPARS_XFEM

!-----------------------------------------------------------
! Brief: Declare the interface for the sparse linear solver
!        Matrix_Solve_LSOE_Sparse operating on CSR K matrices.
!
! Notes:   Interface declaration only; receives CSR aa/ja/ia and a
!          dense right-hand side F.
!-----------------------------------------------------------

module Global_Inter_Matrix_Solve_LSOE_Sparse
INTERFACE
    SUBROUTINE Matrix_Solve_LSOE_Sparse(Key_Indent,Key_LSOE_Sys, c_Key_SLOE,K_CSR_NNZ, K_CSR_aa,K_CSR_ja,K_CSR_ia,F,D,n)
    use Global_Float_Type
    use Global_Common
    use Global_Model
    implicit none
    integer,intent(in)::n,c_Key_SLOE,Key_LSOE_Sys,K_CSR_NNZ
    integer,intent(in)::Key_Indent
    real(kind=FT),intent(in)::K_CSR_aa(1:K_CSR_NNZ)
    integer,intent(in)::K_CSR_ja(1:K_CSR_NNZ)
    integer,intent(in)::K_CSR_ia(1:n+1)
    real(kind=FT),intent(in)::F(n)
    real(kind=FT),intent(out)::D(n)
END SUBROUTINE Matrix_Solve_LSOE_Sparse
END INTERFACE
end module Global_Inter_Matrix_Solve_LSOE_Sparse

!-----------------------------------------------------------
! Brief: Declare the interface for the sparse-version HF residual
!        assembler Cal_HF_Resid_SPARS.
!
! Notes:   Interface declaration only; identical coupling to
!          Cal_HF_Resid but uses CSR storage.
!-----------------------------------------------------------

module Global_Inter_Cal_HF_Resid_SPARS
INTERFACE
    subroutine Cal_HF_Resid_SPARS(ifra,iter,Counter_Iter,R, Total_FD,num_FreeD,num_free_CalP, K_CSR_NNZ,K_CSR_aa, &
    K_CSR_ja,K_CSR_ia, Coupled_Q,F_U,Last_DISP,Last_Last_DISP, freeDOF,freeDOF_HF,Local_freeDOF_HF, &
    Last_CalP_Pres,delta_Time,H,c_S)
    use Global_Float_Type
    use Global_Crack
    use Global_Crack_Common
    use Global_HF
    use Global_Material
    implicit none
    integer, intent(in)::ifra,iter,Counter_Iter,Total_FD,num_FreeD,num_free_CalP,K_CSR_NNZ
    real(kind=FT),intent(in)::Coupled_Q(Total_FD,num_Tol_CalP_Water), &
    F_U(Total_FD),Last_DISP(Total_FD),Last_Last_DISP(Total_FD), Last_CalP_Pres(num_Tol_CalP_Water),delta_Time, &
    H(num_Tol_CalP_Water,num_Tol_CalP_Water), c_S(num_Tol_CalP_Water)
    integer, intent(in)::freeDOF(Total_FD), freeDOF_HF(num_Tol_CalP_Water), Local_freeDOF_HF(num_Tol_CalP_Water)
    real(kind=FT),intent(out)::R(Total_FD+num_Tol_CalP_Water)
    real(kind=FT),intent(in)::K_CSR_aa(1:K_CSR_NNZ)
    integer,intent(in)::K_CSR_ja(1:K_CSR_NNZ)
    integer,intent(in)::K_CSR_ia(1:num_FreeD+1)
END SUBROUTINE Cal_HF_Resid_SPARS
END INTERFACE
end module Global_Inter_Cal_HF_Resid_SPARS

!-----------------------------------------------------------
! Brief: Declare the interface for the sparse NR Jacobian builder
!        Cal_HF_Jacobian_NR_SPARS.
!
! Notes:   Interface declaration only; emits a CSR Jacobian.
!-----------------------------------------------------------

module Global_Inter_Cal_HF_Jacobian_NR_SPARS
INTERFACE
    subroutine Cal_HF_Jacobian_NR_SPARS(Counter,Total_FD, num_FreeD,num_free_CalP, K_CSR_NNZ,K_CSR_aa,K_CSR_ja,K_CSR_ia, &
    NRD_CSR_NNZ_Max, NRD_CSR_NNZ,NRD_CSR_aa,NRD_CSR_ja,NRD_CSR_ia, Coupled_Q,freeDOF,freeDOF_HF, &
    Local_freeDOF_HF,delta_Time,H)
    use Global_Float_Type
    use Global_Crack
    use Global_Crack_Common
    implicit none
    integer, intent(in)::Counter,Total_FD,num_FreeD,num_free_CalP,NRD_CSR_NNZ_Max
    integer, intent(out)::NRD_CSR_NNZ
    integer, intent(in)::K_CSR_NNZ
    real(kind=FT),intent(in)::K_CSR_aa(K_CSR_NNZ)
    integer,intent(in)::K_CSR_ja(K_CSR_NNZ)
    integer,intent(in)::K_CSR_ia(num_FreeD+1)
    real(kind=FT),intent(out)::NRD_CSR_aa(NRD_CSR_NNZ_Max)
    integer,intent(out)::NRD_CSR_ja(NRD_CSR_NNZ_Max)
    integer,intent(out)::NRD_CSR_ia(num_FreeD+num_free_CalP+1)
    real(kind=FT),intent(in):: Coupled_Q(Total_FD,num_Tol_CalP_Water), delta_Time,H(num_Tol_CalP_Water,num_Tol_CalP_Water)
    integer, intent(in)::freeDOF(Total_FD), freeDOF_HF(num_Tol_CalP_Water), Local_freeDOF_HF(num_Tol_CalP_Water)
END SUBROUTINE Cal_HF_Jacobian_NR_SPARS
END INTERFACE
end module Global_Inter_Cal_HF_Jacobian_NR_SPARS

!-----------------------------------------------------------
! Brief: Declare the interface for the contact-reduced residual
!        assembler Cal_Contact_Red_Resid used in 2D HF-Contact NR.
!
! Notes:   Interface declaration only; uses the Schur-complement of
!          the FEM part of the system.
!-----------------------------------------------------------

module Global_Inter_Cal_Contact_Red_Resid
INTERFACE
    subroutine Cal_Contact_Red_Resid( iter,ifra,Counter_Iter,i_NR_P, c_Total_Freedom,c_num_freeDOF,n_freeDOF_FEM,F_U, &
    freeDOF_FEM, U_e,U_all,freeDOF,Part1,Part2,K_ee, K_uu_Inv,K_ue,PC_Gauss_x,PC_Gauss_y,R_red)
    use Global_Float_Type
    use Global_Common
    use Global_Model
    use Global_Crack
    use Global_Crack_Common
    use Global_HF
    use Global_Material
    implicit none
    integer,intent(in)::iter,ifra,Counter_Iter,i_NR_P
    integer,intent(in)::c_Total_Freedom,c_num_freeDOF
    integer,intent(in)::freeDOF_FEM(n_freeDOF_FEM),n_freeDOF_FEM
    real(kind=FT),intent(in):: F_U(c_Total_Freedom),U_e(Enrich_Freedom), U_all(c_Total_Freedom)
    integer, intent(in)::freeDOF(c_Total_Freedom)
    real(kind=FT),intent(in)::Part1(Enrich_Freedom,Enrich_Freedom),Part2(Enrich_Freedom)
    real(kind=FT),intent(in)::K_ee(Enrich_Freedom,Enrich_Freedom)
    real(kind=FT),intent(in)::K_uu_Inv(n_freeDOF_FEM,n_freeDOF_FEM)
    real(kind=FT),intent(in)::K_ue(n_freeDOF_FEM,Enrich_Freedom)
    real(kind=FT),intent(out)::R_red(Enrich_Freedom)
    real(kind=FT),intent(in)::PC_Gauss_x(num_Crack,Max_Num_Cr_CalP-1,2),PC_Gauss_y(num_Crack,Max_Num_Cr_CalP-1,2)
END SUBROUTINE Cal_Contact_Red_Resid
END INTERFACE
end module Global_Inter_Cal_Contact_Red_Resid

!-----------------------------------------------------------
! Brief: Declare the interface for Cal_B_Matrix_Crack_3D which
!        builds the strain-displacement B matrix for a 3D crack
!        element at a given Gauss point.
!
! Notes:   Interface declaration only; returns up to 100 B rows.
!-----------------------------------------------------------

module Global_Inter_Cal_B_Matrix_Crack_3D
INTERFACE
    subroutine Cal_B_Matrix_Crack_3D(kesi,yita,zeta,i_C,i_E,i_G, c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES, tem_B,num_tem_B)
    use Global_Float_Type
    use Global_Crack
    use Global_Crack_Common
    use Global_Model
    use Global_Filename
    use Global_Common
    use Global_Material
    implicit none
    integer,intent(in)::i_C,i_E,i_G
    real(kind=FT),intent(in)::c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
    integer,intent(in)::c_NN(8)
    real(kind=FT),intent(in)::kesi,yita,zeta
    real(kind=FT),intent(out)::tem_B(6,100)
    integer,intent(out)::num_tem_B
END SUBROUTINE Cal_B_Matrix_Crack_3D
END INTERFACE
end module Global_Inter_Cal_B_Matrix_Crack_3D

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D hex shape-function and
!        Jacobian evaluator Cal_N_dNdkesi_J_detJ_3D.
!
! Notes:   Interface declaration only; returns N (3x24), dNdkesi (8x3),
!          J (3x3) and detJ for a single Gauss point.
!-----------------------------------------------------------

module Global_Inter_Cal_N_dNdkesi_J_detJ_3D
INTERFACE
    subroutine Cal_N_dNdkesi_J_detJ_3D(kesi,yita,zeta, X_NODES,Y_NODES,Z_NODES, detJ,J,N,dNdkesi)
    use Global_Float_Type
    implicit none
    real(kind=FT),intent(in)::kesi,yita,zeta,X_NODES(8),Y_NODES(8),Z_NODES(8)
    real(kind=FT),intent(out)::detJ,dNdkesi(8,3)
    real(kind=FT),intent(out):: N(3,24),J(3,3)
END SUBROUTINE Cal_N_dNdkesi_J_detJ_3D
END INTERFACE
end module Global_Inter_Cal_N_dNdkesi_J_detJ_3D

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D XFEM element-DOF locator
!        Location_Element_Stiff_Matrix_3D.
!
! Notes:   Interface declaration only; outputs the FEM+XFEM DOF
!          positions and counts.
!-----------------------------------------------------------

module Global_Inter_Location_Element_Stiff_Matrix_3D
INTERFACE
    SUBROUTINE Location_Element_Stiff_Matrix_3D(i_E,i_C, c_POS_3D_c_Ele,Location_ESM_C_Crack, num_Loc_ESM_C_Crack, &
    Location_ESM_C_Cr_NoFEM, num_Loc_ESM_C_Cr_NoFEM)
    use Global_Float_Type
    use Global_Crack
    use Global_Crack_Common
    use Global_Model
    use Global_Common
    implicit none
    integer,intent(in)::i_E,i_C,c_POS_3D_c_Ele(8)
    integer,intent(out)::Location_ESM_C_Crack(100)
    integer,intent(out)::Location_ESM_C_Cr_NoFEM(100)
    integer,intent(out)::num_Loc_ESM_C_Crack
    integer,intent(out)::num_Loc_ESM_C_Cr_NoFEM
END SUBROUTINE Location_Element_Stiff_Matrix_3D
END INTERFACE
end module Global_Inter_Location_Element_Stiff_Matrix_3D

!-----------------------------------------------------------
! Brief: Declare the interface for Tool_Dis_Point_to_3D_Quad which
!        computes the minimum distance from a point to a 3D quad.
!
! Notes:   Interface declaration only; also returns the closest
!          point and the inside/on flags.
!-----------------------------------------------------------

module Global_Inter_Tool_Cal_Dis_Point_to_3D_Quad
INTERFACE
    subroutine Tool_Dis_Point_to_3D_Quad(Point, Quad_P1,Quad_P2,Quad_P3,Quad_P4, Distance,PER,Yes_PER_in,Yes_PER_on)
    use Global_Float_Type
    use Global_Common
    implicit none
    real(kind=FT),intent(in)::Point(3),Quad_P1(3),Quad_P2(3),Quad_P3(3),Quad_P4(3)
    real(kind=FT),intent(out):: Distance,PER(3)
    logical,intent(out):: Yes_PER_in,Yes_PER_on
END SUBROUTINE Tool_Dis_Point_to_3D_Quad
END INTERFACE
end module Global_Inter_Tool_Cal_Dis_Point_to_3D_Quad

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D EBE contact-Jacobian
!        assembler EBE_Cal_Contact_Jacobian_3D.
!
! Notes:   Interface declaration only; uses OpenMP and the 3D
!          contact Gauss-point penalty stiffnesses.
!-----------------------------------------------------------

module Global_EBE_Cal_Contact_Jacobian_3D
INTERFACE
    subroutine EBE_Cal_Contact_Jacobian_3D(isub,i_NR_P, num_freeD, freeDOF,size_local_0, &
    all_local_0,diag_precon_no_invert_0,diag_precon_no_invert, Kn,Kn_Gauss,Kt1_Gauss,Kt2_Gauss,CT_State_Gauss)

    !-----------------------------
    ! Read public variable module
    !-----------------------------
    use Global_Float_Type
    use Global_Common
    use Global_Model
    use Global_XFEM_Elements
    use Global_Crack_Common
    use Global_Crack_3D
    use Global_HF
    use Global_Material
    use OMP_LIB

    !---------------------------
    ! Variable Type Declaration
    !---------------------------
    implicit none
    integer,intent(in)::isub,i_NR_P
    integer,intent(in)::num_FreeD
    integer,intent(in)::freeDOF(num_FreeD)
    real(kind=FT),intent(in)::Kn
    real(kind=FT),intent(in)::diag_precon_no_invert_0(0:num_FreeD)
    integer,intent(in)::size_local_0(Num_Elem)
    integer,intent(in)::all_local_0(MDOF_3D,Num_Elem)
    real(kind=FT),intent(out)::diag_precon_no_invert(0:num_FreeD)
    real(kind=FT),intent(in)::Kt1_Gauss(num_Crack,Max_Max_N_FluEl_3D)
    real(kind=FT),intent(in)::Kt2_Gauss(num_Crack,Max_Max_N_FluEl_3D)
    real(kind=FT),intent(in)::Kn_Gauss(num_Crack,Max_Max_N_FluEl_3D)
    integer,intent(in)::CT_State_Gauss(num_Crack,Max_Max_N_FluEl_3D)
end subroutine EBE_Cal_Contact_Jacobian_3D
END INTERFACE
end module Global_EBE_Cal_Contact_Jacobian_3D

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D contact residual
!        assembler Cal_Contact_Resid_3D.
!
! Notes:   Interface declaration only; receives the 3D crack-Gauss
!          point coordinates and the current displacement vector.
!-----------------------------------------------------------

module Global_Cal_Contact_Resid_3D
INTERFACE
    subroutine Cal_Contact_Resid_3D(iter,ifra,Counter_Iter,i_NR_P, c_Total_Freedom,c_num_freeDOF,F_U, &
    c_DISP,freeDOF,PC_Gauss_x,PC_Gauss_y,PC_Gauss_z,R_PSI)
    use Global_Float_Type
    use Global_Common
    use Global_Model
    use Global_Crack
    use Global_Crack_Common
    use Global_Crack_3D
    use Global_HF
    use Global_Material
    use OMP_LIB

    ! --------------------------
    ! Variable Type Declaration
    ! --------------------------
    implicit none
    integer, intent(in)::iter,ifra,Counter_Iter,i_NR_P
    integer, intent(in)::c_Total_Freedom,c_num_freeDOF
    real(kind=FT),intent(in)::F_U(c_Total_Freedom),c_DISP(c_Total_Freedom)
    integer, intent(in)::freeDOF(c_Total_Freedom)
    real(kind=FT),intent(in)::PC_Gauss_x(num_Crack,Max_Max_N_FluEl_3D) ,PC_Gauss_y(num_Crack,Max_Max_N_FluEl_3D), &
    PC_Gauss_z(num_Crack,Max_Max_N_FluEl_3D)
    real(kind=FT),intent(out)::R_PSI(c_Total_Freedom)
END subroutine Cal_Contact_Resid_3D
END INTERFACE
end module Global_Cal_Contact_Resid_3D

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D per-Gauss-point contact
!        classifier Cal_Contact_Contact_State_Gauss_3D.
!
! Notes:   Interface declaration only; returns the contact-state
!          array and a global 'any contact?' flag.
!-----------------------------------------------------------

module Global_Cal_Contact_Contact_State_Gauss_3D
INTERFACE
    SUBROUTINE Cal_Contact_Contact_State_Gauss_3D( iter,ifra,Counter_Iter,i_NR_P, c_DISP,Yes_Contact,c_Elem_Conta_Sta, &
    CT_State_Gauss)

    !**********************
    ! Read public variable
    !**********************
    use Global_Float_Type
    use Global_Crack
    use Global_Crack_Common
    use Global_Crack_3D
    use Global_Model
    use Global_Common
    use Global_Contact
    use Global_HF
    use Global_Elem_Area_Vol

    !**********************
    ! Variable Declaration
    !**********************
    implicit none
    integer,intent(in)::iter,ifra,Counter_Iter,i_NR_P
    real(kind=FT),intent(in)::c_DISP(Total_FD)
    logical,intent(out)::Yes_Contact
    integer,intent(out)::CT_State_Gauss(num_Crack,Max_Max_N_FluEl_3D)
    integer,intent(out)::c_Elem_Conta_Sta(Num_Elem,num_Crack)
END SUBROUTINE Cal_Contact_Contact_State_Gauss_3D
END INTERFACE
end module Global_Cal_Contact_Contact_State_Gauss_3D

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D normal/tangential contact
!        force updater Cal_Contact_PN_and_PT_3D.
!
! Notes:   Interface declaration only; updates per-Gauss penalty
!          stiffnesses and the contact state.
!-----------------------------------------------------------

module Global_Cal_Contact_PN_and_PT_3D
INTERFACE
    subroutine Cal_Contact_PN_and_PT_3D(iter,ifra,Counter_Iter, i_NR_P,c_Total_Freedom,c_num_freeDOF, &
    Kn,Kn_Gauss,Kt1_Gauss,Kt2_Gauss, fric_mu,c_DISP, delta_u_a, CT_State_Gauss,c_Elem_Conta_Sta, &
    PC_Gauss_x,PC_Gauss_y,PC_Gauss_z,Num_Indent)

    !-----------------------------
    ! Read public variable module
    !-----------------------------
    use Global_Float_Type
    use Global_Common
    use Global_Model
    use Global_Crack_Common
    use Global_Crack_3D
    use Global_HF
    use Global_Material
    use Global_Field_Problem
    use omp_lib
    !---------------------------
    ! Variable Type Declaration
    !---------------------------
    implicit none
    integer, intent(in)::iter,ifra,Counter_Iter,i_NR_P
    integer, intent(in)::c_Total_Freedom,c_num_freeDOF
    real(kind=FT),intent(in)::c_DISP(c_Total_Freedom)
    real(kind=FT),intent(in)::delta_u_a(c_Total_Freedom)
    real(kind=FT),intent(in)::Kn,fric_mu
    real(kind=FT),intent(inout)::Kn_Gauss(num_Crack,Max_Max_N_FluEl_3D)
    real(kind=FT),intent(inout)::Kt1_Gauss(num_Crack,Max_Max_N_FluEl_3D)
    real(kind=FT),intent(inout)::Kt2_Gauss(num_Crack,Max_Max_N_FluEl_3D)
    integer,intent(inout)::CT_State_Gauss(num_Crack,Max_Max_N_FluEl_3D)
    integer,intent(inout):: c_Elem_Conta_Sta(Num_Elem,num_Crack)
    real(kind=FT),intent(inout):: PC_Gauss_x(num_Crack,Max_Max_N_FluEl_3D), PC_Gauss_y(num_Crack,Max_Max_N_FluEl_3D), &
    PC_Gauss_z(num_Crack,Max_Max_N_FluEl_3D)
    integer Num_Indent
END SUBROUTINE Cal_Contact_PN_and_PT_3D
END INTERFACE
end module Global_Cal_Contact_PN_and_PT_3D

!-----------------------------------------------------------
! Brief: Declare the interface for the contact-iteration convergence
!        test Cal_Contact_Conve_Factor.
!
! Notes:   Interface declaration only; returns the convergence flag
!          and the convergence factor.
!-----------------------------------------------------------

module Global_Cal_Contact_Conve_Factor
INTERFACE
    SUBROUTINE Cal_Contact_Conve_Factor( iter,ifra,Counter_Iter,i_NR_P,Conve_Tolerance, &
    c_Total_Freedom,c_freeDOF,c_num_freeDOF, c_F,R_PSI,Last_R_PSI, delta_U,U,Contact_DISP_0, Yes_Conve,Conve_Factor)
    ! Check for convergence and calculate the convergence factor

    !**********************
    ! Read public variable
    !**********************
    use Global_Float_Type
    use Global_Crack
    use Global_Crack_Common
    use Global_Common
    implicit none
    integer,intent(in)::iter,ifra,Counter_Iter,i_NR_P,c_Total_Freedom,c_num_freeDOF
    real(kind=FT),intent(in)::Conve_Tolerance
    integer, intent(in)::c_freeDOF(c_Total_Freedom)
    real(kind=FT),intent(in)::Contact_DISP_0(c_Total_Freedom)
    real(kind=FT),intent(in)::c_F(c_Total_Freedom),R_PSI(c_Total_Freedom),Last_R_PSI(c_Total_Freedom)
    real(kind=FT),intent(in)::delta_U(c_Total_Freedom),U(c_Total_Freedom)
    logical,intent(out)::Yes_Conve
    real(kind=FT),intent(out)::Conve_Factor
END SUBROUTINE Cal_Contact_Conve_Factor
END INTERFACE
end module Global_Cal_Contact_Conve_Factor

!-----------------------------------------------------------
! Brief: Declare the interface for the explicit shape MATMUL_LP
!        that multiplies two 2D real arrays without invoking BLAS.
!
! Notes:   Interface declaration only; used to avoid the implicit-
!          none / assumed-shape issues with MATMUL on some compilers.
!-----------------------------------------------------------

module Function_MATMUL_LP
INTERFACE
    function MATMUL_LP(A,B) result (Output)
    use Global_Float_Type
    implicit none
    real(kind=FT),dimension(:,:),intent(in):: A,B
    real(kind=FT),dimension(size(A,1),size(B,2))::Output
    integer::m,k,n
END function MATMUL_LP
END INTERFACE
end module Function_MATMUL_LP

!-----------------------------------------------------------
! Brief: Declare the interface for the explicit shape matrix-vector
!        product MVMUL_LP.
!
! Notes:   Interface declaration only; used by the EBE-XFEM PCG
!          solver to avoid temporary allocations.
!-----------------------------------------------------------

module Function_MVMUL_LP
INTERFACE
    function MVMUL_LP(A,V) result (Output)
    use Global_Float_Type
    implicit none
    real(kind=FT),dimension(:,:),intent(in):: A
    real(kind=FT),dimension(:),  intent(in):: V
    real(kind=FT),dimension(size(A,1))::Output
    integer::m,n
END function MVMUL_LP
END INTERFACE
end module Function_MVMUL_LP

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D EBE-XFEM PCG solver
!        EBE_XFEM_PCG_3D_with_K that uses a precomputed K.
!
! Notes:   Interface declaration only; OpenMP-parallel, takes a
!          diagonal preconditioner.
!-----------------------------------------------------------

module Global_EBE_XFEM_PCG_3D_with_K
INTERFACE
    SUBROUTINE EBE_XFEM_PCG_3D_with_K(isub,Lambda,c_cg_tol, max_num_PCG,num_FreeD,freeDOF,F,disp, diag_precon_no_invert)
    !Modified on 2022-06-29.
    use Global_Float_Type
    use Global_Model
    use Global_Filename
    use Global_Common
    use Global_Material
    use Global_Crack
    use Global_Crack_Common
    use Global_Crack_3D
    use Function_MVMUL_LP
    use Global_XFEM_Elements
    use omp_lib
    implicit none
    integer,intent(in)::isub,max_num_PCG,num_FreeD
    real(kind=FT),intent(in)::Lambda,c_cg_tol,F(num_FreeD)
    integer,intent(in)::freeDOF(num_FreeD)
    real(kind=FT),intent(in)::diag_precon_no_invert(0:num_FreeD)
    real(kind=FT),intent(out)::disp(Total_FD)
end SUBROUTINE EBE_XFEM_PCG_3D_with_K
END INTERFACE
end module Global_EBE_XFEM_PCG_3D_with_K

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D EBE contact-state driver
!        EBE_Determine_Contact_State_by_Iteration_3D.
!
! Notes:   Interface declaration only; in-place updates the
!          displacement and preconditions.
!-----------------------------------------------------------

module Global_EBE_Determine_Contact_State_by_Iteration_3D
INTERFACE
    SUBROUTINE EBE_Determine_Contact_State_by_Iteration_3D(isub, c_cg_tol,max_num_PCG,num_FreeD,freeDOF, &
    globalF,DISP,diag_precon_no_invert,Num_Indent)
    use Global_Float_Type
    use Global_Common
    use Global_Filename
    use Global_Model
    use Global_XFEM_Elements
    use Global_Elem_Area_Vol
    use Global_Crack
    use Global_Crack_Common
    use Global_Crack_3D
    use Global_HF
    use Global_Contact
    use Global_Inter_Cal_Contact_Red_Resid
    use omp_lib
    implicit none
    integer,intent(in)::isub,max_num_PCG
    integer,intent(in)::num_FreeD
    integer,intent(in)::freeDOF(num_FreeD)
    real(kind=FT),intent(inout)::DISP(Total_FD)
    real(kind=FT),intent(in)::diag_precon_no_invert(0:num_FreeD)
    real(kind=FT),intent(in)::c_cg_tol
    real(kind=FT),intent(in)::globalF(Total_FD)
    integer Num_Indent
end SUBROUTINE EBE_Determine_Contact_State_by_Iteration_3D
END INTERFACE
end module Global_EBE_Determine_Contact_State_by_Iteration_3D

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D constant-pressure-to-
!        force assembler D3_HF_Const_Pres_to_F_Vector.
!
! Notes:   Interface declaration only; used by the 3D HF driver to
!          build the equivalent nodal-force vector for a fluid
!          pressure prescribed on the crack faces.
!-----------------------------------------------------------

module Global_D3_HF_Const_Pres_to_F_Vector
INTERFACE
    SUBROUTINE D3_HF_Const_Pres_to_F_Vector(isub,c_Pres,num_FreeD, in_Total_FD,in_num_Tol_CalP_Water, freeDOF,in_F_U,F)
    use Global_Float_Type
    use Global_Common
    use Global_Model
    use Global_Elem_Area_Vol
    use Global_Crack
    use Global_Crack_Common
    use Global_Crack_3D
    use Global_HF
    use Function_MVMUL_LP
    implicit none
    integer,intent(in)::isub,num_FreeD,in_num_Tol_CalP_Water,in_Total_FD
    real(kind=FT),intent(in)::c_Pres,in_F_U(in_Total_FD)
    integer,intent(in)::freeDOF(num_FreeD)
    real(kind=FT),intent(out)::F(in_Total_FD)
end SUBROUTINE D3_HF_Const_Pres_to_F_Vector
END INTERFACE
end module Global_D3_HF_Const_Pres_to_F_Vector

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D HF NR driver
!        D3_HF_Get_Pres_by_NR_Method that solves for the fluid
!        pressure under a volume / pressure tolerance.
!
! Notes:   Interface declaration only; takes CSR K and the diagonal
!          preconditioner.
!-----------------------------------------------------------

module Global_D3_HF_Get_Pres_by_NR_Method
INTERFACE
    SUBROUTINE D3_HF_Get_Pres_by_NR_Method(i_WB,i_Stage,i_Prop, isub,i_Time,c_Stage_Q, c_Time,Max_Pres_Steps, &
    num_FreeD,c_Total_FD,c_num_Tol_CalP_Water, freeDOF,c_F_U,F,Lambda, diag_precon_no_invert,DISP, &
    NR_delta_Pres,f_Vol_Tol,Pres_Tol,c_Pres,Output_Pres, K_CSR_NNZ,K_CSR_aa,K_CSR_ja,K_CSR_ia)
    !2022-07-01.
    use Global_Float_Type
    use Global_Common
    use Global_Model
    use Global_Elem_Area_Vol
    use Global_Crack
    use Global_Crack_Common
    use Global_Crack_3D
    use Global_HF
    use Global_XFEM_Elements
    implicit none
    integer,intent(in)::i_WB,i_Stage,i_Prop
    integer,intent(in)::isub,i_Time,Max_Pres_Steps
    integer,intent(in)::num_FreeD,c_Total_FD,c_num_Tol_CalP_Water
    real(kind=FT),intent(in)::c_Stage_Q,c_Time
    integer,intent(in)::freeDOF(num_FreeD)
    real(kind=FT),intent(in)::diag_precon_no_invert(0:num_FreeD)
    real(kind=FT),intent(in)::c_F_U(c_Total_FD)
    !2025-12-02.
    integer,intent(in)::K_CSR_NNZ
    real(kind=FT),intent(in)::K_CSR_aa(K_CSR_NNZ)
    integer,intent(in)::K_CSR_ja(K_CSR_NNZ)
    integer,intent(in)::K_CSR_ia(num_FreeD+1)
    real(kind=FT),intent(out)::F(c_Total_FD)
    real(kind=FT),intent(in)::c_Pres,NR_delta_Pres
    real(kind=FT),intent(in)::Lambda
    real(kind=FT),intent(in)::f_Vol_Tol,Pres_Tol
    real(kind=FT),intent(out)::Output_Pres
    real(kind=FT),intent(out)::DISP(Total_FD)
END SUBROUTINE D3_HF_Get_Pres_by_NR_Method
END INTERFACE
end module Global_D3_HF_Get_Pres_by_NR_Method

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D EBE-XFEM helper
!        EBE_XFEM_PCG_3D_Get_K that builds only the element
!        stiffness matrix and the diagonal preconditioner.
!
! Notes:   Interface declaration only; OpenMP-parallel.
!-----------------------------------------------------------

module Global_EBE_XFEM_PCG_3D_Get_K
INTERFACE
    SUBROUTINE EBE_XFEM_PCG_3D_Get_K(isub,num_FreeD,freeDOF,diag_precon_no_invert)
    ! This program only generates the element stiffness matrix.

    use Global_Float_Type
    use Global_Model
    use Global_Filename
    use Global_Common
    use Global_Material
    use Global_Crack
    use Global_Crack_Common
    use Global_Crack_3D
    use Function_MATMUL_LP
    use Global_XFEM_Elements
    use Global_POST
    use OMP_LIB
    implicit none
    integer,intent(in)::isub,num_FreeD
    integer,intent(in)::freeDOF(num_FreeD)
    real(kind=FT),intent(out)::diag_precon_no_invert(0:num_FreeD)
END SUBROUTINE EBE_XFEM_PCG_3D_Get_K
END INTERFACE
end module Global_EBE_XFEM_PCG_3D_Get_K

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D EBE-XFEM PCG driver
!        Ele_by_Ele_XFEM_PCG_3D that recomputes K on the fly.
!
! Notes:   Interface declaration only; same signature as the
!          precomputed-K version but builds K inside the iteration.
!-----------------------------------------------------------

module Global_Ele_by_Ele_XFEM_PCG_3D
INTERFACE
    SUBROUTINE Ele_by_Ele_XFEM_PCG_3D(isub,Lambda,c_cg_tol, max_num_PCG,num_FreeD,freeDOF,F,disp,diag_precon_no_invert)
    use Global_Float_Type
    use Global_Model
    use Global_Filename
    use Global_Common
    use Global_Material
    use Global_Crack
    use Global_Crack_Common
    use Global_Crack_3D
    use Function_MATMUL_LP
    use Global_XFEM_Elements
    use Global_POST
    use omp_lib
    implicit none
    integer,intent(in)::isub,max_num_PCG,num_FreeD
    real(kind=FT),intent(in)::Lambda,c_cg_tol,F(num_FreeD)
    integer,intent(in)::freeDOF(num_FreeD)
    real(kind=FT),intent(out)::disp(Total_FD)
    real(kind=FT),intent(out)::diag_precon_no_invert(0:num_FreeD)
end SUBROUTINE Ele_by_Ele_XFEM_PCG_3D
END INTERFACE
end module Global_Ele_by_Ele_XFEM_PCG_3D

!-----------------------------------------------------------
! Brief: Declare the interface for the 8-node hexahedral element
!        stiffness assembler Cal_Ele_Stiffness_Matrix_3D_8nodes.
!
! Notes:   Interface declaration only; takes the user-supplied
!          Gauss points and the 6x6 constitutive matrix.
!-----------------------------------------------------------

module Global_Cal_Ele_Stiffness_Matrix_3D_8nodes
INTERFACE
    subroutine Cal_Ele_Stiffness_Matrix_3D_8nodes(i_E,num_Gauss, X_NODES,Y_NODES,Z_NODES, c_D,kesi,yita,zeta,weight,localK)
    use Global_Float_Type
    implicit none
    integer,intent(in)::i_E,num_Gauss
    real(kind=FT),intent(in)::c_D(6,6)
    real(kind=FT),intent(in)::kesi(num_Gauss),yita(num_Gauss),zeta(num_Gauss),weight(num_Gauss)
    real(kind=FT),intent(in)::X_NODES(8),Y_NODES(8),Z_NODES(8)
    real(kind=FT),intent(out)::localK(24,24)

    real(kind=FT) J(3,3),detJ,dNdkesi(8,3),dNdx(8,3),B_FEM(6,24)
    real(kind=FT) Inverse_J(3,3)
    integer i_G,i_N
    real(kind=FT) tem
    real(kind=FT) temp(3,8),Coor(3,8)
    real(kind=FT) one_p_yita,one_m_yita
    real(kind=FT) one_p_zeta,one_m_zeta
    real(kind=FT) one_p_kesi,one_m_kesi
end subroutine Cal_Ele_Stiffness_Matrix_3D_8nodes
END INTERFACE
end module Global_Cal_Ele_Stiffness_Matrix_3D_8nodes

!-----------------------------------------------------------
! Brief: Declare the interface for Cal_Ele_Num_by_Coors_3D which
!        locates the element containing a given (x, y, z) point.
!
! Notes:   Interface declaration only; uses a caller-supplied
!          element-number cache for fast repeated queries.
!-----------------------------------------------------------

module Global_Cal_Ele_Num_by_Coors_3D
INTERFACE
    subroutine Cal_Ele_Num_by_Coors_3D(x,y,z,Ele_Num_Cache,OUT_Elem)
    use Global_Float_Type
    use Global_Model
    implicit none
    real(kind=FT),intent(in):: x,y,z
    integer,intent(out)::OUT_Elem
    integer,intent(inout)::Ele_Num_Cache
end subroutine Cal_Ele_Num_by_Coors_3D
END INTERFACE
end module Global_Cal_Ele_Num_by_Coors_3D

!------------------------------------------------------
! 29. SUBROUTINE D3_Run_a_Fracturing_Step, 2023-03-14.
!------------------------------------------------------
#ifndef github

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D fracturing-step driver
!        D3_Run_a_Fracturing_Step that advances the HF state by
!        one propagation step.
!
! Notes:   Interface declaration only; takes the stage index and
!          a growth flag, dispatches to PCG / contact / HF.
!-----------------------------------------------------------

module Global_D3_Run_a_Fracturing_Step
INTERFACE
    subroutine D3_Run_a_Fracturing_Step(i_Stage,i_Growth)
    use Global_Float_Type
    use Global_Common
    use Global_Filename
    use Global_Model
    use Global_XFEM_Elements
    use Global_Elem_Area_Vol
    use Global_Crack_Common
    use Global_Crack_3D
    use Global_DISP
    use Global_HF
    use Global_Stress
    use Global_Strain
    use Global_POST
    use Global_Ragged_Array_Real_Classs

    ! Read subroutine interface module (activate compiler parameter consistency check)
    use Global_Inter_Matrix_Solve_LSOE
    use Global_Ele_by_Ele_XFEM_PCG_3D
    use Global_EBE_Determine_Contact_State_by_Iteration_3D

    ! Variable Type Declaration
    implicit none
    integer,intent(in)::i_Stage,i_Growth
end subroutine D3_Run_a_Fracturing_Step
END INTERFACE
end module Global_D3_Run_a_Fracturing_Step
#endif

!-----------------------------------------------------------
! Brief: Declare the interface for the point-in-3D-hexahedron test
!        Tool_Yes_Point_in_3D_Hexahedron.
!
! Notes:   Interface declaration only; returns the inside / on
!-----------------------------------------------------------

module Global_INTERFACE_Tool_Yes_Point_in_3D_Hexahedron
INTERFACE
    subroutine Tool_Yes_Point_in_3D_Hexahedron(Point,A,B,C,D,E,F,G,H,Yes_in,Yes_on)

    use Global_Float_Type
    implicit none
    real(kind=FT),intent(in)::Point(3),A(3),B(3),C(3),D(3),E(3),F(3),G(3),H(3)
    logical,intent(out):: Yes_in,Yes_on

end SUBROUTINE Tool_Yes_Point_in_3D_Hexahedron
END INTERFACE
end module Global_INTERFACE_Tool_Yes_Point_in_3D_Hexahedron

!-----------------------------------------------------------
! Brief: Declare the interface for the in-place n-vector normalizer
!        Vector_Normalize.
!
! Notes:   Interface declaration only; modifies the input vector.
!-----------------------------------------------------------

module Global_INTERFACE_Vector_Normalize
INTERFACE
    SUBROUTINE Vector_Normalize(n,Vector)

    use Global_Float_Type
    use Global_Common

    implicit none
    integer,intent(in):: n
    real(kind=FT),intent(inout):: Vector(n)

END SUBROUTINE Vector_Normalize
END INTERFACE
end module Global_INTERFACE_Vector_Normalize

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D rotation estimator
!        Tool_ThetaX_ThetaY_ThetaZ_3D_rotation that returns the
!        rotation matrix and Euler angles between two triads.
!
! Notes:   Interface declaration only; used by the 3D crack-tip
!          baseline construction.
!-----------------------------------------------------------

module Global_INTERFACE_Tool_ThetaX_ThetaY_ThetaZ_3D_rotation
INTERFACE
    subroutine Tool_ThetaX_ThetaY_ThetaZ_3D_rotation(i_V,Vector_Origi_1, &
    Vector_Origi_2,Vector_Origi_3,Vector_Crack_1,Vector_Crack_2, Vector_Crack_3,ThetaX,ThetaY,ThetaZ,T_Matrix)


    use Global_Float_Type
    use Global_Crack_3D

    implicit none
    integer,intent(in)::i_V
    real(kind=FT),intent(in):: Vector_Origi_1(3),Vector_Origi_2(3),Vector_Origi_3(3), &
    Vector_Crack_1(3),Vector_Crack_2(3),Vector_Crack_3(3)
    real(kind=FT),intent(out):: ThetaX,ThetaY,ThetaZ,T_Matrix(3,3)
end SUBROUTINE Tool_ThetaX_ThetaY_ThetaZ_3D_rotation
END INTERFACE
end module Global_INTERFACE_Tool_ThetaX_ThetaY_ThetaZ_3D_rotation

!-----------------------------------------------------------
! Brief: Declare the interface for the point-in-3D-hexahedron test
!        Tool_Yes_Point_in_3D_Hexahedron_with_Tol with a user-
!        supplied tolerance.
!
! Notes:   Interface declaration only; same signature as the
!          tolerance-less version plus a Tol scalar.
!-----------------------------------------------------------

module Global_INTERFACE_Tool_Yes_Point_in_3D_Hexahedron_with_Tol
INTERFACE
    subroutine Tool_Yes_Point_in_3D_Hexahedron_with_Tol(Point,A,B,C,D,E,F,G,H,Tol,Yes_in,Yes_on)

    use Global_Float_Type
    implicit none
    real(kind=FT),intent(in)::Point(3),A(3),B(3),C(3),D(3),E(3),F(3),G(3),H(3),Tol
    logical,intent(out):: Yes_in,Yes_on

end SUBROUTINE Tool_Yes_Point_in_3D_Hexahedron_with_Tol
END INTERFACE
end module Global_INTERFACE_Tool_Yes_Point_in_3D_Hexahedron_with_Tol

!-----------------------------------------------------------
! Brief: Declare the interface for Vector_Location_Int_v2 which
!        returns the index of a given integer in a vector.
!
! Notes:   Interface declaration only; outputs zero on miss.
!-----------------------------------------------------------

module Global_INTERFACE_Vector_Location_Int_v2
INTERFACE
    SUBROUTINE Vector_Location_Int_v2(n,Vector,Variable,location)

    use Global_Float_Type
    implicit none
    integer,intent(in)::n,Vector(n),Variable
    integer,intent(out)::location

end SUBROUTINE Vector_Location_Int_v2
END INTERFACE
end module Global_INTERFACE_Vector_Location_Int_v2

!-----------------------------------------------------------
! Brief: Declare the interface for D3_Get_Signed_Dis_to_Crack_Mesh
!        which returns the signed distance from a point to a 3D
!        crack mesh plus the local perpendicular data.
!
! Notes:   Interface declaration only; also returns the closest
!          point and the crack-side normal vector.
!-----------------------------------------------------------

module Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh
INTERFACE
    SUBROUTINE D3_Get_Signed_Dis_to_Crack_Mesh(Point,i_C,Check_Ball_R, Signed_Dis,Signed_Dis_v2,c_Yes_Node_PER_in_FS, &
    c_PER_Node_to_FS,Yes_Found_Min_Signed_Dis, n_Vector)
    use Global_Float_Type
    use Global_Crack_Common
    use Global_Crack_3D

    implicit none
    integer,intent(in)::i_C
    real(kind=FT),intent(in)::Point(3)
    real(kind=FT),intent(in)::Check_Ball_R
    real(kind=FT),intent(out)::Signed_Dis,c_PER_Node_to_FS(3)
    real(kind=FT),intent(out)::Signed_Dis_v2
    real(kind=FT),intent(out)::n_Vector(3)
    logical,intent(out)::Yes_Found_Min_Signed_Dis
    logical,intent(out)::c_Yes_Node_PER_in_FS

end SUBROUTINE D3_Get_Signed_Dis_to_Crack_Mesh
END INTERFACE
end module Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh

!-----------------------------------------------------------
! Brief: Declare the interface for the in-plane-growth variant of
!        D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth.
!
! Notes:   Interface declaration only; same shape as the generic
!          version but restricted to the in-plane growth rule.
!-----------------------------------------------------------

module Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane
INTERFACE
    SUBROUTINE D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth(Point,i_C, Check_Ball_R,Signed_Dis, &
    c_Yes_Node_PER_in_FS,c_PER_Node_to_FS, Yes_Found_Min_Signed_Dis,n_Vector)
    use Global_Float_Type
    use Global_Crack_Common
    use Global_Crack_3D

    implicit none
    integer,intent(in)::i_C
    real(kind=FT),intent(in)::Point(3)
    real(kind=FT),intent(in)::Check_Ball_R
    real(kind=FT),intent(out)::Signed_Dis,c_PER_Node_to_FS(3)
    real(kind=FT),intent(out)::n_Vector(3)
    logical,intent(out)::c_Yes_Node_PER_in_FS,Yes_Found_Min_Signed_Dis

end SUBROUTINE D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth
END INTERFACE
end module Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane

!-----------------------------------------------------------
! Brief: Declare the interface for the complex-arithmetic dense
!        linear solver Matrix_Solve_LSOE_Complex.
!
! Notes:   Interface declaration only; used by the complex-kernel
!          eigenvalue / frequency-response routines.
!-----------------------------------------------------------

module module_INTERFACE_Matrix_Solve_LSOE_Complex
INTERFACE
    SUBROUTINE Matrix_Solve_LSOE_Complex(Key_Indent,Key_LSOE_Sys,c_Key_SLOE,K,F,D,n)
    ! D = F / K, where K is stored in the form of a full matrix
    use Global_Float_Type
    use Global_Common
    use Global_ITPACK
    use, intrinsic :: ISO_C_BINDING
#ifdef sffortran
#ifndef macos
#ifndef github
    use strumpack
#endif
#endif
#endif  
    implicit none
    integer,intent(in)::n,c_Key_SLOE,Key_LSOE_Sys
    integer,intent(in)::Key_Indent
    Complex*16,intent(in)::K(n,n),F(n)
    Complex*16,intent(out)::D(n)
    integer i,j,NZ_NUM,IERR
    Complex*16 ZR_Complex
end subroutine Matrix_Solve_LSOE_Complex

END INTERFACE
end module module_INTERFACE_Matrix_Solve_LSOE_Complex

!-----------------------------------------------------------
! Brief: Declare the interface for the complex-arithmetic CSR
!        linear solver Matrix_Solve_LSOE_Complex_Sparse.
!
! Notes:   Interface declaration only; takes both real and complex
!          matrix storage representations.
!-----------------------------------------------------------

module module_INTERFACE_Matrix_Solve_LSOE_Complex_Sparse
INTERFACE
    SUBROUTINE Matrix_Solve_LSOE_Complex_Sparse(Key_Indent,Key_LSOE_Sys,c_Key_SLOE,K_CSR_NNZ, K_CSR_aa,K_CSR_ja,K_CSR_ia, &
    Complex_K_vals,Complex_K_rows,Complex_K_cols, F,D,n)
    ! D = F / K, where K is stored in the form of a full matrix
    use Global_Float_Type
    use Global_Common
    use Global_ITPACK
    use, intrinsic :: ISO_C_BINDING
#ifdef sffortran
#ifndef macos
#ifndef github
    use strumpack
#endif
#endif
#endif  
    implicit none
    integer,intent(in)::n,c_Key_SLOE,Key_LSOE_Sys
    integer,intent(in)::K_CSR_NNZ
    integer,intent(in)::Key_Indent
    Complex*16,intent(in)::K_CSR_aa(K_CSR_NNZ)
    integer,intent(in)::K_CSR_ja(K_CSR_NNZ)
    integer,intent(in)::K_CSR_ia(n+1)
    Complex*16 ,intent(in)::Complex_K_vals(K_CSR_NNZ)
    integer,intent(in)::Complex_K_rows(K_CSR_NNZ)  
    integer,intent(in)::Complex_K_cols(K_CSR_NNZ)
    Complex*16,intent(in)::F(n)
    Complex*16,intent(out)::D(n)
end subroutine Matrix_Solve_LSOE_Complex_Sparse
END INTERFACE
end module module_INTERFACE_Matrix_Solve_LSOE_Complex_Sparse

!-----------------------------------------------------------
! Brief: Declare the interface for the complex sparse matrix-vector
!        multiply SPARSKIT2_amux_complex.
!
! Notes:   Interface declaration only; 8-byte integer index and
!          16-byte complex arithmetic.
!-----------------------------------------------------------

module module_INTERFACE_SPARSKIT2_amux_complex
INTERFACE
    subroutine SPARSKIT2_amux_complex(n, x, y, nnz,a,ja,ia) 
    implicit none
    integer(8),intent(in)::nnz
    integer,intent(in)::n
    complex*16 ,intent(in)::x(n),a(nnz) 
    integer,intent(in)::ja(nnz)
    integer(kind=8),intent(in)::ia(n+1)
    complex*16 ,intent(inout)::y(n)
end subroutine SPARSKIT2_amux_complex
END INTERFACE
end module module_INTERFACE_SPARSKIT2_amux_complex

!-----------------------------------------------------------
! Brief: Declare the interface for Vector_belongs_Matrix_Is_Dou
!        which tests if a vector equals any column of a real
!        matrix (within double-precision tolerance).
!
! Notes:   Interface declaration only; returns the column index and
!          a Yes/No flag.
!-----------------------------------------------------------

module module_INTERFACE_Vector_belongs_Matrix_Is_Dou
INTERFACE
    subroutine Vector_belongs_Matrix_Is_Dou(m,n,Matrix,Vector,Location,Yes)   
    use Global_Float_Type
    implicit none
    integer,intent(in)::m,n
    real(kind=FT),intent(in)::Matrix(m,n)
    real(kind=FT),intent(in)::Vector(n)
    integer,intent(out)::Location
    logical,intent(out)::Yes
end subroutine Vector_belongs_Matrix_Is_Dou
END INTERFACE
end module module_INTERFACE_Vector_belongs_Matrix_Is_Dou

!-----------------------------------------------------------
! Brief: Declare the interface for the recursive integer 2-column
!        quick-sort routine Matrix_n_x_2_Quick_Sort_Int.
!
! Notes:   Interface declaration only; sorts by the first column
!          and tracks the original index.
!-----------------------------------------------------------

module module_INTERFACE_Matrix_n_x_2_Quick_Sort_Int
INTERFACE
    recursive subroutine Matrix_n_x_2_Quick_Sort_Int(Matrix, Index, low, high)
    implicit none
    integer, intent(inout) :: Matrix(:, :)
    integer, intent(inout) :: Index(:)
    integer :: low, high, i, j, pivot, tempRow(2), tempIndex
end subroutine Matrix_n_x_2_Quick_Sort_Int
END INTERFACE
end module module_INTERFACE_Matrix_n_x_2_Quick_Sort_Int

!-----------------------------------------------------------
! Brief: Declare the interface for Matrix_n_x_3_Check_Sort which
!        verifies that an n x 3 real matrix is correctly sorted.
!
! Notes:   Interface declaration only; returns a logical flag and
!          an error count.
!-----------------------------------------------------------

module module_INTERFACE_Matrix_n_x_3_Check_Sort
INTERFACE
    subroutine Matrix_n_x_3_Check_Sort(Matrix, n, check_matrix, count_error)
    use Global_Float_Type
    implicit none
    real(kind=FT), intent(in) :: Matrix(n,3)
    integer, intent(in) :: n
    logical, intent(out) :: check_matrix
    integer, intent(out) :: count_error
end subroutine Matrix_n_x_3_Check_Sort
END INTERFACE
end module module_INTERFACE_Matrix_n_x_3_Check_Sort

!-----------------------------------------------------------
! Brief: Declare the interface for Matrix_n_x_3_Check_Sort_Int,
!        the integer version of the sort-checker.
!
! Notes:   Interface declaration only; same signature as the real
!          variant but with integer entries.
!-----------------------------------------------------------

module module_INTERFACE_Matrix_n_x_3_Check_Sort_Int
INTERFACE
    subroutine Matrix_n_x_3_Check_Sort_Int(Matrix, n, check_matrix, count_error)
    implicit none
    integer, intent(in) :: Matrix(n,3)
    integer, intent(in) :: n
    logical, intent(out) :: check_matrix
    integer, intent(out):: count_error
end subroutine Matrix_n_x_3_Check_Sort_Int
END INTERFACE
end module module_INTERFACE_Matrix_n_x_3_Check_Sort_Int

!-----------------------------------------------------------
! Brief: Declare the interface for the SPARSKIT coocsr utility
!        that converts COO storage to CSR storage.
!
! Notes:   Interface declaration only; real*8 values, integer row
!          and column indices.
!-----------------------------------------------------------

module module_INTERFACE_coocsr
INTERFACE
    subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
    real*8 a(*),ao(*),x
    integer ir(*),jc(*),jao(*),iao(*)
END subRoutine coocsr
END INTERFACE
end module module_INTERFACE_coocsr

!-----------------------------------------------------------
! Brief: Declare the interface for Cal_JDomain_Element_of_Point_3D
!        _Simple_HEX8 which collects the J-domain elements around
!        a crack-tip point for 3D IIM integration.
!
! Notes:   Interface declaration only; returns the element list,
!          count, and the corresponding eight corner coordinates.
!-----------------------------------------------------------

module module_Cal_JDomain_Element_of_Point_3D_Simple_HEX8
INTERFACE
    subroutine Cal_JDomain_Element_of_Point_3D_Simple_HEX8(TipCoor, rJ, J_Elem,Max_J_Elem, num_J_Elem, Q_Elem_Nodes)
    use Global_Float_Type
    use Global_Model
    implicit none
    integer,intent(in) :: Max_J_Elem
    real(kind=FT),intent(in) :: TipCoor(3), rJ
    integer,intent(out) :: J_Elem(Max_J_Elem), num_J_Elem
    real(kind=FT),intent(out) :: Q_Elem_Nodes(num_Elem,8)
end subroutine Cal_JDomain_Element_of_Point_3D_Simple_HEX8
END INTERFACE
end module module_Cal_JDomain_Element_of_Point_3D_Simple_HEX8

!-----------------------------------------------------------
! Brief: Declare the interface for the linear-leak-off HF matrix
!        assembler Cal_HF_Matrix_H_Linear.
!
! Notes:   Interface declaration only; takes the last-step crack
!          opening and the current time.
!-----------------------------------------------------------

module module_Cal_HF_Matrix_H_Linear
INTERFACE
    subroutine Cal_HF_Matrix_H_Linear(ifra,Counter,Matrix_H,Last_Cr_CalP_Aper,total_time)
    use Global_Float_Type
    use Global_Crack
    use Global_Crack_Common
    use Global_HF
    implicit none

    integer,intent(in)::ifra,Counter
    real(kind=FT),intent(in)::total_time
    real(kind=FT),intent(in)::Last_Cr_CalP_Aper(Max_Num_Cr,Max_Num_Cr_CalP)
    real(kind=FT),intent(out)::Matrix_H(num_Tol_CalP_Water,num_Tol_CalP_Water)
end subroutine Cal_HF_Matrix_H_Linear
END INTERFACE
end module module_Cal_HF_Matrix_H_Linear

!-----------------------------------------------------------
! Brief: Declare the interface for the second half of the 3D IIM
!        auxiliary-field evaluator Cal_3D_SIFs_IIM_Auxiliary_
!        Fields_Part2 (returns the spatial derivatives).
!
! Notes:   Interface declaration only; mode I/II/III asymptotic
!          fields, driven by (r, theta, G, nu, kappa).
!-----------------------------------------------------------

module module_Cal_3D_SIFs_IIM_Auxiliary_Fields_Part2
INTERFACE
    subroutine Cal_3D_SIFs_IIM_Auxiliary_Fields_Part2(Mode, r_in, theta, G, nu, kappa, &
    Aux_dSig11_dX1, Aux_dSig12_dX2, Aux_dSig13_dX3, Aux_dSig12_dX1, Aux_dSig22_dX2, Aux_dSig23_dX3, &
    Aux_dSig13_dX1, Aux_dSig23_dX2, Aux_dSig33_dX3, Aux_dEps11_dX1, Aux_dEps12_dX1, Aux_dEps13_dX1, &
    Aux_dEps22_dX1, Aux_dEps23_dX1, Aux_dEps33_dX1, Aux_ddu1_dX1dX1, Aux_ddu2_dX1dX1, Aux_ddu3_dX1dX1, &
    Aux_ddu1_dX2dX1, Aux_ddu2_dX2dX1, Aux_ddu3_dX2dX1, Aux_ddu1_dX3dX1, Aux_ddu2_dX3dX1, Aux_ddu3_dX3dX1)
    !-------------------------------------------------------------------------------
    use Global_Float_Type
    integer, intent(in) :: Mode
    real(kind=FT), intent(in) :: r_in, theta, G, nu, kappa
    real(kind=FT), intent(out) :: Aux_dSig11_dX1, Aux_dSig12_dX2, Aux_dSig13_dX3
    real(kind=FT), intent(out) :: Aux_dSig12_dX1, Aux_dSig22_dX2, Aux_dSig23_dX3
    real(kind=FT), intent(out) :: Aux_dSig13_dX1, Aux_dSig23_dX2, Aux_dSig33_dX3
    real(kind=FT), intent(out) :: Aux_dEps11_dX1, Aux_dEps12_dX1, Aux_dEps13_dX1
    real(kind=FT), intent(out) :: Aux_dEps22_dX1, Aux_dEps23_dX1, Aux_dEps33_dX1
    real(kind=FT), intent(out) :: Aux_ddu1_dX1dX1, Aux_ddu2_dX1dX1, Aux_ddu3_dX1dX1
    real(kind=FT), intent(out) :: Aux_ddu1_dX2dX1, Aux_ddu2_dX2dX1, Aux_ddu3_dX2dX1 
    real(kind=FT), intent(out) :: Aux_ddu1_dX3dX1, Aux_ddu2_dX3dX1, Aux_ddu3_dX3dX1
end subroutine Cal_3D_SIFs_IIM_Auxiliary_Fields_Part2
END INTERFACE
end module module_Cal_3D_SIFs_IIM_Auxiliary_Fields_Part2

!-----------------------------------------------------------
! Brief: Declare the interface for the first half of the 3D IIM
!        auxiliary-field evaluator Cal_3D_SIFs_IIM_Auxiliary_
!        Fields_Part1 (returns the asymptotic stress, strain, and
!        displacement gradient).
!
! Notes:   Interface declaration only; same parameter list as Part2
!          but with field arrays as outputs.
!-----------------------------------------------------------

module module_Cal_3D_SIFs_IIM_Auxiliary_Fields_Part1
INTERFACE
    subroutine Cal_3D_SIFs_IIM_Auxiliary_Fields_Part1(Mode, r, theta, G, nu, kappa, Aux_Sig, Aux_Eps, Aux_Grad_u)
    use Global_Float_Type
    integer, intent(in) :: Mode
    real(kind=FT), intent(in) :: r, theta, G, nu, kappa
    real(kind=FT), intent(out) :: Aux_Sig(3,3), Aux_Eps(3,3), Aux_Grad_u(3,3)
end subroutine Cal_3D_SIFs_IIM_Auxiliary_Fields_Part1
END INTERFACE
end module module_Cal_3D_SIFs_IIM_Auxiliary_Fields_Part1

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D IIM auxiliary q-function
!        and its gradient at a Gauss point.
!
! Notes:   Interface declaration only; takes the element ID, nodal
!          q values, the crack-tip frame, and the radial bounds.
!-----------------------------------------------------------

module module_Cal_3D_SIFs_IIM_q_Function
INTERFACE
    subroutine Cal_3D_SIFs_IIM_q_Function(Elem_ID, q_nodal, N, dN_dx, dN_dy, dN_dz, GP_coords, Tip_Coords, e1, e2, e3, &
    R_in, R_out, L_f, q_val, grad_q)
    use Global_Float_Type
    implicit none
    integer, intent(in) :: Elem_ID
    real(kind=FT), intent(in) :: q_nodal(:,:), N(8), dN_dx(8), dN_dy(8), dN_dz(8)
    real(kind=FT), intent(in) :: GP_coords(3), Tip_Coords(3)
    real(kind=FT), intent(in) :: e1(3), e2(3), e3(3)
    real(kind=FT), intent(in) :: R_in, R_out, L_f
    real(kind=FT), intent(out) :: q_val, grad_q(3)
end subroutine Cal_3D_SIFs_IIM_q_Function
END INTERFACE
end module module_Cal_3D_SIFs_IIM_q_Function

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D IIM J-domain element
!        collector Cal_3D_SIFs_IIM_Integration_Domain.
!
! Notes:   Interface declaration only; returns the element list and
!          the q-function value at each element centroid.
!-----------------------------------------------------------

module module_Cal_3D_SIFs_IIM_Integration_Domain
INTERFACE
    subroutine Cal_3D_SIFs_IIM_Integration_Domain(Tip_Point, e1, e2, e3, R_in, R_out, L_f, Num_Elems, Elem_List, q_vals)
    use Global_Float_Type
    use Global_Model

    implicit none
    real(kind=FT), intent(in) :: Tip_Point(3), e1(3), e2(3), e3(3)
    real(kind=FT), intent(in) :: R_in, R_out, L_f
    integer, intent(out) :: Num_Elems
    integer, allocatable, intent(out) :: Elem_List(:)
    real(kind=FT), allocatable, intent(out) :: q_vals(:,:)
end subroutine Cal_3D_SIFs_IIM_Integration_Domain
END INTERFACE
end module module_Cal_3D_SIFs_IIM_Integration_Domain

!-----------------------------------------------------------
! Brief: Declare the interface for the 3-point crack-face Gauss
!        integration-point generator used by the 3D SIF IIM.
!
! Notes:   Interface declaration only; returns point coordinates,
!          weights, and outward normals.
!-----------------------------------------------------------

module module_Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_3_Points
INTERFACE
    subroutine Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_3_Points(i_Crack, Elem_ID, nGP_face, GP_global, GP_w, GP_n_global)
    !Compute integral points of the fluid element at the crack surface.
    !2026-01-28.
    use Global_Float_Type
    use Global_Common
    use Global_Crack_3D
    use Global_Elem_Area_Vol
    implicit none
    integer, intent(in) :: i_Crack, Elem_ID
    integer, intent(out) :: nGP_face
    real(kind=FT), allocatable, intent(out) :: GP_global(:,:), GP_w(:), GP_n_global(:,:)
end subroutine Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_3_Points
END INTERFACE
end module module_Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_3_Points

!-----------------------------------------------------------
! Brief: Declare the interface for the 1-point crack-face Gauss
!        integration-point generator used by the 3D SIF IIM.
!
! Notes:   Interface declaration only; faster but lower order.
!-----------------------------------------------------------

module module_Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_1_Point
INTERFACE
    subroutine Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_1_Point(i_Crack, Elem_ID, nGP_face, GP_global, GP_w, GP_n_global)
    use Global_Float_Type
    use Global_Common
    use Global_Crack_3D
    use Global_Elem_Area_Vol

    implicit none

    integer, intent(in) :: i_Crack, Elem_ID
    integer, intent(out) :: nGP_face
    real(kind=FT), allocatable, intent(out) :: GP_global(:,:), GP_w(:), GP_n_global(:,:)
end subroutine Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_1_Point
END INTERFACE
end module module_Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_1_Point

!-----------------------------------------------------------
! Brief: Declare the (currently commented-out) interface for the
!        7-point crack-face Gauss integration-point generator.
!
! Notes:   Interface declaration only; reserved for higher-order
!          IIM crack-face integrals.
!-----------------------------------------------------------

module module_Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_7_Points
INTERFACE
    !    subroutine Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_7_Points(i_Crack, Elem_ID, nGP_face, GP_global, GP_w, GP_n_global)
    !    implicit none
    !    end subroutine Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_7_Points
END INTERFACE
end module module_Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_7_Points

!-----------------------------------------------------------
! Brief: Declare the interface for the 3D crack-opening calculator
!        Cal_Crack_Point_Aperture_3D used by 3D HF post-processing.
!
! Notes:   Interface declaration only; interpolates the displacement
!          field onto the fluid-element / crack-mesh pair.
!-----------------------------------------------------------

module module_Cal_Crack_Point_Aperture_3D
INTERFACE
    subroutine Cal_Crack_Point_Aperture_3D(c_DISP,Crack_Number,Crack_Point,Crack_Point_Aperture, &
    Fluid_element_number,Fluid_node_number, Crack_mesh_ele_number,Crack_mesh_node_number)
    use Global_Float_Type
    use Global_Common
    use Global_Crack_Common
    use Global_Crack_3D
    use Global_Model
    use Global_Elem_Area_Vol
    use Global_DISP
    use Global_Inter_Tool_Cal_Dis_Point_to_3D_Quad
    use Global_INTERFACE_Tool_ThetaX_ThetaY_ThetaZ_3D_rotation
    implicit none
    real(kind=FT),intent(in)::c_DISP(Total_FD)
    integer,intent(in)::Crack_Number
    real(kind=FT),intent(in)::Crack_Point(3)
    real(kind=FT),intent(out)::Crack_Point_Aperture(3)
    integer,intent(in)::Fluid_element_number,Fluid_node_number,Crack_mesh_ele_number,Crack_mesh_node_number
end subroutine Cal_Crack_Point_Aperture_3D
END INTERFACE
end module module_Cal_Crack_Point_Aperture_3D

!-----------------------------------------------------------
! Brief: Declare the interface for the in-place symmetric matrix
!        diagonalizer Matrix_Diagonalization.
!
! Notes:   Interface declaration only; overwrites the input matrix
!          with its eigen-decomposition (assumed symmetric).
!-----------------------------------------------------------

module module_Matrix_Diagonalization
INTERFACE
    SUBROUTINE Matrix_Diagonalization(n,Matrix)   
    use Global_Float_Type
    implicit none
    integer,       intent(in)    :: n
    real(kind=FT), intent(inout) :: Matrix(:,:)
end subroutine Matrix_Diagonalization
END INTERFACE
end module module_Matrix_Diagonalization

!-----------------------------------------------------------
! Brief: Declare the interface for Local_L2_projection_for_inter_
!        element that re-uses the same signature as the routine
!        of the same name in Src_Main.
!
! Notes:   Interface declaration only; allows the same procedure
!          to be referenced from the dynamic-analysis modules.
!-----------------------------------------------------------

module module_INTERFACE_Local_L2_projection_for_inter_element
INTERFACE
    subroutine Local_L2_projection_for_inter_element( isub, crack_number,num_gp, c_DISP, number_nodes_projected, &
    nodes_to_be_projected, nodes_new_enriched_value)
    ! Local L2 projection for inter-element crack-tip crossing.
    ! Map crack-tip enriched nodal values to Heaviside enriched nodal values
    ! on the specified element.
    ! Assumption: only one tip enrichment function F is used.
    use Global_Float_Type
    use Global_Common
    use Global_Filename
    use Global_Model
    use Global_XFEM_Elements
    use Global_Elem_Area_Vol
    use Global_Crack_Common
    use Global_Crack_3D
    use Global_HF
    use Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh
    use Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane
    use Global_Inter_Cal_N_dNdkesi_J_detJ_3D
    implicit none
    integer, intent(in) :: isub, crack_number, num_gp
    integer, intent(in) :: number_nodes_projected
    integer, intent(in) :: nodes_to_be_projected(number_nodes_projected)
    real(kind=FT), intent(in) :: c_DISP(Total_FD)
    real(kind=FT), intent(out) :: nodes_new_enriched_value(number_nodes_projected*3)
end subroutine Local_L2_projection_for_inter_element
END INTERFACE
end module module_INTERFACE_Local_L2_projection_for_inter_element


!2026-07-11. BUGFIX-2026071111.
module module_INTERFACE_OpenBLAS
    use, intrinsic :: ISO_C_BINDING
    implicit none
    INTERFACE
        SUBROUTINE openblas_set_num_threads(num_threads) &
                   BIND(C, NAME='openblas_set_num_threads')
            import :: C_INT        
            INTEGER(C_INT), VALUE :: num_threads
        END SUBROUTINE openblas_set_num_threads
        FUNCTION openblas_get_num_threads() &
                 BIND(C, NAME='openblas_get_num_threads') RESULT(res)
            import :: C_INT
            INTEGER(C_INT) :: res
        END FUNCTION openblas_get_num_threads
    END INTERFACE
end module module_INTERFACE_OpenBLAS

