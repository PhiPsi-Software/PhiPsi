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
 
module Global_Inter_Cal_HF_Initial_Pressure
INTERFACE
  subroutine Cal_HF_Initial_Pressure(iter,ifra,Counter, &
                                    Yes_Growth,&
                                     Last_Cracks_CalP_Num_ifra,&
                                     Last_Cr_CalP_Pres_ifra,&
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

!----------------
!2. Cal_HF_Resid
!----------------
module Global_Inter_Cal_HF_Resid
INTERFACE
  subroutine Cal_HF_Resid(ifra,iter,Counter_Iter,R,Total_FD,&
                  num_FreeD,num_free_CalP,&
                  globalK,Coupled_Q,F_U,&
                  Last_DISP,Last_Last_DISP,&
                  freeDOF,freeDOF_HF,Local_freeDOF_HF,&
                  Last_CalP_Pres,delta_Time,H,c_S)
  use Global_Float_Type
  use Global_Crack
  use Global_Crack_Common
  use Global_HF
  use Global_Material
  implicit none
  integer, intent(in)::ifra,iter,Counter_Iter,Total_FD,num_FreeD,num_free_CalP
  real(kind=FT),intent(in):: &
                Coupled_Q(Total_FD,num_Tol_CalP_Water),&
                globalK(Total_FD,Total_FD),&
                F_U(Total_FD),Last_DISP(Total_FD),&
                Last_Last_DISP(Total_FD),&
                Last_CalP_Pres(num_Tol_CalP_Water),delta_Time,&
                H(num_Tol_CalP_Water,num_Tol_CalP_Water),&
                c_S(num_Tol_CalP_Water)
  integer, intent(in)::freeDOF(Total_FD),&
                       freeDOF_HF(num_Tol_CalP_Water),&
                       Local_freeDOF_HF(num_Tol_CalP_Water)
  real(kind=FT),intent(out)::R(Total_FD+num_Tol_CalP_Water)
  end subroutine Cal_HF_Resid
END INTERFACE
end module Global_Inter_Cal_HF_Resid

!------------------------
!  3. Cal_HF_Jacobian_NR
!------------------------
module Global_Inter_Cal_HF_Jacobian_NR
INTERFACE
  subroutine Cal_HF_Jacobian_NR(Counter,NR_Deri,Total_FD, &
        num_FreeD,num_free_CalP,globalK,Coupled_Q, &
        freeDOF,freeDOF_HF,Local_freeDOF_HF,delta_Time,H)
  use Global_Float_Type
  use Global_Crack
  use Global_Crack_Common
  implicit none
  integer, intent(in)::Counter,Total_FD,num_FreeD,num_free_CalP
  real(kind=FT),intent(in):: &
                Coupled_Q(Total_FD,num_Tol_CalP_Water), &
                globalK(Total_FD,Total_FD),delta_Time, &
                H(num_Tol_CalP_Water,num_Tol_CalP_Water)
  integer, intent(in)::freeDOF(Total_FD),&
                       freeDOF_HF(num_Tol_CalP_Water),&
                       Local_freeDOF_HF(num_Tol_CalP_Water)
  real(kind=FT),intent(out)::&
    NR_Deri(Total_FD+num_Tol_CalP_Water,&
            Total_FD+num_Tol_CalP_Water)
  end subroutine Cal_HF_Jacobian_NR
END INTERFACE
end module Global_Inter_Cal_HF_Jacobian_NR

!----------------------
! 4. Matrix_Solve_LSOE
!----------------------
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

!------------------------------
! 5. Get_Contact_State_Eles_IN
!------------------------------
module Global_Inter_Get_Contact_State_Eles_IN
INTERFACE
  SUBROUTINE Get_Contact_State_Eles_IN(i_Contact, &
                         c_Cracks_CalP_Aper, &
                         tem_Elem_Conta_Sta, &
                         Yes_Contact, &
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

!--------------------------
! 6. Cal_HF_Line_Searching
!--------------------------
module Global_Inter_Cal_HF_Line_Searching
INTERFACE
  subroutine Cal_HF_Line_Searching(ifra,iter,Counter_Iter, &
      num_ALlDOF,c_Total_FD,num_FreeD,num_free_CalP, &
    Total_freeDOF,c_num_Tol_CalP_Water,num_Total_FD,freeDOF_HF, &
         Local_freeDOF_HF,NR_Deri, &
         delta_x, &
         c_R,Initial_DISP,&
         c_DISP,c_freeDOF,c_CalP_Pres,c_globalK,Coupled_Q,F_U,&
         c_Temp_Cr_CalP_Aper,c_Temp_total_Time,&
         Num_Line_Search,total_Time,x)
  use Global_Float_Type
  use Global_Common
  use Global_Crack
  use Global_Crack_Common
  use Global_HF
  use Global_Material
  implicit none
  integer,intent(in)::ifra,iter,Counter_Iter,c_Total_FD,&
      num_FreeD,num_free_CalP,num_ALlDOF,c_num_Tol_CalP_Water,&
                 num_Total_FD
  integer,intent(in)::Local_freeDOF_HF(c_num_Tol_CalP_Water)
  integer,intent(in)::Total_freeDOF(num_ALlDOF)
  real(kind=FT),intent(in)::c_R(num_ALlDOF)
  real(kind=FT),intent(in)::Coupled_Q(c_Total_FD,c_num_Tol_CalP_Water), &
                c_globalK(c_Total_FD,c_Total_FD),F_U(c_Total_FD)
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

!-----------------------------------------
! 7. Determine_Contact_State_by_Iteration
!-----------------------------------------
module Global_Inter_Determine_Contact_State_by_Iteration
INTERFACE
  SUBROUTINE Determine_Contact_State_by_Iteration(iter,ifra,Counter_Iter,&
               Contact_DISP,c_Total_Freedom,usual_FD,enrich_FD,&
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

!-----------------------------------------
!  8. Assemble_Stiffness_Matrix_SPARS_FEM
!-----------------------------------------
module Global_Inter_Assemble_Stiffness_Matrix_SPARS_FEM
INTERFACE
  SUBROUTINE Assemble_Stiffness_Matrix_SPARS_FEM(isub,freeDOF,num_FreeD, &
                    K_CSR_aa,K_CSR_ja,K_CSR_ia,K_CSR_NNZ_Max,K_CSR_NNZ,&
                    T_Freedom,Total_Num_G_P)
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

!-----------------------------------------
! 9. Assemble_Stiffness_Matrix_SPARS_XFEM
!-----------------------------------------
module Global_Inter_Assemble_Stiffness_Matrix_SPARS_XFEM
INTERFACE
  SUBROUTINE Assemble_Stiffness_Matrix_SPARS_XFEM(isub,     &
                   freeDOF,num_FreeD,&
                   fixedDOF,num_FixedD, &
                   K_CSR_aa,K_CSR_ja,K_CSR_ia,&
                   K_CSR_NNZ_Max,K_CSR_NNZ,&
                   T_Freedom,Total_Num_G_P)
  use Global_Float_Type
  use Global_Crack
  use Global_Crack_Common
  use Global_Model
  use Global_Filename
  use Global_Common
  use Global_Material
  use Global_Contact
#ifndef Silverfrost
  use omp_lib
#endif
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

!------------------------------
! 10. Matrix_Solve_LSOE_Sparse
!------------------------------
module Global_Inter_Matrix_Solve_LSOE_Sparse
INTERFACE
  SUBROUTINE Matrix_Solve_LSOE_Sparse(Key_Indent,Key_LSOE_Sys,&
                      c_Key_SLOE,K_CSR_NNZ,&
                      K_CSR_aa,K_CSR_ja,K_CSR_ia,F,D,n)
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

!------------------------
! 11. Cal_HF_Resid_SPARS
!------------------------
module Global_Inter_Cal_HF_Resid_SPARS
INTERFACE
  subroutine Cal_HF_Resid_SPARS(ifra,iter,Counter_Iter,R, &
             Total_FD,num_FreeD,num_free_CalP,&
                       K_CSR_NNZ,K_CSR_aa,&
                       K_CSR_ja,K_CSR_ia,&
                       Coupled_Q,F_U,Last_DISP,Last_Last_DISP,&
                       freeDOF,freeDOF_HF,Local_freeDOF_HF,&
                       Last_CalP_Pres,delta_Time,H,c_S)
  use Global_Float_Type
  use Global_Crack
  use Global_Crack_Common
  use Global_HF
  use Global_Material
  implicit none
  integer, intent(in)::ifra,iter,Counter_Iter,Total_FD,num_FreeD,num_free_CalP,K_CSR_NNZ
  real(kind=FT),intent(in)::Coupled_Q(Total_FD,num_Tol_CalP_Water), &
                F_U(Total_FD),Last_DISP(Total_FD),Last_Last_DISP(Total_FD),&
                Last_CalP_Pres(num_Tol_CalP_Water),delta_Time,&
                H(num_Tol_CalP_Water,num_Tol_CalP_Water),&
                c_S(num_Tol_CalP_Water)
  integer, intent(in)::freeDOF(Total_FD),&
                       freeDOF_HF(num_Tol_CalP_Water),&
                       Local_freeDOF_HF(num_Tol_CalP_Water)
  real(kind=FT),intent(out)::R(Total_FD+num_Tol_CalP_Water)
  real(kind=FT),intent(in)::K_CSR_aa(1:K_CSR_NNZ)
  integer,intent(in)::K_CSR_ja(1:K_CSR_NNZ)
  integer,intent(in)::K_CSR_ia(1:num_FreeD+1)
  END SUBROUTINE Cal_HF_Resid_SPARS
END INTERFACE
end module Global_Inter_Cal_HF_Resid_SPARS

!------------------------------
! 12. Cal_HF_Jacobian_NR_SPARS
!------------------------------
module Global_Inter_Cal_HF_Jacobian_NR_SPARS
INTERFACE
  subroutine Cal_HF_Jacobian_NR_SPARS(Counter,Total_FD,&
                  num_FreeD,num_free_CalP,&
                  K_CSR_NNZ,K_CSR_aa,K_CSR_ja,K_CSR_ia,&
                  NRD_CSR_NNZ_Max,&
                  NRD_CSR_NNZ,NRD_CSR_aa,NRD_CSR_ja,NRD_CSR_ia,&
                  Coupled_Q,freeDOF,freeDOF_HF,&
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
  real(kind=FT),intent(in):: Coupled_Q(Total_FD,num_Tol_CalP_Water), &
                delta_Time,H(num_Tol_CalP_Water,num_Tol_CalP_Water)
  integer, intent(in)::freeDOF(Total_FD),&
                       freeDOF_HF(num_Tol_CalP_Water),&
                       Local_freeDOF_HF(num_Tol_CalP_Water)
  END SUBROUTINE Cal_HF_Jacobian_NR_SPARS
END INTERFACE
end module Global_Inter_Cal_HF_Jacobian_NR_SPARS

!---------------------------
! 13. Cal_Contact_Red_Resid
!---------------------------
module Global_Inter_Cal_Contact_Red_Resid
INTERFACE
  subroutine Cal_Contact_Red_Resid( &
         iter,ifra,Counter_Iter,i_NR_P, &
         c_Total_Freedom,c_num_freeDOF,n_freeDOF_FEM,F_U, &
         freeDOF_FEM, &
         U_e,U_all,freeDOF,Part1,Part2,K_ee, &
         K_uu_Inv,K_ue,PC_Gauss_x,PC_Gauss_y,R_red)
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
  real(kind=FT),intent(in)::&
                F_U(c_Total_Freedom),U_e(Enrich_Freedom),&
                U_all(c_Total_Freedom)
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

!--------------------------
!14. Cal_B_Matrix_Crack_3D
!--------------------------
module Global_Inter_Cal_B_Matrix_Crack_3D
INTERFACE
  subroutine Cal_B_Matrix_Crack_3D(kesi,yita,zeta,i_C,i_E,i_G,&
                         c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES,&
                         tem_B,num_tem_B)
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

!-----------------------------
! 15. Cal_N_dNdkesi_J_detJ_3D
!-----------------------------
module Global_Inter_Cal_N_dNdkesi_J_detJ_3D
INTERFACE
  subroutine Cal_N_dNdkesi_J_detJ_3D(kesi,yita,zeta,&
                                    X_NODES,Y_NODES,Z_NODES,&
                                    detJ,J,N,dNdkesi)
  use Global_Float_Type
  implicit none
  real(kind=FT),intent(in)::kesi,yita,zeta,X_NODES(8),Y_NODES(8),Z_NODES(8)
  real(kind=FT),intent(out)::detJ,dNdkesi(8,3)
  real(kind=FT),intent(out):: N(3,24),J(3,3)
  END SUBROUTINE Cal_N_dNdkesi_J_detJ_3D
END INTERFACE
end module Global_Inter_Cal_N_dNdkesi_J_detJ_3D

!-------------------------------------
!16. Location_Element_Stiff_Matrix_3D
!-------------------------------------
module Global_Inter_Location_Element_Stiff_Matrix_3D
INTERFACE
  SUBROUTINE Location_Element_Stiff_Matrix_3D(i_E,i_C,&
                           c_POS_3D_c_Ele,Location_ESM_C_Crack,&
                           num_Loc_ESM_C_Crack,&
                           Location_ESM_C_Cr_NoFEM,&
                           num_Loc_ESM_C_Cr_NoFEM)
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

!-------------------------------
! 17. Tool_Dis_Point_to_3D_Quad
!-------------------------------
module Global_Inter_Tool_Cal_Dis_Point_to_3D_Quad
INTERFACE
  subroutine Tool_Dis_Point_to_3D_Quad(Point,&
                         Quad_P1,Quad_P2,Quad_P3,Quad_P4,&
                         Distance,PER,Yes_PER_in,Yes_PER_on)
  use Global_Float_Type
  use Global_Common
  implicit none
  real(kind=FT),intent(in)::Point(3),Quad_P1(3),Quad_P2(3),Quad_P3(3),Quad_P4(3)
  real(kind=FT),intent(out):: Distance,PER(3)
  logical,intent(out):: Yes_PER_in,Yes_PER_on
  END SUBROUTINE Tool_Dis_Point_to_3D_Quad
END INTERFACE
end module Global_Inter_Tool_Cal_Dis_Point_to_3D_Quad

!---------------------------------
! 18. EBE_Cal_Contact_Jacobian_3D
!---------------------------------
module Global_EBE_Cal_Contact_Jacobian_3D
INTERFACE
subroutine EBE_Cal_Contact_Jacobian_3D(isub,i_NR_P, num_freeD, &
freeDOF,size_local_0, &
all_local_0,diag_precon_no_invert_0,diag_precon_no_invert,&
Kn,Kn_Gauss,Kt1_Gauss,Kt2_Gauss,CT_State_Gauss)

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
#ifndef Silverfrost
use OMP_LIB
#endif

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

!--------------------------------------
!18.1 Cal_Contact_Resid_3D  2022-08-11
!--------------------------------------
module Global_Cal_Contact_Resid_3D
INTERFACE
  subroutine Cal_Contact_Resid_3D(iter,ifra,Counter_Iter,i_NR_P,      &
    c_Total_Freedom,c_num_freeDOF,F_U,               &
    c_DISP,freeDOF,PC_Gauss_x,PC_Gauss_y,PC_Gauss_z,R_PSI)
  use Global_Float_Type
  use Global_Common
  use Global_Model
  use Global_Crack
  use Global_Crack_Common
  use Global_Crack_3D
  use Global_HF
  use Global_Material
#ifndef Silverfrost
  use OMP_LIB
#endif  

! --------------------------
! Variable Type Declaration
! --------------------------
  implicit none
  integer, intent(in)::iter,ifra,Counter_Iter,i_NR_P
  integer, intent(in)::c_Total_Freedom,c_num_freeDOF
  real(kind=FT),intent(in)::F_U(c_Total_Freedom),c_DISP(c_Total_Freedom)
  integer, intent(in)::freeDOF(c_Total_Freedom)
  real(kind=FT),intent(in)::PC_Gauss_x(num_Crack,Max_Max_N_FluEl_3D)    &
          ,PC_Gauss_y(num_Crack,Max_Max_N_FluEl_3D),    &
        PC_Gauss_z(num_Crack,Max_Max_N_FluEl_3D)
  real(kind=FT),intent(out)::R_PSI(c_Total_Freedom)
  END subroutine Cal_Contact_Resid_3D
END INTERFACE
end module Global_Cal_Contact_Resid_3D

!------------------------------------------------------------
!18.2 Global_Cal_Contact_Contact_State_Gauss_3D.  2022-08-12
!------------------------------------------------------------
module Global_Cal_Contact_Contact_State_Gauss_3D
INTERFACE
  SUBROUTINE Cal_Contact_Contact_State_Gauss_3D( &
                  iter,ifra,Counter_Iter,i_NR_P, &
                  c_DISP,Yes_Contact,c_Elem_Conta_Sta, &
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


!----------------------------------------------------
!  18.3 Global_Cal_Contact_PN_and_PT_3D.  2022-08-12
!----------------------------------------------------
module Global_Cal_Contact_PN_and_PT_3D
INTERFACE
  subroutine Cal_Contact_PN_and_PT_3D(iter,ifra,Counter_Iter,  &
                   i_NR_P,c_Total_Freedom,c_num_freeDOF, &
                           Kn,Kn_Gauss,Kt1_Gauss,Kt2_Gauss, &
                           fric_mu,c_DISP, &
                           delta_u_a, &
                           CT_State_Gauss,c_Elem_Conta_Sta, &
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
#ifndef Silverfrost
  use omp_lib
#endif
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
  real(kind=FT),intent(inout):: &
                PC_Gauss_x(num_Crack,Max_Max_N_FluEl_3D),&
                PC_Gauss_y(num_Crack,Max_Max_N_FluEl_3D),&
                PC_Gauss_z(num_Crack,Max_Max_N_FluEl_3D)
  integer Num_Indent
  END SUBROUTINE Cal_Contact_PN_and_PT_3D
END INTERFACE
end module Global_Cal_Contact_PN_and_PT_3D

!----------------------------------------------------
!  18.4 Global_Cal_Contact_Conve_Factor.  2022-08-12
!----------------------------------------------------
module Global_Cal_Contact_Conve_Factor
INTERFACE
  SUBROUTINE Cal_Contact_Conve_Factor(&
                iter,ifra,Counter_Iter,i_NR_P,Conve_Tolerance,&
                c_Total_Freedom,c_freeDOF,c_num_freeDOF,&
                c_F,R_PSI,Last_R_PSI,&
                delta_U,U,Contact_DISP_0,&
                Yes_Conve,Conve_Factor)
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

!------------------------------------
! 19. Function MATMUL_LP, 2021-02-14
!------------------------------------
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

!-----------------------------------
! 20. Function MVMUL_LP, 2022-07-22
!-----------------------------------
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

!------------------------------------------------------
!   21. SUBROUTINE EBE_XFEM_PCG_3D_with_K, 2022-08-12.
!------------------------------------------------------
module Global_EBE_XFEM_PCG_3D_with_K
INTERFACE
  SUBROUTINE EBE_XFEM_PCG_3D_with_K(isub,Lambda,c_cg_tol,&
               max_num_PCG,num_FreeD,freeDOF,F,disp,  &
               diag_precon_no_invert)
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
#ifndef Silverfrost
  use omp_lib
#endif
  implicit none
  integer,intent(in)::isub,max_num_PCG,num_FreeD
  real(kind=FT),intent(in)::Lambda,c_cg_tol,F(num_FreeD)
  integer,intent(in)::freeDOF(num_FreeD)
  real(kind=FT),intent(in)::diag_precon_no_invert(0:num_FreeD)
  real(kind=FT),intent(out)::disp(Total_FD)
  end SUBROUTINE EBE_XFEM_PCG_3D_with_K
END INTERFACE
end module Global_EBE_XFEM_PCG_3D_with_K

!-------------------------------------------------------------------------
! 22. SUBROUTINE EBE_Determine_Contact_State_by_Iteration_3D, 2022-08-12.
!-------------------------------------------------------------------------
module Global_EBE_Determine_Contact_State_by_Iteration_3D
INTERFACE
  SUBROUTINE EBE_Determine_Contact_State_by_Iteration_3D(isub,&
               c_cg_tol,max_num_PCG,num_FreeD,freeDOF, &
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
#ifndef Silverfrost
  use omp_lib
#endif
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

!-------------------------------------------------------------------------
! 23. SUBROUTINE EBE_Determine_Contact_State_by_Iteration_3D, 2022-08-12.
!-------------------------------------------------------------------------
module Global_D3_HF_Const_Pres_to_F_Vector
INTERFACE
  SUBROUTINE D3_HF_Const_Pres_to_F_Vector(isub,c_Pres,num_FreeD,&
                           in_Total_FD,in_num_Tol_CalP_Water,&
                           freeDOF,in_F_U,F)
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

!---------------------------------------------------------
! 24. SUBROUTINE D3_HF_Get_Pres_by_NR_Method, 2022-08-12.
!---------------------------------------------------------
module Global_D3_HF_Get_Pres_by_NR_Method
INTERFACE
  SUBROUTINE D3_HF_Get_Pres_by_NR_Method(i_WB,i_Stage,i_Prop, &
            isub,i_Time,c_Stage_Q,&
            c_Time,Max_Pres_Steps,         &
            num_FreeD,c_Total_FD,c_num_Tol_CalP_Water,&
            freeDOF,c_F_U,F,Lambda,   &
            diag_precon_no_invert,DISP,         &
            NR_delta_Pres,f_Vol_Tol,Pres_Tol,c_Pres,Output_Pres,&
            K_CSR_NNZ,K_CSR_aa,K_CSR_ja,K_CSR_ia)
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

!---------------------------------------------------
! 25. SUBROUTINE EBE_XFEM_PCG_3D_Get_K, 2022-08-12.
!---------------------------------------------------
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
#ifndef Silverfrost
  use OMP_LIB
#endif
  implicit none
  integer,intent(in)::isub,num_FreeD
  integer,intent(in)::freeDOF(num_FreeD)
  real(kind=FT),intent(out)::diag_precon_no_invert(0:num_FreeD)
  END SUBROUTINE EBE_XFEM_PCG_3D_Get_K
END INTERFACE
end module Global_EBE_XFEM_PCG_3D_Get_K

!---------------------------------------------------
!26. SUBROUTINE Ele_by_Ele_XFEM_PCG_3D, 2022-08-12.
!---------------------------------------------------
module Global_Ele_by_Ele_XFEM_PCG_3D
INTERFACE
  SUBROUTINE Ele_by_Ele_XFEM_PCG_3D(isub,Lambda,c_cg_tol, &
     max_num_PCG,num_FreeD,freeDOF,F,disp,diag_precon_no_invert)
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
#ifndef Silverfrost  
  use omp_lib
#endif
  implicit none
  integer,intent(in)::isub,max_num_PCG,num_FreeD
  real(kind=FT),intent(in)::Lambda,c_cg_tol,F(num_FreeD)
  integer,intent(in)::freeDOF(num_FreeD)
  real(kind=FT),intent(out)::disp(Total_FD)
  real(kind=FT),intent(out)::diag_precon_no_invert(0:num_FreeD)
  end SUBROUTINE Ele_by_Ele_XFEM_PCG_3D
END INTERFACE
end module Global_Ele_by_Ele_XFEM_PCG_3D

!-----------------------------------------------------------------------
! 27. SUBROUTINE Global_Cal_Ele_Stiffness_Matrix_3D_8nodes, 2022-09-06.
!-----------------------------------------------------------------------
module Global_Cal_Ele_Stiffness_Matrix_3D_8nodes
INTERFACE
subroutine Cal_Ele_Stiffness_Matrix_3D_8nodes(i_E,num_Gauss, &
          X_NODES,Y_NODES,Z_NODES, &
          c_D,kesi,yita,zeta,weight,localK)
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
!28. SUBROUTINE Global_Cal_Ele_Num_by_Coors_3D, 2022-09-24.
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
!-------------------------------------------------------------
! 30. SUBROUTINE Tool_Yes_Point_in_3D_Hexahedron, 2023-08-14.
!-------------------------------------------------------------
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

!----------------------------------------------
! 31. SUBROUTINE Vector_Normalize, 2023-08-14.
!----------------------------------------------
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

!-------------------------------------------------------------------
! 32. SUBROUTINE Tool_ThetaX_ThetaY_ThetaZ_3D_rotation, 2023-08-14.
!-------------------------------------------------------------------
module Global_INTERFACE_Tool_ThetaX_ThetaY_ThetaZ_3D_rotation
INTERFACE
subroutine Tool_ThetaX_ThetaY_ThetaZ_3D_rotation(i_V,Vector_Origi_1, &
                Vector_Origi_2,Vector_Origi_3,Vector_Crack_1,Vector_Crack_2,&
                Vector_Crack_3,ThetaX,ThetaY,ThetaZ,T_Matrix)


use Global_Float_Type
use Global_Crack_3D

implicit none
integer,intent(in)::i_V
real(kind=FT),intent(in)::  &
   Vector_Origi_1(3),Vector_Origi_2(3),Vector_Origi_3(3),&
   Vector_Crack_1(3),Vector_Crack_2(3),Vector_Crack_3(3)
real(kind=FT),intent(out):: ThetaX,ThetaY,ThetaZ,T_Matrix(3,3)
end SUBROUTINE Tool_ThetaX_ThetaY_ThetaZ_3D_rotation
END INTERFACE
end module Global_INTERFACE_Tool_ThetaX_ThetaY_ThetaZ_3D_rotation

!----------------------------------------------------------------------
! 33. SUBROUTINE Tool_Yes_Point_in_3D_Hexahedron_with_Tol, 2023-08-14.
!----------------------------------------------------------------------
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

!----------------------------------------------------
! 34. SUBROUTINE Vector_Location_Int_v2, 2023-08-14.
!----------------------------------------------------
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


!-------------------------------------------------------------
! 35. SUBROUTINE D3_Get_Signed_Dis_to_Crack_Mesh, 2023-08-14.
!-------------------------------------------------------------
module Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh
INTERFACE
SUBROUTINE D3_Get_Signed_Dis_to_Crack_Mesh(Point,i_C,Check_Ball_R,  &
                  Signed_Dis,Signed_Dis_v2,c_Yes_Node_PER_in_FS,  &
                  c_PER_Node_to_FS,Yes_Found_Min_Signed_Dis,  &
                  n_Vector)
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

!--------------------------------------------------------------------------------
! 36. SUBROUTINE D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth, 2023-08-14.
!--------------------------------------------------------------------------------
module Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane
INTERFACE
SUBROUTINE D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth(Point,i_C,&
                               Check_Ball_R,Signed_Dis,                 &
                               c_Yes_Node_PER_in_FS,c_PER_Node_to_FS,   &
                               Yes_Found_Min_Signed_Dis,n_Vector)
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

!----------------------------------------------------
! 37. Matrix_Solve_LSOE_Complex Subroutine INTERFACE
!----------------------------------------------------
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
! 38. Matrix_Solve_LSOE_Complex_Sparse Subroutine INTERFACE
!-----------------------------------------------------------
module module_INTERFACE_Matrix_Solve_LSOE_Complex_Sparse
    INTERFACE
    SUBROUTINE Matrix_Solve_LSOE_Complex_Sparse(Key_Indent,Key_LSOE_Sys,c_Key_SLOE,K_CSR_NNZ,&
                     K_CSR_aa,K_CSR_ja,K_CSR_ia,  &
                     Complex_K_vals,Complex_K_rows,Complex_K_cols,   &
                     F,D,n)   
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

!-------------------------------------------------
! 39. SPARSKIT2_amux_complex subroutine INTERFACE
!-------------------------------------------------
#ifndef Silverfrost
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
#endif 

!-------------------------------------------------------
! 40. Vector_belongs_Matrix_Is_Dou Subroutine INTERFACE
!-------------------------------------------------------
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

!-------------------------------------------------------------------
! 41. Matrix_n_x_2_Quick_Sort_Int Subroutine INTERFACE. 2024-11-07.
!-------------------------------------------------------------------
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
  
!---------------------------------------------------------------
! 42. Matrix_n_x_3_Check_Sort Subroutine INTERFACE. 2024-11-07.
!---------------------------------------------------------------
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

!-------------------------------------------------------------------
! 44. Matrix_n_x_3_Check_Sort_Int Subroutine INTERFACE. 2024-11-07.
!-------------------------------------------------------------------
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

!----------------------------------------------
! 45. coocsr subroutine INTERFACE. 2024-11-07.
!----------------------------------------------
module module_INTERFACE_coocsr
    INTERFACE
    subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
    real*8 a(*),ao(*),x
    integer ir(*),jc(*),jao(*),iao(*)
    END subRoutine coocsr
    END INTERFACE
end module module_INTERFACE_coocsr
