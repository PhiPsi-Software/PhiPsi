!-----------------------------------------------------------
! Brief: Compute 2D quasi-static energy quantities for energy
!        balance verification.
!
!        Elastic strain energy:  Pi_e = integral (1/2)*eps:D:eps dOmega
!        Fracture energy:       E_f  = G_c * L   (Griffith, 2D)
!        External work:         W_ext = accumulated F^T * delta_U
!        Residual:              R    = W_ext - Pi_e - E_f
!        Normalized residual:   r    = |R| / |W_ext| * 100%
!
! Parameters:
!   Input:  isub             - current load substep index
!           DISP(:)          - global displacement vector
!           globalF(:)       - global external force vector
!
!   Output: Elastic_Strain_Energy     (module Global_Common)
!           Fracture_Energy           (module Global_Common)
!           External_Work             (module Global_Common)
!           Residual_Energy           (module Global_Common)
!           Normalized_Residual_Energy (module Global_Common)
!
! Notes:   Elastic strain energy is computed by Gauss-point
!          integration.  Standard FEM elements use the exact
!          B-matrix (identical to Cal_Ele_Stress_N4).  Enriched
!          XFEM elements use the stored Gauss-point stress with
!          the compliance matrix S and an element-area-based
!          volume approximation (exact in total).
!          External work is accumulated incrementally across
!          substeps using a SAVEd previous-displacement vector.
!          Crack length is computed from the Crack_Coor polyline.
!-----------------------------------------------------------

subroutine Cal_Energy_2D(isub, DISP, globalF)

!-------------------------------
! Read public variable modules
!-------------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Material
use Global_Crack
use Global_Crack_Common
use Global_Stress
use omp_lib

implicit none

integer, intent(in) :: isub
real(kind=FT), intent(in) :: DISP(Total_FD), globalF(Total_FD)

integer :: i_E, i_G, i_C, i_pt, mat_num, c_NN(4), n_pts
real(kind=FT) :: c_D(3,3), c_S(3,3), c_thick, U_e(8)
real(kind=FT) :: c_X_NODES(4), c_Y_NODES(4)
real(kind=FT) :: c_stress(3), strain_energy_density, vol_per_gp
integer :: gp_idx
real(kind=FT), allocatable :: kesi_array(:), yita_array(:), weight_array(:)
real(kind=FT) :: kesi, yita, weight, detJ
real(kind=FT) :: X_c(4), Y_c(4), JM(4,4)
real(kind=FT) :: aa, bb, cc, dd, Nkesi(4), Nyita(4)
real(kind=FT) :: B1(3,2), B2(3,2), B3(3,2), B4(3,2), ToTal_B(3,8)
real(kind=FT) :: epsilon_vec(3), sigma_vec(3), energy_contrib
real(kind=FT) :: local_Elastic_Strain_Energy
integer :: max_GP, num_GP_this_elem

real(kind=FT) :: crack_length, Gc_val, E_eff_val
real(kind=FT) :: crack_x, crack_y
integer :: OUT_Elem
real(kind=FT) :: local_Fracture_Energy

real(kind=FT), allocatable, save :: Previous_DISP(:)
real(kind=FT), allocatable, save :: Previous_globalF(:)
real(kind=FT), save :: Accumulated_External_Work = 0.0_FT
logical, save :: First_Call = .True.
real(kind=FT) :: delta_W, delta_U, avg_F
integer :: i_dof

local_Elastic_Strain_Energy = ZR
local_Fracture_Energy       = ZR


max_GP = Num_Gauss_P_FEM
do i_E = 1, Num_Elem
    if (num_GP_Elem(i_E) > max_GP) max_GP = num_GP_Elem(i_E)
end do

allocate(kesi_array(max_GP), yita_array(max_GP), weight_array(max_GP))

call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM, kesi_array, yita_array, weight_array)

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP& PRIVATE(i_E, i_G, mat_num, c_D, c_thick, c_NN, c_X_NODES, c_Y_NODES, U_e, &
!$OMP&         num_GP_this_elem, kesi, yita, weight, X_c, Y_c, JM, aa, bb, cc, dd, &
!$OMP&         Nkesi, Nyita, detJ, B1, B2, B3, B4, ToTal_B, epsilon_vec, sigma_vec, &
!$OMP&         c_S, gp_idx, c_stress, strain_energy_density, vol_per_gp, energy_contrib) &
!$OMP& REDUCTION(+:local_Elastic_Strain_Energy)
do i_E = 1, Num_Elem

    if (Elem_Break(i_E)) cycle

    mat_num = Elem_Mat(i_E)
    c_D     = D(mat_num, 1:3, 1:3)
    c_thick = thick(mat_num)

    c_NN      = G_NN(:, i_E)
    c_X_NODES = G_X_NODES(:, i_E)
    c_Y_NODES = G_Y_NODES(:, i_E)

    U_e(1) = DISP(c_NN(1)*2-1)
    U_e(2) = DISP(c_NN(1)*2)
    U_e(3) = DISP(c_NN(2)*2-1)
    U_e(4) = DISP(c_NN(2)*2)
    U_e(5) = DISP(c_NN(3)*2-1)
    U_e(6) = DISP(c_NN(3)*2)
    U_e(7) = DISP(c_NN(4)*2-1)
    U_e(8) = DISP(c_NN(4)*2)

    num_GP_this_elem = num_GP_Elem(i_E)

    do i_G = 1, num_GP_this_elem

        if (num_GP_this_elem <= Num_Gauss_P_FEM) then
            kesi   = kesi_array(i_G)
            yita   = yita_array(i_G)
            weight = weight_array(i_G)

            X_c = c_X_NODES
            Y_c = c_Y_NODES

            JM(1,1) = ZR
            JM(1,2) = ONE - yita
            JM(1,3) = yita - kesi
            JM(1,4) = kesi - ONE
            JM(2,1) = yita - ONE
            JM(2,2) = ZR
            JM(2,3) = ONE + kesi
            JM(2,4) = -kesi - yita
            JM(3,1) = kesi - yita
            JM(3,2) = -kesi - ONE
            JM(3,3) = ZR
            JM(3,4) = yita + ONE
            JM(4,1) = ONE - kesi
            JM(4,2) = yita + kesi
            JM(4,3) = -yita - ONE
            JM(4,4) = ZR

            Nkesi(1)  = ( yita - ONE) / FOU
            Nkesi(2)  = (-yita + ONE) / FOU
            Nkesi(3)  = ( yita + ONE) / FOU
            Nkesi(4)  = (-yita - ONE) / FOU
            Nyita(1)  = ( kesi - ONE) / FOU
            Nyita(2)  = (-kesi - ONE) / FOU
            Nyita(3)  = ( kesi + ONE) / FOU
            Nyita(4)  = (-kesi + ONE) / FOU

            aa = (Y_c(1)*(kesi-ONE)   + Y_c(2)*(-ONE-kesi) + &
                  Y_c(3)*(ONE+kesi)   + Y_c(4)*(ONE-kesi))  / FOU
            bb = (Y_c(1)*(yita-ONE)   + Y_c(2)*(ONE-yita)  + &
                  Y_c(3)*(ONE+yita)   + Y_c(4)*(-ONE-yita)) / FOU
            cc = (X_c(1)*(yita-ONE)   + X_c(2)*(ONE-yita)  + &
                  X_c(3)*(ONE+yita)   + X_c(4)*(-ONE-yita)) / FOU
            dd = (X_c(1)*(kesi-ONE)   + X_c(2)*(-ONE-kesi) + &
                  X_c(3)*(ONE+kesi)   + X_c(4)*(ONE-kesi))  / FOU

            detJ = dot_product(MATMUL(X_c, JM), Y_c) / EIG

            B1(1,1) = aa*Nkesi(1) - bb*Nyita(1)
            B1(1,2) = ZR
            B1(2,1) = ZR
            B1(2,2) = cc*Nyita(1) - dd*Nkesi(1)
            B1(3,1) = cc*Nyita(1) - dd*Nkesi(1)
            B1(3,2) = aa*Nkesi(1) - bb*Nyita(1)

            B2(1,1) = aa*Nkesi(2) - bb*Nyita(2)
            B2(1,2) = ZR
            B2(2,1) = ZR
            B2(2,2) = cc*Nyita(2) - dd*Nkesi(2)
            B2(3,1) = cc*Nyita(2) - dd*Nkesi(2)
            B2(3,2) = aa*Nkesi(2) - bb*Nyita(2)

            B3(1,1) = aa*Nkesi(3) - bb*Nyita(3)
            B3(1,2) = ZR
            B3(2,1) = ZR
            B3(2,2) = cc*Nyita(3) - dd*Nkesi(3)
            B3(3,1) = cc*Nyita(3) - dd*Nkesi(3)
            B3(3,2) = aa*Nkesi(3) - bb*Nyita(3)

            B4(1,1) = aa*Nkesi(4) - bb*Nyita(4)
            B4(1,2) = ZR
            B4(2,1) = ZR
            B4(2,2) = cc*Nyita(4) - dd*Nkesi(4)
            B4(3,1) = cc*Nyita(4) - dd*Nkesi(4)
            B4(3,2) = aa*Nkesi(4) - bb*Nyita(4)

            ToTal_B(1:3, 1:2) = B1
            ToTal_B(1:3, 3:4) = B2
            ToTal_B(1:3, 5:6) = B3
            ToTal_B(1:3, 7:8) = B4

            epsilon_vec = MATMUL(ToTal_B, U_e) / detJ
            sigma_vec   = MATMUL(c_D, epsilon_vec)

            energy_contrib = HLF * dot_product(epsilon_vec, sigma_vec) * &
                             detJ * weight * c_thick
            local_Elastic_Strain_Energy = local_Elastic_Strain_Energy + energy_contrib

        else
            gp_idx = Ele_GP_Start_Num(i_E) + i_G - 1
            c_stress(1) = Stress_xx_Gauss(gp_idx)
            c_stress(2) = Stress_yy_Gauss(gp_idx)
            c_stress(3) = Stress_xy_Gauss(gp_idx)
            c_S = S(mat_num, 1:3, 1:3)
            strain_energy_density = HLF * dot_product(c_stress, MATMUL(c_S, c_stress))
            vol_per_gp = Elem_Area(i_E) * c_thick / real(num_GP_this_elem, FT)
            energy_contrib = strain_energy_density * vol_per_gp
            local_Elastic_Strain_Energy = local_Elastic_Strain_Energy + energy_contrib

        end if

    end do
end do
!$OMP END PARALLEL DO

deallocate(kesi_array, yita_array, weight_array)

Elastic_Strain_Energy = local_Elastic_Strain_Energy

if (num_Crack > 0) then
    do i_C = 1, num_Crack
        n_pts = Each_Cr_Poi_Num(i_C)
        if (n_pts < 2) cycle
        
        crack_length = ZR

        do i_pt = 1, n_pts - 1
            crack_length = crack_length + &
                sqrt((Crack_Coor(i_C, i_pt+1, 1) - Crack_Coor(i_C, i_pt, 1))**2 + &
                     (Crack_Coor(i_C, i_pt+1, 2) - Crack_Coor(i_C, i_pt, 2))**2)
        end do
        
        if(crack_length>=(Cracks_Growth_Length_Last_Step(i_C)-Tol_5)) then
            crack_length = crack_length - Cracks_Growth_Length_Last_Step(i_C)
        endif
        
        
        crack_x = Crack_Coor(i_C, 1, 1)
        crack_y = Crack_Coor(i_C, 1, 2)
        call Cal_Ele_Num_by_Coors(crack_x, crack_y, OUT_Elem)
        if (OUT_Elem > 0 .and. OUT_Elem <= Num_Elem) then
            mat_num = Elem_Mat(OUT_Elem)
        else
            mat_num = 1
        end if

        if (Key_Type_2D == 2) then
            E_eff_val = E(mat_num, 1) / (ONE - v(mat_num, 1)**2)
        else
            E_eff_val = E(mat_num, 1)
        end if

        Gc_val = KIc(mat_num, 1)**2 / E_eff_val
        
        if(crack_length >= (Initial_Cracks_Length(i_C)-Tol_5)) then
            crack_length = crack_length - Initial_Cracks_Length(i_C)
        endif

        local_Fracture_Energy = local_Fracture_Energy + &
                                 Gc_val * crack_length * thick(mat_num)
    end do
end if

Fracture_Energy = local_Fracture_Energy

if (First_Call) then
    allocate(Previous_DISP(Total_FD))
    allocate(Previous_globalF(Total_FD))
    Previous_DISP(:)    = ZR
    Previous_globalF(:) = ZR
    Accumulated_External_Work = ZR
    First_Call = .False.
end if

if (size(Previous_DISP) /= Total_FD) then
    block
        real(kind=FT), allocatable :: tmp(:)
        integer :: old_sz
        old_sz = size(Previous_DISP)
        allocate(tmp(Total_FD))
        tmp(:) = ZR
        tmp(1:min(old_sz, Total_FD)) = Previous_DISP(1:min(old_sz, Total_FD))
        call move_alloc(tmp, Previous_DISP)
        allocate(tmp(Total_FD))
        tmp(:) = ZR
        tmp(1:min(old_sz, Total_FD)) = Previous_globalF(1:min(old_sz, Total_FD))
        call move_alloc(tmp, Previous_globalF)
    end block
end if

delta_W = ZR
do i_dof = 1, Total_FD
    delta_U = DISP(i_dof) - Previous_DISP(i_dof)
    avg_F   = HLF * (Previous_globalF(i_dof) + globalF(i_dof))
    delta_W = delta_W + avg_F * delta_U
end do
Accumulated_External_Work = Accumulated_External_Work + delta_W

External_Work = Accumulated_External_Work

do i_dof = 1, Total_FD
    Previous_DISP(i_dof)    = DISP(i_dof)
    Previous_globalF(i_dof) = globalF(i_dof)
end do

Residual_Energy = External_Work - Elastic_Strain_Energy - Fracture_Energy

if (abs(External_Work) > Tol_30) then
    Normalized_Residual_Energy = abs(Residual_Energy) / abs(External_Work) * 100.0_FT
else
    Normalized_Residual_Energy = ZR
end if

end subroutine Cal_Energy_2D
