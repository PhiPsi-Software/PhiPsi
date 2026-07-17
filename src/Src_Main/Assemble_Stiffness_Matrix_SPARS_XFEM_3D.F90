!-----------------------------------------------------------
! Brief: Assemble the 3D XFEM stiffness matrix in CSR
!        sparse storage, including crack enrichments.
!
! Parameters:
!   Input:  isub, freeDOF, fixedDOF, num_FreeD, num_FixedD,
!           K_CSR_NNZ_Max, T_Freedom
!   Output: K_CSR_aa, K_CSR_ja, K_CSR_ia - CSR storage
!           K_CSR_NNZ - actual nonzero count
!           Total_Num_G_P - total Gauss points
!
! Notes:   Supports multiple cracks, composite materials,
!          sub-cube and multi-crack integration schemes.
!-----------------------------------------------------------

SUBROUTINE Assemble_Stiffness_Matrix_SPARS_XFEM_3D(isub, freeDOF,num_FreeD, fixedDOF,num_FixedD, &
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
use module_INTERFACE_coocsr
use Global_POST,only:Key_Save_Nothing,Key_Post_Elements_Gauss_Num
use Global_Crack_3D,only:c_POS_3D,Enriched_Node_Type_3D,Elem_num_Related_Cracks

implicit none
integer,intent(in)::isub,num_FreeD,T_Freedom
integer,intent(in)::K_CSR_NNZ_Max
integer,intent(in)::freeDOF(1:num_FreeD)
integer,intent(in)::num_FixedD
integer,intent(in)::fixedDOF(1:num_FixedD)
integer,intent(out)::Total_Num_G_P
real(kind=FT),intent(out)::K_CSR_aa(K_CSR_NNZ_Max)
integer,intent(out)::K_CSR_ja(K_CSR_NNZ_Max)
integer,intent(out)::K_CSR_ia(num_FreeD+1)
integer,intent(out)::K_CSR_NNZ
integer i_C,i_E,i_G
real(kind=FT) c_D(6,6)
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
integer c_NN(8)
real(kind=FT) kesi_Enr(Num_Gau_Points_3D),yita_Enr(Num_Gau_Points_3D), &
zeta_Enr(Num_Gau_Points_3D), weight_Enr(Num_Gau_Points_3D)
real(kind=FT) kesi_No_Enr(Num_Gauss_P_FEM_3D),yita_No_Enr(Num_Gauss_P_FEM_3D), &
zeta_No_Enr(Num_Gauss_P_FEM_3D),weight_No_Enr(Num_Gauss_P_FEM_3D)
real(kind=FT),ALLOCATABLE::kesi_Enr_Cubes(:),yita_Enr_Cubes(:),zeta_Enr_Cubes(:),weight_Enr_Cubes(:)
real(kind=FT)    kesi_Enr_MC(Num_Gau_Points_3D_MC), yita_Enr_MC(Num_Gau_Points_3D_MC), &
zeta_Enr_MC(Num_Gau_Points_3D_MC), weight_Enr_MC(Num_Gau_Points_3D_MC)
real(kind=FT) ,ALLOCATABLE::kesi(:),yita(:),zeta(:),weight(:)
integer:: Location_ESM(MDOF_3D)
integer::Location_ESM_C_Crack(MDOF_3D),Location_ESM_C_Cr_NoFEM(MDOF_3D)
integer num_Loc_ESM_C_Cr_NoFEM
real(kind=FT) B(6,MDOF_3D),tem_B(6,MDOF_3D)
integer num_B,num_tem_B
integer num_Loc_ESM,num_Loc_ESM_C_Crack
integer c_Num_Gauss_Point
real(kind=FT) detJ
real(kind=FT) Rot_c_D_Comp(6,6),c_D_Comp(6,6),Volume_Ratio
real(kind=FT) T_Matrix(6,6),TT_Matrix(6,6)
integer mat_num
real(kind=FT) localK(MDOF_3D,MDOF_3D)
integer c_POS_3D_c_Ele(8)
integer c_row,c_col
integer real_row,real_col
integer i,j
integer k_each_nnz_max
integer, parameter :: integer_1 = selected_int_kind(1)
integer i_row,i_col

! ==== low-memory COO storage (encoded key + value) ====
integer(kind=8),allocatable :: key_arr(:)
real(kind=FT),  allocatable :: val_arr(:)
integer(kind=8)             :: NF_8
integer                     :: Temp_K_Each_Count
integer                     :: M
real(kind=FT)               :: accum
integer,allocatable:: K_rows(:), K_cols(:)
real(kind=FT),allocatable:: K_vals(:)


! ==== OpenMP related variables ====
integer,allocatable :: num_K_Each_Elem(:)
integer,allocatable :: elem_K_offset(:)
integer :: l_count

1112 FORMAT(5X,'Real density ratio of K is ',F12.4,'%')

!***************************************************************
! Step 0
!***************************************************************
K_CSR_aa(1:K_CSR_NNZ_Max)   = ZR
K_CSR_ja(1:K_CSR_NNZ_Max)   = 0
K_CSR_ia(1:num_FreeD+1)   = 0
K_CSR_NNZ = 0
Total_Num_G_P  = 0

call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,kesi_No_Enr,yita_No_Enr,zeta_No_Enr,weight_No_Enr)

if (Key_Integral_Sol  == 2) then
    call Cal_Gauss_Points_3D_8nodes(Num_Gau_Points_3D,kesi_Enr,yita_Enr,zeta_Enr,weight_Enr)
    call Cal_Gauss_Points_3D_8nodes(Num_Gau_Points_3D_MC,kesi_Enr_MC,yita_Enr_MC,zeta_Enr_MC,weight_Enr_MC)
elseif (Key_Integral_Sol  == 3) then
    allocate(kesi_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    allocate(yita_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    allocate(zeta_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    allocate(weight_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    call Cal_Gauss_Points_3D_for_SUBCUBES(Num_Sub_3D_Cubes,Num_Gau_Points_3D_Cube,kesi_Enr_Cubes, &
    yita_Enr_Cubes,zeta_Enr_Cubes,weight_Enr_Cubes)
endif

allocate(num_K_Each_Elem(Num_Elem))
allocate(elem_K_offset(Num_Elem))
num_K_Each_Elem(1:Num_Elem) = 0

!********************************************************
! Step 1: count nonzeros per element  (PARALLEL)
!********************************************************
EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gau_Points_3D)= .False.

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(i_E,i_C,c_NN,c_Num_Gauss_Point,Location_ESM,num_Loc_ESM, &
!$OMP         Location_ESM_C_Crack,num_Loc_ESM_C_Crack, &
!$OMP         Location_ESM_C_Cr_NoFEM,num_Loc_ESM_C_Cr_NoFEM, &
!$OMP         c_POS_3D_c_Ele,i_row,i_col,c_row,c_col,l_count) &
!$OMP REDUCTION(+:Total_Num_G_P) &
!$OMP SCHEDULE(DYNAMIC)
do i_E = 1,Num_Elem
    c_NN    = G_NN(:,i_E)
    l_count = 0
    select case(Key_Integral_Sol)
    case(2)
        if (num_Crack/= 0 .and. (maxval(Enriched_Node_Type_3D(c_NN,1:num_Crack)).ne.0))then
            c_Num_Gauss_Point            = Num_Gau_Points_3D
            Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
        else
            c_Num_Gauss_Point            = Num_Gauss_P_FEM_3D
            Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
        end if
        if(Elem_num_Related_Cracks(i_E)>=2) then
            c_Num_Gauss_Point               = Num_Gau_Points_3D_MC
            Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
        endif
    case(3)
        if (num_Crack/= 0 .and. (maxval(Enriched_Node_Type_3D(c_NN,1:num_Crack)).ne.0))then
            c_Num_Gauss_Point            = Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube
            Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
        else
            c_Num_Gauss_Point            = Num_Gauss_P_FEM_3D
            Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
        end if
    case(4)
        !To Be Done.
    end select

    if(Key_Save_Nothing==0) then
        if(Key_Post_Elements_Gauss_Num==1) then
            Elements_Gauss_Num(i_E)      =  c_Num_Gauss_Point
        endif
    endif

    num_GP_Elem(i_E) = c_Num_Gauss_Point

    Location_ESM(1:MDOF_3D)   = 0
    num_Loc_ESM           = 0
    if(num_Crack/=0)then
        do i_C =1,num_Crack
            c_POS_3D_c_Ele(1:8) = c_POS_3D(c_NN,i_C)
            if(i_C >1 .and. sum(c_POS_3D_c_Ele(1:8))==0) cycle
            call Location_Element_Stiff_Matrix_3D(i_E,i_C,c_POS_3D_c_Ele(1:8),Location_ESM_C_Crack, &
            num_Loc_ESM_C_Crack,Location_ESM_C_Cr_NoFEM,num_Loc_ESM_C_Cr_NoFEM)
            Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
            num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
        end do
    endif

    do i_row = 1,num_Loc_ESM
        do i_col = 1,num_Loc_ESM
            c_row = Location_ESM(i_row)
            c_col = Location_ESM(i_col)
            if(Flag_FreeDOF(c_row)==1  .and. Flag_FreeDOF(c_col)==1) then
                l_count =  l_count + 1
            endif
        end do
    end do
    num_K_Each_Elem(i_E) = l_count
end do
!$OMP END PARALLEL DO

! ---- Prefix sum: compute per-element write offset (serial, cheap) ----
elem_K_offset(1) = 0
do i_E = 2,Num_Elem
    elem_K_offset(i_E) = elem_K_offset(i_E-1) + num_K_Each_Elem(i_E-1)
end do
Temp_K_Each_Count = elem_K_offset(Num_Elem) + num_K_Each_Elem(Num_Elem)

!**************************************************************
! Step 2: assemble nonzeros of stiffness matrix  (PARALLEL)
!   -> store as encoded int64 key + real value (low memory)
!**************************************************************
K_Each_NNZ_Max = Temp_K_Each_Count
NF_8 = int(num_FreeD,8)

allocate(key_arr(max(Temp_K_Each_Count,1)))
allocate(val_arr(max(Temp_K_Each_Count,1)))

EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gau_Points_3D)= .False.


!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(i_E,mat_num,localK,c_D,Volume_Ratio,c_D_Comp,T_Matrix,TT_Matrix, &
!$OMP         Rot_c_D_Comp,c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES, &
!$OMP         kesi,yita,zeta,weight,c_Num_Gauss_Point, &
!$OMP         Location_ESM,num_Loc_ESM,i_C,c_POS_3D_c_Ele, &
!$OMP         Location_ESM_C_Crack,num_Loc_ESM_C_Crack, &
!$OMP         Location_ESM_C_Cr_NoFEM,num_Loc_ESM_C_Cr_NoFEM, &
!$OMP         i_G,B,num_B,detJ,tem_B,num_tem_B, &
!$OMP         i_row,i_col,c_row,c_col,real_row,real_col,l_count) &
!$OMP SCHEDULE(DYNAMIC)
do i_E = 1,Num_Elem
    l_count = elem_K_offset(i_E)
    mat_num = Elem_Mat(i_E)
    localK(1:MDOF_3D,1:MDOF_3D) = ZR
    c_D(1:6,1:6)     = D(Elem_Mat(i_E),1:6,1:6)
    if (Material_Type(mat_num)==5)then
        Volume_Ratio = Material_Para_Added(mat_num,10)
        c_D_Comp = D_Comp(mat_num,1:6,1:6)
        T_Matrix = Ele_ComMat_RotMatrix(i_E,1:6,1:6)
        TT_Matrix= TRANSPOSE(T_Matrix)
        Rot_c_D_Comp = MATMUL(TT_Matrix,c_D_Comp)
        Rot_c_D_Comp = MATMUL(Rot_c_D_Comp,T_Matrix)
        c_D =(ONE-Volume_Ratio)*c_D + Volume_Ratio*Rot_c_D_Comp
    endif
    c_NN    = G_NN(1:8,i_E)
    c_X_NODES = G_X_NODES(1:8,i_E)
    c_Y_NODES = G_Y_NODES(1:8,i_E)
    c_Z_NODES = G_Z_NODES(1:8,i_E)

    select case(Key_Integral_Sol)
    case(2)
        if (num_Crack/= 0 .and. (maxval(Enriched_Node_Type_3D(c_NN,1:num_Crack)).ne.0))then
            if(allocated(kesi)) deallocate(kesi)
            if(allocated(yita)) deallocate(yita)
            if(allocated(zeta)) deallocate(zeta)
            if(allocated(weight)) deallocate(weight)
            allocate(kesi(Num_Gau_Points_3D))
            allocate(yita(Num_Gau_Points_3D))
            allocate(zeta(Num_Gau_Points_3D))
            allocate(weight(Num_Gau_Points_3D))
            kesi(1:Num_Gau_Points_3D)    = kesi_Enr
            yita(1:Num_Gau_Points_3D)    = yita_Enr
            zeta(1:Num_Gau_Points_3D)    = zeta_Enr
            weight(1:Num_Gau_Points_3D)  = weight_Enr
            c_Num_Gauss_Point            = Num_Gau_Points_3D
        else
            if(allocated(kesi)) deallocate(kesi)
            if(allocated(yita)) deallocate(yita)
            if(allocated(zeta)) deallocate(zeta)
            if(allocated(weight)) deallocate(weight)
            allocate(kesi(Num_Gauss_P_FEM_3D))
            allocate(yita(Num_Gauss_P_FEM_3D))
            allocate(zeta(Num_Gauss_P_FEM_3D))
            allocate(weight(Num_Gauss_P_FEM_3D))
            kesi(1:Num_Gauss_P_FEM_3D)   = kesi_No_Enr
            yita(1:Num_Gauss_P_FEM_3D)   = yita_No_Enr
            zeta(1:Num_Gauss_P_FEM_3D)   = zeta_No_Enr
            weight(1:Num_Gauss_P_FEM_3D) = weight_No_Enr
            c_Num_Gauss_Point            = Num_Gauss_P_FEM_3D
        end if
        if(Elem_num_Related_Cracks(i_E)>=2) then
            if(allocated(kesi)) deallocate(kesi)
            if(allocated(yita)) deallocate(yita)
            if(allocated(zeta)) deallocate(zeta)
            if(allocated(weight)) deallocate(weight)
            allocate(kesi(Num_Gau_Points_3D_MC))
            allocate(yita(Num_Gau_Points_3D_MC))
            allocate(zeta(Num_Gau_Points_3D_MC))
            allocate(weight(Num_Gau_Points_3D_MC))
            kesi(1:Num_Gau_Points_3D_MC)    = kesi_Enr_MC
            yita(1:Num_Gau_Points_3D_MC)    = yita_Enr_MC
            zeta(1:Num_Gau_Points_3D_MC)    = zeta_Enr_MC
            weight(1:Num_Gau_Points_3D_MC)  = weight_Enr_MC
            c_Num_Gauss_Point               = Num_Gau_Points_3D_MC
        endif
    case(3)
        if (num_Crack/= 0 .and. (maxval(Enriched_Node_Type_3D(c_NN,1:num_Crack)).ne.0))then
            if(allocated(kesi)) deallocate(kesi)
            if(allocated(yita)) deallocate(yita)
            if(allocated(zeta)) deallocate(zeta)
            if(allocated(weight)) deallocate(weight)
            allocate(kesi(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
            allocate(yita(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
            allocate(zeta(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
            allocate(weight(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
            kesi(1:Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube)    = kesi_Enr_Cubes
            yita(1:Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube)    = yita_Enr_Cubes
            zeta(1:Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube)    = zeta_Enr_Cubes
            weight(1:Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube)  = weight_Enr_Cubes
            c_Num_Gauss_Point            = Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube
        else
            if(allocated(kesi)) deallocate(kesi)
            if(allocated(yita)) deallocate(yita)
            if(allocated(zeta)) deallocate(zeta)
            if(allocated(weight)) deallocate(weight)
            allocate(kesi(Num_Gauss_P_FEM_3D))
            allocate(yita(Num_Gauss_P_FEM_3D))
            allocate(zeta(Num_Gauss_P_FEM_3D))
            allocate(weight(Num_Gauss_P_FEM_3D))
            kesi(1:Num_Gauss_P_FEM_3D)   = kesi_No_Enr
            yita(1:Num_Gauss_P_FEM_3D)   = yita_No_Enr
            zeta(1:Num_Gauss_P_FEM_3D)   = zeta_No_Enr
            weight(1:Num_Gauss_P_FEM_3D) = weight_No_Enr
            c_Num_Gauss_Point            = Num_Gauss_P_FEM_3D
        end if
    case(4)
        !To Be Done.
    end select

    num_GP_Elem(i_E) = c_Num_Gauss_Point

    Location_ESM(1:MDOF_3D)   = 0
    num_Loc_ESM           = 0
    if(num_Crack/=0)then
        do i_C =1,num_Crack
            c_POS_3D_c_Ele(1:8) = c_POS_3D(c_NN,i_C)
            if(i_C >1 .and. sum(c_POS_3D_c_Ele(1:8))==0) cycle
            call Location_Element_Stiff_Matrix_3D(i_E,i_C,c_POS_3D_c_Ele(1:8),Location_ESM_C_Crack, &
            num_Loc_ESM_C_Crack,Location_ESM_C_Cr_NoFEM,num_Loc_ESM_C_Cr_NoFEM)
            Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
            num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
        end do
    endif

    do i_G = 1,c_Num_Gauss_Point
        B(1:6,1:MDOF_3D) = ZR
        num_B = 0
        call Cal_detJ_3D(kesi(i_G),yita(i_G),zeta(i_G),c_X_NODES,c_Y_NODES,c_Z_NODES,detJ)
        if(num_Crack/=0)then
            do i_C =1,num_Crack
                c_POS_3D_c_Ele(1:8) = c_POS_3D(c_NN,i_C)
                if(i_C >1 .and. sum(c_POS_3D_c_Ele(1:8))==0) cycle
                call Cal_B_Matrix_Crack_3D(kesi(i_G),yita(i_G),zeta(i_G),i_C,i_E,i_G, &
                c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES,tem_B,num_tem_B)
                B(1:6,num_B+1:num_B+num_tem_B) = tem_B(1:6,1:num_tem_B)
                num_B = num_B + num_tem_B
            end do
        endif
        localK(1:num_B,1:num_B) = localK(1:num_B,1:num_B) + weight(i_G)*detJ* &
        MATMUL(MATMUL(transpose(B(1:6,1:num_B)),c_D),B(1:6,1:num_B))
    end do

    do i_row = 1,num_Loc_ESM
        do i_col = 1,num_Loc_ESM
            c_row = Location_ESM(i_row)
            c_col = Location_ESM(i_col)
            if(Flag_FreeDOF(c_row)==1  .and. Flag_FreeDOF(c_col)==1) then
                real_row =  Location_FreeDOF(c_row)
                real_col =  Location_FreeDOF(c_col)
                l_count =  l_count + 1
                ! encode (row,col) into a single int64 key -> saves memory & keeps exact indices
                key_arr(l_count) = int(real_row-1,8)*NF_8 + int(real_col-1,8)
                val_arr(l_count) = localK(i_row,i_col)
            endif
        end do
    end do
end do
!$OMP END PARALLEL DO

deallocate(num_K_Each_Elem, elem_K_offset)

!*******************************************************
! Step 3: sort COO triplets by encoded key (in place)
!*******************************************************
call quicksort_kv(Temp_K_Each_Count, key_arr, val_arr)

!*******************************************************
! Step 4: merge duplicate (row,col) entries in place,
!         dropping exact zeros
!*******************************************************
M = 0
i = 1
do while (i <= Temp_K_Each_Count)
    accum = val_arr(i)
    j = i
    do while (j < Temp_K_Each_Count)
        if (key_arr(j+1) /= key_arr(i)) exit
        j = j + 1
        accum = accum + val_arr(j)
    end do
    if (accum /= ZR) then
        M = M + 1
        key_arr(M) = key_arr(i)
        val_arr(M) = accum
    end if
    i = j + 1
end do
K_CSR_NNZ = M

!*******************************************************
! Step 5: bounds check + decode key back to (row,col)
!*******************************************************
if(K_CSR_NNZ > K_CSR_NNZ_Max) then
    print *,'    ERROR :: K_CSR_NNZ > K_CSR_NNZ_Max in Assemble_Stiffness_Matrix_SPARS_XFEM.f90!'
    print *,"         K_CSR_NNZ:     ",K_CSR_NNZ
    print *,"         K_CSR_NNZ_Max: ",K_CSR_NNZ_Max
    call Warning_Message('S',Keywords_Blank)
endif

allocate(K_rows(K_CSR_NNZ), K_cols(K_CSR_NNZ), K_vals(K_CSR_NNZ))

do i = 1,K_CSR_NNZ
    K_rows(i) = int(key_arr(i)/NF_8) + 1
    K_cols(i) = int(mod(key_arr(i),NF_8)) + 1
    K_vals(i) = val_arr(i)
end do

deallocate(key_arr,val_arr)

write(*,1112) DBLE(K_CSR_NNZ)/T_Freedom/T_Freedom*100.0D0

!************************
!           Step 6
!************************
call coocsr(num_FreeD,K_CSR_NNZ,K_vals(1:K_CSR_NNZ),K_rows(1:K_CSR_NNZ),K_cols(1:K_CSR_NNZ), &
K_CSR_aa(1:K_CSR_NNZ),K_CSR_ja(1:K_CSR_NNZ),K_CSR_ia(1:num_FreeD+1))

deallocate(K_rows,K_cols,K_vals)

if (Key_CS_Natural_Crack==1) then
    print *, '    Error :: Key_CS_Natural_Crack=1 not avaliable yet when *Key_K_Sparse=1!'
    print *, '    Error :: In Assemble_Stiffness_Matrix_SPARS_XFEM_3D.F90!'
    call Warning_Message('S',Keywords_Blank)
endif

RETURN

contains

    !-----------------------------------------------------------
    ! Iterative quicksort (median-of-three, insertion sort for
    ! small ranges) sorting int64 keys and carrying real values.
    ! O(n log n) time, O(log n) auxiliary stack -> almost no
    ! extra memory.
    !-----------------------------------------------------------
    subroutine quicksort_kv(n, a, b)
        integer,        intent(in)    :: n
        integer(kind=8),intent(inout) :: a(n)
        real(kind=FT),  intent(inout) :: b(n)
        integer, parameter :: MSTACK = 64, MINSIZE = 16
        integer :: lstk(MSTACK), rstk(MSTACK), sp
        integer :: l, r, ii, jj, mid
        integer(kind=8) :: pivot, ta
        real(kind=FT)   :: tb

        if (n < 2) return
        sp = 1; lstk(1) = 1; rstk(1) = n

        do while (sp >= 1)
            l = lstk(sp); r = rstk(sp); sp = sp - 1

            do while (r - l > MINSIZE)
                mid = (l + r)/2
                ! median-of-three: order a(l) <= a(mid) <= a(r)
                if (a(mid) < a(l)) then
                    ta=a(mid);a(mid)=a(l);a(l)=ta; tb=b(mid);b(mid)=b(l);b(l)=tb
                end if
                if (a(r) < a(l)) then
                    ta=a(r);a(r)=a(l);a(l)=ta; tb=b(r);b(r)=b(l);b(l)=tb
                end if
                if (a(r) < a(mid)) then
                    ta=a(r);a(r)=a(mid);a(mid)=ta; tb=b(r);b(r)=b(mid);b(mid)=tb
                end if
                pivot = a(mid)
                ! hide pivot at position r-1
                ta=a(mid);a(mid)=a(r-1);a(r-1)=ta; tb=b(mid);b(mid)=b(r-1);b(r-1)=tb

                ii = l; jj = r - 1
                do
                    do
                        ii = ii + 1
                        if (a(ii) >= pivot) exit
                    end do
                    do
                        jj = jj - 1
                        if (a(jj) <= pivot) exit
                    end do
                    if (ii >= jj) exit
                    ta=a(ii);a(ii)=a(jj);a(jj)=ta; tb=b(ii);b(ii)=b(jj);b(jj)=tb
                end do
                ! restore pivot
                ta=a(ii);a(ii)=a(r-1);a(r-1)=ta; tb=b(ii);b(ii)=b(r-1);b(r-1)=tb

                ! recurse into smaller side, loop on larger side
                if (ii - l > r - ii) then
                    sp = sp + 1; lstk(sp) = l;      rstk(sp) = ii - 1; l = ii + 1
                else
                    sp = sp + 1; lstk(sp) = ii + 1; rstk(sp) = r;      r = ii - 1
                end if
            end do

            ! insertion sort for the small final range
            do ii = l + 1, r
                ta = a(ii); tb = b(ii); jj = ii - 1
                do while (jj >= l)
                    if (a(jj) <= ta) exit
                    a(jj+1) = a(jj); b(jj+1) = b(jj); jj = jj - 1
                end do
                a(jj+1) = ta; b(jj+1) = tb
            end do
        end do
    end subroutine quicksort_kv

END SUBROUTINE Assemble_Stiffness_Matrix_SPARS_XFEM_3D