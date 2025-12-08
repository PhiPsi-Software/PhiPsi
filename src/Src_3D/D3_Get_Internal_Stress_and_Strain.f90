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
 
SUBROUTINE D3_Get_Internal_Stress_and_Strain
! Obtain internal force (stress at each conventional Gauss point and 
! corresponding strain): obtained according to user-defined settings or input files.
! Store variables such as InSitu_Strs_Gaus_xx(num_Elem, Num_Gauss_P_FEM_3D).
! Ref: XSite Theory and Verification Examples_Description_of_Formulation
!      _With_Validation_Problems_Rev1_Red.pdf, Equation (79) and Equation (78).
!
!2022-06-03.
!



!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Stress
use Global_Material
use Global_Filename
use omp_lib

!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
integer i_E
real(kind=FT) Sigma_A,Sigma_B,Sigma_C
real(kind=FT) n_A(3),n_B(3),n_C(3),c_Vector_6(6)
real(kind=FT) c_Matrix_6_3(6,3)
integer mat_num,i_G
real(kind=FT) c_S(6,6)
real(kind=FT) c_El_x,c_El_y,c_El_z
integer c_Range_Num
logical c_Yes_In
real(kind=FT) c_D(6,6)
real(kind=FT) c_T_Alpha,c_TStress(6)
integer c_NN(8)
real(kind=FT),ALLOCATABLE::istn_File_Data(:,:)
LOGICAL ALIVE
character*200 temp_name
integer Tool_Count_Lines
logical Flag_Blank
real(kind=FT) c_Ele_Biot,c_Ele_PoreP


!----------------------------
! Variable Memory Allocation
!----------------------------
if (allocated(InSitu_Strs_Gaus_xx).eqv. .false.) then
  ALLOCATE(InSitu_Strs_Gaus_xx(Num_Elem,Num_Gauss_P_FEM_3D))
endif
if (allocated(InSitu_Strs_Gaus_yy).eqv. .false.) then
  ALLOCATE(InSitu_Strs_Gaus_yy(Num_Elem,Num_Gauss_P_FEM_3D))
endif
if (allocated(InSitu_Strs_Gaus_zz).eqv. .false.) then
  ALLOCATE(InSitu_Strs_Gaus_zz(Num_Elem,Num_Gauss_P_FEM_3D))
endif
if (allocated(InSitu_Strs_Gaus_xy).eqv. .false.) then
  ALLOCATE(InSitu_Strs_Gaus_xy(Num_Elem,Num_Gauss_P_FEM_3D))
endif
if (allocated(InSitu_Strs_Gaus_yz).eqv. .false.) then
  ALLOCATE(InSitu_Strs_Gaus_yz(Num_Elem,Num_Gauss_P_FEM_3D))
endif
if (allocated(InSitu_Strs_Gaus_xz).eqv. .false.) then
  ALLOCATE(InSitu_Strs_Gaus_xz(Num_Elem,Num_Gauss_P_FEM_3D))
endif

InSitu_Strs_Gaus_xx(1:Num_Elem,1:Num_Gauss_P_FEM_3D) = ZR
InSitu_Strs_Gaus_yy(1:Num_Elem,1:Num_Gauss_P_FEM_3D) = ZR
InSitu_Strs_Gaus_zz(1:Num_Elem,1:Num_Gauss_P_FEM_3D) = ZR
InSitu_Strs_Gaus_xy(1:Num_Elem,1:Num_Gauss_P_FEM_3D) = ZR
InSitu_Strs_Gaus_yz(1:Num_Elem,1:Num_Gauss_P_FEM_3D) = ZR
InSitu_Strs_Gaus_xz(1:Num_Elem,1:Num_Gauss_P_FEM_3D) = ZR


!-----------------------------------------------------------------
!
! If XA calls the Fortran library, then jump to NEWFTU2023032402.
!
!-----------------------------------------------------------------
if(Key_Cpp_Call_Fortran_Lib == 1) then
    goto 600
endif

!---------------------------------------------------------
!
! If it is for quickly defining the initial stress field.
!
!---------------------------------------------------------
if(Key_Read_Initial_Node_Stress_File==0)then
    n_A = InSitu_S1_nv_3D
    n_B = InSitu_S2_nv_3D
    n_C = InSitu_S3_nv_3D
    c_Matrix_6_3(1,1:3) = [n_A(1)*n_A(1),n_B(1)*n_B(1),n_C(1)*n_C(1)]
    c_Matrix_6_3(2,1:3) = [n_A(2)*n_A(2),n_B(2)*n_B(2),n_C(2)*n_C(2)]
    c_Matrix_6_3(3,1:3) = [n_A(3)*n_A(3),n_B(3)*n_B(3),n_C(3)*n_C(3)]
    c_Matrix_6_3(4,1:3) = [n_A(1)*n_A(2),n_B(1)*n_B(2),n_C(1)*n_C(2)]
    c_Matrix_6_3(5,1:3) = [n_A(2)*n_A(3),n_B(2)*n_B(3),n_C(2)*n_C(3)]
    c_Matrix_6_3(6,1:3) = [n_A(1)*n_A(3),n_B(1)*n_B(3),n_C(1)*n_C(3)]

    !$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_E,c_El_x,c_El_y,c_El_z, &
    !$OMP              c_Yes_In,c_Range_Num,Sigma_A,Sigma_B,Sigma_C, &
    !$OMP              c_Vector_6)       &
    !$OMP            SCHEDULE(static)
    ! element cycle.
    do i_E = 1,Num_Elem

      !BUGFIX2022070901.
      Sigma_A = InSitu_S1_3D
      Sigma_B = InSitu_S2_3D
      Sigma_C = InSitu_S3_3D

      ! The element's x, y, z coordinates
      c_El_x = Elem_Centroid(i_E,1)
      c_El_y = Elem_Centroid(i_E,2)
      c_El_z = Elem_Centroid(i_E,3)

      !********************************************************************
      ! Rapid definition of non-uniform stress field (x, y, z directions).
      ! NEWFTU2022070601.
      !********************************************************************
      ! The stress in the X direction changes with the Z coordinate
      if(Key_Nonuniform_InSitu_X_with_Z==1)then
          ! Find the Z-coordinate corresponding to the centroid of the element
          call Vector_Range_Number_of_Value_Dou(c_El_z,100,  &
                 InSitu_Sx_3D_Seg_Loca_X_with_Z(1:100),c_Yes_In, &
                 c_Range_Num)
          ! If it cannot be found, an error message will pop up.
          if(c_Yes_In .eqv. .False.)then
              print *, '    Error :: in D3_Get_Internal_Stress_and_Strain!'
              call Warning_Message('S',Keywords_Blank)
          endif
          ! True stress in the x-direction
          Sigma_A = InSitu_Sx_3D_Seg_Strs_X_with_Z(c_Range_Num)
      endif
      ! The stress in the X direction changes with the Y coordinate
      if(Key_Nonuniform_InSitu_X_with_Y==1)then
          ! Find the Y-coordinate corresponding to the centroid of the element
          call Vector_Range_Number_of_Value_Dou(c_El_y,100, &
                InSitu_Sx_3D_Seg_Loca_X_with_Y(1:100),c_Yes_In, &
                c_Range_Num)
          ! If it cannot be found, an error message will pop up.
          if(c_Yes_In .eqv. .False.)then
              print *, '    Error :: in D3_Get_Internal_Stress_and_Strain!'
              call Warning_Message('S',Keywords_Blank)
          endif
          ! True stress in the x-direction
          Sigma_A = InSitu_Sx_3D_Seg_Strs_X_with_Y(c_Range_Num)
      endif

      ! The Y-direction stress varies with the Z coordinate.
      if(Key_Nonuniform_InSitu_Y_with_Z==1)then
          ! Find the Z-coordinate corresponding to the centroid of the element
          call Vector_Range_Number_of_Value_Dou(c_El_z,100, &
                 InSitu_Sy_3D_Seg_Loca_Y_with_Z(1:100),c_Yes_In, &
                 c_Range_Num)
          ! If it cannot be found, an error message will pop up.
          if(c_Yes_In .eqv. .False.)then
              print *, '    Error :: in D3_Get_Internal_Stress_and_Strain!'
              call Warning_Message('S',Keywords_Blank)
          endif
          ! True Y-direction stress
          Sigma_B = InSitu_Sy_3D_Seg_Strs_Y_with_Z(c_Range_Num)
      endif
      ! The Y-direction stress varies with the X coordinate.
      if(Key_Nonuniform_InSitu_Y_with_X==1)then
          ! Find the Y-coordinate corresponding to the centroid of the element
          call Vector_Range_Number_of_Value_Dou(c_El_x,100,&
                 InSitu_Sy_3D_Seg_Loca_Y_with_X(1:100),c_Yes_In,&
                 c_Range_Num)
          ! If it cannot be found, an error message will pop up.
          if(c_Yes_In .eqv. .False.)then
              print *, '    Error :: in D3_Get_Internal_Stress_and_Strain!'
              call Warning_Message('S',Keywords_Blank)
          endif
          ! True Y-direction stress
          Sigma_B = InSitu_Sy_3D_Seg_Strs_Y_with_X(c_Range_Num)
      endif

      ! The stress in the Z direction changes with the X coordinate.
      if(Key_Nonuniform_InSitu_Z_with_X==1)then
          ! Find the Z-coordinate corresponding to the centroid of the element
          call Vector_Range_Number_of_Value_Dou(c_El_x,100, &
                 InSitu_Sz_3D_Seg_Loca_Z_with_X(1:100),c_Yes_In,&
                 c_Range_Num)
          ! If it cannot be found, an error message will pop up.
          if(c_Yes_In .eqv. .False.)then
              print *, '    Error :: in D3_Get_Internal_Stress_and_Strain!'
              call Warning_Message('S',Keywords_Blank)
          endif
          ! True x-direction stress
          Sigma_C = InSitu_Sz_3D_Seg_Strs_Z_with_X(c_Range_Num)
      endif
      ! The stress in the Z direction changes with the Y coordinate
      if(Key_Nonuniform_InSitu_Z_with_Y==1)then
          ! Find the Y-coordinate corresponding to the centroid of the element
          call Vector_Range_Number_of_Value_Dou(c_El_y,100, &
                 InSitu_Sz_3D_Seg_Loca_Z_with_Y(1:100),c_Yes_In, &
                 c_Range_Num)
          ! If it cannot be found, an error message will pop up.
          if(c_Yes_In .eqv. .False.)then
              print *, '    Error :: in D3_Get_Internal_Stress_and_Strain!'
              call Warning_Message('S',Keywords_Blank)
          endif
          ! True stress in the x direction
          Sigma_C = InSitu_Sz_3D_Seg_Strs_Z_with_Y(c_Range_Num)
      endif



      c_Vector_6(1:6) = MATMUL(c_Matrix_6_3,[Sigma_A,Sigma_B,Sigma_C])


      InSitu_Strs_Gaus_xx(i_E,1:Num_Gauss_P_FEM_3D) = c_Vector_6(1)
      InSitu_Strs_Gaus_yy(i_E,1:Num_Gauss_P_FEM_3D) = c_Vector_6(2)
      InSitu_Strs_Gaus_zz(i_E,1:Num_Gauss_P_FEM_3D) = c_Vector_6(3)
      InSitu_Strs_Gaus_xy(i_E,1:Num_Gauss_P_FEM_3D) = c_Vector_6(4)
      InSitu_Strs_Gaus_yz(i_E,1:Num_Gauss_P_FEM_3D) = c_Vector_6(5)
      InSitu_Strs_Gaus_xz(i_E,1:Num_Gauss_P_FEM_3D) = c_Vector_6(6)
    enddo
    !$omp end parallel do
endif

!-------------------------------------------------------------------
!
! If reading node stress from *.istn. 2022-12-21. NEWFTU2022122101.
!
!-------------------------------------------------------------------
if(Key_Read_Initial_Node_Stress_File==1)then
    temp_name = trim(trim(Full_Pathname)//'.istn')
    inquire(file=temp_name, exist=alive)
    if(alive.EQV..FALSE.)then
        print *, "    ERROR-2022122101 :: Can not find *.istn files!"
        call Warning_Message('S',Keywords_Blank)
    else
        if (Tool_Count_Lines(temp_name) /=Num_Node) then
          print *, "    ERROR-2022122102 :: illegal *.istn files!"
          print *, "                        Data line number /= Num_Node!"
          call Warning_Message('S',Keywords_Blank)
        endif
        ALLOCATE(istn_File_Data(Num_Node,6))
        Call Tool_Read_File(temp_name,"istn",Num_Node,6,istn_File_Data,Flag_Blank)
    endif

    ! element cycle.
    !$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_E,c_NN,c_Vector_6) &
    !$OMP            SCHEDULE(static)
    do i_E = 1,Num_Elem
        ! Node of the element.
        c_NN    = G_NN(1:8,i_E)
        ! The average initial stress of the 8 nodes of the current element.
        c_Vector_6(1) = sum(istn_File_Data(c_NN,1))/dble(8)
        c_Vector_6(2) = sum(istn_File_Data(c_NN,2))/dble(8)
        c_Vector_6(3) = sum(istn_File_Data(c_NN,3))/dble(8)
        c_Vector_6(4) = sum(istn_File_Data(c_NN,4))/dble(8)
        c_Vector_6(5) = sum(istn_File_Data(c_NN,5))/dble(8)
        c_Vector_6(6) = sum(istn_File_Data(c_NN,6))/dble(8)

        ! Note, InSitu_Strs_Gaus_xx considers compression as positive, 
        ! while the stress in *.istn files considers tension as positive,
        ! so a 'negative sign' needs to be added here.
        InSitu_Strs_Gaus_xx(i_E,1:Num_Gauss_P_FEM_3D) = -c_Vector_6(1)
        InSitu_Strs_Gaus_yy(i_E,1:Num_Gauss_P_FEM_3D) = -c_Vector_6(2)
        InSitu_Strs_Gaus_zz(i_E,1:Num_Gauss_P_FEM_3D) = -c_Vector_6(3)
        InSitu_Strs_Gaus_xy(i_E,1:Num_Gauss_P_FEM_3D) = -c_Vector_6(4)
        InSitu_Strs_Gaus_yz(i_E,1:Num_Gauss_P_FEM_3D) = -c_Vector_6(5)
        InSitu_Strs_Gaus_xz(i_E,1:Num_Gauss_P_FEM_3D) = -c_Vector_6(6)
    enddo
    !$omp end parallel do
    if(allocated(istn_File_Data)) deallocate(istn_File_Data)
endif

600 continue

!----------------------------------------------------
!
! If XA calls the Fortran library. NEWFTU2023032402.
!
!----------------------------------------------------
if(Key_Cpp_Call_Fortran_Lib == 1) then
    !$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_E,&
    !$OMP              Sigma_A,Sigma_B,Sigma_C, n_A,n_B,n_C,&
    !$OMP              c_Matrix_6_3,c_Vector_6)       &
    !$OMP            SCHEDULE(static)
    ! element cycle.
    do i_E = 1,Num_Elem
      Sigma_A = XA_Ele_InSitu_S1_S2_S3(i_E,1)
      Sigma_B = XA_Ele_InSitu_S1_S2_S3(i_E,2)
      Sigma_C = XA_Ele_InSitu_S1_S2_S3(i_E,3)
      
      n_A = XA_Ele_InSitu_S1_Vector(i_E,1:3)
      n_B = XA_Ele_InSitu_S2_Vector(i_E,1:3)
      n_C = XA_Ele_InSitu_S3_Vector(i_E,1:3)
      
      c_Matrix_6_3(1,1:3) = [n_A(1)*n_A(1),n_B(1)*n_B(1),n_C(1)*n_C(1)]
      c_Matrix_6_3(2,1:3) = [n_A(2)*n_A(2),n_B(2)*n_B(2),n_C(2)*n_C(2)]
      c_Matrix_6_3(3,1:3) = [n_A(3)*n_A(3),n_B(3)*n_B(3),n_C(3)*n_C(3)]
      c_Matrix_6_3(4,1:3) = [n_A(1)*n_A(2),n_B(1)*n_B(2),n_C(1)*n_C(2)]
      c_Matrix_6_3(5,1:3) = [n_A(2)*n_A(3),n_B(2)*n_B(3),n_C(2)*n_C(3)]
      c_Matrix_6_3(6,1:3) = [n_A(1)*n_A(3),n_B(1)*n_B(3),n_C(1)*n_C(3)]

      c_Vector_6(1:6) = MATMUL(c_Matrix_6_3,[Sigma_A,Sigma_B,Sigma_C])

      InSitu_Strs_Gaus_xx(i_E,1:Num_Gauss_P_FEM_3D) = c_Vector_6(1)
      InSitu_Strs_Gaus_yy(i_E,1:Num_Gauss_P_FEM_3D) = c_Vector_6(2)
      InSitu_Strs_Gaus_zz(i_E,1:Num_Gauss_P_FEM_3D) = c_Vector_6(3)
      InSitu_Strs_Gaus_xy(i_E,1:Num_Gauss_P_FEM_3D) = c_Vector_6(4)
      InSitu_Strs_Gaus_yz(i_E,1:Num_Gauss_P_FEM_3D) = c_Vector_6(5)
      InSitu_Strs_Gaus_xz(i_E,1:Num_Gauss_P_FEM_3D) = c_Vector_6(6)
    enddo
    !$omp end parallel do
endif

!---------------------------------------------------------------
!
! 3D problem thermal stress load, 2022-10-03. NEWFTU2022100301.
!
!---------------------------------------------------------------
if(Key_Thermal_Stress==1)then
  ! Cycle of each element
  !$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_E,c_D,c_T_Alpha,c_TStress) &
  !$OMP            SCHEDULE(static)
  do i_E = 1,Num_Elem
      if(Key_XA/=2) then
          c_D       = D(Elem_Mat(i_E),1:6,1:6)
          c_T_Alpha = T_Alpha(Elem_Mat(i_E))
      !Key_XA==2. IMPROV2023031905.
      else
          c_D       = Elem_D_XA(i_E,1:6,1:6)
          c_T_Alpha = Elem_TEC_XA(i_E)
      endif
      
      !c_TStress = c_T_Alpha*Thermal_Str_Temper(Elem_Mat(i_E))*MATMUL(c_D,[ONE,ONE,ONE,ZR,ZR,ZR])
      c_TStress = c_T_Alpha*Elem_T_for_Stress(i_E)*MATMUL(c_D,[ONE,ONE,ONE,ZR,ZR,ZR])
      ! Updated ground stress: add thermal stress. Note that compression is positive here.
      InSitu_Strs_Gaus_xx(i_E,1:Num_Gauss_P_FEM_3D) = InSitu_Strs_Gaus_xx(i_E,1:Num_Gauss_P_FEM_3D)+ c_TStress(1)
      InSitu_Strs_Gaus_yy(i_E,1:Num_Gauss_P_FEM_3D) = InSitu_Strs_Gaus_yy(i_E,1:Num_Gauss_P_FEM_3D)+ c_TStress(2)
      InSitu_Strs_Gaus_zz(i_E,1:Num_Gauss_P_FEM_3D) = InSitu_Strs_Gaus_zz(i_E,1:Num_Gauss_P_FEM_3D)+ c_TStress(3)
      InSitu_Strs_Gaus_xy(i_E,1:Num_Gauss_P_FEM_3D) = InSitu_Strs_Gaus_xy(i_E,1:Num_Gauss_P_FEM_3D)+ c_TStress(4)
      InSitu_Strs_Gaus_yz(i_E,1:Num_Gauss_P_FEM_3D) = InSitu_Strs_Gaus_yz(i_E,1:Num_Gauss_P_FEM_3D)+ c_TStress(5)
      InSitu_Strs_Gaus_xz(i_E,1:Num_Gauss_P_FEM_3D) = InSitu_Strs_Gaus_xz(i_E,1:Num_Gauss_P_FEM_3D)+ c_TStress(6)
  enddo
  !$omp end parallel do
endif

!---------------------------------------------------------
!
! 3D Problem Pore Pressure, 2022-10-03. NEWFTU2022100302.
!
!---------------------------------------------------------
if(Key_PoreP==1)then
  ! Cycle of each element
  !$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_E,c_Ele_Biot,c_Ele_PoreP) &
  !$OMP            SCHEDULE(static)
  do i_E = 1,Num_Elem
      ! Updated ground stress: subtract pore pressure. Note that compression is positive here.
      ! IMPROV2023031902 (replacing Initial_PoreP with element values, considering Biot coefficient).
      c_Ele_Biot  = Elem_Biots(i_E) 
      c_Ele_PoreP = Elem_Current_PoreP(i_E)
      InSitu_Strs_Gaus_xx(i_E,1:Num_Gauss_P_FEM_3D) = InSitu_Strs_Gaus_xx(i_E,1:Num_Gauss_P_FEM_3D)- c_Ele_Biot*c_Ele_PoreP
      InSitu_Strs_Gaus_yy(i_E,1:Num_Gauss_P_FEM_3D) = InSitu_Strs_Gaus_yy(i_E,1:Num_Gauss_P_FEM_3D)- c_Ele_Biot*c_Ele_PoreP
      InSitu_Strs_Gaus_zz(i_E,1:Num_Gauss_P_FEM_3D) = InSitu_Strs_Gaus_zz(i_E,1:Num_Gauss_P_FEM_3D)- c_Ele_Biot*c_Ele_PoreP
      InSitu_Strs_Gaus_xy(i_E,1:Num_Gauss_P_FEM_3D) = InSitu_Strs_Gaus_xy(i_E,1:Num_Gauss_P_FEM_3D)- c_Ele_Biot*c_Ele_PoreP
      InSitu_Strs_Gaus_yz(i_E,1:Num_Gauss_P_FEM_3D) = InSitu_Strs_Gaus_yz(i_E,1:Num_Gauss_P_FEM_3D)- c_Ele_Biot*c_Ele_PoreP
      InSitu_Strs_Gaus_xz(i_E,1:Num_Gauss_P_FEM_3D) = InSitu_Strs_Gaus_xz(i_E,1:Num_Gauss_P_FEM_3D)- c_Ele_Biot*c_Ele_PoreP
  enddo
  !$omp end parallel do
endif

!--------------------------------------------
! Obtain the initial strain at Gauss points.
!--------------------------------------------
! Ref: XSite Theory and Validation Examples_Description Of Formulation With Validation
! Problems_Rev1_Red.pdf
! Equation (78)
if (allocated(InSitu_Strain_Gaus_xx).eqv. .false.) then
  ALLOCATE(InSitu_Strain_Gaus_xx(Num_Elem,Num_Gauss_P_FEM_3D))
endif
if (allocated(InSitu_Strain_Gaus_yy).eqv. .false.) then
  ALLOCATE(InSitu_Strain_Gaus_yy(Num_Elem,Num_Gauss_P_FEM_3D))
endif
if (allocated(InSitu_Strain_Gaus_zz).eqv. .false.) then
  ALLOCATE(InSitu_Strain_Gaus_zz(Num_Elem,Num_Gauss_P_FEM_3D))
endif
if (allocated(InSitu_Strain_Gaus_xy).eqv. .false.) then
  ALLOCATE(InSitu_Strain_Gaus_xy(Num_Elem,Num_Gauss_P_FEM_3D))
endif
if (allocated(InSitu_Strain_Gaus_yz).eqv. .false.) then
  ALLOCATE(InSitu_Strain_Gaus_yz(Num_Elem,Num_Gauss_P_FEM_3D))
endif
if (allocated(InSitu_Strain_Gaus_xz).eqv. .false.) then
  ALLOCATE(InSitu_Strain_Gaus_xz(Num_Elem,Num_Gauss_P_FEM_3D))
endif

InSitu_Strain_Gaus_xx(1:Num_Elem,1:Num_Gauss_P_FEM_3D) = ZR
InSitu_Strain_Gaus_yy(1:Num_Elem,1:Num_Gauss_P_FEM_3D) = ZR
InSitu_Strain_Gaus_zz(1:Num_Elem,1:Num_Gauss_P_FEM_3D) = ZR
InSitu_Strain_Gaus_xy(1:Num_Elem,1:Num_Gauss_P_FEM_3D) = ZR
InSitu_Strain_Gaus_yz(1:Num_Elem,1:Num_Gauss_P_FEM_3D) = ZR
InSitu_Strain_Gaus_xz(1:Num_Elem,1:Num_Gauss_P_FEM_3D) = ZR

!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_E,mat_num,c_S,i_G, &
!$OMP               c_Vector_6)        &
!$OMP             SCHEDULE(static)
do i_E = 1,Num_Elem
  if(Key_XA/=2) then
      mat_num = Elem_Mat(i_E)
      c_S = S(mat_num,1:6,1:6)
  !Key_XA==2. IMPROV2023031905.
  else
      ! S matrix
      call Matrix_Inverse(Elem_D_XA(i_E,1:6,1:6),c_S,6)      
  endif

  do i_G = 1,Num_Gauss_P_FEM_3D
      c_Vector_6(1:6) = MATMUL(c_S, &
                           [InSitu_Strs_Gaus_xx(i_E,i_G),&
                            InSitu_Strs_Gaus_yy(i_E,i_G),&
                            InSitu_Strs_Gaus_zz(i_E,i_G),&
                            InSitu_Strs_Gaus_xy(i_E,i_G),&
                            InSitu_Strs_Gaus_yz(i_E,i_G),&
                            InSitu_Strs_Gaus_xz(i_E,i_G)])
    InSitu_Strain_Gaus_xx(i_E,i_G) = c_Vector_6(1)
    InSitu_Strain_Gaus_yy(i_E,i_G) = c_Vector_6(2)
    InSitu_Strain_Gaus_zz(i_E,i_G) = c_Vector_6(3)
    InSitu_Strain_Gaus_xy(i_E,i_G) = c_Vector_6(4)
    InSitu_Strain_Gaus_yz(i_E,i_G) = c_Vector_6(5)
    InSitu_Strain_Gaus_xz(i_E,i_G) = c_Vector_6(6)
  end do
enddo
!$omp end parallel do

RETURN
END SUBROUTINE D3_Get_Internal_Stress_and_Strain
