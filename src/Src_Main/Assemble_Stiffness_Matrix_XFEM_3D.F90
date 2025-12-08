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
 
SUBROUTINE Assemble_Stiffness_Matrix_XFEM_3D(isub,c_globalK,c_Total_Freedom,Total_Num_G_P)
! Assembled stiffness matrix (considering cracks, XFEM).

!**********************
! Read public variable
!**********************
use Global_Float_Type
use Global_Crack_Common
use Global_Crack_3D
use Global_Model
use Global_Filename
use Global_Common
use Global_Material
use Global_HF
use Global_Contact
use Global_Post
use CMD_Progress
use Global_XFEM_Elements
use omp_lib

!----------------------------------------------------------------------------------
! Read subroutine interface module (activate compiler parameter consistency check)
!----------------------------------------------------------------------------------
use Global_Inter_Cal_B_Matrix_Crack_3D
use Global_Inter_Location_Element_Stiff_Matrix_3D

!**********************
! Variable Declaration
!**********************
implicit none
integer,intent(in)::isub,c_Total_Freedom
integer,intent(out)::Total_Num_G_P
real(kind=FT) ,intent(out)::c_globalK(c_Total_Freedom,c_Total_Freedom)
integer i_C,i_E,i_G
real(kind=FT) c_D(6,6)
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
integer c_NN(8)  
real(kind=FT) kesi_Enr(Num_Gau_Points_3D),yita_Enr(Num_Gau_Points_3D),&
              zeta_Enr(Num_Gau_Points_3D), weight_Enr(Num_Gau_Points_3D)
real(kind=FT) kesi_No_Enr(Num_Gauss_P_FEM_3D),yita_No_Enr(Num_Gauss_P_FEM_3D),&
              zeta_No_Enr(Num_Gauss_P_FEM_3D),weight_No_Enr(Num_Gauss_P_FEM_3D)
real(kind=FT),ALLOCATABLE::kesi_Enr_Cubes(:),yita_Enr_Cubes(:),zeta_Enr_Cubes(:),weight_Enr_Cubes(:)
real(kind=FT)    kesi_Enr_MC(Num_Gau_Points_3D_MC), yita_Enr_MC(Num_Gau_Points_3D_MC),& 
                 zeta_Enr_MC(Num_Gau_Points_3D_MC), weight_Enr_MC(Num_Gau_Points_3D_MC)
real(kind=FT) ,ALLOCATABLE::kesi(:),yita(:),zeta(:),weight(:)
integer:: Location_ESM(MDOF_3D)
integer::Location_ESM_C_Crack(MDOF_3D),Location_ESM_C_Cr_NoFEM(MDOF_3D)
integer num_Loc_ESM_C_Cr_NoFEM
real(kind=FT) B(6,MDOF_3D),tem_B(6,MDOF_3D)
integer num_B,num_tem_B
integer num_Loc_ESM,num_Loc_ESM_C_Crack
integer c_Num_Gauss_Point
real(kind=FT) detJ,c_N(3,24)
real(kind=FT) c_G_x,c_G_y,c_G_z
real(kind=FT) Rot_c_D_Comp(6,6),c_D_Comp(6,6),Volume_Ratio
real(kind=FT) T_Matrix(6,6),TT_Matrix(6,6)  
integer mat_num
real(kind=FT) localK(MDOF_3D,MDOF_3D)
integer k
type( CLS_CMD_Progress ) ::Progress
integer Yes_Node_Dealed(num_Node)
integer num_t,num_h,num_j,cnt,tem,i_f
integer c_LocXFEM(3)
integer i,j
real(kind=FT) c_Beta(3),c_Alpha,c_Penalty_CS
integer i_Node
character(200) c_File_name_1
integer c_POS_3D_c_Ele(8)


! Progress Bar Settings
#ifndef Silverfrost
call Progress % Set( N = Num_Elem , L = 25 )
Progress % Prefix = "     K assembling process:  "
Progress % M = "#"
Progress % O = "."
#endif

!************************************************************************************************
! Initialize the stiffness matrix and calculate the local coordinates and weights of the element
! Gauss points
!************************************************************************************************
c_globalK(1:c_Total_Freedom,1:c_Total_Freedom) = ZR
Total_Num_G_P = 0

call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,kesi_No_Enr,yita_No_Enr,zeta_No_Enr,weight_No_Enr)

! If it is a fixed-point calculation.
if (Key_Integral_Sol  == 2) then    
    call Cal_Gauss_Points_3D_8nodes(Num_Gau_Points_3D,kesi_Enr,yita_Enr,zeta_Enr,weight_Enr)
    call Cal_Gauss_Points_3D_8nodes(Num_Gau_Points_3D_MC,kesi_Enr_MC,yita_Enr_MC,zeta_Enr_MC,weight_Enr_MC)
! If it is 3D fixed block integration. 2022-07-29. NEWFTU2022072901.
elseif (Key_Integral_Sol  == 3) then     
    allocate(kesi_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    allocate(yita_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    allocate(zeta_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    allocate(weight_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    call Cal_Gauss_Points_3D_for_SUBCUBES(Num_Sub_3D_Cubes,Num_Gau_Points_3D_Cube,kesi_Enr_Cubes,&
                                      yita_Enr_Cubes,zeta_Enr_Cubes,weight_Enr_Cubes)                                
endif
      
!**************************************************************
! Calculate the stiffness matrix, looping through each element
!**************************************************************
EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gau_Points_3D)= .False.

do i_E = 1,Num_Elem
      mat_num = Elem_Mat(i_E)
      localK(1:MDOF_3D,1:MDOF_3D) = ZR
      c_D(1:6,1:6)     = D(Elem_Mat(i_E),1:6,1:6)     
      !If the mat is composite material.
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

      ! Starting number of Gauss points for each element
      Ele_GP_Start_Num(i_E) = Total_Num_G_P + 1
      
      select case(Key_Integral_Sol)
      !////////////////
      ! Fixed Points /
      !////////////////
      case(2)
        ! For enhancement elements
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
              Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
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
              Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
        end if
                    
        ! For Integration Scheme 2: If the element is related to multiple cracks, more Gauss points are
        ! used. 2022-07-16. NEWFTU2022071604.
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
              Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
        endif
      !///////////////////////////
      ! Fixed-block integration /
      !///////////////////////////
      !2022-07-29. NEWFTU2022072901.
      case(3)
        ! For enhancement elements
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
              Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
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
              Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
        end if
        
      !////////////////////////////////////////////////
      ! Integration by parts (temporarily abandoned) /
      !////////////////////////////////////////////////
      !2022-07-27.
      case(4)
        !To Be Done.
      end select      
      
      
      
      ! Save the number of Gauss integration points for each element to a global variable for
      ! post-processing. 2022-07-16.
      if(Key_Save_Nothing==0) then
          if(Key_Post_Elements_Gauss_Num==1) then
              Elements_Gauss_Num(i_E)      =  c_Num_Gauss_Point
          endif
      endif
      
      
      ! Save the number of Gauss points for each element to a global variable
      num_GP_Elem(i_E) = c_Num_Gauss_Point
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !Decide the location of each element stiffness 
      !matrix in the global stiffness matrix.   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Location_ESM(1:MDOF_3D)   = 0
      num_Loc_ESM           = 0
      if(num_Crack/=0)then
        do i_C =1,num_Crack
        
          c_POS_3D_c_Ele(1:8) = c_POS_3D(c_NN,i_C)
          if(i_C >1 .and. sum(c_POS_3D_c_Ele(1:8))==0) cycle
          
          call Location_Element_Stiff_Matrix_3D(i_E,i_C,c_POS_3D_c_Ele(1:8),Location_ESM_C_Crack, &
                                                num_Loc_ESM_C_Crack,Location_ESM_C_Cr_NoFEM,num_Loc_ESM_C_Cr_NoFEM)
          ! Includes FEM degrees of freedom
          Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
          num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
        end do
      endif
      
      do i_G = 1,c_Num_Gauss_Point
          B(1:6,1:MDOF_3D) = ZR
          num_B = 0
          
          !call Cal_N_3D(kesi(i_G),yita(i_G),zeta(i_G),c_N)
          call Cal_detJ_3D(kesi(i_G),yita(i_G),zeta(i_G),c_X_NODES,c_Y_NODES,c_Z_NODES,detJ)  
          
          !Calculate the B Matrix, Loop through each crack.
          if(num_Crack/=0)then
            do i_C =1,num_Crack 
            
              c_POS_3D_c_Ele(1:8) = c_POS_3D(c_NN,i_C)
              if(i_C >1 .and. sum(c_POS_3D_c_Ele(1:8))==0) cycle
              
              call Cal_B_Matrix_Crack_3D(kesi(i_G),yita(i_G),zeta(i_G),i_C,i_E,i_G,  & 
                             c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES,tem_B,num_tem_B)
              !tem_B=24
              !num_tem_B(1:)
              B(1:6,num_B+1:num_B+num_tem_B) = tem_B(1:6,1:num_tem_B)
              num_B = num_B + num_tem_B
            end do
          endif
          localK(1:num_B,1:num_B) = localK(1:num_B,1:num_B) + weight(i_G)*detJ*   &
                           MATMUL(MATMUL(transpose(B(1:6,1:num_B)),c_D),B(1:6,1:num_B))
      end do

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !Assemble the global stiffness matrix. 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !num_Loc_ESM =ttt
      c_globalK(Location_ESM(1:num_Loc_ESM),Location_ESM(1:num_Loc_ESM)) = &
            c_globalK(Location_ESM(1:num_Loc_ESM),Location_ESM(1:num_Loc_ESM)) + localK(1:num_B,1:num_B)
  
  ! Percentage of Feedback Completed
#ifndef Silverfrost
  if(i_E==1.or.i_E==5.or.i_E==10.or.i_E==20.or.i_E==30)then
      call Progress % Put(i_E, CMD_PROGRESS_ABSOLUTE )
  endif  
  if(Mod(i_E,50)==0 .and. (i_E/=Num_Elem))then
      call Progress % Put(i_E, CMD_PROGRESS_ABSOLUTE )
  endif    
  if(i_E==Num_Elem)then
      call Progress % Put(i_E, CMD_PROGRESS_ABSOLUTE )
  endif  
#endif

end do  

!----------------------------------------------------------------------------------------------------
! Penalty function method for handling compression-shear fracture surfaces, 2022-05-12.
! NEWFTU2022051201.
! Ref: \theory_documents\034 Efficient Algorithm for 3D Compression-Shear Fracture Treatment Based
! on Penalty Function Method - Faculty Visit - 2022-05-11.pdf
! Ref: \theory_documents\033 Three-Point Constraint Penalty Function Method - Equations 9.8-9.9-9.10
! - 2022-05-11.pdf
!
! Note, the types of enhanced nodes are not distinguished here, which may cause issues. For example,
! there can be multiple tip-enriched elements. Refer to EBE_XFEM_PCG_3D_Get_K.f90. 2023-08-27.
!
!----------------------------------------------------------------------------------------------------
! Penalty_CS = 1.0D20    !Penalty parameter (global variable)
if (Key_CS_Natural_Crack==1) then
    Yes_Node_Dealed(1:num_Node) =0
    c_Alpha = ZR
    !=====OPTION 1=====
    ! c_Penalty_CS = maxval(abs(c_globalK(:,:))) * 1.0D4 !Penalty parameter, according to Equation 2-86
    ! in my book
    
    !=====OPTION 2=====
    ! Directly given. 2022-12-23. IMPROV2022122301.
    c_Penalty_CS =  Penalty_CS_Natural_Crack  
    
    ! Each fracture loop
    do i_C=1,num_Crack
    ! If the crack is marked as a shear crack
    if(Key_CS_Crack(i_C)==1)then
    ! Loop through each node
    do i_Node =1,Num_Node
        !if(i_Node==2199 .or.i_Node==2200) then     
        ! If it's an enhanced node
        if (Enriched_Node_Type_3D(i_Node,i_C)/=0)then
        ! if (Enriched_Node_Type_3D(i_Node,i_C)==2) then ! if it is a Heaviside enriched node
          ! Obtain the numbers (global numbers) corresponding to the 3 degrees of freedom of the enhanced node
          c_LocXFEM(1) = 3*c_POS_3D(i_Node,i_C) - 2
          c_LocXFEM(2) = 3*c_POS_3D(i_Node,i_C) - 1
          c_LocXFEM(3) = 3*c_POS_3D(i_Node,i_C)    
          
          ! If the current augmented node is not processed by the penalty function
          if (Yes_Node_Dealed(i_Node)==0) then
            Yes_Node_Dealed(i_Node)=1
            !---------------------------------------------------------------------------------------
            ! Obtain the outward normal vector of the fracture surface corresponding to the current
            ! enhanced node
            !---------------------------------------------------------------------------------------
            !c_Beta(1:3)=Enriched_Node_Crack_n_Vector_3D(i_Node,i_C,1:3)
            c_Beta(1:3)=Enriched_Node_Crack_n_Vector_3D(i_Node)%row(i_C,1:3)
            
            !-----------------------------------------------------------------------------
            ! Penalty Function Correction (Corrected Load Array)
            ! Ref: \theory_documents\033 Three-Point Constraint Penalty Function Method -
            ! Formula 9.8-9.9-9.10-2022-05-11.pdf
            !-----------------------------------------------------------------------------
            if(Key_Penalty_CS_Method==1) then
                do i=1,3
                    do j=1,3                
                    c_globalK(c_LocXFEM(i),c_LocXFEM(j))= c_globalK(c_LocXFEM(i),c_LocXFEM(j)) + &
                               c_Penalty_CS*c_Beta(i)*c_Beta(j)         
                    enddo
                enddo 
            elseif(Key_Penalty_CS_Method==2) then
                do i=3,3
                    do j=3,3                
                    c_globalK(c_LocXFEM(i),c_LocXFEM(j))= c_globalK(c_LocXFEM(i),c_LocXFEM(j)) + &
                               c_Penalty_CS*c_Beta(i)*c_Beta(j)         
                    enddo
                enddo 
            endif
          endif
        endif
        !endif
    enddo
    endif
    enddo
endif
   
RETURN
END SUBROUTINE Assemble_Stiffness_Matrix_XFEM_3D
