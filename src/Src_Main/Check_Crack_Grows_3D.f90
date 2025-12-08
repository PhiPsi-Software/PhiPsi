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
 
subroutine Check_Crack_Grows_3D(ifra,iter,Yes_Grow)
! Check whether the crack has propagated according to the weighted average maximum principal stress
! criterion, and determine the crack propagation direction as well as the new discrete crack nodes
! and elements.

!***************************
! SECTION 0: Public Module.
!***************************
use Global_Float_Type
use Global_Crack_Common
use Global_Crack_3D
use Global_DISP
use Global_Elem_Area_Vol
use Global_Model
use Global_HF
use Global_Common
use Global_Material
use Global_Filename
use Global_POST
use omp_lib
use Global_Ragged_Array_Real_Classs
use Global_Cal_Ele_Num_by_Coors_3D  

!**********************************
! SECTION 1: Variable Declaration.
!**********************************
implicit none
integer,intent(in)::ifra,iter
logical,intent(out)::Yes_Grow(Max_Num_Cr_3D)
integer i_C,j_C,i_Tip,Tip_Elem,mat_num,Old_num_Poi
real(kind=FT) pi_180,pi_2
integer Num_CrMesh_Outlines,i_Out_Node,c_Mesh_Node
integer in_Ele_Num
real(kind=FT) c_x,c_y,c_z
integer node_of_Ele(8)
integer all_possi_Ele(200),all_possi_Ele_Uniqued(200)
integer i_Node,total_Possi_Ele,c_num
integer Uniqued_num_Ele,c_Elem,i_Ele
real(kind=FT) c_Sxx_GP,c_Syy_GP,c_Szz_GP
real(kind=FT) c_Sxy_GP,c_Syz_GP,c_Sxz_GP    
real(kind=FT) c_Kesi,c_Yita,c_Zeta
real(kind=FT) Kesi(num_Gauss_S1),Yita(num_Gauss_S1),Zeta(num_Gauss_S1),Gauss_Weight(num_Gauss_S1)
integer i_G
real(kind=FT) R_Factor,c_G_x,c_G_y,c_G_z
real(kind=FT) detJ,c_N(3,24)
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
real(kind=FT) Tool_Function_2Point_Dis_3D
real(kind=FT) c_DIS,Check_R
real(kind=FT) Ave_Sxx_GP,Ave_Syy_GP,Ave_Szz_GP
real(kind=FT) Ave_Sxy_GP,Ave_Syz_GP,Ave_Sxz_GP  
real(kind=FT) S_1,S_2,S_3
real(kind=FT) Vector_S1(3),Vector_S2(3),Vector_S3(3)      
real(kind=FT) All_weight,weight,a_Weight
real(kind=FT) St_Critical
real(kind=FT) poss_Vector(5000,3)
integer num_poss_Vector,i_Find,i_N
real(kind=FT) MAX_DIS,delta_L,c_P(3),Final_Point(5000,3) 
integer new_node_1 ,new_node_2  ,new_node_3
integer old_node_num,old_ele_num
integer in_Elem_num,i_outline
integer c_P1,c_P2
real(kind=FT) c_P1_Coor(3),c_P2_Coor(3),MID_P(3)
real(kind=FT) L_Factor,Check_L 
integer old_Ele,old_Ele_N1,old_Ele_N2,old_Ele_N3
integer i_theta,num_divison
real(kind=FT) a(3),b(3),theta,norm_a,norm_b
integer i_Crack_Ele
integer Crack_Node1,Crack_Node2,Crack_Node3
real(kind=FT) Point1(3),Point2(3),Point3(3)
real(kind=FT) Point_Factor,centroid_x,centroid_y,centroid_z
real(kind=FT) new_Point(3)
real(kind=FT) Line_AB(2,3), new_Line_AB(2,3)
integer cc_Elem
real(kind=FT) Vector_V1(3),Vector_V2(3),Vector_Normal(3)      
integer i_Crack_Node
integer c_count_Ele
real(kind=FT) c_All_Nor_Vector(3)      
real(kind=FT),ALLOCATABLE:: S1_all(:),S1_Points(:,:)
real(kind=FT),ALLOCATABLE:: St_Critical_all(:)
real(kind=FT),ALLOCATABLE:: Prop_Vector_all(:,:)
real(kind=FT) c_KI_3D,c_KII_3D,c_KIII_3D
integer num_of_check,i_Check_Theta 
real(kind=FT),ALLOCATABLE::Check_Theta(:)
real(kind=FT),ALLOCATABLE::Spe_Pri_Stress(:)      
real(kind=FT) tem_part1,tem_part2,Stress_Theta,Tao_Theta
real(kind=FT) tem_root_2pir,c_Theta 
integer i_Theta_min
real(kind=FT) r_for_c_Theta_finding
real(kind=FT) c_tem_part1,c_tem_part2,c_Stress_Theta,c_Tao_Theta,c_Sai
integer inner_Node
real(kind=FT) c_Vec_x(3)
real(kind=FT) cc_Part1,cc_Part2,cc_Part3
real(kind=FT) c_KIc 
real(kind=FT) Rot_Matrix(3,3)
real(kind=FT) c_poss_Vector(5000,3)
real(kind=FT) n_x,n_y,n_z
integer num_Cr_Nodes
real(kind=FT) c_x_old,c_y_old,c_z_old
real(kind=FT) c_x_stress,c_y_stress,c_z_stress
real(kind=FT),ALLOCATABLE::Sai_Cr(:),Ld_Cr(:),Ks_Cr(:),D_Cr(:) 
real(kind=FT) c_Sai_w,c_Ld_w,c_Ks_w,c_D_w
integer i_FluidEl
real(kind=FT) kw(2,2),fw(2)
real(kind=FT),ALLOCATABLE::Global_kw(:,:),Global_fw(:),Global_w(:)
real(kind=FT),ALLOCATABLE::KI_eq(:),Theta_All(:)
real(kind=FT) Sai_Cr_1,Sai_Cr_2
integer c_Mesh_Node_Next
real(kind=FT) Penalty_Value
real(kind=FT) c_x_check,c_y_check,c_z_check
integer c_OUT_Elem
real(kind=FT) cc_kesi,cc_yita,cc_zeta
real(kind=FT) delta_L_Stoped,c_Vector_x(3)
character(200) c_File_name_test
integer in_Ele_Num_new
real(kind=FT) Big_Value,Small_Value,Eq_a,Eq_b
real(kind=FT) Growth_Factor
real(kind=FT) c_Sum,num_c_Count,c_Ave_S_1
character(5) temp
integer Pre_Mesh_Node1,Pre_Mesh_Node2,All_Out
integer Nex_Mesh_Node1,Nex_Mesh_Node2
integer Elem_Pre1,Elem_Pre2,Elem_Nex1,Elem_Nex2
real(kind=FT) c_x_Pre1,c_y_Pre1,c_z_Pre1,c_x_Pre2,c_y_Pre2
real(kind=FT) c_z_Pre2,c_x_Nex1,c_y_Nex1,c_z_Nex1,c_x_Nex2,c_y_Nex2,c_z_Nex2
real(kind=FT) c_Prop_Vector(3)
integer c_OUT_Elem_Old(5000)
real(kind=FT) c_s_x,c_s_y,c_s_z
integer i_ball_theta,i_ball_phi
real(kind=FT) c_ball_theta,c_ball_phi
integer i_v,n_Sigma
logical Flag_Propagation_this_Crack
real(kind=FT) All_Growth_Factor(5000)
integer Max_GF_iC,Max_GF_iV
real(kind=FT) Schollm_Max_Theta_in_pi
real(kind=FT),ALLOCATABLE::Old_Vertex_Points(:,:)
logical,ALLOCATABLE::Flag_Ver_Grow(:)
integer,ALLOCATABLE::New_Ver_Node_Num(:)
real(kind=FT) c_Grow_Dis
integer c_Grow_Count
logical Curr_Flag,Next_Flag
real(kind=FT) Cros_Product_Vector_Old(3)
real(kind=FT) Cros_Product_Vector_New(3)
real(kind=FT) c_angle,c_Growth_Factor
integer Ele_Num_Cache
integer num_Crack_Do,Crack_List(num_Crack),c_Crack
integer c_Outlines_num
integer c_P1_i,c_P1_j,c_P2_i,c_P2_j
integer j_outline
real(kind=FT) c_P1_Coor_i(3),c_P1_Coor_j(3),c_P2_Coor_i(3),c_P2_Coor_j(3)
logical c_Yes_Inter
real(kind=FT) c_InterSection_P(3)
integer c_Max_N_Node_3D
real(kind=FT),ALLOCATABLE::Points_to_Smooth(:,:),Points_Smoothed(:,:)
real(kind=FT),ALLOCATABLE::Pre_Fronts(:,:)
integer Natural_Crack,Num_List_Elements
integer c_Ele_location
logical c_Yes_Ele_In
integer Num_Suspended_Cracks

!************************************************************************************************
!SECTION 2: 
! Open the file to store crack tip stress calculation points (and the calculation sphere radius;
! vspt file) and Gaussian points within the calculation sphere (vspg file)
! Maximum principal stress S1 at crack tip Vertex (cvps file, 2022-04-24)
! Equivalent stress intensity factor Keq of crack tip Vertex (cvke file, 2022-04-24)
!************************************************************************************************
write(temp,'(I5)') iter
if(Key_Save_Nothing  /= 1) then 
    if(CFCP==2)then
        c_File_name_test   = trim(Full_Pathname)//'.vspt_'//ADJUSTL(temp)  
        open(101,file=c_File_name_test,status='unknown')  
        c_File_name_test   = trim(Full_Pathname)//'.vspg_'//ADJUSTL(temp)  
        open(102,file=c_File_name_test,status='unknown')       
        c_File_name_test   = trim(Full_Pathname)//'.cvps_'//ADJUSTL(temp)  
        open(103,file=c_File_name_test,status='unknown')     
        c_File_name_test   = trim(Full_Pathname)//'.cvpr_'//ADJUSTL(temp)
        open(104,file=c_File_name_test,status='unknown')               
    elseif(CFCP==3)then
        c_File_name_test   = trim(Full_Pathname)//'.cvke_'//ADJUSTL(temp)  
        open(104,file=c_File_name_test,status='unknown')             
    endif
endif      

!*********************************
! SECTION 4: Interactive Display.
!*********************************
Num_Suspended_Cracks = count(Crack_Type_Status_3D(1:num_Crack,3)==0)
if(Num_Suspended_Cracks>=2) then
    print *, '    Number of suspended cracks:',  Num_Suspended_Cracks
endif
print *,'    Checking propagation of each crack front vertex...'  
1031 FORMAT(5X,'-- S1 of vertex ',I5,' of crack',I3, ' is', F16.6,' MPa')  
1032 FORMAT(5X,'-- Twist angle of vertex ',I5,' of crack',I3, ' is', F9.3,' degress')   
1035 FORMAT(5X,'-- Max, min, average S1 at vertexs of crack',I4,' are',3E16.5,' MPa')     
1131 FORMAT(5X,'-- S1 of vertex ',I5,' of crack',I3, ' is', F16.6,' MPa;','   Growth factor:',F6.2)       
1061 FORMAT(5X,'-- Warning :: vertex ',I5,' of crack ',I4, ' is going beyong the model, crack ',I3, ' stop propagting!')       
1071 FORMAT(5X,'Warning :: vertex ',I6,' of crack ',I4, ' is going outside the model!')   
1062 FORMAT(5X,'-- Attention :: crack ',I4,' will not grow!')      

 
!*********************************
! SECTION 5: Data Initialization.
!*********************************
Yes_Grow(1:Max_Num_Cr_3D) = .False.
!Crack3D_Vector_S1(1:Max_Num_Cr_3D,1:Max_N_Node_3D,1:3) = ZR      
!tem_save_1(1:100,1:100,1:3) = ZR


!**********************************
! SECTION 6: Variable Preparation.
!**********************************
pi_180    = pi/Con_180
pi_2      = TWO*pi
! Point_Factor = 0.2D0 ! Stress calculation point movement factor (moves outward from the crack
! centroid; key parameter, do not modify lightly) (default value: 0.2);
                                  ! Has been replaced by the global variable S1_Point_Factor, 2021-02-07
! R_Factor = 1.0D0                ! Search sphere radius factor (default is 1.0)
! L_Factor = 1.5D0 !Boundary splitting factor, 1.5 (replaced by global variable Prp_Bisec_Factor_3D,
! 2021-02-07)

Check_R   = Factor_Ave_R*Ave_Elem_L_Enrich
Check_L   = Prp_Bisec_Factor_3D*Ave_Elem_L_Enrich
delta_L   = Factor_Propagation*Ave_Elem_L_Enrich

! If local encryption is used, the cell feature length before local encryption is adopted, because
! crack propagation still occurs according to the scale of the larger grid before encryption
! (2021-08-22).
if (Key_Local_Mesh_Refine>=1) then
  Check_L = Prp_Bisec_Factor_3D*Ave_Elem_L_Enrich_Unlocalrefined
  delta_L = Factor_Propagation*Ave_Elem_L_Enrich_Unlocalrefined
endif

! If a fixed propagation step length is specified through the keyword *Propagation_Length
! (2021-08-03), it overrides the original delta_L.
if (Propagation_Length> Tol_10)then   
  delta_L = Propagation_Length
  !Check Propagation_Length, 2021-08-16
  if (delta_L < 0.5D0*Ave_Elem_L_Enrich)then
      print *, '    ERROR :: *Propagation_Length is too small!'
      print *, '             Please check *Propagation_Length!'     
      print *, '             Error in Check_Crack_Grows_3D.f'   
      call Warning_Message('S',Keywords_Blank)
  endif
  if (delta_L > 10.0D0*Ave_Elem_L_Enrich)then
      print *, '    ERROR :: *Propagation_Length is too large!'
      print *, '             Please check *Propagation_Length!'     
      print *, '             Error in Check_Crack_Grows_3D.f'   
      call Warning_Message('S',Keywords_Blank)
  endif         
endif  
    


! Given the virtual step length without expanding vertices (provided differently according to the
! crack surface update algorithm. Optimized on 2022-07-12. NEWFTU2022071201).
if (Key_3D_Cr_Update_Scheme == 1)then
  ! delta_L_Stoped = 0.01D0*Ave_Elem_L_Enrich ! Crack propagation step length (virtual step length
  ! when not extending the vertex).
  ! delta_L_Stoped = 0.3D0*Ave_Elem_L_Enrich ! Crack propagation step (virtual step without extending
  ! vertices).
  delta_L_Stoped   = 0.0001D0*Ave_Elem_L_Enrich
elseif(Key_3D_Cr_Update_Scheme == 2)then
  delta_L_Stoped   = ZR         
endif      

! Gaussian points (only calculating the conventional FEM Gaussian points), used to compute the
! weighted average of the maximum principal stress at the Gaussian points within the stress sphere
! at the crack front.
call Cal_Gauss_Points_3D_8nodes(num_Gauss_S1,Kesi,Yita,Zeta,Gauss_Weight)  

! Max_Growth_Factor = -Con_Big_15 !For testing only
      
      
!****************************************************
!                                                  *
!                                                  *
! SECTION 7: Crack Cycle.                          *
!                                                  *
!                                                  *
!****************************************************
do i_C =1,num_Crack 
    ! Current cracks discrete fracture surface triangle node count c_Max_N_Node_3D.
    c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
    Crack3D_Vector_S1(i_C)%row(1:c_Max_N_Node_3D,1:3) = ZR   
    ! If the current crack is not allowed to propagate, skip it directly.
    !if(Cracks_Allow_Propa(i_C) == 0)then  
    if(Cracks_Allow_Propa(i_C) == 0 .or. Crack_Type_Status_3D(i_C,3)==0)then
      goto 100
    endif
    ! Temporarily mark this crack as not meeting the crack propagation criteria, 2022-05-08.
    Flag_Propagation_this_Crack = .False.
    ! Mark that the crack did not propagate in the previous step for subsequent analysis.
    Crack_Type_Status_3D(i_C,5) = 0   
    Num_CrMesh_Outlines = Crack3D_Meshed_Outline_num(i_C)
    ! Variable Initialization
    c_OUT_Elem_Old(1:5000)    = 1
    Final_Point(1:5000,1:3)   = ZR 
    ! The boundary point is equal to the original boundary point (to prevent some boundary points from
    ! not meeting the expansion conditions and being expanded to the origin), 2019-10-04
    do i_Out_Node = 1,Num_CrMesh_Outlines
      c_Mesh_Node = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
      ! Deviate from the previous point.
      c_x_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1) 
      c_y_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,2) 
      c_z_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,3)
      Final_Point(i_Out_Node,1)=c_x_old
      Final_Point(i_Out_Node,2)=c_y_old
      Final_Point(i_Out_Node,3)=c_z_old
    enddo
    
    ! Select according to the crack propagation criteria.
    select case(CFCP)
    
    !/////////////////////////////////////////////////////////////////////////
    !/                                                            /
    !/                                                            /
    !/                                                            /
    ! / CFCP=2, Weighted Average Maximum Principal Tensile Stress Criterion /
    !/                                                            /
    !/                                                            /
    !/                                                            /
    !/////////////////////////////////////////////////////////////////////////
    case(2)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! STEP 1, Loop through the discrete points along the crack surface boundary to calculate the
    ! magnitude and direction of the principal stress.
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(allocated(S1_all))deallocate(S1_all)
    ALLOCATE(S1_all(Num_CrMesh_Outlines))
    if(allocated(Prop_Vector_all))deallocate(Prop_Vector_all)
    ALLOCATE(Prop_Vector_all(Num_CrMesh_Outlines,3))
    if(allocated(S1_Points))deallocate(S1_Points)
    ALLOCATE(S1_Points(Num_CrMesh_Outlines,3))
    if(allocated(St_Critical_all))deallocate(St_Critical_all)
    ALLOCATE(St_Critical_all(Num_CrMesh_Outlines))
    S1_all(1:Num_CrMesh_Outlines) = -Con_Big_20
    St_Critical_all(1:Num_CrMesh_Outlines)   = ZR
    Prop_Vector_all(1:Num_CrMesh_Outlines,1:3) =ZR
    S1_Points(1:Num_CrMesh_Outlines,1:3)=ZR
    
    do i_Out_Node = 1,Num_CrMesh_Outlines
      Ele_Num_Cache = 1
      c_Mesh_Node = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
      in_Ele_Num  = Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(c_Mesh_Node)
      ! Deviate from the previous point.
      c_x_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1) 
      c_y_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,2) 
      c_z_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,3)
      ! Calculate the points after deviation.
      !....................................................................
      ! option 1, deviating outward from the centroid of the crack surface
      !....................................................................
      !Line_AB(1,1:3) = Crack3D_Centroid(i_C,1:3)
      !Line_AB(2,1:3) = [c_x_old,c_y_old,c_z_old]
      !call Tool_Shorten_or_Extend_Line_3D(Line_AB,Point_Factor*Ave_Elem_L,'B',new_Line_AB,new_Point)
 
      !.............................................................................................
      ! option 2, deviation of the x-axis in the local coordinate system at the boundary node along
      ! the crack surface
      !.............................................................................................
      c_Vec_x = Crack3D_Meshed_Vertex_x_Vector(i_C)%row(i_Out_Node,1:3)
      Line_AB(1,1:3) = [c_x_old,c_y_old,c_z_old]
      Line_AB(2,1:3) = [c_x_old+5.0D0*Ave_Elem_L*c_Vec_x(1), & 
                        c_y_old+5.0D0*Ave_Elem_L*c_Vec_x(2), & 
                        c_z_old+5.0D0*Ave_Elem_L*c_Vec_x(3)]
      call Tool_Shorten_or_Extend_Line_3D(Line_AB,-S1_Point_Factor*Ave_Elem_L,'A',new_Line_AB,new_Point)
      ! Coordinates of points used for calculating the stress field
      c_x_stress =new_Point(1)
      c_y_stress =new_Point(2)
      c_z_stress =new_Point(3)
      
      ! Save the coordinates of the S1 stress calculation points for post-processing
      S1_Points(i_Out_Node,1:3) = new_Point(1:3)
      ! Element number of the point used to calculate the stress field in in_Ele_Num_new, 2020-03-11
      call Cal_Ele_Num_by_Coors_3D(c_x_stress,c_y_stress,c_z_stress,Ele_Num_Cache,in_Ele_Num_new)
      ! If the discrete point is inside the model
      if (in_Ele_Num_new /=0) then
          node_of_Ele(1:8) = Elem_Node(in_Ele_Num_new,1:8)
          !ooooooooooooooooooooooooooooooooooooooooooooo
          ! Weighted - Standard Algorithm in Literature
          !ooooooooooooooooooooooooooooooooooooooooooooo
          if(Key_Ave_Stress == 2)then
              !TO-BE-DONE
          !ooooooooooooooooooooooooooooooooooooooooo
          ! Algorithms in Dr. Shi's Weighted Thesis
          !ooooooooooooooooooooooooooooooooooooooooo
          elseif(Key_Ave_Stress == 1)then
              !#####################################################################################
              ! Find all the elements connected to the nodes of the element containing the discrete
              ! point
              !#####################################################################################
              all_possi_Ele(1:200) = 0
              total_Possi_Ele = 0
              do i_Node =1,8
                  c_num = num_Node_Elements(node_of_Ele(i_Node))
                  !all_possi_Ele(total_Possi_Ele+1:total_Possi_Ele+c_num)= Node_Elements(node_of_Ele(i_Node),1:c_num)
                  all_possi_Ele(total_Possi_Ele+1:total_Possi_Ele+c_num)= Node_Elements_3D(node_of_Ele(i_Node))%row(1:c_num)
                  total_Possi_Ele = total_Possi_Ele + c_num  
              enddo
              ! Remove duplicate cells
              call Vector_Unique_Int(total_Possi_Ele,total_Possi_Ele,all_possi_Ele(1:total_Possi_Ele), & 
                   all_possi_Ele_Uniqued(1:total_Possi_Ele),Uniqued_num_Ele)  
              !##################################################
              ! Loop through each Gaussian point of each element
              !##################################################
              !count_Gauss = 0
              All_weight = ZR
              Ave_Sxx_GP  = ZR
              Ave_Syy_GP  = ZR
              Ave_Szz_GP  = ZR
              Ave_Sxy_GP  = ZR
              Ave_Syz_GP  = ZR
              Ave_Sxz_GP  = ZR    
              do i_Ele = 1,Uniqued_num_Ele
                c_Elem = all_possi_Ele_Uniqued(i_Ele)
                c_X_NODES = G_X_NODES(1:8,c_Elem)
                c_Y_NODES = G_Y_NODES(1:8,c_Elem)  
                c_Z_NODES = G_Z_NODES(1:8,c_Elem)                     
                do i_G = 1,num_Gauss_S1
                  ! Calculate the real coordinates of this Gaussian point.
                  call Cal_N_3D(Kesi(i_G),Yita(i_G),Zeta(i_G),c_N)    
                  c_G_x = DOT_PRODUCT(c_N(1,1:24:3),c_X_NODES(1:8))
                  c_G_y = DOT_PRODUCT(c_N(1,1:24:3),c_Y_NODES(1:8))   
                  c_G_z = DOT_PRODUCT(c_N(1,1:24:3),c_Z_NODES(1:8))                   
                  ! Calculate the distance from the Gauss points to the discrete points on the fracture surface.
                  c_DIS = Tool_Function_2Point_Dis_3D([c_x_stress,c_y_stress,c_z_stress],[c_G_x,c_G_y,c_G_z])
                  if(c_DIS <= Check_R)then
                    ! Store the Gauss points inside the computational sphere (for Matlab post-processing display only).
                    !Notice: this writing operation has no effect on CPU comsumption.
                    if(Key_Save_Nothing  /= 1) write(102, '(3E20.12)') c_G_x,c_G_y,c_G_z                         
                    !count_Gauss = count_Gauss +1
                    weight = (ONE-((c_DIS/Check_R)**a_Ave_Shi))**3
                    !weight = 1.0
                    All_weight = All_weight + weight
                    ! Calculate the stress at this Gauss point.
                    call Cal_Any_Point_Str_KesiYita_3D(c_Elem,0,1,1,Kesi(i_G),Yita(i_G),Zeta(i_G),i_G,DISP, & 
                          c_Sxx_GP,c_Syy_GP,c_Szz_GP,c_Sxy_GP,c_Syz_GP,c_Sxz_GP)
                    Ave_Sxx_GP  =   Ave_Sxx_GP + c_Sxx_GP*weight
                    Ave_Syy_GP  =   Ave_Syy_GP + c_Syy_GP*weight
                    Ave_Szz_GP  =   Ave_Szz_GP + c_Szz_GP*weight
                    Ave_Sxy_GP  =   Ave_Sxy_GP + c_Sxy_GP*weight
                    Ave_Syz_GP  =   Ave_Syz_GP + c_Syz_GP*weight
                    Ave_Sxz_GP  =   Ave_Sxz_GP + c_Sxz_GP*weight
                  end if               
                enddo
              enddo
             
              Ave_Sxx_GP = Ave_Sxx_GP/All_weight
              Ave_Syy_GP = Ave_Syy_GP/All_weight
              Ave_Szz_GP = Ave_Szz_GP/All_weight
              Ave_Syz_GP = Ave_Syz_GP/All_weight
              Ave_Sxz_GP = Ave_Sxz_GP/All_weight
              Ave_Sxy_GP = Ave_Sxy_GP/All_weight
          !ooooooooooooooooooooooooooooooooooooooooooooooooooooooo
          ! Unweighted, stress is calculated using a single point
          ! It doesn't work well and is not recommended for use.
          !ooooooooooooooooooooooooooooooooooooooooooooooooooooooo
          elseif(Key_Ave_Stress == 3)  then
              ! element number of the single point
              call Cal_Ele_Num_by_Coors_3D(c_x_stress,c_y_stress,c_z_stress,Ele_Num_Cache,cc_Elem)
              if (cc_Elem <=0)then
                  print *,'   Warning :: illegal cc_Elem!'
                  print *,'              in Check_Crack_Grows_3D.f!'                          
              endif
              ! Single-point local coordinates
              call Cal_KesiYita_by_Coor_3D([c_x_stress,c_y_stress,c_z_stress],cc_Elem,cc_kesi,cc_yita,cc_zeta)
              ! Calculate the stress at this Gauss point.
              i_G = 1
              call Cal_Any_Point_Str_KesiYita_3D(cc_Elem,0,1,1,     &
                      cc_kesi,cc_yita,cc_zeta,i_G,DISP,  & 
                      Ave_Sxx_GP,Ave_Syy_GP,Ave_Szz_GP,Ave_Sxy_GP,Ave_Syz_GP,Ave_Sxz_GP)
          !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
          ! Do not use weighting; instead, use several points within the sphere and then take the
          ! average
          ! (2021-08-24)
          ! It doesn't work well and is not recommended for use.
          !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
          elseif(Key_Ave_Stress == 4)  then
              All_weight  = ZR
              Ave_Sxx_GP  = ZR
              Ave_Syy_GP  = ZR
              Ave_Szz_GP  = ZR
              Ave_Sxy_GP  = ZR
              Ave_Syz_GP  = ZR
              Ave_Sxz_GP  = ZR    
              do i_ball_theta = 1,10
                  c_ball_theta = pi/10.0D0*dble((i_ball_theta-1))
                  do i_ball_phi = 1,10
                      c_ball_phi=two*pi/10.0D0*dble((i_ball_phi-1))
                      c_s_x = c_x_stress + Check_R*0.4D0*sin(c_ball_theta)*cos(c_ball_phi)
                      c_s_y = c_y_stress + Check_R*0.4D0*sin(c_ball_theta)*sin(c_ball_phi)
                      c_s_z = c_z_stress + Check_R*0.4D0*cos(c_ball_theta)
                      ! Store the Gauss points inside the computational sphere (for Matlab post-processing display only).
                      !Notice: this writing operation has no effect on CPU comsumption.      
                      if(Key_Save_Nothing  /= 1) write(102, '(3E20.12)') c_s_x,c_s_y,c_s_z
                      ! element number of the single point
                      call Cal_Ele_Num_by_Coors_3D(c_s_x,c_s_y,c_s_z,Ele_Num_Cache,cc_Elem)
                      if (cc_Elem <=0)then
                          print *,'    Error :: illegal cc_Elem!'
                          print *,'      in Check_Crack_Grows_3D.f!' 
                          !call Warning_Message('S',Keywords_Blank)  
                          cycle
                      endif
                      ! Single-point local coordinates
                      call Cal_KesiYita_by_Coor_3D([c_s_x,c_s_y,c_s_z],cc_Elem,cc_kesi,cc_yita,cc_zeta)
                      ! Calculate the stress at this Gauss point.
                      i_G = 1
                      call Cal_Any_Point_Str_KesiYita_3D(cc_Elem,0,  &
                              1,1,cc_kesi,cc_yita,cc_zeta,i_G,DISP, &
                               c_Sxx_GP,c_Syy_GP,c_Szz_GP,c_Sxy_GP,c_Syz_GP,c_Sxz_GP) 
                      weight = ONE
                      All_weight = All_weight + weight
                      Ave_Sxx_GP  =   Ave_Sxx_GP + c_Sxx_GP*weight
                      Ave_Syy_GP  =   Ave_Syy_GP + c_Syy_GP*weight
                      Ave_Szz_GP  =   Ave_Szz_GP + c_Szz_GP*weight
                      Ave_Sxy_GP  =   Ave_Sxy_GP + c_Sxy_GP*weight
                      Ave_Syz_GP  =   Ave_Syz_GP + c_Syz_GP*weight
                      Ave_Sxz_GP  =   Ave_Sxz_GP + c_Sxz_GP*weight   
                  enddo                  
              enddo
              Ave_Sxx_GP = Ave_Sxx_GP/All_weight
              Ave_Syy_GP = Ave_Syy_GP/All_weight
              Ave_Szz_GP = Ave_Szz_GP/All_weight
              Ave_Syz_GP = Ave_Syz_GP/All_weight
              Ave_Sxz_GP = Ave_Sxz_GP/All_weight
              Ave_Sxy_GP = Ave_Sxy_GP/All_weight                  
          endif
          ! Calculate the principal stresses and their directions.
          call Tool_Principal_Stresses_3D(Ave_Sxx_GP,Ave_Syy_GP,Ave_Szz_GP,Ave_Syz_GP,Ave_Sxz_GP,Ave_Sxy_GP, &
                       S_1,S_2,S_3,Vector_S1,Vector_S2,Vector_S3)    
 
          Crack3D_Vector_S1(i_C)%row(c_Mesh_Node,1:3)  = Vector_S1(1:3)     

          
          ! Related to Sinopec.
          if(Key_XA/=2) then
              if(in_Ele_Num ==0) then
                  mat_num = 1   
              else
                  mat_num = Elem_Mat(in_Ele_Num)
              endif
              St_Critical  = Material_Para(mat_num,5)
          !Key_XA==2. IMPROV2023031905.
          else
              if(in_Ele_Num ==0) then
                  St_Critical = Elem_St_XA(1) 
              else
                  St_Critical = Elem_St_XA(in_Ele_Num)
              endif
          endif
          
          St_Critical_all(i_Out_Node) = St_Critical
          ! If any one vertex meets the stress fracture criterion, mark the crack as an extended crack.
          if(S_1>=St_Critical)then
              Flag_Propagation_this_Crack =.True.
              Crack_Type_Status_3D(i_C,5) = 1
          endif
          S1_all(i_Out_Node) = S_1
          
          
          ! Calculate the update direction of possible crack surface discrete nodes, Save to
          ! Prop_Vector_all(i_Out_Node, 1:3)
          ! Discretize the circle and identify the possible directions in which discrete crack edge nodes may
          ! extend.
          ! References: http://blog.sina.com.cn/s/blog_622fbc040102wt9o.html or theoretical_documents
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !  OPTION 1: Brute force search algorithm
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !num_divison = 200
          !num_divison = 360
          !num_divison = 720
          num_divison = 1440
          num_poss_Vector = num_divison
          do i_theta = 1,num_divison
              theta = (i_theta-1)*TWO*pi/num_divison
              call Vector_Cross_Product_3(Vector_S1,[ONE,ZR,ZR],a)                 
              if(sum(abs(a))<=Tol_20) then
                call Vector_Cross_Product_3(Vector_S1,[ZR,ONE,ZR],a)   
              endif
              if(sum(abs(a))<=Tol_20) then
                call Vector_Cross_Product_3(Vector_S1,[ZR,ZR,ONE],a)   
              endif
              ! b = cross(n, a)! Calculate the b vector
              call Vector_Cross_Product_3(Vector_S1,a,b)   
              ! a = a / norm(a)! Normalize the vector a
              ! b = b / norm(b)! Normalizing the b vector
              call Vector_Norm2(3,a,norm_a)   
              call Vector_Norm2(3,b,norm_b)  
              a=a/norm_a
              b=b/norm_b
              poss_Vector(i_theta,1)=a(1)*cos(theta)+b(1)*sin(theta)
              poss_Vector(i_theta,2)=a(2)*cos(theta)+b(2)*sin(theta)
              poss_Vector(i_theta,3)=a(3)*cos(theta)+b(3)*sin(theta)
          enddo
          ! The finally selected point is the one with the smallest sum of distances to all discrete points of
          ! the cracks (i.e., expanding outward).
          MAX_DIS = ZR
          do i_Find=1,num_poss_Vector
            c_P(1) = c_x_old + delta_L*poss_Vector(i_Find,1)
            c_P(2) = c_y_old + delta_L*poss_Vector(i_Find,2)
            c_P(3) = c_z_old + delta_L*poss_Vector(i_Find,3)
            c_DIS  = ZR
            do i_N = 1,Crack3D_Meshed_Node_num(i_C)
                c_DIS = c_DIS+Tool_Function_2Point_Dis_3D(c_P,Crack3D_Meshed_Node(i_C)%row(i_N,1:3))
            enddo
            if (c_DIS>=MAX_DIS)then
                MAX_DIS = c_DIS
                !Final_Point(i_Out_Node,1:3) = c_P(1:3)
                Prop_Vector_all(i_Out_Node,1:3) =  poss_Vector(i_Find,1:3)
            endif
          enddo
          !~~~~~~~~~~~~~~~~~~~~~~~~
          !  OPTION 2: To Be Done.
          !~~~~~~~~~~~~~~~~~~~~~~~~
          !To Be Done.
          
      endif
    end do
     
    ! If it is not a simplified post-processing, save part of the data.
    if (Key_Simple_Post/=1 .and. Key_Save_Nothing  /= 1 ) then
      ! Store the coordinates of the crack tip stress calculation points (for Matlab post-processing
      ! display only)
      do i_Out_Node = 1,Num_CrMesh_Outlines
        write(101, '(4E20.12)') S1_Points(i_Out_Node,1:3),Check_R
      enddo
      ! Save the principal stress S1_all(i_Out_Node) of each vertex of the current crack for Matlab
      ! post-processing, 2022-04-24. NEWFTU2022042401.
      print *,'    Saving cvpr file for crack ',i_C,'...'
      write(104,'(50000E20.12)') (S1_all(i_v),i_v=1,Num_CrMesh_Outlines)        
    endif
    ! Denoising. 2022-04-25. NEWFTU2022042501.
    if (Key_Denoise_Vertex_Value>=1)then
        print *,'    Denoising vertex values of S1...'   
        n_Sigma =2
        call Tool_Denoise_Data(S1_all(1:Num_CrMesh_Outlines),Num_CrMesh_Outlines, &
                    Key_Denoise_Vertex_Value,n_Sigma,.True.)
    endif
    
    ! Data Smoothing Processing. 2022-04-25. NEWFTU2022042501.
    if (Key_Smooth_Vertex_Value>=1)then
        print *,'    Smoothing vertex values of S1...'
        call Tool_Smooth_Data(S1_all(1:Num_CrMesh_Outlines),Num_CrMesh_Outlines, &
                 Key_Smooth_Vertex_Value,Smooth_Vertex_n, .True.)
      ! Data Smoothing (Secondary Processing). 2022-05-25. NEWFTU2022052501.
      if (Key_Smooth_Vertex_Value2>=1)then
          print *,'    Smoothing vertex values of S1...'
          call Tool_Smooth_Data(S1_all(1:Num_CrMesh_Outlines),Num_CrMesh_Outlines, &
                Key_Smooth_Vertex_Value2,Smooth_Vertex_n2, .True.)
      endif           
    endif
    
    ! Save the principal stress S1_all(i_Out_Node) of each vertex of the current crack for Matlab
    ! post-processing, 2022-04-24. NEWFTU2022042401.
    if(Key_Save_Nothing  /= 1) then 
        print *,'    Saving cvps file for crack ',i_C,'...'
        write(103, '(50000E20.12)') (S1_all(i_v),i_v=1,Num_CrMesh_Outlines)
    endif
 
    ! If none of the vertices of the current crack satisfy the stress crack propagation criterion, then
    ! move on to the next crack.
    if(Flag_Propagation_this_Crack .eqv. .False.)then
        write(*,1062) i_C
        goto 100
    endif
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !STEP 2: Determine the propagated fracture front vertex.
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !Loop over i_Out_Node.
    Ele_Num_Cache = 1
    do i_Out_Node = 1,Num_CrMesh_Outlines
      S_1 = S1_all(i_Out_Node)
      St_Critical = St_Critical_all(i_Out_Node)
      !write(*,1031) i_Out_Node,i_C,S_1/1.0D6
      c_Mesh_Node = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
      ! Pre_Mesh_Node1 and Pre_Mesh_Node2
      All_Out = Num_CrMesh_Outlines
      if(i_Out_Node == 1)then
        Pre_Mesh_Node1=Crack3D_Meshed_Outline(i_C)%row(All_Out,1)
        Pre_Mesh_Node2=Crack3D_Meshed_Outline(i_C)%row(All_Out-1,1)
      elseif(i_Out_Node == 2)then
        Pre_Mesh_Node1=Crack3D_Meshed_Outline(i_C)%row(1,1)
        Pre_Mesh_Node2=Crack3D_Meshed_Outline(i_C)%row(All_Out,1)
      else
        Pre_Mesh_Node1=Crack3D_Meshed_Outline(i_C)%row(i_Out_Node-1,1)
        Pre_Mesh_Node2=Crack3D_Meshed_Outline(i_C)%row(i_Out_Node-2,1)
      endif
      if(i_Out_Node == All_Out)then
        Nex_Mesh_Node1=Crack3D_Meshed_Outline(i_C)%row(1,1)
        Nex_Mesh_Node2=Crack3D_Meshed_Outline(i_C)%row(2,1)
      elseif(i_Out_Node == All_Out-1)then
        Nex_Mesh_Node1=Crack3D_Meshed_Outline(i_C)%row(All_Out,1)
        Nex_Mesh_Node2=Crack3D_Meshed_Outline(i_C)%row(1,1)
      else
        Nex_Mesh_Node1=Crack3D_Meshed_Outline(i_C)%row(i_Out_Node+1,1)
        Nex_Mesh_Node2=Crack3D_Meshed_Outline(i_C)%row(i_Out_Node+2,1)
      endif
      
      in_Ele_Num  = Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(c_Mesh_Node)
      ! Deviation from the previous point
      c_x_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1) 
      c_y_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,2) 
      c_z_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,3)
      
      
      c_OUT_Elem_Old(i_Out_Node) = in_Ele_Num
      if(in_Ele_Num ==0) then
          ! Determine whether Pre_Mesh_Node1 and Pre_Mesh_Node2 are inside the model
          c_x_Pre1 = Crack3D_Meshed_Node(i_C)%row(Pre_Mesh_Node1,1) 
          c_y_Pre1 = Crack3D_Meshed_Node(i_C)%row(Pre_Mesh_Node1,2) 
          c_z_Pre1 = Crack3D_Meshed_Node(i_C)%row(Pre_Mesh_Node1,3) 
          c_x_Pre2 = Crack3D_Meshed_Node(i_C)%row(Pre_Mesh_Node2,1) 
          c_y_Pre2 = Crack3D_Meshed_Node(i_C)%row(Pre_Mesh_Node2,2) 
          c_z_Pre2 = Crack3D_Meshed_Node(i_C)%row(Pre_Mesh_Node2,3) 
          c_x_Nex1 = Crack3D_Meshed_Node(i_C)%row(Nex_Mesh_Node1,1) 
          c_y_Nex1 = Crack3D_Meshed_Node(i_C)%row(Nex_Mesh_Node1,2) 
          c_z_Nex1 = Crack3D_Meshed_Node(i_C)%row(Nex_Mesh_Node1,3) 
          c_x_Nex2 = Crack3D_Meshed_Node(i_C)%row(Nex_Mesh_Node2,1) 
          c_y_Nex2 = Crack3D_Meshed_Node(i_C)%row(Nex_Mesh_Node2,2) 
          c_z_Nex2 = Crack3D_Meshed_Node(i_C)%row(Nex_Mesh_Node2,3) 
          call Cal_Ele_Num_by_Coors_3D(c_x_Pre1,c_y_Pre1,c_z_Pre1,Ele_Num_Cache,Elem_Pre1)
          call Cal_Ele_Num_by_Coors_3D(c_x_Pre2,c_y_Pre2,c_z_Pre2,Ele_Num_Cache,Elem_Pre2)
          call Cal_Ele_Num_by_Coors_3D(c_x_Nex1,c_y_Nex1,c_z_Nex1,Ele_Num_Cache,Elem_Nex1)
          call Cal_Ele_Num_by_Coors_3D(c_x_Nex2,c_y_Nex2,c_z_Nex2,Ele_Num_Cache,Elem_Nex2)     
          
          ! If both adjacent boundary points are 0
          if(Elem_Pre1 ==0 .and. Elem_Nex1==0)then
              c_Vector_x=Crack3D_Meshed_Vertex_x_Vector(i_C)%row(i_Out_Node,1:3)
              c_P(1)=c_x_old+delta_L_Stoped*c_Vector_x(1)
              c_P(2)=c_y_old+delta_L_Stoped*c_Vector_x(2)
              c_P(3)=c_z_old+delta_L_Stoped*c_Vector_x(3)    
              Final_Point(i_Out_Node,1:3) = c_P(1:3)
              cycle 
          endif
          ! If the previous adjacent boundary point is inside the model, then the vector in the extension
          ! direction is equal to the extension direction vector of the previous adjacent boundary point.
          if(Elem_Pre1 /=0)then
              if (i_Out_Node==1) then
                  c_Prop_Vector=Prop_Vector_all(Num_CrMesh_Outlines,1:3) 
              else
                  c_Prop_Vector = Prop_Vector_all(i_Out_Node-1,1:3) 
              endif
              c_P(1)=c_x_old+delta_L*c_Prop_Vector(1)
              c_P(2)=c_y_old+delta_L*c_Prop_Vector(2)
              c_P(3)=c_z_old+delta_L*c_Prop_Vector(3)    
              Final_Point(i_Out_Node,1:3) = c_P(1:3)
              cycle 
          endif
          ! If the next adjacent boundary point is inside the model, then the vector in the extension
          ! direction is equal to the extension direction vector of the next adjacent boundary point.
          if(Elem_Nex1 /=0)then
              if (i_Out_Node==Num_CrMesh_Outlines) then
                  c_Prop_Vector = Prop_Vector_all(1,1:3) 
              else
                  c_Prop_Vector = Prop_Vector_all(i_Out_Node+1,1:3) 
              endif
              c_P(1)=c_x_old+delta_L*c_Prop_Vector(1)
              c_P(2)=c_y_old+delta_L*c_Prop_Vector(2)
              c_P(3)=c_z_old+delta_L*c_Prop_Vector(3)    
              Final_Point(i_Out_Node,1:3) = c_P(1:3)
              cycle 
          endif              
      endif
      
      ! If the element containing this vertex is a Junction enhanced element, then extend with a very
      ! small step. BUGFIX2022052701.
      if (Elem_Type_3D(in_Ele_Num,i_C)==4) then
        c_P(1)=c_x_old+delta_L_Stoped*Prop_Vector_all(i_Out_Node,1)
        c_P(2)=c_y_old+delta_L_Stoped*Prop_Vector_all(i_Out_Node,2)
        c_P(3)=c_z_old+delta_L_Stoped*Prop_Vector_all(i_Out_Node,3) 
        goto 737
      endif      
      
      ! If the maximum principal stress exceeds the tensile strength
      if(S_1>=St_Critical)then
        Yes_Grow(i_C) = .TRUE.
        !Every fracture front vertex has the same propagation increament.
        if(Key_Adj_Prp_Step_3D==0)then
          c_P(1) = c_x_old + delta_L*Prop_Vector_all(i_Out_Node,1)
          c_P(2) = c_y_old + delta_L*Prop_Vector_all(i_Out_Node,2)
          c_P(3) = c_z_old + delta_L*Prop_Vector_all(i_Out_Node,3)
        !Fracture front vertex has different propagation increament.
        elseif(Key_Adj_Prp_Step_3D==1)then
          if(S_1/St_Critical<ZR) then
              Growth_Factor = ZR
          else
              Growth_Factor = (S_1/St_Critical)**Prp_3D_Factor_m
          endif
          if(Growth_Factor>=Adj_Prp_Step_3D_Max)then
              Growth_Factor=Adj_Prp_Step_3D_Max
          endif
          write(*,1131) i_Out_Node,i_C,S_1/1.0D6,Growth_Factor
          c_P(1) = c_x_old + delta_L*Growth_Factor*Prop_Vector_all(i_Out_Node,1)
          c_P(2) = c_y_old + delta_L*Growth_Factor*Prop_Vector_all(i_Out_Node,2)
          c_P(3) = c_z_old + delta_L*Growth_Factor*Prop_Vector_all(i_Out_Node,3)
        endif
      else
        c_P(1)=c_x_old+delta_L_Stoped*Prop_Vector_all(i_Out_Node,1)
        c_P(2)=c_y_old+delta_L_Stoped*Prop_Vector_all(i_Out_Node,2)
        c_P(3)=c_z_old+delta_L_Stoped*Prop_Vector_all(i_Out_Node,3) 
      endif  
      
737     continue    
      
      ! If this vertex is a Segmentation suppression vertex (2021-08-20)
      if (key_front_segmentation==1)then
        if (Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,4)==0)then
         c_P(1)=c_x_old+delta_L_Stoped*Prop_Vector_all(i_Out_Node,1)
         c_P(2)=c_y_old+delta_L_Stoped*Prop_Vector_all(i_Out_Node,2)
         c_P(3)=c_z_old+delta_L_Stoped*Prop_Vector_all(i_Out_Node,3)  
         ! Potential suppression point coordinates (in the subsequent D3_Prepare_Front_Segmentation, the
         ! suppression points will be determined based on the coordinates of these points)
         ! It is used to set Crack3D_Meshed_Outline(i_C, i_Out_Node, 4) = 0
         num_Suspended_Point = num_Suspended_Point +1
         Suspended_Points(num_Suspended_Point,1:3)=c_P
        endif
      endif
      
      Final_Point(i_Out_Node,1:3) = c_P
      !---------------------------------------------------------------------------------------------
      ! If it is in-plane growth (Key_InPlane_Growth = 1), then correct the Final_Point back to the
      ! plane of the initial crack.
      ! Project onto the plane, mark the coordinates of the foot of the perpendicular as the new
      ! Final_Point, 2022-04-18, NEWFTU2022041801
      !---------------------------------------------------------------------------------------------
      if(Key_InPlane_Growth== 1) then
        call D3_In_Plane_Growth(Final_Point(i_Out_Node,1:3),i_C)
      endif          
    enddo
    
    
    !FOR TEST ONLY!
    !Final_Point(1:i_Out_Node,3) =10.0d-2
    
    
    ! Display the maximum, minimum, and average principal stresses on the screen
    if(mod(i_Out_Node,1)==0)then
       write(*,1035) i_C, maxval(S1_all,mask=(S1_all>-Con_Big_15))/Cn_M,&
                          minval(S1_all,mask=(S1_all>-Con_Big_15))/Cn_M,&
                          sum(S1_all,mask=(S1_all>-Con_Big_15))/Num_CrMesh_Outlines/Cn_M
    endif
    
    if (allocated(S1_all)) DEALLOCATE(S1_all) 
    if (allocated(Prop_Vector_all)) DEALLOCATE(Prop_Vector_all)
    if (allocated(S1_Points)) DEALLOCATE(S1_Points)  

    !//////////////////////////////////////////////////////////////
    !/                                                            /
    !/                                                            /
    !/                                                            /
    !/                                                            /
    !/     CFCP=3, Schollmanns criterion.                         /
    !/     References: Three-dimensional crack growth with        /
    !/              hp-generalized finite element and             /
    !/              face offsetting methods.                      /
    !/                                                            /
    !/                                                            /
    !/                                                            /
    !//////////////////////////////////////////////////////////////
    case(3)        
        ! Set the Theta angle density for searching special principal stress
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
        
        r_for_c_Theta_finding = 1.0D0*Ave_Elem_L
        tem_root_2pir = sqrt(TWO*pi*r_for_c_Theta_finding)
        
        Ele_Num_Cache = 1
        do i_Out_Node = 1,Num_CrMesh_Outlines
          c_Mesh_Node = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
          c_Mesh_Node_Next =Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,2)
          in_Ele_Num  = Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(c_Mesh_Node)
          
          ! A element number of 0 indicates that the fracture discrete point is outside the model. 2022-06-11.
          if(in_Ele_Num ==0 ) then
              !cycle  !Start the next do i_Out_Node loop.
              cycle
          endif
          
      
          ! Related to Sinopec.
          if(Key_XA/=2) then
              mat_num = Elem_Mat(in_Ele_Num)
              c_KIc   = Material_Para(mat_num,6)
          !Key_XA==2. IMPROV2023031905.
          else
              c_KIc   = Elem_KIc_XA(in_Ele_Num)
          endif
          
          ! Special handling for Key_NaCr_Active_Scheme_3D = 3. NEWFTU2023011204. 2023-01-12.
          if(Key_NaCr_Active_Scheme_3D == 3)then
              ! If the current crack was caused by the activation of a natural fracture.
              if (Crack_Type_Status_3D(i_C,6) >=1) then
                  Natural_Crack = Crack_Type_Status_3D(i_C,6)
                  ! Check whether the element containing the crack tip is in the list of elements with natural
                  ! fractures.
                  Num_List_Elements = NaCr3D_Status(Natural_Crack,2)
                  ! Using Vector_Location_Int(n, Vector, Variable, location, Yes_In)
                  call Vector_Location_Int(Num_List_Elements,&
                       Na_Crack3D_Ele_List(Natural_Crack)%row(1:Num_List_Elements),in_Ele_Num,&
                       c_Ele_location,c_Yes_Ele_In)   
                  ! If the element where the crack tip is located belongs to the list of elements containing natural
                  ! fractures, then the fracture toughness is taken as the fracture toughness of the natural fracture.
                  if(c_Yes_Ele_In .eqv. .True.)then
                      if (Key_Cpp_Call_Fortran_Lib/=1) then
                          c_KIc  = KIc_NaCr(Natural_Crack) 
                      elseif(Key_Cpp_Call_Fortran_Lib==1)then
                          c_KIc  = Na_Crack3D_St(Natural_Crack) 
                      endif
                  endif
              endif
          endif
          
          ! Coordinates of the crack edge nodes.
          c_x  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1) 
          c_y  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,2) 
          c_z  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,3) 
          
          !****************************************************************************************
          ! Check whether the boundary point is outside the model; if it is, the boundary point is
          ! not allowed to expand (expand in very small steps).
          !****************************************************************************************
          call Cal_Ele_Num_by_Coors_3D(c_x,c_y,c_z,Ele_Num_Cache,c_OUT_Elem)
          if(c_OUT_Elem ==0) then
              cycle
          endif
      
          ! Calculate Ld_Cr for the subsequent calculation of the torsion vector w
          Ld_Cr(i_Out_Node) = Tool_Function_2Point_Dis_3D(Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1:3), &
                 Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node_Next,1:3))
  
          ! KS (used to calculate the twist of the crack surface w)
          Ks_Cr(i_Out_Node) = ONE/Ld_Cr(i_Out_Node)
          ! D (used to calculate the twisting amount w of the crack surface)
          D_Cr(i_Out_Node) = Ld_Cr(i_Out_Node)**3
          
          ! Stress Intensity Factor
          !c_KI_3D   = KI_3D(i_C,i_Out_Node)
          !c_KII_3D  = KII_3D(i_C,i_Out_Node)
          !c_KIII_3D = KIII_3D(i_C,i_Out_Node)
          c_KI_3D   = KI_3D(i_C)%row(i_Out_Node)
          c_KII_3D  = KII_3D(i_C)%row(i_Out_Node)
          c_KIII_3D = KIII_3D(i_C)%row(i_Out_Node)              
          
          !.......................................................................
          ! Find c_Theta such that the special principal stress (S1) is maximized
          ! Literature: Three-dimensional crack growth with hp-generalized
          ! Eq. (14) of the finite element and face offsetting methods
          !.......................................................................
          do i_Check_Theta =1,num_of_check+1
            !////////////////////////////////////////////////////////////
            ! OPTION 1: Allowed crack propagation angle is (-pi/2, pi/2)
            !////////////////////////////////////////////////////////////
            !Check_Theta(i_Check_Theta)=-pi/TWO+pi/num_of_check*(i_Check_Theta-1)

            !////////////////////////////////////////////////////////
            ! OPTION 2: Allowed crack propagation angle is (-pi, pi)
            !////////////////////////////////////////////////////////
            !Check_Theta(i_Check_Theta)=-pi+TWO*pi/num_of_check*(i_Check_Theta-1) 

            !/////////////////////////////////////////////////////////////////////////
            ! OPTION 3: The allowable crack propagation angle is (-Schollm_Max_Theta,
            ! Schollm_Max_Theta)
            ! 2022-07-10. The default value of Schollm_Max_Theta is 55 degrees.
            !          NEWFTU2022071001.
            !/////////////////////////////////////////////////////////////////////////
            Schollm_Max_Theta_in_pi = Schollm_Max_Theta*pi/Con_180
            Check_Theta(i_Check_Theta)=-Schollm_Max_Theta_in_pi+TWO*Schollm_Max_Theta_in_pi/num_of_check*(i_Check_Theta-1) 
 
            
            tem_part1 = THR*cos(Check_Theta(i_Check_Theta)/TWO) +cos(THR*Check_Theta(i_Check_Theta)/TWO)
            tem_part2 = THR*sin(Check_Theta(i_Check_Theta)/TWO) +THR*sin(THR*Check_Theta(i_Check_Theta)/TWO)     
            Stress_Theta = c_KI_3D/FOU/tem_root_2pir*tem_part1 -c_KII_3D/FOU/tem_root_2pir*tem_part2
            Tao_Theta=c_KIII_3D*cos(Check_Theta(i_Check_Theta)/TWO)/tem_root_2pir
            Spe_Pri_Stress(i_Check_Theta) = Stress_Theta/TWO + ZP5*sqrt(Stress_Theta**2 + FOU*Tao_Theta**2)
          enddo 
          i_Theta_min = maxloc(Spe_Pri_Stress(1:num_of_check+1),1)
          Theta_All(i_Out_Node) = Check_Theta(i_Theta_min)

          c_Theta = Theta_All(i_Out_Node)
          !............................................................................
          ! Calculate c_Sai when Key_Schollmann_Twist = 1. Used for twist calculation.
          !............................................................................
          if (Key_Schollmann_Twist == 1)then
              c_tem_part1=THR*cos(c_Theta/TWO)+ cos(THR*c_Theta/TWO)
              c_tem_part2=THR*sin(c_Theta/TWO)+THR*sin(THR*c_Theta/TWO)     
              c_Stress_Theta=c_KI_3D/FOU/tem_root_2pir*c_tem_part1 -c_KII_3D/FOU/tem_root_2pir*c_tem_part2          
              c_Tao_Theta = c_KIII_3D*cos(c_Theta/TWO)/tem_root_2pir          
              c_Sai = ZP5*atan(TWO*c_Tao_Theta/c_Stress_Theta)
              Sai_Cr(i_Out_Node) = c_Sai
          endif
          
          !........................................
          ! Calculate the equivalent stress KI_eq.
          !........................................
          cc_Part1 = c_KI_3D*(cos(c_Theta/TWO))**2
          cc_Part2 = THR/TWO*c_KII_3D*sin(c_Theta)
          !cc_Part3 =sqrt((cc_Part1 -cc_Part2)**2 + FOU*c_KIII_3D*2)
          cc_Part3 =sqrt((cc_Part1 -cc_Part2)**2 + FOU*c_KIII_3D**2)
          KI_eq(i_Out_Node)  =ZP5*cos(c_Theta/TWO)*(cc_Part1-cc_Part2+cc_Part3)
          !KI_eq(i_Out_Node) = KI_eq_3D(i_C,i_Out_Node)    !2022-06-29
          
          !......................................
          ! Mark whether the crack has extended.
          !......................................
          ! The equivalent stress intensity factor must be greater than the fracture toughness in order to
          ! propagate.!NEWFTU2022070201.
          if(Key_CFCP_3_Type   == 1) then  
              if(KI_eq(i_Out_Node)>=c_KIc)then
                  Flag_Propagation_this_Crack = .True.
                  Crack_Type_Status_3D(i_C,5) = 1
              endif
          ! Cracks can propagate as long as the equivalent stress intensity factor is greater than 0. Ref:
          ! Tang_2019_Analysis of stress interference among_Eq.16.
          elseif(Key_CFCP_3_Type   == 2) then
              if(KI_eq(i_Out_Node)>=Zr+Tol_3)then
                  Flag_Propagation_this_Crack = .True.
                  Crack_Type_Status_3D(i_C,5) = 1
              endif
          endif
        end do
        
        ! If none of the current crack's vertices meet the stress intensity factor crack propagation
        ! criterion, then proceed to the next crack.
        if(Flag_Propagation_this_Crack .eqv. .False.)then
            write(*,1062) i_C
            goto 100
        endif            
        
        !.......................................................................................
        ! Calculate the crack surface distortion vector w, theoretical Pereira_2010_Generalized
        ! finite element
        ! Methods for three-dimensional crack growth simulations, Equation 3.12.
        !.......................................................................................
        if (Key_Schollmann_Twist == 1)then
          allocate(Global_kw(Num_CrMesh_Outlines,Num_CrMesh_Outlines))
          allocate(Global_fw(Num_CrMesh_Outlines))       
          allocate(Global_w(Num_CrMesh_Outlines))    
          Global_kw(1:Num_CrMesh_Outlines,1:Num_CrMesh_Outlines)=ZR
          Global_fw(1:Num_CrMesh_Outlines) = ZR
          Global_w(1:Num_CrMesh_Outlines)  = ZR
          ! Fluid Element Cycle
          do i_FluidEl = 1,Num_CrMesh_Outlines
            Sai_Cr_1 = Sai_Cr(i_FluidEl)
            if(i_FluidEl<Num_CrMesh_Outlines)then
                Sai_Cr_2 = Sai_Cr(i_FluidEl+1)
            else
                Sai_Cr_2 = Sai_Cr(1)
            endif
            c_Ld_w  = Ld_Cr(i_FluidEl)
            c_Ks_w  = Ks_Cr(i_FluidEl)
            c_D_w   = D_Cr(i_FluidEl)
            kw(1,1) = (c_Ks_w*13.0D0*c_Ld_w**4 + 420.0D0*c_D_w)/(35.0D0*c_Ld_w**3)
            kw(1,2) = THR/TWO*(c_Ks_w*THR*c_Ld_w**4-280.0D0*c_D_w)/(35.0D0*c_Ld_w**3)
            kw(2,1) = kw(1,2)
            kw(2,2) = kw(1,1)
            fw(1) = -((c_Ks_w*11.0D0*c_Ld_w**4 +1260.0D0*c_D_w)*Sai_Cr_1 - &
              ONE/TWO*(c_Ks_w*13.0D0*c_Ld_w**4-2520.0D0*c_D_w)*Sai_Cr_2)/(210.0D0*c_Ld_w**2)
            fw(2) =  -(ONE/TWO*(c_Ks_w*13.0D0*c_Ld_w**4 -2520.0D0*c_D_w)   &
                  *Sai_Cr_1-(c_Ks_w*11.0D0*c_Ld_w**4 +1260.0D0*c_D_w)*Sai_Cr_2)/(210.0D0*c_Ld_w**2)  
            Global_fw(i_FluidEl:i_FluidEl+1)  =  Global_fw(i_FluidEl:i_FluidEl+1) + fw(1:2)
            Global_kw(i_FluidEl,i_FluidEl)    =  Global_kw(i_FluidEl,i_FluidEl)     + kw(1,1)     
            Global_kw(i_FluidEl,i_FluidEl+1)  =  Global_kw(i_FluidEl,i_FluidEl+1)   + kw(1,2)         
            Global_kw(i_FluidEl+1,i_FluidEl)  =  Global_kw(i_FluidEl+1,i_FluidEl)   + kw(2,1)     
            Global_kw(i_FluidEl+1,i_FluidEl+1)=  Global_kw(i_FluidEl+1,i_FluidEl+1) + kw(2,2)        
          enddo
          ! Penalty function constraint, the first and last w values are equal (Equation 2-87 in my book)
          Penalty_Value = 1.0D4*maxval(abs(Global_kw))
          
          Global_kw(1,1) = Global_kw(1,1) + Penalty_Value
          Global_kw(1,Num_CrMesh_Outlines) =Global_kw(1,Num_CrMesh_Outlines) - Penalty_Value
          Global_kw(Num_CrMesh_Outlines,1) =Global_kw(Num_CrMesh_Outlines,1) - Penalty_Value     
          Global_kw(Num_CrMesh_Outlines,Num_CrMesh_Outlines) =Global_kw(Num_CrMesh_Outlines,Num_CrMesh_Outlines)+ Penalty_Value          
          ! Solve the system of linear equations and calculate the w vector
          call Matrix_Solve_LSOE(-1,1,Key_SLOE,Global_kw,Global_fw,Global_w,Num_CrMesh_Outlines)    
          deallocate(Global_kw)
          deallocate(Global_fw)     
        endif

        !......................................................
        ! Calculate the crack front Growth Factor. 2022-07-14.
        !......................................................
        !BUGFIX2023021401.
        if(Num_CrMesh_Outlines>5000) then
            print *,'    ERROR-2023021401 :: Num_CrMesh_Outlines > 5000 in Check_Crack_Grows_3D.f!' 
            print *,'                        Try to increase size of vector All_Growth_Factor!' 
            call Warning_Message('S',Keywords_Blank)  
        endif
        
        All_Growth_Factor(1:Num_CrMesh_Outlines) = ZR
        
        do i_Out_Node =1,Num_CrMesh_Outlines 
          ! BUGFIX2023021201. Added the following lines, major bug.
          c_Mesh_Node = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
          
          in_Ele_Num =Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(c_Mesh_Node)
          
          !2023-08-16.
          if(in_Ele_Num ==0 ) then
              All_Growth_Factor(i_Out_Node) = ZR
              cycle
          endif 
 
          ! Not Neo-O.
          if(Key_XA/=2) then
              mat_num = Elem_Mat(in_Ele_Num)
              c_KIc  = Material_Para(mat_num,6)
          !Key_XA==2. IMPROV2023031905.
          else
              c_KIc  = Elem_KIc_XA(in_Ele_Num)
          endif
          
          ! Special handling for Key_NaCr_Active_Scheme_3D = 3. NEWFTU2023011204. 2023-01-12.
          if(Key_NaCr_Active_Scheme_3D == 3)then
              ! If the current crack was caused by the activation of a natural fracture.
              if (Crack_Type_Status_3D(i_C,6) >=1) then
                  Natural_Crack = Crack_Type_Status_3D(i_C,6)
                  ! Check whether the element containing the crack tip is in the list of elements with natural
                  ! fractures.
                  Num_List_Elements = NaCr3D_Status(Natural_Crack,2)
                  ! Using Vector_Location_Int(n, Vector, Variable, location, Yes_In)
                  call Vector_Location_Int(Num_List_Elements,&
                       Na_Crack3D_Ele_List(Natural_Crack)%row(1:Num_List_Elements),in_Ele_Num,&
                       c_Ele_location,c_Yes_Ele_In)   
                  !2023-05-16. BUGFIX2023051601.
                  if (c_Yes_Ele_In) then
                      ! If the element where the crack tip is located belongs to the list of elements containing natural
                      ! fractures, the fracture toughness is taken as that of the natural fractures.
                      if (Key_Cpp_Call_Fortran_Lib/=1) then
                          c_KIc  = KIc_NaCr(Natural_Crack) 
                      elseif(Key_Cpp_Call_Fortran_Lib==1)then
                          c_KIc  = Na_Crack3D_St(Natural_Crack) 
                      endif
                  endif
              endif
          endif
          
          !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
          ! Crack propagation step OPTION 1: Only propagate when exceeding fracture toughness.
          ! NEWFTU2022070201.
          !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
          if(Key_CFCP_3_Type   == 1) then
            ! Check whether the equivalent stress KI_eq is greater than the fracture toughness
            if (KI_eq(i_Out_Node) >=c_KIc) then
              c_Growth_Factor = (KI_eq(i_Out_Node)/c_KIc)**Prp_3D_Factor_m
              if(c_Growth_Factor>=Adj_Prp_Step_3D_Max)then
                  c_Growth_Factor=Adj_Prp_Step_3D_Max
              endif  
              Yes_Grow(i_C) = .TRUE.
            else
              c_Growth_Factor     = ZR            
            endif
          !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
          ! Crack propagation step length OPTION 2: Extend as long as it is greater than 0.
          ! 2022-07-02.
          ! NEWFTU2022070201.
          !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
          elseif(Key_CFCP_3_Type   == 2) then
            if (KI_eq(i_Out_Node) >=(ZR+Tol_3)) then
                c_Growth_Factor = KI_eq(i_Out_Node)/c_KIc
                if(c_Growth_Factor>=Adj_Prp_Step_3D_Max)then
                    c_Growth_Factor=Adj_Prp_Step_3D_Max
                endif      
                Yes_Grow(i_C) = .TRUE.
            else
                c_Growth_Factor = ZR
            endif
          endif
          
          ! Interactive display.
          
          ! If this vertex is a Segmentation suppression vertex (2021-08-20)
          if (key_front_segmentation==1)then
            if (Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,4)==0)then
              c_Growth_Factor = ZR
            endif
          endif                 
          All_Growth_Factor(i_Out_Node) = c_Growth_Factor
        enddo
        
        !.........................................................................
        ! Perform denoising and smoothing on the Growth Factor. NEWFTU2022071402.
        !.........................................................................
        ! Denoising.
        if (Key_Denoise_GF_Value>=1)then
            n_Sigma =2
            call Tool_Denoise_Data(All_Growth_Factor(1:Num_CrMesh_Outlines), &
                           Num_CrMesh_Outlines,Key_Denoise_Vertex_Value,n_Sigma,.True.)
        endif
        ! Data smoothing.
        if (Key_Smooth_GF_Value>=1)then
            call Tool_Smooth_Data(All_Growth_Factor(1:Num_CrMesh_Outlines), Num_CrMesh_Outlines, &
                        Key_Smooth_GF_Value,Smooth_GF_n,.True.)
            ! Data smoothing (secondary processing).
            if (Key_Smooth_GF_Value2>=1)then
              call Tool_Smooth_Data(All_Growth_Factor(1:Num_CrMesh_Outlines),Num_CrMesh_Outlines, &
                        Key_Smooth_GF_Value2,Smooth_GF_n2,.True.)
            endif  
        endif  
        
        !...........................................................................................
        ! Calculate the crack front propagation vector, theoretical Pereira_2010_Generalized finite
        ! element
        ! methods for three-dimensional crack growth simulations, P122, Figure 3.15
        !...........................................................................................
        do i_Out_Node =1,Num_CrMesh_Outlines
          c_Mesh_Node =Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
          in_Ele_Num  =Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(c_Mesh_Node)
          ! A element number of 0 indicates that the fracture discrete point is outside the model. 2022-06-11.
          if(in_Ele_Num ==0 ) then
              ! If cracks are allowed to extend outside the model. 2023-08-16.
              if(Key_Allow_3D_Outside_Crack==1) then
                  !do nothing. 2023-08-16.
              ! If cracks are not allowed to extend outside the model.
              else
                  ! Key_Stop_Outside_Crack =1, once the crack extends outside the model, it is marked as no longer
                  ! expanding. IMPROV2022102301.
                  if(Key_Stop_Outside_Crack ==1) then
                      print *,'    Warn-2022102311 :: in_Ele_Num =0 in Check_Crack_Grows_3D.f!' 
                      print *,'                       Maybe caused by outside crack vertex.' 
                      print *,'                       Crack ',i_C,'stops propagation!'
                      !call Warning_Message('S',Keywords_Blank)  
                      Crack_Type_Status_3D(i_C,3) = 0
                      goto 100
                  elseif(Key_Stop_Outside_Crack ==0) then 
                      print *,'    ERROR-2022061101 :: in_Ele_Num =0 in Check_Crack_Grows_3D.f!' 
                      print *,'                        Maybe caused by outside crack vertex.' 
                      print *,'                        Try *Key_Stop_Outside_Crack=1.'
                      call Warning_Message('S',Keywords_Blank)  
                  endif

              endif
          endif              
          c_x  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1)
          c_y  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,2) 
          c_z  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,3)     
          
          c_Theta = Theta_All(i_Out_Node)
 
          
          !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
          ! Consider the effect of the twisted vector w (Key_Schollmann_Twist==1)
          !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
          if(Key_Schollmann_Twist==1)then
              ! c_Theta = c_Theta   Global_w(i_Out_Node)/delta_L   !Angle equals arc length divided by radius
              c_Theta = c_Theta + TWO*asin(Global_w(i_Out_Node)/delta_L/TWO)   
 
              WRITE(*,1032)  i_Out_Node,i_C,Global_w(i_Out_Node)/delta_L*180/pi               
          endif

          !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ! Calculate the extension direction vector poss_Vector(1:3) of the discrete nodes at the crack
          ! boundary
          ! The theory is as follows: (1) In the local coordinate system xy plane at the discrete nodes of the
          ! crack boundary, the x-axis
          ! Direction vector, rotating about the y-axis from the local coordinate origin by c_Theta
          ! The angle obtained is to rotate the original local coordinate system around the z-axis by an angle
          ! of c_Theta.
          ! (2) The rotation algorithm can be found in 3D Math Primer for Graphics and Game
          ! Development_2011_Introduction to Three-Dimensional Geometry P144
          !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&       
          
          n_x = Crack3D_Meshed_Vertex_z_Vector(i_C)%row(i_Out_Node,1)
          n_y = Crack3D_Meshed_Vertex_z_Vector(i_C)%row(i_Out_Node,2)
          n_z = Crack3D_Meshed_Vertex_z_Vector(i_C)%row(i_Out_Node,3)
          Rot_Matrix(1,1)=n_x**2*(ONE-cos(c_Theta))+cos(c_Theta)
          Rot_Matrix(1,2)=n_x*n_y*(ONE- cos(c_Theta))+ n_z*sin(c_Theta)
          Rot_Matrix(1,3)= n_x*n_z*(ONE- cos(c_Theta))- n_y*sin(c_Theta)   
          Rot_Matrix(2,1)= n_x*n_y*(ONE- cos(c_Theta))- n_z*sin(c_Theta) 
          Rot_Matrix(2,2)=n_y**2*(ONE-cos(c_Theta))+cos(c_Theta)     
          Rot_Matrix(2,3)= n_y*n_z*(ONE- cos(c_Theta))+ n_x*sin(c_Theta)   
          Rot_Matrix(3,1)= n_x*n_z*(ONE- cos(c_Theta))+ n_y*sin(c_Theta)      
          Rot_Matrix(3,2)= n_y*n_z*(ONE- cos(c_Theta))- n_x*sin(c_Theta)       
          Rot_Matrix(3,3)=n_z**2*(ONE-cos(c_Theta))+cos(c_Theta)      
          c_poss_Vector(i_Out_Node,1:3) = MATMUL(Crack3D_Meshed_Vertex_x_Vector(i_C)%row(i_Out_Node,1:3),Rot_Matrix(1:3,1:3))
          if(sum(abs(c_poss_Vector(i_Out_Node,1:3)))<=Tol_10)then
              print *,'    Error :: c_poss_Vector!'
              print *,'             in Check_Crack_Grows_3D.f'
              print *,'             i_C,i_Out_Node',i_C,i_Out_Node
              call Warning_Message('S',Keywords_Blank)                      
          endif
        end do

        !............................................................................................
        ! Calculate the points after the crack front expansion, theoretical Pereira_2010_Generalized
        ! finite element
        ! methods for three-dimensional crack growth simulations, P122, Figure 3.15
        !   Modified on 2021-09-04.
        !   BUGFIX2021090401
        !............................................................................................
        Ele_Num_Cache = 1
        
        if(allocated(Pre_Fronts)) deallocate(Pre_Fronts)
        allocate(Pre_Fronts(Num_CrMesh_Outlines,3))
        
        do i_Out_Node =1,Num_CrMesh_Outlines
          c_Mesh_Node=Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
          in_Ele_Num =Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(c_Mesh_Node)
          c_x  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1)
          c_y  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,2) 
          c_z  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,3)     
          
          Pre_Fronts(i_Out_Node,1:3) = [c_x,c_y,c_z]
          
          ! For 3D edge cracks
          All_Out = Num_CrMesh_Outlines
          if(i_Out_Node == 1)then
            Pre_Mesh_Node1=Crack3D_Meshed_Outline(i_C)%row(All_Out,1)
            Pre_Mesh_Node2=Crack3D_Meshed_Outline(i_C)%row(All_Out-1,1)
          elseif(i_Out_Node == 2)then
            Pre_Mesh_Node1=Crack3D_Meshed_Outline(i_C)%row(1,1)
            Pre_Mesh_Node2=Crack3D_Meshed_Outline(i_C)%row(All_Out,1)
          else
            Pre_Mesh_Node1=Crack3D_Meshed_Outline(i_C)%row(i_Out_Node-1,1)
            Pre_Mesh_Node2=Crack3D_Meshed_Outline(i_C)%row(i_Out_Node-2,1)
          endif
          if(i_Out_Node == All_Out)then
            Nex_Mesh_Node1=Crack3D_Meshed_Outline(i_C)%row(1,1)
            Nex_Mesh_Node2=Crack3D_Meshed_Outline(i_C)%row(2,1)
          elseif(i_Out_Node == All_Out-1)then
            Nex_Mesh_Node1=Crack3D_Meshed_Outline(i_C)%row(All_Out,1)
            Nex_Mesh_Node2=Crack3D_Meshed_Outline(i_C)%row(1,1)
          else
            Nex_Mesh_Node1=Crack3D_Meshed_Outline(i_C)%row(i_Out_Node+1,1)
            Nex_Mesh_Node2=Crack3D_Meshed_Outline(i_C)%row(i_Out_Node+2,1)
          endif
      
          !****************************************************************************************
          ! Check whether the boundary point is outside the model.
          ! OPTION 1: If so, this boundary point is not allowed to extend (only a very small step
          ! size for extension)
          ! OPTION 2: If so, then further determine whether the left and right sides of this
          ! boundary point
          ! If the boundary point is also outside the model, then extending this boundary point is
          ! not allowed.
          ! Otherwise, extend according to the extension direction of its adjacent boundary points
          ! (2021-04-20)
          !****************************************************************************************
          call Cal_Ele_Num_by_Coors_3D(c_x,c_y,c_z,Ele_Num_Cache,c_OUT_Elem)
          c_OUT_Elem_Old(i_Out_Node) = c_OUT_Elem
          if(c_OUT_Elem ==0) then
              !~~~~~~~~~~~~~~~~~     OPTION 1   ~~~~~~~~~~~~~~~~~~~~~~~~
              !call Cal_Ele_Num_by_Coors_3D(c_x,c_y,c_z,c_OUT_Elem)
              !c_OUT_Elem_Old(i_Out_Node) = c_OUT_Elem
              !if(c_OUT_Elem ==0) then
              !    !cycle  !Start the next do i_Out_Node loop.
              ! c_Vector_x = Crack3D_Meshed_Vertex_x_Vector(i_C, i_Out_Node, 1:3)! Local x-axis vector of the 3D
              ! crack boundary point after discretization
              !    c_P(1)=c_x+delta_L_Stoped*c_Vector_x(1)
              !    c_P(2)=c_y+delta_L_Stoped*c_Vector_x(2)
              !    c_P(3)=c_z+delta_L_Stoped*c_Vector_x(3)    
              !    Final_Point(i_Out_Node,1:3) = c_P(1:3)
              !    cycle
              !endif
              !~~~~~~~~~~~~~~~~~     OPTION 2   ~~~~~~~~~~~~~~~~~~~~~~~~
              !Added on 2021-04-20.
              ! Determine whether Pre_Mesh_Node1 and Pre_Mesh_Node2 are inside the model
              c_x_Pre1=Crack3D_Meshed_Node(i_C)%row(Pre_Mesh_Node1,1) 
              c_y_Pre1=Crack3D_Meshed_Node(i_C)%row(Pre_Mesh_Node1,2) 
              c_z_Pre1=Crack3D_Meshed_Node(i_C)%row(Pre_Mesh_Node1,3) 
              c_x_Pre2=Crack3D_Meshed_Node(i_C)%row(Pre_Mesh_Node2,1) 
              c_y_Pre2=Crack3D_Meshed_Node(i_C)%row(Pre_Mesh_Node2,2) 
              c_z_Pre2=Crack3D_Meshed_Node(i_C)%row(Pre_Mesh_Node2,3) 
              c_x_Nex1=Crack3D_Meshed_Node(i_C)%row(Nex_Mesh_Node1,1) 
              c_y_Nex1=Crack3D_Meshed_Node(i_C)%row(Nex_Mesh_Node1,2) 
              c_z_Nex1=Crack3D_Meshed_Node(i_C)%row(Nex_Mesh_Node1,3) 
              c_x_Nex2=Crack3D_Meshed_Node(i_C)%row(Nex_Mesh_Node2,1) 
              c_y_Nex2=Crack3D_Meshed_Node(i_C)%row(Nex_Mesh_Node2,2) 
              c_z_Nex2=Crack3D_Meshed_Node(i_C)%row(Nex_Mesh_Node2,3) 
              call Cal_Ele_Num_by_Coors_3D(c_x_Pre1,c_y_Pre1,c_z_Pre1,Ele_Num_Cache,Elem_Pre1)
              call Cal_Ele_Num_by_Coors_3D(c_x_Pre2,c_y_Pre2,c_z_Pre2,Ele_Num_Cache,Elem_Pre2)
              call Cal_Ele_Num_by_Coors_3D(c_x_Nex1,c_y_Nex1,c_z_Nex1,Ele_Num_Cache,Elem_Nex1)
              call Cal_Ele_Num_by_Coors_3D(c_x_Nex2,c_y_Nex2,c_z_Nex2,Ele_Num_Cache,Elem_Nex2)     
              
              ! If both adjacent boundary points are 0
              if(Elem_Pre1 ==0 .and. Elem_Nex1==0)then
                  c_Vector_x=Crack3D_Meshed_Vertex_x_Vector(i_C)%row(i_Out_Node,1:3)
                  c_P(1)=c_x+delta_L_Stoped*c_Vector_x(1)
                  c_P(2)=c_y+delta_L_Stoped*c_Vector_x(2)
                  c_P(3)=c_z+delta_L_Stoped*c_Vector_x(3)    
                  Final_Point(i_Out_Node,1:3) = c_P(1:3)
                  cycle 
              endif
              ! If the previous adjacent boundary point is inside the model, then the vector in the extension
              ! direction is equal to the extension direction vector of the previous adjacent boundary point.
              if(Elem_Pre1 /=0)then
                  if (i_Out_Node==1) then
                      c_Prop_Vector=c_poss_Vector(Num_CrMesh_Outlines,1:3) 
                  else
                     c_Prop_Vector=c_poss_Vector(i_Out_Node-1,1:3) 
                  endif
                  c_P(1)=c_x+delta_L*c_Prop_Vector(1)
                  c_P(2)=c_y+delta_L*c_Prop_Vector(2)
                  c_P(3)=c_z+delta_L*c_Prop_Vector(3)    
                  Final_Point(i_Out_Node,1:3) = c_P(1:3)
                  cycle 
              endif
              ! If the next adjacent boundary point is inside the model, then the vector in the extension
              ! direction is equal to the extension direction vector of the next adjacent boundary point.
              if(Elem_Nex1 /=0)then
                  if (i_Out_Node==Num_CrMesh_Outlines) then
                      c_Prop_Vector = c_poss_Vector(1,1:3) 
                  else
                     c_Prop_Vector=c_poss_Vector(i_Out_Node+1,1:3) 
                  endif
                  c_P(1)=c_x+delta_L*c_Prop_Vector(1)
                  c_P(2)=c_y+delta_L*c_Prop_Vector(2)
                  c_P(3)=c_z+delta_L*c_Prop_Vector(3)    
                  Final_Point(i_Out_Node,1:3) = c_P(1:3)
                  cycle 
              endif              
          endif
          
          
          !*************************
          ! New vertex coordinates.
          !*************************
          Final_Point(i_Out_Node,1)=c_x+delta_L*All_Growth_Factor(i_Out_Node)*c_poss_Vector(i_Out_Node,1)
          Final_Point(i_Out_Node,2)=c_y+delta_L*All_Growth_Factor(i_Out_Node)*c_poss_Vector(i_Out_Node,2)
          Final_Point(i_Out_Node,3)=c_z+delta_L*All_Growth_Factor(i_Out_Node)*c_poss_Vector(i_Out_Node,3)   
          
          ! If this vertex is a Segmentation suppression vertex (2021-08-20)
          if (key_front_segmentation==1)then
            if (Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,4)==0)then
             ! Potential suppression point coordinates (in the subsequent D3_Prepare_Front_Segmentation, the
             ! suppression points will be determined based on the coordinates of these points)
             ! This is used to set Crack3D_Meshed_Outline(i_C, i_Out_Node, 4) = 0
             num_Suspended_Point = num_Suspended_Point +1
             Suspended_Points(num_Suspended_Point,1:3)=Final_Point(i_Out_Node,1:3)
            endif
          endif    
          
          !------------------------------------------------------------------------------------------
          ! If it is in-plane growth (Key_InPlane_Growth = 1), then correct the Final_Point back to
          ! the plane of the initial crack.
          ! Project onto the plane, mark the coordinates of the foot of the perpendicular as the new
          ! Final_Point, 2022-04-18, NEWFTU2022041801
          !------------------------------------------------------------------------------------------
          if(Key_InPlane_Growth== 1) then
            call D3_In_Plane_Growth(Final_Point(i_Out_Node,1:3),i_C)
          endif
        end do
        
        
        if(Key_Schollmann_Twist==1)then
            if(allocated(Global_w))deallocate(Global_w)    
        endif
        ! Clear other temporary variables from memory
        if(allocated(Check_Theta))DEALLOCATE(Check_Theta)
        if(allocated(Spe_Pri_Stress))DEALLOCATE(Spe_Pri_Stress)
        if(allocated(KI_eq))DEALLOCATE(KI_eq)
        if(allocated(theta_all))DEALLOCATE(theta_all)
    end select
    
    !------------------------------------------------------------------------------
    ! If the crack does not extend, jump to statement 100; bug fixed on 2019-08-03
    !------------------------------------------------------------------------------
    if(Yes_Grow(i_C) .eqv. .False.)then
        goto 100
    else
        print *,'    Crack ',i_C,' propagated!'
    endif
    
    !------------------------------------------------------------------------------------------
    ! Taubin Smoothing. 2022-12-18. IMPROV2022122002.
    ! Taubin algorithm. Key_Smooth_Front=6. IMPROV2022121001.
    !  Ref: \theory_documents\042 Mesh Smoothing-2022-12-10.pdf，P15-P25.
    ! Only applicable to the leading edge of closed cracks.
    ! It should be ensured that the apex of the front edge of all crack surfaces is within the
    ! model.
    ! Applicable to various crack propagation criteria (CFCP).
    !------------------------------------------------------------------------------------------
    if(Key_Smooth_Front==6)then
        ! Smoothing is performed only when the number of vertices is 12 or more.
        if(Num_CrMesh_Outlines>=12)then
            allocate(Points_to_Smooth(Num_CrMesh_Outlines,3))
            allocate(Points_Smoothed(Num_CrMesh_Outlines,3))
            Points_to_Smooth(1:Num_CrMesh_Outlines,1:3) = Final_Point(1:Num_CrMesh_Outlines,1:3)
            ! Taubin smoothing.
            call Tool_Fit_3D_Points_Taubin(i_C,Num_CrMesh_Outlines, &
                            Points_to_Smooth(1:Num_CrMesh_Outlines,1:3), &
                            Points_Smoothed(1:Num_CrMesh_Outlines,1:3))   
                            
            ! Prevent the vertex at the front edge of the crack from retreating back into the original crack
            ! front polygon after optimization. For theoretical details, see PhiPsi Record V1, P91.
            ! IMPROV2022122001.
            ! For out-of-plane extension, see PhiPsi Record V1, p.133 for details.
            call Tool_Fit_3D_Points_Check_Rollback(i_C,Num_CrMesh_Outlines,&
                             Pre_Fronts(1:Num_CrMesh_Outlines,1:3),&
                             Points_to_Smooth(1:Num_CrMesh_Outlines,1:3),&
                             Points_Smoothed(1:Num_CrMesh_Outlines,1:3))
            
            Final_Point(1:Num_CrMesh_Outlines,1:3) = Points_Smoothed(1:Num_CrMesh_Outlines,1:3)
            deallocate(Points_to_Smooth)
            deallocate(Points_Smoothed)
        endif
    endif
    
    !------------------------------------------------------------!
    !                                                            !
    !                                                            !
    !                                                            !
    ! Generate new discrete fracture elements and nodes.
    !                                                            !
    ! Note: Depending on the crack surface update algorithm Key_3D_Cr_Update_Scheme!
    ! Use different processing algorithms. NEWFTU2022071201. !
    !                                                            !
    !                                                            !
    !------------------------------------------------------------!
    old_node_num = Crack3D_Meshed_Node_num(i_C)
    old_ele_num  = Crack3D_Meshed_Ele_num(i_C)
    
    ! Handling the connection of data at the beginning and end.
    select case(Key_3D_Cr_Update_Scheme)
    !//////////////////////////////////////////////////////////////////////////////////////////////
    !/                                                            /
    !/                                                            /
    !/     Key_3D_Cr_Update_Scheme=1                              /
    ! / Regardless of whether it expands or not, always extend by one step; for crack tips that do
    ! not expand, /
    ! /     Add a tiny amount. Original algorithm.                                /
    !/                                                            /
    !//////////////////////////////////////////////////////////////////////////////////////////////
    case(1)
      Ele_Num_Cache = 1
      do i_Out_Node = 1,Num_CrMesh_Outlines          
        !********************************************************************************************
        ! Check whether the boundary point is outside the model; if it is, the boundary point is not
        ! allowed to expand (the expansion should be very small).
        !********************************************************************************************
        call Cal_Ele_Num_by_Coors_3D(Final_Point(i_Out_Node,1),Final_Point(i_Out_Node,2),Final_Point(i_Out_Node,3),&
                                     Ele_Num_Cache,c_OUT_Elem)
        ! It was originally inside the model, but now it has expanded outside the model, so a warning is
        ! issued.
        if(c_OUT_Elem ==0 .and. c_OUT_Elem_Old(i_Out_Node) /=0) then
            print *,'    Warning::Vertex ',i_Out_Node,' is going out of model.'
            !----------------------------------------------------------------------------------------
            ! If Key_Stop_Outside_Crack=1, when the leading edge vertex of the crack extends outside
            ! the crack,
            ! Then stop expanding. NEWFTU2022100202. 2022-10-02.
            !----------------------------------------------------------------------------------------
            if (Key_Stop_Outside_Crack ==1) then
                ! It was originally inside the model, but now it has expanded outside the model, so a warning is
                ! issued.
                if(c_OUT_Elem ==0 .and. c_OUT_Elem_Old(i_Out_Node) /=0) then
                    print *,'    Warning :: crack ',i_C,'will stop propagation!'
                    Crack_Type_Status_3D(i_C,3) = 0
                    goto 100
                endif  
            endif
        endif            
        ! Check if the array is out of bounds. 2022-10-04. IMPROV2022100403.
        c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
        if((old_node_num+i_Out_Node) > c_Max_N_Node_3D) then
          ! Expand memory. NEWFTU2022110501.
          call D3_Allocate_Crack_Memory(i_C,1,1)
          !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
        endif  
        
        ! Add discrete node
        Crack3D_Meshed_Node(i_C)%row(old_node_num+i_Out_Node,1:3)= Final_Point(i_Out_Node,1:3) 
        ! Update the number of discrete nodes
        Crack3D_Meshed_Node_num(i_C) =Crack3D_Meshed_Node_num(i_C)  +1     
        
        c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
        if(Crack3D_Meshed_Ele_num(i_C)> c_Max_N_Node_3D)then
            !call Warning_Message('S',Keywords_Blank) 
            ! Expand memory. NEWFTU2022110501.
            call D3_Allocate_Crack_Memory(i_C,1,1)
            !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
        endif 
        
        ! Element number of the 3D crack node after discretization
        call Cal_Ele_Num_by_Coors_3D(Final_Point(i_Out_Node,1),Final_Point(i_Out_Node,2), &
                                     Final_Point(i_Out_Node,3),Ele_Num_Cache,in_Elem_num)               
        Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(old_node_num+i_Out_Node)=in_Elem_num
        ! Local coordinates of the elements containing the 3D crack nodes after discretization
        if (in_Elem_num>0)then
            call Cal_KesiYita_by_Coor_3D(Final_Point(i_Out_Node,1:3),in_Elem_num,c_Kesi,c_Yita,c_Zeta)
            Cr3D_Meshed_Node_in_Ele_Local(i_C)%row(old_node_num+i_Out_Node,1:3) =[c_Kesi,c_Yita,c_Zeta] 
        endif    
        ! Update the number of discrete crack elements
        Crack3D_Meshed_Ele_num(i_C)=Crack3D_Meshed_Ele_num(i_C)+1      

        ! Check the number of discrete nodes.
        c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
        if(Crack3D_Meshed_Ele_num(i_C)> c_Max_N_Node_3D)then
            !call Warning_Message('S',Keywords_Blank) 
            ! Expand memory. NEWFTU2022110501.
            call D3_Allocate_Crack_Memory(i_C,1,1)
            !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
        endif  
        
        ! Node 1 of the new element
        new_node_1 =  Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
        ! Node 2 of the new element
        new_node_2 = old_node_num+i_Out_Node    
        ! Node 3 of the new element (the next node after the original outer boundary)
        new_node_3 = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,2)       
        Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),1) = new_node_1   
        Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),2) = new_node_2 
        Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),3) = new_node_3   
        
        ! Update the number of discrete crack elements
        Crack3D_Meshed_Ele_num(i_C)=Crack3D_Meshed_Ele_num(i_C)+1   
        
        ! Check the number of discrete nodes.
        c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
        if(Crack3D_Meshed_Ele_num(i_C)> c_Max_N_Node_3D)then
            !call Warning_Message('S',Keywords_Blank) 
            ! Expand memory. NEWFTU2022110501.
            call D3_Allocate_Crack_Memory(i_C,1,1)
            !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
        endif  
        ! Node 1 of the new element
        new_node_1 =  old_node_num+i_Out_Node  
        ! Node 2 of the new element (requires categorized handling)
        if (i_Out_Node < Num_CrMesh_Outlines)then
            new_node_2 = old_node_num+i_Out_Node+1   
        elseif(i_Out_Node == Num_CrMesh_Outlines)then
            new_node_2 = old_node_num+1 
        endif
        ! Node 3 of the new element (the next node after the original outer boundary)
        new_node_3 =Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,2) 
        ! Three nodes of the new element
        Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),1) = new_node_1   
        Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),2) = new_node_2 
        Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),3) = new_node_3
      enddo
    !/////////////////////////////////////////////////////////////////////////////////////////
    !/                                                            
    !/                                                            
    !/     Key_3D_Cr_Update_Scheme=2 (default)                       
    !/ New algorithm, non-propagating crack tips do not extend, propagating crack tips update
    !coordinates with very small steps 
    !/     See my PhiPsi development notebook, page 38.                      
    !/     2022-07-12.                                            
    !/                                                            
    !/////////////////////////////////////////////////////////////////////////////////////////
    case(2)
      !..............................................
      ! Get the coordinates of the old crack vertex.
      !..............................................
      if (allocated(Old_Vertex_Points))deallocate(Old_Vertex_Points)
      allocate(Old_Vertex_Points(Num_CrMesh_Outlines,3))
      do i_Out_Node = 1,Num_CrMesh_Outlines 
        !--------------------------------------------------------------------------------------------
        ! If Key_Stop_Outside_Crack=1, when the leading edge vertex of the crack extends outside the
        ! crack,
        ! Then stop expanding. NEWFTU2022100202. 2022-10-02.
        !--------------------------------------------------------------------------------------------
        if (Key_Stop_Outside_Crack ==1) then
            ! Check whether the boundary point is outside the model.
            call Cal_Ele_Num_by_Coors_3D(Final_Point(i_Out_Node,1),Final_Point(i_Out_Node,2),Final_Point(i_Out_Node,3),&
                                         Ele_Num_Cache,c_OUT_Elem)
            ! It was originally inside the model, but now it has expanded outside the model, so a warning is
            ! issued.
            if(c_OUT_Elem ==0 .and. c_OUT_Elem_Old(i_Out_Node) /=0) then
                write (*,1071) i_Out_Node,i_C
                print *,'    Warning :: crack ',i_C,'will stop propagation!'
                Crack_Type_Status_3D(i_C,3) = 0
                goto 100
            endif  
        endif
        c_Mesh_Node = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
        Old_Vertex_Points(i_Out_Node,1:3) = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1:3)
      enddo

      !.....................................................................................
      ! Determine whether each point should be expanded (points that expand very little are
      ! considered not to expand here (only update their position)).
      !.....................................................................................
      if (allocated(Flag_Ver_Grow))deallocate(Flag_Ver_Grow)
      allocate(Flag_Ver_Grow(Num_CrMesh_Outlines+1))
      Flag_Ver_Grow(1:Num_CrMesh_Outlines+1) = .False.
      do i_Out_Node = 1,Num_CrMesh_Outlines   
        ! The extension length of the current vertex.
        c_Grow_Dis = Tool_Function_2Point_Dis_3D(Old_Vertex_Points(i_Out_Node,1:3),Final_Point(i_Out_Node,1:3))
        if(c_Grow_Dis>=0.2D0*Ave_Elem_L_Enrich)then
        !if(c_Grow_Dis>=0.1D0*Ave_Elem_L_Enrich)then !2023-05-16.
            Flag_Ver_Grow(i_Out_Node) = .True.
        endif
      enddo
      Flag_Ver_Grow(Num_CrMesh_Outlines+1)  = Flag_Ver_Grow(1) 
      
      !..........................................................................................
      ! For points where no expansion has occurred (including points with very small expansion),
      ! directly update their coordinates.
      ! without adding points to the discrete fracture surfaces.
      !..........................................................................................
      Ele_Num_Cache = 1
      do i_Out_Node = 1,Num_CrMesh_Outlines   
        if (Flag_Ver_Grow(i_Out_Node).eqv. .False.)then 
          c_Mesh_Node = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
          ! Update its coordinates.
          Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1:3) =  Final_Point(i_Out_Node,1:3)          
          ! Update the element number.
          call Cal_Ele_Num_by_Coors_3D(Final_Point(i_Out_Node,1),Final_Point(i_Out_Node,2), &
                                       Final_Point(i_Out_Node,3),Ele_Num_Cache,in_Elem_num)               
          Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(c_Mesh_Node)=in_Elem_num
          ! Update the local coordinates of the element number.
          if (in_Elem_num>0)then
              call Cal_KesiYita_by_Coor_3D(Final_Point(i_Out_Node,1:3),in_Elem_num,c_Kesi,c_Yita,c_Zeta)
              Cr3D_Meshed_Node_in_Ele_Local(i_C)%row(c_Mesh_Node,1:3) =[c_Kesi,c_Yita,c_Zeta] 
          endif
        endif
      enddo
      
      !.........................................................................................
      ! For the points where the extension occurs, increase the number of discrete nodes on the
      ! crack surface and retain other information.
      ! In addition, it also stores the crack surface discrete node numbers New_Ver_Node_Num
      ! corresponding to each old vertex.
      !.........................................................................................
      if (allocated(New_Ver_Node_Num))deallocate(New_Ver_Node_Num)
      allocate(New_Ver_Node_Num(Num_CrMesh_Outlines+1))
      New_Ver_Node_Num(1:Num_CrMesh_Outlines+1) = 0
      c_Grow_Count = 0
      Ele_Num_Cache = 1
      do i_Out_Node = 1,Num_CrMesh_Outlines   
        if (Flag_Ver_Grow(i_Out_Node).eqv. .True.)then
            c_Grow_Count = c_Grow_Count +1
            ! The crack surface discrete node numbers of the new vertices corresponding to each old vertex.
            New_Ver_Node_Num(i_Out_Node)= old_node_num+c_Grow_Count

            ! Update the number of discrete nodes
            Crack3D_Meshed_Node_num(i_C) =Crack3D_Meshed_Node_num(i_C)  +1   
            ! Check the number of discrete nodes.
            c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
            !if(Crack3D_Meshed_Ele_num(i_C) > c_Max_N_Node_3D)then
            if((old_node_num+c_Grow_Count) > c_Max_N_Node_3D)then
              !call Warning_Message('S',Keywords_Blank) 
              ! Expand memory. NEWFTU2022110501.
              call D3_Allocate_Crack_Memory(i_C,1,1)
              !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
            endif  
            
            ! Add discrete node
            Crack3D_Meshed_Node(i_C)%row(old_node_num+c_Grow_Count,1:3)= Final_Point(i_Out_Node,1:3) 
            
            ! Element number of the 3D crack node after discretization
            call Cal_Ele_Num_by_Coors_3D(Final_Point(i_Out_Node,1),Final_Point(i_Out_Node,2), &
                                         Final_Point(i_Out_Node,3),Ele_Num_Cache,in_Elem_num)               
            Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(old_node_num+c_Grow_Count)= in_Elem_num
            ! Local coordinates of the elements containing the 3D crack nodes after discretization
            if (in_Elem_num>0)then
                call Cal_KesiYita_by_Coor_3D(Final_Point(i_Out_Node,1:3),in_Elem_num,c_Kesi,c_Yita,c_Zeta)
                Cr3D_Meshed_Node_in_Ele_Local(i_C)%row(old_node_num+c_Grow_Count,1:3) =[c_Kesi,c_Yita,c_Zeta] 
            endif    
        endif
      enddo
      New_Ver_Node_Num(Num_CrMesh_Outlines+1)  = New_Ver_Node_Num(1)
      
      !....................................
      ! Handle according to the situation.
      !....................................
      do i_Out_Node = 1,Num_CrMesh_Outlines   
        Curr_Flag = Flag_Ver_Grow(i_Out_Node)
        Next_Flag = Flag_Ver_Grow(i_Out_Node+1)
        
        !////////////////////////////////////////////////////////////////////////////////////////////
        ! CASE 1: The current vertex is expanded, and the next vertex is also expanded, adding 2 new
        ! elements
        !////////////////////////////////////////////////////////////////////////////////////////////
        if(Curr_Flag .and. Next_Flag)then
          !oooooooooooooooooooooo
          ! Processing element 1
          !oooooooooooooooooooooo
          ! Update the number of discrete crack elements
          Crack3D_Meshed_Ele_num(i_C)=Crack3D_Meshed_Ele_num(i_C)+1 
          ! Check whether the number of discrete nodes is out of bounds. 2022-10-04. IMPROV2022100403.
          c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
          if(Crack3D_Meshed_Ele_num(i_C)> c_Max_N_Node_3D)then
              !call Warning_Message('S',Keywords_Blank) 
              ! Expand memory. NEWFTU2022110501.
              call D3_Allocate_Crack_Memory(i_C,1,1)
              !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
          endif            
          ! Node 1 of the newly added element.
          new_node_1 = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
          ! Node 2 of the new element: the crack surface discrete node numbers New_Ver_Node_Num corresponding
          ! to each old vertex's new vertex.
          new_node_2 = New_Ver_Node_Num(i_Out_Node)  
          ! Node 3 of the new element (the next node after the original outer boundary)
          new_node_3 = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,2)     
          ! Update the three nodes of the element.
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),1) = new_node_1   
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),2) = new_node_2 
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),3) = new_node_3 

          !ooooooooooooooooooooooooooooooo
          ! Processing the second element
          !ooooooooooooooooooooooooooooooo
          ! Update the number of discrete crack elements
          Crack3D_Meshed_Ele_num(i_C)=Crack3D_Meshed_Ele_num(i_C)+1 
          ! Check whether the number of discrete nodes is out of bounds. 2022-10-04. IMPROV2022100403.
          c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
          if(Crack3D_Meshed_Ele_num(i_C)> c_Max_N_Node_3D)then
              ! Expand memory. NEWFTU2022110501.
              call D3_Allocate_Crack_Memory(i_C,1,1)
              !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
          endif           
          ! Node 1 of the new element
          new_node_1 = New_Ver_Node_Num(i_Out_Node)  
          ! Node 2 of the new element
          new_node_2 = New_Ver_Node_Num(i_Out_Node+1)  
          ! Node 3 of the new element (the next node after the original outer boundary)
          new_node_3 =Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,2) 
          ! Three nodes of the new element
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),1) = new_node_1   
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),2) = new_node_2 
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),3) = new_node_3
          !if(i_C==2)then
          !endif            
        !/////////////////////////////////////////////////////////////////////////////
        ! CASE 2: Current vertex expands, next vertex does not expand, add 1 new cell
        !/////////////////////////////////////////////////////////////////////////////
        elseif(Curr_Flag .and.(Next_Flag.eqv..False.)) then
          ! Update the number of discrete crack elements
          Crack3D_Meshed_Ele_num(i_C)=Crack3D_Meshed_Ele_num(i_C)+1    
          ! Check whether the number of discrete nodes is out of bounds. 2022-10-04. IMPROV2022100403.
          c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
          if(Crack3D_Meshed_Ele_num(i_C)> c_Max_N_Node_3D)then
              !call Warning_Message('S',Keywords_Blank) 
              ! Expand memory. NEWFTU2022110501.
              call D3_Allocate_Crack_Memory(i_C,1,1)
              !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
          endif           
          ! Node 1 of the newly added element.
          new_node_1 = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
          ! Node 2 of the new element: the crack surface discrete node numbers New_Ver_Node_Num corresponding
          ! to each old vertex's new vertex.
          new_node_2 = New_Ver_Node_Num(i_Out_Node)  
          ! Node 3 of the new element (the next node after the original outer boundary)
          new_node_3 = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,2)     
          ! Update the three nodes of the element.
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),1) = new_node_1   
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),2) = new_node_2 
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),3) = new_node_3 
          !if(i_C==2)then
          !endif  
        !////////////////////////////////////////////////////////////////////////////////////
        ! CASE 3: The current vertex does not expand, the next vertex expands, add 1 element
        !////////////////////////////////////////////////////////////////////////////////////
        elseif((Curr_Flag.eqv..False.) .and.Next_Flag) then
          ! Update the number of discrete crack elements
          Crack3D_Meshed_Ele_num(i_C)=Crack3D_Meshed_Ele_num(i_C)+1 
          ! Check whether the number of discrete nodes is out of bounds. 2022-10-04. IMPROV2022100403.
          c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
          if(Crack3D_Meshed_Ele_num(i_C)> c_Max_N_Node_3D)then
              !call Warning_Message('S',Keywords_Blank) 
              ! Expand memory. NEWFTU2022110501.
              call D3_Allocate_Crack_Memory(i_C,1,1)
              !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
          endif           
          ! Nodes 1-3 of the new element.
          new_node_1 = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
          new_node_2 = New_Ver_Node_Num(i_Out_Node+1)  
          new_node_3 = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,2)    
          ! Update the three nodes of the element.
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),1) = new_node_1   
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),2) = new_node_2 
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),3) = new_node_3 
          !if(i_C==2)then
          !endif  
        !/////////////////////////////////////////////////////////////////////////////////////////
        ! CASE 4: The current vertex is not expanded, and the next vertex is also not expanded; 0
        ! new elements are added.
        !/////////////////////////////////////////////////////////////////////////////////////////
        elseif((Curr_Flag.eqv..False.) .and.(Next_Flag.eqv..False.)) then
          !DO NOTHING   
        endif
        
      enddo
      
      
    end select
    
    !--------------------------------------------------------------------------------------------
    ! Obtain the average of the cross products of adjacent boundary lines for the old outline in
    ! clockwise and counterclockwise directions. 2022-07-13.
    !  BUGFIX2022071303.
    !--------------------------------------------------------------------------------------------
    call D3_Get_Crack_Mesh_Outline_Clockwise(i_C,Crack3D_Meshed_Outline_num(i_C), &
            Crack3D_Meshed_Outline(i_C)%row(1:Crack3D_Meshed_Outline_num(i_C),1:2),   &
            Cros_Product_Vector_Old)
    
    !---------------------------------------------------------------------------------------------
    ! Determine the outer boundary of the fracture surface discretization elements, and store the
    ! results in the following two global variables:
    ! integer Crack3D_Meshed_Outline(Max_Num_Cr, Max_N_Node_3D, 3) !3D crack outer boundary after
    ! discretization (data 3 corresponds to the element number)
    ! integer Crack3D_Meshed_Outline_num(Max_Num_Cr) !Number of 3D crack boundary lines after
    ! discretization
    !---------------------------------------------------------------------------------------------
    call D3_Find_Crack_Mesh_Outline(i_C,Crack3D_Meshed_Ele_num(i_C) ,Crack3D_Meshed_Node_num(i_C), &
                                    Crack3D_Meshed_Ele(i_C)%row(1:Crack3D_Meshed_Ele_num(i_C),1:3))   
 
    !-------------------------------------------------------------------------------------------
    ! Check whether the new outline is consistent with the old one in the clockwise and
    ! counterclockwise directions. 2022-07-13.
    !  BUGFIX2022071303.
    ! Note: The outer boundary outline can be either clockwise or counterclockwise. If they are
    ! inconsistent, reverse it.
    !-------------------------------------------------------------------------------------------
    call D3_Get_Crack_Mesh_Outline_Clockwise(i_C,Crack3D_Meshed_Outline_num(i_C), &
             Crack3D_Meshed_Outline(i_C)%row(1:Crack3D_Meshed_Outline_num(i_C),1:2),   &
             Cros_Product_Vector_New)
    call Tool_Angle_of_Vectors_a_and_b_3D(Cros_Product_Vector_Old,Cros_Product_Vector_New,c_angle,2)
    
    
    ! If the angle is greater than 90 degrees, flip the outline.
    if (c_angle>pi/TWO) then
        call D3_Flip_Crack_Mesh_Outline(i_C)
    endif
     
    !-----------------------------------------------------------------------
    ! Check the length of each outer boundary; if it is too long, split it.
    !-----------------------------------------------------------------------
    Ele_Num_Cache = 1
    do i_outline =1,Crack3D_Meshed_Outline_num(i_C)
      
      c_P1 = Crack3D_Meshed_Outline(i_C)%row(i_outline,1) 
      c_P2 = Crack3D_Meshed_Outline(i_C)%row(i_outline,2) 
      c_P1_Coor = Crack3D_Meshed_Node(i_C)%row(c_P1,1:3)
      c_P2_Coor = Crack3D_Meshed_Node(i_C)%row(c_P2,1:3)
      c_DIS = Tool_Function_2Point_Dis_3D(c_P1_Coor,c_P2_Coor)
      if(c_DIS >= Check_L)then
        old_Ele = Crack3D_Meshed_Outline(i_C)%row(i_outline,3)
        old_Ele_N1 = Crack3D_Meshed_Ele(i_C)%row(old_Ele,1)
        old_Ele_N2 = Crack3D_Meshed_Ele(i_C)%row(old_Ele,2)
        old_Ele_N3 = Crack3D_Meshed_Ele(i_C)%row(old_Ele,3) 
        ! Determine the non-boundary nodes of the corresponding discrete elements
        if((old_Ele_N1 .ne. c_P1) .and. (old_Ele_N1 .ne. c_P2))then
            inner_Node =  old_Ele_N1 
        endif
        if((old_Ele_N2 .ne. c_P1) .and. (old_Ele_N2 .ne. c_P2))then
            inner_Node =  old_Ele_N2 
            ! IMPROV2022100301. Do not perform binary division in abnormal situations.
            print *,'    WARNING-2022071302 :: unforeseen case in Check_Crack_Grows_3D.f'    
            print *,'                          This segment will not be divided!'  
            cycle
        endif        
        if((old_Ele_N3 .ne. c_P1) .and. (old_Ele_N3 .ne. c_P2))then
            inner_Node =  old_Ele_N3 
        endif 
        
        MID_P(1:3) = (c_P1_Coor(1:3) + c_P2_Coor(1:3))/TWO
        ! Slightly move the newly added split points outward to prevent the new division edges from being
        ! collinear (collinearity is unfavorable for determining the local coordinate system of edge nodes).
        ! Angle Bisector Algorithm in D3_Crack_Vertex_Local_Coor_Sys
        Line_AB(1,1:3) =Crack3D_Meshed_Node(i_C)%row(inner_Node,1:3)
        Line_AB(2,1:3) =  MID_P(1:3)
        call Tool_Shorten_or_Extend_Line_3D(Line_AB,0.02D0*Ave_Elem_L,'B',new_Line_AB,MID_P)
        ! Add Node Number
        Crack3D_Meshed_Node_num(i_C)=Crack3D_Meshed_Node_num(i_C)+1 
        
        ! Check whether the number of discrete nodes is out of bounds. 2022-10-04. IMPROV2022100403.
        c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
        if(Crack3D_Meshed_Node_num(i_C)> c_Max_N_Node_3D)then
            !call Warning_Message('S',Keywords_Blank) 
            ! Expand memory. NEWFTU2022110501.
            call D3_Allocate_Crack_Memory(i_C,1,1)
            !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
        endif    
        ! Coordinates of the new node
        Crack3D_Meshed_Node(i_C)%row(Crack3D_Meshed_Node_num(i_C),1:3) =  MID_P(1:3)  
        ! element number of the newly added crack node
        call Cal_Ele_Num_by_Coors_3D(MID_P(1),MID_P(2),MID_P(3),Ele_Num_Cache,in_Elem_num)               
        Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(Crack3D_Meshed_Node_num(i_C)) = in_Elem_num
        ! Local coordinates of the elements containing the 3D crack nodes after discretization
        if (in_Elem_num>0)then
            call Cal_KesiYita_by_Coor_3D(MID_P(1:3),in_Elem_num,c_Kesi,c_Yita,c_Zeta)
                     Cr3D_Meshed_Node_in_Ele_Local(i_C)%row(Crack3D_Meshed_Node_num(i_C),1:3)=[c_Kesi,c_Yita,c_Zeta] 
        endif
        
        ! Original element modification node number (Add new element 1)
        Crack3D_Meshed_Ele(i_C)%row(old_Ele,1) = old_Ele_N1
        Crack3D_Meshed_Ele(i_C)%row(old_Ele,2) = Crack3D_Meshed_Node_num(i_C)
        Crack3D_Meshed_Ele(i_C)%row(old_Ele,3) = old_Ele_N3
        
        ! new element 2: Two situations should be considered. BUGFIX2022071301. Ref: My PhiPsi Development
        ! Notebook V1-P40.
        if (inner_Node ==  old_Ele_N3) then
          Crack3D_Meshed_Ele_num(i_C) =  Crack3D_Meshed_Ele_num(i_C) +1  
          ! Check whether the number of discrete nodes is out of bounds. 2022-10-04. IMPROV2022100403.
          c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
          if(Crack3D_Meshed_Ele_num(i_C)> c_Max_N_Node_3D)then
              !call Warning_Message('S',Keywords_Blank) 
              ! Expand memory. NEWFTU2022110501.
              call D3_Allocate_Crack_Memory(i_C,1,1)
              !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
          endif  
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),1) =Crack3D_Meshed_Node_num(i_C)
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),2) =old_Ele_N2
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),3) =old_Ele_N3    
        elseif(inner_Node ==  old_Ele_N1)then
          Crack3D_Meshed_Ele_num(i_C) =  Crack3D_Meshed_Ele_num(i_C) +1  
          ! Check whether the number of discrete nodes is out of bounds. 2022-10-04. IMPROV2022100403.
          c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
          if(Crack3D_Meshed_Ele_num(i_C)> c_Max_N_Node_3D)then
              ! Expand memory. BUGFIX2022122001.
              call D3_Allocate_Crack_Memory(i_C,1,1)
          endif            
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),1) =old_Ele_N2
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),2) =Crack3D_Meshed_Node_num(i_C)
          Crack3D_Meshed_Ele(i_C)%row(Crack3D_Meshed_Ele_num(i_C),3) =old_Ele_N1   
        elseif(inner_Node ==  old_Ele_N2)then 
              !call Warning_Message('S',Keywords_Blank) 
              ! For details, see IMPROV2022100301. 2022-10-03.
        endif
 
      endif
    enddo
    
    
    !---------------------------------------------------------------------------------------------
    ! Reconfirm the outer boundary of the discrete elements on the crack surface, and store the
    ! results in the following two global variables:
    ! integer Crack3D_Meshed_Outline(Max_Num_Cr, Max_N_Node_3D, 3) !3D crack outer boundary after
    ! discretization (data 3 corresponds to the element number)
    ! integer Crack3D_Meshed_Outline_num(Max_Num_Cr) !Number of 3D crack boundary lines after
    ! discretization
    !---------------------------------------------------------------------------------------------
    call D3_Find_Crack_Mesh_Outline(i_C,Crack3D_Meshed_Ele_num(i_C), &
        Crack3D_Meshed_Node_num(i_C),Crack3D_Meshed_Ele(i_C)%row(1:Crack3D_Meshed_Ele_num(i_C),1:3))    
 
    !-----------------------------------------------------------------------------------
    ! Check whether the new outline is consistent with the old one in the clockwise and
    ! counterclockwise directions. 2022-07-13.
    ! Note: The outer boundary Outline can be either clockwise or counterclockwise.
    ! BUGFIX2022071303.
    !-----------------------------------------------------------------------------------
    call D3_Get_Crack_Mesh_Outline_Clockwise(i_C,Crack3D_Meshed_Outline_num(i_C),      &
             Crack3D_Meshed_Outline(i_C)%row(1:Crack3D_Meshed_Outline_num(i_C),1:2),   &
             Cros_Product_Vector_New)
    call Tool_Angle_of_Vectors_a_and_b_3D(Cros_Product_Vector_Old,Cros_Product_Vector_New,c_angle,2)
    
    
    ! If the angle is greater than 90 degrees, flip the outline.
    if (c_angle>pi/TWO) then
        call D3_Flip_Crack_Mesh_Outline(i_C)
    endif

    100 continue   

end do
  
  
      
!**************************************************************************************
! SECTION 8: For in-plane expansion, check whether the outline intersects; if it does,
! Then it is marked as a non-propagating crack. 2022-10-29.
!                    IMPROV2022102901.
!**************************************************************************************
1531 FORMAT(5X,'INFO -- Crack, i_outline, j_outline: ',I5,I5,I5)  
! If it is not an in-plane extension, skip this section.
if(Key_InPlane_Growth/= 1)then
    goto 300
endif
! Extract the crack numbers that need to be processed.
num_Crack_Do = 0 
Crack_List(1:num_Crack) = 0
do i_C =1,num_Crack  
    if(Cracks_Allow_Propa(i_C) == 0 .or. Crack_Type_Status_3D(i_C,3)==0)then    
        cycle
    endif 
    ! If the number of boundary lines outside the crack surface is less than 20, no processing will be
    ! done.
    c_Outlines_num = Crack3D_Meshed_Outline_num(i_C)
    if(c_Outlines_num<=20) then
        cycle
    endif
    
    num_Crack_Do = num_Crack_Do  + 1
    Crack_List(num_Crack_Do) = i_C
enddo
! Crack cycling.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,c_Crack,c_Outlines_num,&
!$OMP            i_outline,c_P1_i,c_P2_i,                            &
!$OMP            c_P1_Coor_i,c_P2_Coor_i,j_outline,c_P1_j,c_P2_j,    &
!$OMP            c_P1_Coor_j,c_P2_Coor_j,c_Yes_Inter,c_InterSection_P)   
do i_C =1,num_Crack_Do   
    c_Crack = Crack_List(i_C)
    c_Outlines_num = Crack3D_Meshed_Outline_num(c_Crack)
    ! Crack outline i loop.
    do i_outline =1,c_Outlines_num
        c_P1_i = Crack3D_Meshed_Outline(c_Crack)%row(i_outline,1) 
        c_P2_i = Crack3D_Meshed_Outline(c_Crack)%row(i_outline,2) 
        c_P1_Coor_i = Crack3D_Meshed_Node(c_Crack)%row(c_P1_i,1:3)
        c_P2_Coor_i = Crack3D_Meshed_Node(c_Crack)%row(c_P2_i,1:3)
        
        ! Crack outline j loop.
        do j_outline = (i_outline+2),(c_Outlines_num-1)
            c_P1_j = Crack3D_Meshed_Outline(c_Crack)%row(j_outline,1) 
            c_P2_j = Crack3D_Meshed_Outline(c_Crack)%row(j_outline,2) 
            c_P1_Coor_j = Crack3D_Meshed_Node(c_Crack)%row(c_P1_j,1:3)
            c_P2_Coor_j = Crack3D_Meshed_Node(c_Crack)%row(c_P2_j,1:3)
            
            ! Check if they intersect.
            call Tool_Intersection_of_AB_and_CD_3D(&
                 c_P1_Coor_i,c_P2_Coor_i,c_P1_Coor_j,c_P2_Coor_j,c_Yes_Inter,c_InterSection_P) 
                 
            ! If they intersect, mark the crack as a non-propagable crack and exit the Crack loop.
            if(c_Yes_Inter) then
                Crack_Type_Status_3D(c_Crack,3) = 0
                print *,'Warning-2022102901 :: crossed outlines found!'
                write(*,1531) c_Crack,i_outline,j_outline
                print *,'    Crack ',c_Crack,' will stop propagation!'
                ! Exit the Crack loop.
                goto 288
            endif
        enddo
    enddo
    288 continue
enddo
!$OMP END PARALLEL DO

300 continue      
            
!*************************************************************************************
! SECTION 9: Calculation of the characteristics of discrete crack elements, stored in
! Crack3D_Meshed_Ele_Attri(i_Crack_Ele,1:5)
!*************************************************************************************
! NEWFTU2022093001. OpenMP Parallelization.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,i_Crack_Ele,Crack_Node1,Crack_Node2,Crack_Node3,&
!$OMP          Point1,Point2,Point3,Vector_V1,Vector_V2,Vector_Normal)  
do i_C =1,num_Crack    
  do i_Crack_Ele =1,Crack3D_Meshed_Ele_num(i_C)
      Crack_Node1 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,1)
      Crack_Node2 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,2)
      Crack_Node3 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,3)
      Point1(1:3) =Crack3D_Meshed_Node(i_C)%row(Crack_Node1,1:3)
      Point2(1:3) =Crack3D_Meshed_Node(i_C)%row(Crack_Node2,1:3)
      Point3(1:3) =Crack3D_Meshed_Node(i_C)%row(Crack_Node3,1:3)
      ! Perimeter of a discrete element (spatial triangle)
      Crack3D_Meshed_Ele_Attri(i_C)%row(i_Crack_Ele,1)  = Tool_Function_2Point_Dis_3D(Point1,Point2) +  &
                                                          Tool_Function_2Point_Dis_3D(Point2,Point3) +  &
                                                          Tool_Function_2Point_Dis_3D(Point1,Point3)
      ! Area of discrete elements (spatial triangles)
      call Tool_Area_Tri_3D(Point1,Point2,Point3,Crack3D_Meshed_Ele_Attri(i_C)%row(i_Crack_Ele,2))
      ! Outer normal vector Crack3D_Meshed_Ele_Nor_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3)
      Vector_V1 = Point2- Point1
      Vector_V2 = Point3- Point1
      call Vector_Normalize(3,Vector_V1)  
      call Vector_Normalize(3,Vector_V2)  
      call Vector_Cross_Product_3(Vector_V1,Vector_V2,Vector_Normal) 
      call Vector_Normalize(3,Vector_Normal)  
      ! Check whether Vector_Normal is reasonable. 2022-07-10.
      if(sum(abs(Vector_Normal))<=Tol_10)then
          print *,'    Error-2023081601 :: illegal Vector_Normal!'
          print *,'                        in Check_Crack_Grows_3D.f90'
          print *,'                        i_C,i_Crack_Ele',i_C,i_Crack_Ele
          print *,'Point1:',Point1
          print *,'Point2:',Point2
          print *,'Point3:',Point3
          print *,'Crack_Node1:',Crack_Node1
          print *,'Crack_Node2:',Crack_Node2
          print *,'Crack_Node3:',Crack_Node3
          print *,'Vector_V1:',Vector_V1
          print *,'Vector_V2:',Vector_V2
          call Warning_Message('S',Keywords_Blank)  
      endif
      Crack3D_Meshed_Ele_Nor_Vector(i_C)%row(i_Crack_Ele,1:3)= Vector_Normal
  enddo
enddo
!$OMP END PARALLEL DO


!******************************************************************************
! SECTION 10: Outward Normal Vectors of 3D Fracture Nodes After Discretization
!   Crack3D_Meshed_Node_Nor_Vector(Max_Num_Cr_3D,Max_N_Node_3D,3) 
!******************************************************************************
! NEWFTU2022093001. OpenMP Parallelization.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,i_Crack_Node,c_count_Ele,c_All_Nor_Vector,i_Crack_Ele,&
!$OMP               Crack_Node1,Crack_Node2,Crack_Node3,c_Max_N_Node_3D)   
do i_C =1,num_Crack   
  c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
  Crack3D_Meshed_Node_Nor_Vector(i_C)%row(1:c_Max_N_Node_3D,1:3)   = ZR
  do i_Crack_Node =1,Crack3D_Meshed_Node_num(i_C) 
      c_count_Ele = 0
      c_All_Nor_Vector = [ZR,ZR,ZR]
      do i_Crack_Ele =1,Crack3D_Meshed_Ele_num(i_C)
          Crack_Node1=Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,1)
          Crack_Node2=Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,2)
          Crack_Node3=Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,3)
          if(i_Crack_Node == Crack_Node1 .or.i_Crack_Node == Crack_Node2 .or.i_Crack_Node == Crack_Node3 )then
              c_All_Nor_Vector  = c_All_Nor_Vector +Crack3D_Meshed_Ele_Nor_Vector(i_C)%row(i_Crack_Ele,1:3)
              c_count_Ele = c_count_Ele +1
          endif
      enddo
      Crack3D_Meshed_Node_Nor_Vector(i_C)%row(i_Crack_Node,1:3)  = c_All_Nor_Vector/c_count_Ele
      call Vector_Normalize(3,Crack3D_Meshed_Node_Nor_Vector(i_C)%row(i_Crack_Node,1:3))
      
      if(sum(abs(c_All_Nor_Vector/c_count_Ele))<=Tol_10)then
          !IMPROV2022092801.
          ! For some interactive overlapping elements (\changelog.src\2022-09-28-01.png),
          ! Normal directions may cancel each other out, causing c_All_Nor_Vector to be 0.
           
          do i_Crack_Ele =1,Crack3D_Meshed_Ele_num(i_C)
              Crack_Node1=Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,1)
              Crack_Node2=Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,2)
              Crack_Node3=Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,3)
              if(i_Crack_Node == Crack_Node1 .or.i_Crack_Node == Crack_Node2 .or.i_Crack_Node == Crack_Node3 )then
                  Crack3D_Meshed_Node_Nor_Vector(i_C)%row(i_Crack_Node,1:3)  = &
                                    Crack3D_Meshed_Ele_Nor_Vector(i_C)%row(i_Crack_Ele,1:3)
                  exit
              endif
          enddo
          
          ! 2022-10-28. Marked as non-propagating crack. IMPROV2022102801.
          if (Crack_Type_Status_3D(i_C,3) /= 0) then  
              print *,'    WARNING-2022092801 :: possible overlaped crack front!'
              print *,'    Info -- i_C,i_Crack_Node:',i_C,i_Crack_Node 
              
              Crack_Type_Status_3D(i_C,3) = 0
              print *,'    This crack will stop propagation!'
          endif
      endif
      
  enddo
enddo
!$OMP END PARALLEL DO
  
!***********************************************************************
! SECTION 11: Centroid of the 3D fracture surface after discretization.
!***********************************************************************
! NEWFTU2022093001. OpenMP Parallelization.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,num_Cr_Nodes,centroid_x,centroid_y,centroid_z)   
do i_C =1,num_Crack  
    num_Cr_Nodes = Crack3D_Meshed_Node_num(i_C)
    centroid_x=sum(Crack3D_Meshed_Node(i_C)%row(1:num_Cr_Nodes,1))/dble(num_Cr_Nodes)
    centroid_y=sum(Crack3D_Meshed_Node(i_C)%row(1:num_Cr_Nodes,2))/dble(num_Cr_Nodes)
    centroid_z=sum(Crack3D_Meshed_Node(i_C)%row(1:num_Cr_Nodes,3))/dble(num_Cr_Nodes)
    Crack3D_Centroid(i_C,1) = centroid_x
    Crack3D_Centroid(i_C,2) = centroid_y
    Crack3D_Centroid(i_C,3) = centroid_z
enddo
!$OMP END PARALLEL DO

  
!***************************************
! SECTION 12: Close the relevant files.
!***************************************
if(Key_Save_Nothing  /= 1) then 
    if(CFCP==2)then
        close(101)
        close(102)
        close(103)
        close(104)
    elseif(CFCP==3)then
        close(104)
    endif          
endif
      
RETURN
END SUBROUTINE Check_Crack_Grows_3D