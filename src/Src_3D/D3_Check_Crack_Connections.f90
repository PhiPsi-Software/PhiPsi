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
 
SUBROUTINE D3_Check_Crack_Connections(i_WB,i_Stage,i_Prop,isub)
! Determine the connectivity of 3D cracks and update the crack types based on the connectivity.
!2022-06-10.

! The global variable Crack_Type_Status_3D(i_C,10) is used to indicate the type and status of
! cracks.
! Column 1 (Fracture Type): =1, HF Fracture; =2, Natural Fracture; =3, Post-Fracturing Hydraulic
! Fracture
! Note: Natural fractures and propped hydraulic fractures may potentially turn into HF fractures.
! Column 2 (Fracture Status): =1, HF fracturing not completed; =2, HF fracturing completed
! Column 3 (Can the crack continue to propagate): =1, yes; =0, no
! Column 4 (Whether the fracture has obtained a fluid node): =1, Yes; =0, No
! Column 5 (Did the crack propagate in the previous step?): =1, Yes; =0, No
! Column 6 (if activated from natural fractures secondarily, it corresponds to the number of the
! respective natural fracture)
!   
! Global variable NaCr3D_Status(:,:), stores the activation status and other state variables for
! each 3D natural fracture. 2023-01-09.
! Column 1 stores the active status;
! Column 2 records the number of solid cells containing the natural fracture.

!**********************************
! Read the public variable module.
!**********************************
use Global_Float_Type
use Global_Crack_Common
use Global_Crack_3D
use Global_Model
use Global_Common
use Global_DISP

!***********************
! Variable declaration.
!***********************
implicit none
integer,intent(in)::i_WB,i_Stage,i_Prop,isub
integer i_C,j_C
logical Crack_Coor_Overlap_Status(num_Crack,num_Crack)
integer i_Cr_Node1,i_Cr_Node2,i_Cr_Node3
integer j_Cr_Node1,j_Cr_Node2,j_Cr_Node3
integer j_Crack_Ele
real(kind=FT) c_Tri_1(3,3),c_Tri_2(3,3)
real(kind=FT) c_Inter_Points(10,3)
logical c_Logical_Inter
integer c_num_Inters
integer i_Crack_Ele
logical c_Logical_Parallel
logical Yes_Overlap,Yes_Connect
real(kind=FT) c_HF_Coor_Ranges(3,2),c_NF_Coor_Ranges(3,2)
logical,ALLOCATABLE::NF_Flag_Active(:)
integer,ALLOCATABLE::NF_Active_HF(:)
real(kind=FT),ALLOCATABLE::NF_Average_Inter_Points(:,:)
integer num_New_Cracks
real(kind=FT) New_Crack_Center(3),n_vector_NF(3)
integer Ele_Num_Cache,c_OUT_Elem
real(kind=FT) NF_Diameter
integer c_num_Element,c_List(10000)
real(kind=FT) HF_NF_Inter_Point(3)
real(kind=FT) c_Kesi,c_Yita,c_Zeta
real(kind=FT) c_Sxx,c_Syy,c_Szz,c_Sxy,c_Syz,c_Sxz
real(kind=FT) Normal_Stress

1001 FORMAT(5X,'Activated natural fractures include NF',I5)   
                                        
print *,"    Counting cracks' connections......"

! Perform different calculations based on the Key_NaCr_Active_Scheme_3D.
select case(Key_NaCr_Active_Scheme_3D)

!////////////////////////////////////////////////////
!                                                   /
!                                                   /
!                                                   /
! Key_NaCr_Active_Scheme_3D = 1                     /
!                                                   /
!                                                   /
! All natural fractures are activated from the      /
! initial moment,and all natural fractures are equal/
! For actual natural fracture sizes, after          /
! communication with HF, the fractures are filled   /
! with fracturing fluid. Default.                   /
!                                                   /
!                                                   /
!                                                   /
!                                                   /
!////////////////////////////////////////////////////
case(1)
    !*************************
    ! Variable Initialization
    !*************************
    ! 3D Crack Intersection State
    if (allocated(Cracks_3D_Inter_Status)) deallocate(Cracks_3D_Inter_Status)
    allocate(Cracks_3D_Inter_Status(num_Crack,num_Crack))
    Cracks_3D_Inter_Status(1:num_Crack,1:num_Crack) = 0

    !###############################################################
    ! Obtain the coordinate range of each discrete fracture surface
    ! and store it in Crack_Coor_Range(num_Crack, 3, 2)
    !###############################################################
    !call D3_Get_Cracks_Coor_Ranges(Crack_Coor_Range)
    call D3_Get_Cracks_Coor_Ranges


    !###############################################################
    ! Determine whether the coordinate ranges between the various
    ! cracks overlap.
    ! Store into Crack_Coor_Overlap_Status(1:num_Crack,1:num_Crack)
    !###############################################################
    !call D3_Get_Cracks_Coor_Overlap_Status(Crack_Coor_Range,Crack_Coor_Overlap_Status)
    call D3_Get_Cracks_Coor_Overlap_Status(Crack_Coor_Overlap_Status)


    !####################################################
    ! Determine whether there is an intersection through 
    ! the intersection points of the discrete triangles 
    ! on the fracture surface
    !####################################################
    do i_C=1,num_Crack
        do j_C=1,num_Crack
        ! If there is an overlap in the coordinate ranges of the two cracks.
        if(Crack_Coor_Overlap_Status(i_C,j_C) .eqv. .True.) then
            ! num_i_j_Inter = 0 !The total number of intersection points of these two cracks
            do i_Crack_Ele =1,Crack3D_Meshed_Ele_num(i_C)
              i_Cr_Node1 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,1)
              i_Cr_Node2 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,2)
              i_Cr_Node3 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,3)
              c_Tri_1(1,1:3)=Crack3D_Meshed_Node(i_C)%row(i_Cr_Node1,1:3)
              c_Tri_1(2,1:3)=Crack3D_Meshed_Node(i_C)%row(i_Cr_Node2,1:3)
              c_Tri_1(3,1:3)=Crack3D_Meshed_Node(i_C)%row(i_Cr_Node3,1:3)    
              
              ! Discrete triangular loop of j_C cracks
              do j_Crack_Ele =1,Crack3D_Meshed_Ele_num(j_C)
                j_Cr_Node1 = Crack3D_Meshed_Ele(j_C)%row(j_Crack_Ele,1)
                j_Cr_Node2 = Crack3D_Meshed_Ele(j_C)%row(j_Crack_Ele,2)
                j_Cr_Node3 = Crack3D_Meshed_Ele(j_C)%row(j_Crack_Ele,3)
                c_Tri_2(1,1:3)=Crack3D_Meshed_Node(j_C)%row(j_Cr_Node1,1:3)
                c_Tri_2(2,1:3)=Crack3D_Meshed_Node(j_C)%row(j_Cr_Node2,1:3)
                c_Tri_2(3,1:3)=Crack3D_Meshed_Node(j_C)%row(j_Cr_Node3,1:3)                
                ! Check whether the discrete fracture surface triangles of the two cracks intersect.
                call Tool_Intersections_of_Two_Triangles_3D(c_Tri_1,c_Tri_2,c_Logical_Inter,c_num_Inters, &
                    c_Inter_Points,c_Logical_Parallel)  
                ! If there is an intersection.
                if(c_Logical_Inter .eqv. .True.)then        
                    Cracks_3D_Inter_Status(i_C,j_C) = 1
                    goto 777
                endif
              enddo
            enddo
        endif
        enddo
    777   continue
    enddo


    !#####################################################################
    ! Update the fracture type based on the fracture intersection status: 
    ! originally a natural fracture, updated to an HF fracture after 
    ! communication with HF.
    !#####################################################################
    1101 FORMAT(5X,'Natural crack ',I4,' connected by HF crack ',I4)      
    do i_C=1,num_Crack
        do j_C=1,num_Crack
        ! If i_C or j_C was originally a natural fracture
        if(Crack_Type_Status_3D(i_C,1)==2 .or. Crack_Type_Status_3D(j_C,1)==2) then
              ! If i_C and j_C intersect
              if(Cracks_3D_Inter_Status(i_C,j_C) == 1 .or. Cracks_3D_Inter_Status(j_C,i_C) == 1) then
                
                ! If i_C is an HF crack and j_C is an NF. 2022-09-10. BUGFIX2022091001.
                if(Crack_Type_Status_3D(i_C,1)==1 .and. Crack_Type_Status_3D(j_C,1)==2) then
                    Crack_Type_Status_3D(j_C,1)=1
                    Crack_Type_Status_3D(j_C,2)=1
                    Crack_Type_Status_3D(j_C,3)=1
                    Crack_Type_Status_3D(j_C,4)=0
                    Key_CS_Crack(j_C) = 0
                    write(*,1101) j_C,i_C
                endif 
                
                ! If j_C is an HF crack and i_C is an NF. 2022-09-10. BUGFIX2022091001.
                if(Crack_Type_Status_3D(j_C,1)==1 .and. Crack_Type_Status_3D(i_C,1)==2) then
                    Crack_Type_Status_3D(i_C,1)=1
                    Crack_Type_Status_3D(i_C,2)=1
                    Crack_Type_Status_3D(i_C,3)=1
                    Crack_Type_Status_3D(i_C,4)=0
                    Crack_Type_Status_3D(i_C,6)=j_C
                    Key_CS_Crack(i_C) = 0
                    write(*,1101) i_C,j_C
                endif                
              endif
        endif
        enddo
    enddo
    
!////////////////////////////////////////////////////
!                                                   /
!                                                   /
!                                                   /
! Key_NaCr_Active_Scheme_3D = 2                     /
!                                                   /
!                                                   /
! Activated only after communication with HF; after /
! communication with HF, the fractures are filled   /
! with fracturing fluid.                            /
! Only supports rectangular or polygonal natural    /
! fractures.                                        /
!                                                   /
! NEWFTU2023010901. 2023-01-09                      /
!                                                   /
!                                                   /
!                                                   /
!////////////////////////////////////////////////////
case(2)
    !*********************
    ! Temporary variable.
    !*********************
    allocate(NF_Flag_Active(num_Rand_Na_Crack))
    allocate(NF_Active_HF(num_Rand_Na_Crack))
    NF_Flag_Active(1:num_Rand_Na_Crack) = .False.

    !**********************************************
    ! Loop through each crack to check whether the 
    ! cracks intersect with natural fractures.
    ! If there is, then mark NF_Flag_Active.
    !**********************************************
    !BUGFIX2023092301.
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,c_HF_Coor_Ranges,j_C,c_NF_Coor_Ranges, &
    !$OMP                                     Yes_Overlap,Yes_Connect,c_num_Inters,c_Inter_Points)  
    do i_C=1,num_Crack
        ! If it's an HF crack. BUGFIX2023021903.
        if(Crack_Type_Status_3D(i_C,1) == 1) then 
            !...............................................................
            ! Obtain the coordinate range of the HF crack c_HF_Coor_Ranges.
            !...............................................................
            call D3_Get_Crack_Coor_Ranges(i_C,c_HF_Coor_Ranges)
            ! Natural fracture circulation
            do j_C=1,num_Rand_Na_Crack
                ! If the natural fracture is not activated
                if(NaCr3D_Status(j_C,1) /=1)then
                    ! Obtain the coordinate range of natural fractures.
                    call D3_Get_NF_Coor_Ranges(j_C,c_NF_Coor_Ranges)
                    
                    ! Check if there is an overlap between the HF and NF coordinate ranges.
                    call D3_Check_HF_NF_Overlap_Status(c_HF_Coor_Ranges,c_NF_Coor_Ranges,Yes_Overlap)
                    
                    ! If the HF and NF coordinate ranges do not overlap, proceed to the next loop.
                    if (Yes_Overlap .eqv. .False.) then
                        CYCLE
                    endif
            
                    ! Check whether hydraulic fractures and natural fractures intersect.
                    call D3_Check_hf_nf_connect(i_C,j_C,Yes_Connect,c_num_Inters,c_Inter_Points)
                    
                    ! If HF and NF do not intersect, proceed to the next loop.
                    if (Yes_Connect .eqv. .False.) then
                        CYCLE
                    endif
                    
                    ! Mark the natural fracture as an activated natural fracture.
                    NF_Flag_Active(j_C) = .True.
                    
                    !HF crack number. 2023-09-22. BUGFIX2023092201.
                    NF_Active_HF(j_C)   = i_C
                endif
            enddo
        endif
    enddo
    !$OMP END PARALLEL DO    
 
    !***********************
    ! Information exchange.
    !***********************
    1111 FORMAT(5X,'Number of newly activated natural fractures: ',I6)  
    num_New_Cracks = count(NF_Flag_Active(1:num_Rand_Na_Crack))
    if(num_New_Cracks >= 1) then
        print *,'    Attention :: New natural fractures are activated!'
        write(*,1111) num_New_Cracks
    endif
    ! If the number of activated natural fractures is less than 5, list them one by one. 2024-02-22.
    if(num_New_Cracks<=5) then
        do j_C=1,num_Rand_Na_Crack
            if (NF_Flag_Active(j_C)) write (*,1001) j_C
        enddo
    endif
    
    !*******************************************************************
    ! Update the number of cracks, crack coordinates, and crack status.
    !*******************************************************************
    do i_C=1,num_Rand_Na_Crack
        if(NF_Flag_Active(i_C).eqv. .True.) then
            ! Update the number of cracks
            num_Crack            = num_Crack + 1
            ! Update crack coordinates
            Each_Cr_Poi_Num(num_Crack)   = Each_NaCr3D_Poi_Num(i_C)
            Crack3D_Coor(num_Crack,1:Each_NaCr3D_Poi_Num(i_C),1:3) = Na_Crack3D_Coor(i_C,1:Each_NaCr3D_Poi_Num(i_C),1:3)
            ! Update crack status data
            NaCr3D_Status(i_C,1) = 1
            Crack_Type_Status_3D(num_Crack,1) = 1
            Crack_Type_Status_3D(num_Crack,2) = 1
            Crack_Type_Status_3D(num_Crack,3) = 1
            Crack_Type_Status_3D(num_Crack,4) = 0
            Crack_Type_Status_3D(num_Crack,6) = i_C
            Key_CS_Crack(num_Crack)           = 0
            
            !Apply pressure to this natural crack. 2023-09-22. BUGFIX2023092201.
            Crack_Pressure(num_Crack) = Crack_Pressure(NF_Active_HF(i_C))
            
        endif
    enddo
    

   !****************************
   ! Clear temporary variables.
   !****************************
   if(allocated(NF_Flag_Active)) deallocate(NF_Flag_Active)


   
!////////////////////////////////////////////////////
!                                                   /
!                                                   /
! Key_NaCr_Active_Scheme_3D = 3                     /
!                                                   /
!                                                   /
!                                                   /
! Activated only after communication with HF; after /
! communication with HF, only some cracks open,     /
! along the                                         /
! The surface where natural cracks are located      /
! expands, and the fracture toughness of elements   /
! through which the natural cracks pass is low.     /
! Only supports rectangular or polygonal natural    /
! fractures.                                        /
!                                                   /
! NEWFTU2023011101. 2023-01-11.                     /
!                                                   /
!////////////////////////////////////////////////////
case(3)
    !*********************
    ! Temporary variable.
    !*********************
    allocate(NF_Flag_Active(num_Rand_Na_Crack))
    NF_Flag_Active(1:num_Rand_Na_Crack) = .False.
    allocate(NF_Active_HF(num_Rand_Na_Crack))
    allocate(NF_Average_Inter_Points(num_Rand_Na_Crack,3))
    NF_Average_Inter_Points =  ZR
    
    !***************************************************************************************
    ! Loop through each crack to check whether the cracks intersect with natural fractures.
    ! If there is, then mark NF_Flag_Active.
    !***************************************************************************************
    !BUGFIX2023092301.
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,c_HF_Coor_Ranges,j_C,c_NF_Coor_Ranges, &
    !$OMP                                     Yes_Overlap,Yes_Connect,c_num_Inters,c_Inter_Points,&
    !$OMP                                     HF_NF_Inter_Point,Ele_Num_Cache,c_OUT_Elem,&
    !$OMP                                     c_Kesi,c_Yita,c_Zeta,c_Sxx,c_Syy,c_Szz,c_Sxy,c_Syz,c_Sxz,&
    !$OMP                                     n_vector_NF,Normal_Stress)  
    do i_C=1,num_Crack
        ! If it's an HF crack. 2023-02-19. BUGFIX2023021903.
        if(Crack_Type_Status_3D(i_C,1) == 1) then 
            !...............................................................
            ! Obtain the coordinate range of the HF crack c_HF_Coor_Ranges.
            !...............................................................
            call D3_Get_Crack_Coor_Ranges(i_C,c_HF_Coor_Ranges)
            ! Natural Fracture Circulation
            do j_C=1,num_Rand_Na_Crack
                ! If the natural fracture is not activated
                if(NaCr3D_Status(j_C,1) /=1)then
                    ! Obtain the coordinate range of natural fractures.
                    call D3_Get_NF_Coor_Ranges(j_C,c_NF_Coor_Ranges)
                    
                    ! Check if there is an overlap between the HF and NF coordinate ranges.
                    call D3_Check_HF_NF_Overlap_Status(c_HF_Coor_Ranges,c_NF_Coor_Ranges,Yes_Overlap)
                    
                    ! If the HF and NF coordinate ranges do not intersect, proceed to the next loop.
                    if (Yes_Overlap .eqv. .False.) then
                        CYCLE
                    endif
            
                    ! Check whether hydraulic fractures and natural fractures intersect.
                    call D3_Check_hf_nf_connect(i_C,j_C,Yes_Connect,c_num_Inters,c_Inter_Points)
                    
                    ! If HF and NF do not intersect, proceed to the next loop.
                    if (Yes_Connect .eqv. .False.) then
                        CYCLE
                    endif
                    
                    ! Obtain the average coordinates of the intersection. 2024-02-22. IMPROV2024022201.
                    HF_NF_Inter_Point(1) = sum(c_Inter_Points(1:c_num_Inters,1))/c_num_Inters
                    HF_NF_Inter_Point(2) = sum(c_Inter_Points(1:c_num_Inters,2))/c_num_Inters
                    HF_NF_Inter_Point(3) = sum(c_Inter_Points(1:c_num_Inters,3))/c_num_Inters
                    
                    
                    ! Obtain the element number where the intersection is located and store it in c_OUT_Elem.
                    ! 2024-02-22. IMPROV2024022201.
                    Ele_Num_Cache = 1
                    call Cal_Ele_Num_by_Coors_3D(HF_NF_Inter_Point(1),HF_NF_Inter_Point(2),HF_NF_Inter_Point(3),&
                                                 Ele_Num_Cache,c_OUT_Elem)
                    
                    ! Obtain the coordinates in the local coordinate system of the intersection. 2024-02-22.
                    ! IMPROV2024022201.
                    call Cal_KesiYita_by_Coor_3D(HF_NF_Inter_Point,c_OUT_Elem,c_Kesi,c_Yita,c_Zeta)
                    ! Check the stress state of natural fractures and calculate the stress state at the intersection of
                    ! hydraulic fractures and natural fractures (HF_NF_Inter_Point). 2024-02-22.
                    
                    ! Calculate the stress state at the intersection. 2024-02-22. IMPROV2024022201.
                    call Cal_Any_Point_Str_KesiYita_3D(c_OUT_Elem,0,1,1,c_Kesi,c_Yita,c_Zeta,1,DISP, &
                                                       c_Sxx,c_Syy,c_Szz,c_Sxy,c_Syz,c_Sxz)
                  
                    ! Calculate the outward normal vector of natural fractures. 2024-02-22. IMPROV2024022201.
                    call Tool_Normal_vector_of_3D_Tri(Na_Crack3D_Coor(j_C,1,1:3),&
                                                      Na_Crack3D_Coor(j_C,2,1:3),&
                                                      Na_Crack3D_Coor(j_C,3,1:3),n_vector_NF(1:3))
                    
                    ! Calculate the normal stress on natural fracture surfaces. 2024-02-22. IMPROV2024022201.
                    call Tool_Get_Normal_Stress_on_Plane_by_Stress_Tensor([c_Sxx,c_Syy,c_Szz,c_Sxy,c_Syz,c_Sxz], &
                                                                          n_vector_NF,Normal_Stress)
                    
                    !IMPROV2024022201.
                    ! Check whether the water pressure of the HF can open the natural fracture intersection points.
                    ! 2024-02-22.
                    ! If the net pressure exceeds the tensile strength of the natural fracture, the natural fracture
                    ! will open and be activated.
                    ! That is: HF water pressure - normal stress on natural fractures (positive value) > tensile
                    ! strength of the natural fractures.
                    ! Explanation: In fact, the stress calculated by Cal_Any_Point_Str_KesiYita_3D has already taken the
                    ! effect of HF into account.
                    ! However, since NF and HF did not communicate, the impact of natural fractures is difficult to
                    ! fully reflect. Therefore, in order to be more consistent with
                    ! In practice, we consider the water pressure applied by HF as the actual water pressure within NF,
                    ! which is then used to calculate the net water pressure. 2024-02-22.
                    if ((Crack_Pressure(i_C) - (-Normal_Stress)) < St_NaCr(j_C)) then
                        cycle
                    else
                        !do nothing.
                        ! The natural fracture has been activated.
                    endif
                    
                    ! Mark the natural fracture as an activated natural fracture.
                    NF_Flag_Active(j_C) = .True.
                    
                    !HF crack number. 2023-09-22. BUGFIX2023092201.
                    NF_Active_HF(j_C)   = i_C
                    
                    ! Average coordinates of the HF and NF intersection. 2023-01-11.
                    NF_Average_Inter_Points(j_C,1:3) = HF_NF_Inter_Point(1:3)
                endif
            enddo
        endif
    enddo
    !$OMP END PARALLEL DO    
 
    !***********************
    ! Information exchange.
    !***********************
    2111 FORMAT(5X,'Number of newly activated natural fractures: ',I6)  
    num_New_Cracks = count(NF_Flag_Active(1:num_Rand_Na_Crack))
    if(num_New_Cracks >= 1) then
        print *,'    Attention :: New natural fractures are activated!'
        write(*,2111) num_New_Cracks
    endif
    ! If the number of activated natural fractures is less than 5, list them one by one. 2024-02-22.
    if(num_New_Cracks<=5) then
        do j_C=1,num_Rand_Na_Crack
            if (NF_Flag_Active(j_C)) write (*,1001) j_C
        enddo
    endif

    !*******************************************************************
    ! Update the number of cracks, crack coordinates, and crack status.
    !*******************************************************************
    do i_C=1,num_Rand_Na_Crack
        if(NF_Flag_Active(i_C).eqv. .True.) then
            !------------------------------
            ! Update the number of cracks.
            !------------------------------
            num_Crack  = num_Crack + 1
            
            !---------------------------
            !
            ! Update crack coordinates.
            !
            !---------------------------
            ! STEP 1: Center coordinates of the new fracture: obtained by 
            ! averaging the coordinates of the intersection points between
            ! hydraulic fractures and natural fractures (there may be more
            ! than one).
            New_Crack_Center = NF_Average_Inter_Points(i_C,1:3)
            
            ! STEP 2: Obtain the cell number where the intersection is located.
            Ele_Num_Cache = 1
            call Cal_Ele_Num_by_Coors_3D(New_Crack_Center(1),New_Crack_Center(2),New_Crack_Center(3),Ele_Num_Cache,c_OUT_Elem)
            if(c_OUT_Elem<=0)then
                print *, "    Error-2023011101 :: illegal c_OUT_Elem!" 
                print *, "                        In D3_Check_Crack_Connections.f90."  
                print *, "                        Point Coordinate:",New_Crack_Center(1:3)
                call Warning_Message('S',Keywords_Blank) 
            endif
            
            ! STEP 3: Calculate the outward normal vector of the natural fracture.
            call Tool_Normal_vector_of_3D_Tri(Na_Crack3D_Coor(i_C,1,1:3),&
                                              Na_Crack3D_Coor(i_C,2,1:3),&
                                              Na_Crack3D_Coor(i_C,3,1:3),n_vector_NF(1:3))
                                              
            ! STEP 4: Calculate the diameter of the generated circular crack (Size_Factor_of_Active_NaCr is the
            ! size factor).
            NF_Diameter = (Elem_Vol(c_OUT_Elem)**(1.0D0/3.0D0)) * Size_Factor_of_Active_NaCr
            
            ! STEP 5: Initial circular crack.
            Crack3D_Cir_Coor(num_Crack,1:7) = [New_Crack_Center(1),New_Crack_Center(2),New_Crack_Center(3),&
                                               n_vector_NF(1),n_vector_NF(2),n_vector_NF(3),NF_Diameter/TWO]  
                                               
            ! STEP 6: Obtain the list of solid elements containing natural cracks and store it in
            ! Na_Crack3D_Ele_List, with the number of solid elements stored in NaCr3D_Status(i_C,2).
            ! NEWFTU2023011202.
            call D3_Get_Element_List_of_Natural_Crack(i_C,c_num_Element,c_List(1:10000))
            if(.not. allocated(Na_Crack3D_Ele_List(i_C)%row))then
                allocate(Na_Crack3D_Ele_List(i_C)%row(c_num_Element))
                Na_Crack3D_Ele_List(i_C)%row(1:c_num_Element) = c_List(1:c_num_Element)
            endif
            NaCr3D_Status(i_C,2) = c_num_Element
            
            !--------------------------
            ! Update crack status data
            !--------------------------
            NaCr3D_Status(i_C,1) = 1
            Crack_Type_Status_3D(num_Crack,1)=1
            Crack_Type_Status_3D(num_Crack,2)=1
            Crack_Type_Status_3D(num_Crack,3)=1
            Crack_Type_Status_3D(num_Crack,4)=0
            Crack_Type_Status_3D(num_Crack,6)=i_C
            Key_CS_Crack(num_Crack)          =0
            
            !Apply pressure to this natural crack. 2023-09-22. BUGFIX2023092201.
            Crack_Pressure(num_Crack) = Crack_Pressure(NF_Active_HF(i_C))
        endif
    enddo
    

   !****************************
   ! Clear temporary variables.
   !****************************
   if(allocated(NF_Flag_Active)) deallocate(NF_Flag_Active)
   if(allocated(NF_Average_Inter_Points)) deallocate(NF_Average_Inter_Points)


endselect

RETURN
END SUBROUTINE D3_Check_Crack_Connections
