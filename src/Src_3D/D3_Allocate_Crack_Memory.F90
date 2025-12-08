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
 
SUBROUTINE D3_Allocate_Crack_Memory(i_C,Value_Type,Key_Extend)
! Allocate memory or expand memory for crack-related variables. 2022-09-03.
!
! Value_Type=1, related to discrete fracture surfaces. Max_N_Node_3D.
! Value_Type=2, related to fluid elements. Max_N_FluEl_3D.
! Value_Type=3, related to fluid nodes. Max_N_CalP_3D.
!
! Added Key_Extend: 1.5 times expanded space. 2022-11-05. IMPROV2022110502.
!
! Memory Expansion Example:
! Expand the dimension of Ragged_Array(3)%row
! allocate(tem_Ragged_Array(10))
! tem_Ragged_Array(1:size(Ragged_Array(3)%row)) = Ragged_Array(3)%row
! deallocate(Ragged_Array(3)%row)
! Note: The dimension of Ragged_Array(3)%row has been doubled, and tem_Ragged_Array is automatically
! deallocated
! call move_alloc(tem_Ragged_Array, Ragged_Array(3)%row)    

!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common,only: KEYWORDS_BLANK
use Global_Crack_3D
use Global_Post

implicit none
integer,intent(in):: i_C,Value_Type,Key_Extend
integer c_Max_N_Node_3D,c_Max_N_Node_3D_New
integer c_Max_N_FluEl_3D,c_Max_N_FluEl_3D_New
integer c_Max_N_CalP_3D,c_Max_N_CalP_3D_New
integer c_Max_ele_num_CalP,c_Max_ele_num_CalP_New
real(kind=FT) Scale_Factor
real(kind=FT),ALLOCATABLE::tem_vector_1(:)    
real(kind=FT),ALLOCATABLE::tem_vector_2(:,:)     
real(kind=FT),ALLOCATABLE::tem_vector_3(:,:,:)     
integer,ALLOCATABLE::tem_Int_vector_1(:) 
integer,ALLOCATABLE::tem_Int_vector_2(:,:)     
integer,ALLOCATABLE::tem_Int_vector_3(:,:,:)    

Scale_Factor = 1.5D0

!-------------------------------------------------------------------------------
! Allocate memory for different jagged array variables according to Value_Type.
!-------------------------------------------------------------------------------
select case(Value_Type)      

!/////////////////////////////////////////////////////
!                                                   /
!Value_Type=1,related to discrete fracture surfaces./
!                                                   /
!/////////////////////////////////////////////////////
case(1)
    !||||||||||||||||||||||||||||
    !                          |
    ! If it is a new variable. |
    !                          |
    !||||||||||||||||||||||||||||
    if(Key_Extend==0) then
        if(.not. allocated(Crack3D_Meshed_Node(i_C)%row))then
            allocate(Crack3D_Meshed_Node(i_C)%row(Max_N_Node_3D(i_C),3))
            Crack3D_Meshed_Node(i_C)%row(1:Max_N_Node_3D(i_C),1:3) = ZR
        endif
        if(.not. allocated(Crack3D_Meshed_Ele(i_C)%row))then
            allocate(Crack3D_Meshed_Ele(i_C)%row(Max_N_Node_3D(i_C),3))
            Crack3D_Meshed_Ele(i_C)%row(1:Max_N_Node_3D(i_C),1:3)  = 0
        endif 
        if(.not. allocated(Crack3D_Meshed_Node_Value(i_C)%row))then
            allocate(Crack3D_Meshed_Node_Value(i_C)%row(Max_N_Node_3D(i_C),3))
            Crack3D_Meshed_Node_Value(i_C)%row(1:Max_N_Node_3D(i_C),1:3) = ZR
        endif     
        if(.not. allocated(Cr3D_Meshed_Node_in_Ele_Num(i_C)%row))then
            allocate(Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(Max_N_Node_3D(i_C)))
            Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(1:Max_N_Node_3D(i_C)) = 0
        endif         
        if(.not. allocated(Cr3D_Meshed_Node_in_Ele_Local(i_C)%row))then
            allocate(Cr3D_Meshed_Node_in_Ele_Local(i_C)%row(Max_N_Node_3D(i_C),3))
            Cr3D_Meshed_Node_in_Ele_Local(i_C)%row(1:Max_N_Node_3D(i_C),1:3) = ZR
        endif      
        if(.not. allocated(Crack3D_Meshed_Ele_Attri(i_C)%row))then
            allocate(Crack3D_Meshed_Ele_Attri(i_C)%row(Max_N_Node_3D(i_C),5))
            Crack3D_Meshed_Ele_Attri(i_C)%row(1:Max_N_Node_3D(i_C),1:5) =  ZR
        endif    
        
        if(.not. allocated(Crack3D_Meshed_Ele_Nor_Vector(i_C)%row))then
            allocate(Crack3D_Meshed_Ele_Nor_Vector(i_C)%row(Max_N_Node_3D(i_C),3))
            Crack3D_Meshed_Ele_Nor_Vector(i_C)%row(1:Max_N_Node_3D(i_C),1:3) = ZR
        endif   
        if(.not. allocated(Crack3D_Meshed_Node_Nor_Vector(i_C)%row))then
            allocate(Crack3D_Meshed_Node_Nor_Vector(i_C)%row(Max_N_Node_3D(i_C),3))
            Crack3D_Meshed_Node_Nor_Vector(i_C)%row(1:Max_N_Node_3D(i_C),1:3) = ZR
        endif  
        if(.not. allocated(Crack3D_Meshed_Outline(i_C)%row))then
            allocate(Crack3D_Meshed_Outline(i_C)%row(Max_N_Node_3D(i_C),4))
            Crack3D_Meshed_Outline(i_C)%row(1:Max_N_Node_3D(i_C),1:4) = 0
        endif  

        if(.not. allocated(Crack3D_Meshed_Vertex_x_Vector(i_C)%row))then
            allocate(Crack3D_Meshed_Vertex_x_Vector(i_C)%row(Max_N_Node_3D(i_C),3))
            Crack3D_Meshed_Vertex_x_Vector(i_C)%row(1:Max_N_Node_3D(i_C),1:3) = ZR
        endif  
        if(.not. allocated(Crack3D_Meshed_Vertex_y_Vector(i_C)%row))then
            allocate(Crack3D_Meshed_Vertex_y_Vector(i_C)%row(Max_N_Node_3D(i_C),3))
            Crack3D_Meshed_Vertex_y_Vector(i_C)%row(1:Max_N_Node_3D(i_C),1:3) = ZR
        endif   
        if(.not. allocated(Crack3D_Meshed_Vertex_z_Vector(i_C)%row))then
            allocate(Crack3D_Meshed_Vertex_z_Vector(i_C)%row(Max_N_Node_3D(i_C),3))
            Crack3D_Meshed_Vertex_z_Vector(i_C)%row(1:Max_N_Node_3D(i_C),1:3) = ZR
        endif   
        if(.not. allocated(Crack3D_Meshed_Vertex_T_Matrx(i_C)%row))then
            allocate(Crack3D_Meshed_Vertex_T_Matrx(i_C)%row(Max_N_Node_3D(i_C),3,3))
            Crack3D_Meshed_Vertex_T_Matrx(i_C)%row(1:Max_N_Node_3D(i_C),1:3,1:3) = ZR
        endif  
        if(.not. allocated(Crack3D_Meshed_Vertex_T_Theta(i_C)%row))then
            allocate(Crack3D_Meshed_Vertex_T_Theta(i_C)%row(Max_N_Node_3D(i_C),3))
            Crack3D_Meshed_Vertex_T_Theta(i_C)%row(1:Max_N_Node_3D(i_C),1:3) = ZR
        endif  
        if(.not. allocated(Crack3D_Vector_S1(i_C)%row))then
            allocate(Crack3D_Vector_S1(i_C)%row(Max_N_Node_3D(i_C),3))
            Crack3D_Vector_S1(i_C)%row(1:Max_N_Node_3D(i_C),1:3) = ZR
        endif
        
    endif
    !||||||||||||||||||||||||||||||||||||||||||||||||
    !                                              |
    ! If it is an expanded array. IMPROV2022110502.|
    !                                              |
    !||||||||||||||||||||||||||||||||||||||||||||||||
    if(Key_Extend==1) then
#ifndef Silverfrost
        c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
        c_Max_N_Node_3D_New = int(dble(c_Max_N_Node_3D)*Scale_Factor)
        
        if(c_Max_N_Node_3D_New > Max_N_Node_3D(i_C)) then
            Max_N_Node_3D(i_C) = c_Max_N_Node_3D_New
        endif
        
        allocate(tem_vector_2(c_Max_N_Node_3D_New,3))  
        tem_vector_2(1:c_Max_N_Node_3D,1:3) = Crack3D_Meshed_Node(i_C)%row
        deallocate(Crack3D_Meshed_Node(i_C)%row)
        ! Note: Crack3D_Meshed_Node(i_C)%row has been expanded, tem_vector_2 is automatically deleted
        call move_alloc(tem_vector_2,Crack3D_Meshed_Node(i_C)%row)    
    
        allocate(tem_Int_vector_2(c_Max_N_Node_3D_New,3))  
        tem_Int_vector_2(1:c_Max_N_Node_3D,1:3) = Crack3D_Meshed_Ele(i_C)%row
        deallocate(Crack3D_Meshed_Ele(i_C)%row)
        call move_alloc(tem_Int_vector_2,Crack3D_Meshed_Ele(i_C)%row)    
   
        allocate(tem_vector_2(c_Max_N_Node_3D_New,3))  
        tem_vector_2(1:c_Max_N_Node_3D,1:3) = Crack3D_Meshed_Node_Value(i_C)%row
        deallocate(Crack3D_Meshed_Node_Value(i_C)%row)
        call move_alloc(tem_vector_2,Crack3D_Meshed_Node_Value(i_C)%row)    
             
        allocate(tem_Int_vector_1(c_Max_N_Node_3D_New))  
        tem_Int_vector_1(1:c_Max_N_Node_3D) = Cr3D_Meshed_Node_in_Ele_Num(i_C)%row
        deallocate(Cr3D_Meshed_Node_in_Ele_Num(i_C)%row)
        call move_alloc(tem_Int_vector_1,Cr3D_Meshed_Node_in_Ele_Num(i_C)%row) 

        allocate(tem_vector_2(c_Max_N_Node_3D_New,3))  
        tem_vector_2(1:c_Max_N_Node_3D,1:3) = Cr3D_Meshed_Node_in_Ele_Local(i_C)%row
        deallocate(Cr3D_Meshed_Node_in_Ele_Local(i_C)%row)
        call move_alloc(tem_vector_2,Cr3D_Meshed_Node_in_Ele_Local(i_C)%row)           
         
        allocate(tem_vector_2(c_Max_N_Node_3D_New,5))  
        tem_vector_2(1:c_Max_N_Node_3D,1:5) = Crack3D_Meshed_Ele_Attri(i_C)%row
        deallocate(Crack3D_Meshed_Ele_Attri(i_C)%row)
        call move_alloc(tem_vector_2,Crack3D_Meshed_Ele_Attri(i_C)%row)     
         
        allocate(tem_vector_2(c_Max_N_Node_3D_New,3))  
        tem_vector_2(1:c_Max_N_Node_3D,1:3) = Crack3D_Meshed_Ele_Nor_Vector(i_C)%row
        deallocate(Crack3D_Meshed_Ele_Nor_Vector(i_C)%row)
        call move_alloc(tem_vector_2,Crack3D_Meshed_Ele_Nor_Vector(i_C)%row)          
         
        allocate(tem_vector_2(c_Max_N_Node_3D_New,3))  
        tem_vector_2(1:c_Max_N_Node_3D,1:3) = Crack3D_Meshed_Node_Nor_Vector(i_C)%row
        deallocate(Crack3D_Meshed_Node_Nor_Vector(i_C)%row)
        call move_alloc(tem_vector_2,Crack3D_Meshed_Node_Nor_Vector(i_C)%row)    
        
        !BUGFIX2022112301. 2022-11-23.
        allocate(tem_Int_vector_2(c_Max_N_Node_3D_New,4))  
        tem_Int_vector_2(1:c_Max_N_Node_3D,1:4) = Crack3D_Meshed_Outline(i_C)%row
        deallocate(Crack3D_Meshed_Outline(i_C)%row)
        call move_alloc(tem_Int_vector_2,Crack3D_Meshed_Outline(i_C)%row)    
        
        allocate(tem_vector_2(c_Max_N_Node_3D_New,3))  
        tem_vector_2(1:c_Max_N_Node_3D,1:3) = Crack3D_Meshed_Vertex_x_Vector(i_C)%row
        deallocate(Crack3D_Meshed_Vertex_x_Vector(i_C)%row)
        call move_alloc(tem_vector_2,Crack3D_Meshed_Vertex_x_Vector(i_C)%row)  
        
        allocate(tem_vector_2(c_Max_N_Node_3D_New,3))  
        tem_vector_2(1:c_Max_N_Node_3D,1:3) = Crack3D_Meshed_Vertex_y_Vector(i_C)%row
        deallocate(Crack3D_Meshed_Vertex_y_Vector(i_C)%row)
        call move_alloc(tem_vector_2,Crack3D_Meshed_Vertex_y_Vector(i_C)%row)  
        
        allocate(tem_vector_2(c_Max_N_Node_3D_New,3))  
        tem_vector_2(1:c_Max_N_Node_3D,1:3) = Crack3D_Meshed_Vertex_z_Vector(i_C)%row
        deallocate(Crack3D_Meshed_Vertex_z_Vector(i_C)%row)
        call move_alloc(tem_vector_2,Crack3D_Meshed_Vertex_z_Vector(i_C)%row)  
         
        allocate(tem_vector_3(c_Max_N_Node_3D_New,3,3))  
        tem_vector_3(1:c_Max_N_Node_3D,1:3,1:3) = Crack3D_Meshed_Vertex_T_Matrx(i_C)%row
        deallocate(Crack3D_Meshed_Vertex_T_Matrx(i_C)%row)
        call move_alloc(tem_vector_3,Crack3D_Meshed_Vertex_T_Matrx(i_C)%row)  
        
        allocate(tem_vector_2(c_Max_N_Node_3D_New,3))  
        tem_vector_2(1:c_Max_N_Node_3D,1:3) = Crack3D_Meshed_Vertex_T_Theta(i_C)%row
        deallocate(Crack3D_Meshed_Vertex_T_Theta(i_C)%row)
        call move_alloc(tem_vector_2,Crack3D_Meshed_Vertex_T_Theta(i_C)%row)  

        allocate(tem_vector_2(c_Max_N_Node_3D_New,3))  
        tem_vector_2(1:c_Max_N_Node_3D,1:3) = Crack3D_Vector_S1(i_C)%row
        deallocate(Crack3D_Vector_S1(i_C)%row)
        call move_alloc(tem_vector_2,Crack3D_Vector_S1(i_C)%row)  
#endif

#ifdef Silverfrost  
        print *, '    Error :: Key_Extend=1 not valid for Silverfrost compiler!'
        print *, '             move_alloc function in Silverfrost is faulty, code is commented out!'
        print *, '             In D3_Allocate_Crack_Memory.F90.'
        call Warning_Message('S',Keywords_Blank)
#endif
    endif    
    
!///////////////////////////////////////////////
!                                             /
! Value_Type=2, related to fluid elements.    /
!                                             /
!///////////////////////////////////////////////
case(2)
    !|||||||||||||||||||||||||||||||||
    !                          |
    ! If it is a new variable.      |
    !                          |
    !|||||||||||||||||||||||||||||||||
    if(Key_Extend==0) then
        ! Fluid element. IMPROV2022100401.
        if(.not. allocated(Cracks_FluidEle_CalP_3D(i_C)%row))then
            !IMPROV2023081303.
            allocate(Cracks_FluidEle_CalP_3D(i_C)%row(Max_N_FluEl_3D(i_C),Max_ele_num_CalP)) 
            Cracks_FluidEle_CalP_3D(i_C)%row(1:Max_N_FluEl_3D(i_C),1:Max_ele_num_CalP) = 0 
        endif
        if(.not. allocated(Cracks_FluidEle_Glo_CalP_3D(i_C)%row))then
            !IMPROV2023081303.
            allocate(Cracks_FluidEle_Glo_CalP_3D(i_C)%row(Max_N_FluEl_3D(i_C),Max_ele_num_CalP)) 
            Cracks_FluidEle_Glo_CalP_3D(i_C)%row(1:Max_N_FluEl_3D(i_C),1:Max_ele_num_CalP) = 0  
        endif
        if(.not. allocated(Cracks_FluidEle_Aper_3D(i_C)%row))then
            allocate(Cracks_FluidEle_Aper_3D(i_C)%row(Max_N_FluEl_3D(i_C)))
            Cracks_FluidEle_Aper_3D(i_C)%row(1:Max_N_FluEl_3D(i_C)) = ZR
        endif
        if(.not. allocated(Cracks_FluidEle_num_CalP_3D(i_C)%row))then
            allocate(Cracks_FluidEle_num_CalP_3D(i_C)%row(Max_N_FluEl_3D(i_C)))
            Cracks_FluidEle_num_CalP_3D(i_C)%row(1:Max_N_FluEl_3D(i_C)) = 0
        endif
        
        if(.not. allocated(Cracks_FluidEle_Area_3D(i_C)%row))then
            allocate(Cracks_FluidEle_Area_3D(i_C)%row(Max_N_FluEl_3D(i_C)))
            Cracks_FluidEle_Area_3D(i_C)%row(1:Max_N_FluEl_3D(i_C)) = ZR
        endif
        if(.not. allocated(Cracks_FluidEle_EleNum_3D(i_C)%row))then
            allocate(Cracks_FluidEle_EleNum_3D(i_C)%row(Max_N_FluEl_3D(i_C)))
            Cracks_FluidEle_EleNum_3D(i_C)%row(1:Max_N_FluEl_3D(i_C)) = 0
        endif
        if(.not. allocated(Cracks_FluidEle_Centroid_3D(i_C)%row))then
            allocate(Cracks_FluidEle_Centroid_3D(i_C)%row(Max_N_FluEl_3D(i_C),3))
            Cracks_FluidEle_Centroid_3D(i_C)%row(1:Max_N_FluEl_3D(i_C),1:3) = ZR
        endif
        if(.not. allocated(Cracks_FluidEle_Vector_3D(i_C)%row))then
            allocate(Cracks_FluidEle_Vector_3D(i_C)%row(Max_N_FluEl_3D(i_C),3))
            Cracks_FluidEle_Vector_3D(i_C)%row(1:Max_N_FluEl_3D(i_C),1:3) = ZR
        endif
        if(.not. allocated(Cracks_FluidEle_LCS_x_3D(i_C)%row))then
            allocate(Cracks_FluidEle_LCS_x_3D(i_C)%row(Max_N_FluEl_3D(i_C),3))
            Cracks_FluidEle_LCS_x_3D(i_C)%row(1:Max_N_FluEl_3D(i_C),1:3) = ZR
        endif
        if(.not. allocated(Cracks_FluidEle_LCS_y_3D(i_C)%row))then
            allocate(Cracks_FluidEle_LCS_y_3D(i_C)%row(Max_N_FluEl_3D(i_C),3))
            Cracks_FluidEle_LCS_y_3D(i_C)%row(1:Max_N_FluEl_3D(i_C),1:3) = ZR
        endif
        if(.not. allocated(Cracks_FluidEle_LCS_z_3D(i_C)%row))then
            allocate(Cracks_FluidEle_LCS_z_3D(i_C)%row(Max_N_FluEl_3D(i_C),3))
            Cracks_FluidEle_LCS_z_3D(i_C)%row(1:Max_N_FluEl_3D(i_C),1:3) = ZR
        endif
        if(.not. allocated(Cracks_FluidEle_LCS_T_3D(i_C)%row))then
            allocate(Cracks_FluidEle_LCS_T_3D(i_C)%row(Max_N_FluEl_3D(i_C),3,3))
            Cracks_FluidEle_LCS_T_3D(i_C)%row(1:Max_N_FluEl_3D(i_C),1:3,1:3) = ZR
        endif    
    endif
    !||||||||||||||||||||||||||||||||||||||||||||||||
    !                                              |
    ! If it is an expanded array. IMPROV2022110503.|
    !                                              |
    !||||||||||||||||||||||||||||||||||||||||||||||||
    if(Key_Extend==1) then
#ifndef Silverfrost
        c_Max_N_FluEl_3D = size(Cracks_FluidEle_CalP_3D(i_C)%row,1)
        c_Max_N_FluEl_3D_New = int(dble(c_Max_N_FluEl_3D)*Scale_Factor)
        
        ! Update Max_N_FluEl_3D.
        if(c_Max_N_FluEl_3D_New > Max_N_FluEl_3D(i_C)) then
            Max_N_FluEl_3D(i_C) = c_Max_N_FluEl_3D_New
        endif

        !IMPROV2023081303.
        allocate(tem_Int_vector_2(c_Max_N_FluEl_3D_New,Max_ele_num_CalP))  
        tem_Int_vector_2(1:c_Max_N_FluEl_3D,1:Max_ele_num_CalP) = Cracks_FluidEle_CalP_3D(i_C)%row
        deallocate(Cracks_FluidEle_CalP_3D(i_C)%row)
        call move_alloc(tem_Int_vector_2,Cracks_FluidEle_CalP_3D(i_C)%row)   

        allocate(tem_Int_vector_2(c_Max_N_FluEl_3D_New,Max_ele_num_CalP))  
        tem_Int_vector_2(1:c_Max_N_FluEl_3D,1:Max_ele_num_CalP) = Cracks_FluidEle_Glo_CalP_3D(i_C)%row
        deallocate(Cracks_FluidEle_Glo_CalP_3D(i_C)%row)
        call move_alloc(tem_Int_vector_2,Cracks_FluidEle_Glo_CalP_3D(i_C)%row)           
        
        !
        allocate(tem_vector_1(c_Max_N_FluEl_3D_New))  
        tem_vector_1(1:c_Max_N_FluEl_3D) = Cracks_FluidEle_Aper_3D(i_C)%row
        deallocate(Cracks_FluidEle_Aper_3D(i_C)%row)
        call move_alloc(tem_vector_1,Cracks_FluidEle_Aper_3D(i_C)%row) 
        
        allocate(tem_Int_vector_1(c_Max_N_FluEl_3D_New))  
        tem_Int_vector_1(1:c_Max_N_FluEl_3D) = Cracks_FluidEle_num_CalP_3D(i_C)%row
        deallocate(Cracks_FluidEle_num_CalP_3D(i_C)%row)
        call move_alloc(tem_Int_vector_1,Cracks_FluidEle_num_CalP_3D(i_C)%row) 
        
        allocate(tem_vector_1(c_Max_N_FluEl_3D_New))  
        tem_vector_1(1:c_Max_N_FluEl_3D) = Cracks_FluidEle_Area_3D(i_C)%row
        deallocate(Cracks_FluidEle_Area_3D(i_C)%row)
        call move_alloc(tem_vector_1,Cracks_FluidEle_Area_3D(i_C)%row) 
        
        allocate(tem_Int_vector_1(c_Max_N_FluEl_3D_New))  
        tem_Int_vector_1(1:c_Max_N_FluEl_3D) = Cracks_FluidEle_EleNum_3D(i_C)%row
        deallocate(Cracks_FluidEle_EleNum_3D(i_C)%row)
        call move_alloc(tem_Int_vector_1,Cracks_FluidEle_EleNum_3D(i_C)%row) 
        
        allocate(tem_vector_2(c_Max_N_FluEl_3D_New,3))  
        tem_vector_2(1:c_Max_N_FluEl_3D,1:3) = Cracks_FluidEle_Centroid_3D(i_C)%row
        deallocate(Cracks_FluidEle_Centroid_3D(i_C)%row)
        call move_alloc(tem_vector_2,Cracks_FluidEle_Centroid_3D(i_C)%row)  
        
        allocate(tem_vector_2(c_Max_N_FluEl_3D_New,3))  
        tem_vector_2(1:c_Max_N_FluEl_3D,1:3) = Cracks_FluidEle_Vector_3D(i_C)%row
        deallocate(Cracks_FluidEle_Vector_3D(i_C)%row)
        call move_alloc(tem_vector_2,Cracks_FluidEle_Vector_3D(i_C)%row) 
        
        allocate(tem_vector_2(c_Max_N_FluEl_3D_New,3))  
        tem_vector_2(1:c_Max_N_FluEl_3D,1:3) = Cracks_FluidEle_LCS_x_3D(i_C)%row
        deallocate(Cracks_FluidEle_LCS_x_3D(i_C)%row)
        call move_alloc(tem_vector_2,Cracks_FluidEle_LCS_x_3D(i_C)%row) 
        
        allocate(tem_vector_2(c_Max_N_FluEl_3D_New,3))  
        tem_vector_2(1:c_Max_N_FluEl_3D,1:3) = Cracks_FluidEle_LCS_y_3D(i_C)%row
        deallocate(Cracks_FluidEle_LCS_y_3D(i_C)%row)
        call move_alloc(tem_vector_2,Cracks_FluidEle_LCS_y_3D(i_C)%row) 
        
        allocate(tem_vector_2(c_Max_N_FluEl_3D_New,3))  
        tem_vector_2(1:c_Max_N_FluEl_3D,1:3) = Cracks_FluidEle_LCS_z_3D(i_C)%row
        deallocate(Cracks_FluidEle_LCS_z_3D(i_C)%row)
        call move_alloc(tem_vector_2,Cracks_FluidEle_LCS_z_3D(i_C)%row) 
        
        allocate(tem_vector_3(c_Max_N_FluEl_3D_New,3,3))  
        tem_vector_3(1:c_Max_N_FluEl_3D,1:3,1:3) = Cracks_FluidEle_LCS_T_3D(i_C)%row
        deallocate(Cracks_FluidEle_LCS_T_3D(i_C)%row)
        call move_alloc(tem_vector_3,Cracks_FluidEle_LCS_T_3D(i_C)%row)  
#endif

#ifdef Silverfrost  
        print *, '    Error :: Key_Extend=1 not valid for Silverfrost compiler!'
        print *, '             move_alloc function in Silverfrost is faulty, code is commented out!'
        print *, '             In D3_Allocate_Crack_Memory.F90.'
        call Warning_Message('S',Keywords_Blank)
#endif
    endif
    
!///////////////////////////////////////////////
!                                             /
! Value_Type=3, related to fluid nodes.       /
!                                             /
!///////////////////////////////////////////////
case(3)
    !||||||||||||||||||||||||||||
    !                          |
    ! If it is a new variable. |
    !                          |
    !||||||||||||||||||||||||||||
    if(Key_Extend==0) then
        if(.not. allocated(Cracks_CalP_Coors_3D(i_C)%row))then
            allocate(Cracks_CalP_Coors_3D(i_C)%row(Max_N_CalP_3D(i_C),3))
            Cracks_CalP_Coors_3D(i_C)%row(1:Max_N_CalP_3D(i_C),1:3) = ZR
        endif
        if(.not. allocated(Cracks_CalP_Orient_3D(i_C)%row))then
            allocate(Cracks_CalP_Orient_3D(i_C)%row(Max_N_CalP_3D(i_C),3))
            Cracks_CalP_Orient_3D(i_C)%row(1:Max_N_CalP_3D(i_C),1:3) = ZR
        endif
        if(.not. allocated(Cracks_CalP_MeshedEl_3D(i_C)%row))then
            allocate(Cracks_CalP_MeshedEl_3D(i_C)%row(Max_N_CalP_3D(i_C)))
            Cracks_CalP_MeshedEl_3D(i_C)%row(1:Max_N_CalP_3D(i_C)) = 0
        endif
        if(.not. allocated(Cracks_CalP_Elem_3D(i_C)%row))then
            allocate(Cracks_CalP_Elem_3D(i_C)%row(Max_N_CalP_3D(i_C),2))
            Cracks_CalP_Elem_3D(i_C)%row(1:Max_N_CalP_3D(i_C),1:2) = 0
        endif
        if(.not. allocated(Cracks_CalP_Aper_3D(i_C)%row))then
            allocate(Cracks_CalP_Aper_3D(i_C)%row(Max_N_CalP_3D(i_C)))
            Cracks_CalP_Aper_3D(i_C)%row(1:Max_N_CalP_3D(i_C)) = ZR
        endif  
        if(.not. allocated(Cracks_CalP_UpDis_3D(i_C)%row))then
            allocate(Cracks_CalP_UpDis_3D(i_C)%row(Max_N_CalP_3D(i_C),3))
            Cracks_CalP_UpDis_3D(i_C)%row(1:Max_N_CalP_3D(i_C),1:3) = ZR
        endif  
        if(.not. allocated(Cracks_CalP_LowDis_3D(i_C)%row))then
            allocate(Cracks_CalP_LowDis_3D(i_C)%row(Max_N_CalP_3D(i_C),3))
            Cracks_CalP_LowDis_3D(i_C)%row(1:Max_N_CalP_3D(i_C),1:3) = ZR
        endif  
        if(.not. allocated(Cracks_CalP_Pres_3D(i_C)%row))then
            allocate(Cracks_CalP_Pres_3D(i_C)%row(Max_N_CalP_3D(i_C)))
            Cracks_CalP_Pres_3D(i_C)%row(1:Max_N_CalP_3D(i_C)) = ZR
        endif  
        if(.not. allocated(Cracks_CalP_Tractions_3D(i_C)%row))then
            allocate(Cracks_CalP_Tractions_3D(i_C)%row(Max_N_CalP_3D(i_C),3))
            Cracks_CalP_Tractions_3D(i_C)%row(1:Max_N_CalP_3D(i_C),1:3) = ZR
        endif  
        if(.not. allocated(Cracks_CalP_Pgra_3D(i_C)%row))then
            allocate(Cracks_CalP_Pgra_3D(i_C)%row(Max_N_CalP_3D(i_C)))
            Cracks_CalP_Pgra_3D(i_C)%row(1:Max_N_CalP_3D(i_C)) = ZR
        endif
        if(.not. allocated(Cracks_CalP_Velo_3D(i_C)%row))then
            allocate(Cracks_CalP_Velo_3D(i_C)%row(Max_N_CalP_3D(i_C)))
            Cracks_CalP_Velo_3D(i_C)%row(1:Max_N_CalP_3D(i_C)) = ZR
        endif  
        if(.not. allocated(Cracks_CalP_Quan_3D(i_C)%row))then
            allocate(Cracks_CalP_Quan_3D(i_C)%row(Max_N_CalP_3D(i_C)))
            Cracks_CalP_Quan_3D(i_C)%row(1:Max_N_CalP_3D(i_C)) = ZR
        endif  
        if(.not. allocated(Cracks_CalP_Conc_3D(i_C)%row))then
            allocate(Cracks_CalP_Conc_3D(i_C)%row(Max_N_CalP_3D(i_C)))
            Cracks_CalP_Conc_3D(i_C)%row(1:Max_N_CalP_3D(i_C)) = ZR
        endif   
        if(.not. allocated(Cracks_CalP_Remo_Strs_3D(i_C)%row))then
            allocate(Cracks_CalP_Remo_Strs_3D(i_C)%row(Max_N_CalP_3D(i_C)))
            Cracks_CalP_Remo_Strs_3D(i_C)%row(1:Max_N_CalP_3D(i_C)) = ZR
        endif  
    endif
    !||||||||||||||||||||||||||||||||||||||||||||||||
    !                                              |
    ! If it is an expanded array. IMPROV2022110504.|
    !                                              |
    !||||||||||||||||||||||||||||||||||||||||||||||||
    if(Key_Extend==1) then
#ifndef Silverfrost
        c_Max_N_CalP_3D = size(Cracks_CalP_Coors_3D(i_C)%row,1)
        c_Max_N_CalP_3D_New = int(dble(c_Max_N_CalP_3D)*Scale_Factor)
        
        ! Update Max_N_CalP_3D
        if(c_Max_N_CalP_3D_New > Max_N_CalP_3D(i_C)) then
            Max_N_CalP_3D(i_C) = c_Max_N_CalP_3D_New
        endif
        
        allocate(tem_vector_2(c_Max_N_CalP_3D_New,3))  
        tem_vector_2(1:c_Max_N_CalP_3D,1:3) = Cracks_CalP_Coors_3D(i_C)%row
        deallocate(Cracks_CalP_Coors_3D(i_C)%row)
        call move_alloc(tem_vector_2,Cracks_CalP_Coors_3D(i_C)%row)  
        
        allocate(tem_vector_2(c_Max_N_CalP_3D_New,3))  
        tem_vector_2(1:c_Max_N_CalP_3D,1:3) = Cracks_CalP_Orient_3D(i_C)%row
        deallocate(Cracks_CalP_Orient_3D(i_C)%row)
        call move_alloc(tem_vector_2,Cracks_CalP_Orient_3D(i_C)%row) 
        
        allocate(tem_Int_vector_1(c_Max_N_CalP_3D_New))  
        tem_Int_vector_1(1:c_Max_N_CalP_3D) = Cracks_CalP_MeshedEl_3D(i_C)%row
        deallocate(Cracks_CalP_MeshedEl_3D(i_C)%row)
        call move_alloc(tem_Int_vector_1,Cracks_CalP_MeshedEl_3D(i_C)%row) 
        
        allocate(tem_Int_vector_2(c_Max_N_CalP_3D_New,2))  
        tem_Int_vector_2(1:c_Max_N_CalP_3D,1:2) = Cracks_CalP_Elem_3D(i_C)%row
        deallocate(Cracks_CalP_Elem_3D(i_C)%row)
        call move_alloc(tem_Int_vector_2,Cracks_CalP_Elem_3D(i_C)%row)  
          
        allocate(tem_vector_1(c_Max_N_CalP_3D_New))  
        tem_vector_1(1:c_Max_N_CalP_3D) = Cracks_CalP_Aper_3D(i_C)%row
        deallocate(Cracks_CalP_Aper_3D(i_C)%row)
        call move_alloc(tem_vector_1,Cracks_CalP_Aper_3D(i_C)%row) 
         
        allocate(tem_vector_2(c_Max_N_CalP_3D_New,3))  
        tem_vector_2(1:c_Max_N_CalP_3D,1:3) = Cracks_CalP_UpDis_3D(i_C)%row
        deallocate(Cracks_CalP_UpDis_3D(i_C)%row)
        call move_alloc(tem_vector_2,Cracks_CalP_UpDis_3D(i_C)%row)  
         
        allocate(tem_vector_2(c_Max_N_CalP_3D_New,3))  
        tem_vector_2(1:c_Max_N_CalP_3D,1:3) = Cracks_CalP_LowDis_3D(i_C)%row
        deallocate(Cracks_CalP_LowDis_3D(i_C)%row)
        call move_alloc(tem_vector_2,Cracks_CalP_LowDis_3D(i_C)%row)  
        
        allocate(tem_vector_1(c_Max_N_CalP_3D_New))  
        tem_vector_1(1:c_Max_N_CalP_3D) = Cracks_CalP_Pres_3D(i_C)%row
        deallocate(Cracks_CalP_Pres_3D(i_C)%row)
        call move_alloc(tem_vector_1,Cracks_CalP_Pres_3D(i_C)%row) 
        
        allocate(tem_vector_2(c_Max_N_CalP_3D_New,3))  
        tem_vector_2(1:c_Max_N_CalP_3D,1:3) = Cracks_CalP_Tractions_3D(i_C)%row
        deallocate(Cracks_CalP_Tractions_3D(i_C)%row)
        call move_alloc(tem_vector_2,Cracks_CalP_Tractions_3D(i_C)%row) 
        
        allocate(tem_vector_1(c_Max_N_CalP_3D_New))  
        tem_vector_1(1:c_Max_N_CalP_3D) = Cracks_CalP_Pgra_3D(i_C)%row
        deallocate(Cracks_CalP_Pgra_3D(i_C)%row)
        call move_alloc(tem_vector_1,Cracks_CalP_Pgra_3D(i_C)%row) 
         
        allocate(tem_vector_1(c_Max_N_CalP_3D_New))  
        tem_vector_1(1:c_Max_N_CalP_3D) = Cracks_CalP_Velo_3D(i_C)%row
        deallocate(Cracks_CalP_Velo_3D(i_C)%row)
        call move_alloc(tem_vector_1,Cracks_CalP_Velo_3D(i_C)%row) 
         
        allocate(tem_vector_1(c_Max_N_CalP_3D_New))  
        tem_vector_1(1:c_Max_N_CalP_3D) = Cracks_CalP_Quan_3D(i_C)%row
        deallocate(Cracks_CalP_Quan_3D(i_C)%row)
        call move_alloc(tem_vector_1,Cracks_CalP_Quan_3D(i_C)%row) 
           
        allocate(tem_vector_1(c_Max_N_CalP_3D_New))  
        tem_vector_1(1:c_Max_N_CalP_3D) = Cracks_CalP_Conc_3D(i_C)%row
        deallocate(Cracks_CalP_Conc_3D(i_C)%row)
        call move_alloc(tem_vector_1,Cracks_CalP_Conc_3D(i_C)%row) 
         
        allocate(tem_vector_1(c_Max_N_CalP_3D_New))  
        tem_vector_1(1:c_Max_N_CalP_3D) = Cracks_CalP_Remo_Strs_3D(i_C)%row
        deallocate(Cracks_CalP_Remo_Strs_3D(i_C)%row)
        call move_alloc(tem_vector_1,Cracks_CalP_Remo_Strs_3D(i_C)%row) 
#endif
#ifdef Silverfrost  
        print *, '    Error :: Key_Extend=1 not valid for Silverfrost compiler!'
        print *, '             move_alloc function in Silverfrost is faulty, code is commented out!'
        print *, '             In D3_Allocate_Crack_Memory.F90.'
        call Warning_Message('S',Keywords_Blank)
#endif
    endif
!//////////////////////////
!Value_Type=4, 2023-08-13.
!//////////////////////////
case(4)

end select
      
!=====================================================================
! Updated the maximum values of several global variables. 2023-08-13. 
!=====================================================================
Max_Max_N_FluEl_3D = maxval(Max_N_FluEl_3D)  
Max_Max_N_CalP_3D  = maxval(Max_N_CalP_3D) 
Max_Max_N_Node_3D  = maxval(Max_N_Node_3D) 

RETURN
END SUBROUTINE D3_Allocate_Crack_Memory
