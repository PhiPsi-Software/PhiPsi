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
 
SUBROUTINE D3_Offset_Initial_Crack_to_Regular_n_Vector(i_C,c_Crack_Type)
! Adjust the normal direction of the fracture surface to the nearest regular direction (smallest
! angle).
! 2022-08-02.
!
!
! c_Crack_Type = 1  ! Rectangular crack surface or polygonal crack surface
! c_Crack_Type = 2  !Circular crack surface
! c_Crack_Type = 3  !Elliptical crack surface (To be done)

!#############################
! Read public variable module
!#############################
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Crack_Common
use Global_Crack_3D
use Global_Filename

!###########################
! Variable Type Declaration
!###########################
implicit none
integer,intent(in)::i_C,c_Crack_Type
integer num_Cr_Poi
real(kind=FT) CR_Center_Old(3)
real(kind=FT) c_Elem_Centroid(3)
real(kind=FT) Center_Offset_Vector(3)
integer OUT_Elem,i_Cr_Poi
integer num_Angle_Division_Theta,num_Angle_Division_Phi
integer i_Theta,i_Phi
real(kind=FT),allocatable::Tem_n_vectors(:,:)
real(kind=FT),allocatable::Uniqued_Matrix(:,:)
real(kind=FT),allocatable::Possible_n_vectors(:,:)
integer num_Possible_n_vectors,Tem_num_n_vectors
real(kind=FT) c_Theta,c_Phi,c_x,c_y,c_z
integer c_vector_Count
integer Uniqued_m
integer,allocatable::Uni_Mat_Count(:)
real(kind=FT) Nor_Vector_Old(3)
real(kind=FT) angle_1,angle_2,angle_3
integer i
integer i_Posi
real(kind=FT),allocatable::Tem_Angles(:)
real(kind=FT) Expected_n_Vector(3),Expected_Angle
integer c_min_Location
real(kind=FT) k_vector(3),cross_pro_old_new(3)
real(kind=FT) c_Norm_2_cross_pro_old_new
real(kind=FT) c_Point(3),c_Point_New(3),tem_Cross_Pro_Vector(3)
real(kind=FT) tem_Dot_Product,c_V(3)
integer i_Point
real(kind=FT) n_Vector_Old(3)


!#################################################
!#                                              #
!#  Obtain possible normal vectors: OPTION 2    #
!#  Directly given according to the resolution  #
!#  Ref: My Notebook V1-P61~P62.                #
!#                                              #
!#################################################
! 90 degrees (6 pieces)
if (Adjust_Cracks_Resolution==6) then
     num_Possible_n_vectors = 6
     allocate(Possible_n_vectors(6,3))
     Possible_n_vectors(1,1:3) = [  ONE,   ZR,    ZR]
     Possible_n_vectors(2,1:3) = [ -ONE,   ZR,    ZR]
     Possible_n_vectors(3,1:3) = [  ZR,   ONE,    ZR]
     Possible_n_vectors(4,1:3) = [  ZR,  -ONE,    ZR]
     Possible_n_vectors(5,1:3) = [  ZR,    ZR,   ONE]
     Possible_n_vectors(6,1:3) = [  ZR,    ZR,  -ONE]
! 90 degrees (6 pieces) 45-degree cuts (12 pieces)
elseif(Adjust_Cracks_Resolution==18) then
     num_Possible_n_vectors = 18
     allocate(Possible_n_vectors(18,3))
     Possible_n_vectors(1,1:3) = [  ONE,   ZR,    ZR]
     Possible_n_vectors(2,1:3) = [ -ONE,   ZR,    ZR]
     Possible_n_vectors(3,1:3) = [  ZR,   ONE,    ZR]
     Possible_n_vectors(4,1:3) = [  ZR,  -ONE,    ZR]
     Possible_n_vectors(5,1:3) = [  ZR,    ZR,   ONE]
     Possible_n_vectors(6,1:3) = [  ZR,    ZR,  -ONE]
     
     Possible_n_vectors(7,1:3) = [  sqrt(TWO)/TWO,    sqrt(TWO)/TWO,               ZR]     
     Possible_n_vectors(8,1:3) = [  sqrt(TWO)/TWO,   -sqrt(TWO)/TWO,               ZR]  
     Possible_n_vectors(9,1:3) = [ -sqrt(TWO)/TWO,    sqrt(TWO)/TWO,               ZR]  
     Possible_n_vectors(10,1:3)= [ -sqrt(TWO)/TWO,   -sqrt(TWO)/TWO,               ZR]  
     Possible_n_vectors(11,1:3)= [  sqrt(TWO)/TWO,               ZR,    sqrt(TWO)/TWO]  
     Possible_n_vectors(12,1:3)= [  sqrt(TWO)/TWO,               ZR,   -sqrt(TWO)/TWO]  
     Possible_n_vectors(13,1:3)= [ -sqrt(TWO)/TWO,               ZR,    sqrt(TWO)/TWO]  
     Possible_n_vectors(14,1:3)= [ -sqrt(TWO)/TWO,               ZR,   -sqrt(TWO)/TWO]  
     Possible_n_vectors(15,1:3)= [             ZR,    sqrt(TWO)/TWO,    sqrt(TWO)/TWO]  
     Possible_n_vectors(16,1:3)= [             ZR,    sqrt(TWO)/TWO,   -sqrt(TWO)/TWO]  
     Possible_n_vectors(17,1:3)= [             ZR,   -sqrt(TWO)/TWO,    sqrt(TWO)/TWO]  
     Possible_n_vectors(18,1:3)= [             ZR,   -sqrt(TWO)/TWO,   -sqrt(TWO)/TWO]  
! 90 degrees (6 pieces) 45-degree cuts (12 pieces) Regular octahedron (8 pieces)
elseif(Adjust_Cracks_Resolution==26) then
     num_Possible_n_vectors = 26
     allocate(Possible_n_vectors(26,3))
     Possible_n_vectors(1,1:3) = [  ONE,   ZR,    ZR]
     Possible_n_vectors(2,1:3) = [ -ONE,   ZR,    ZR]
     Possible_n_vectors(3,1:3) = [  ZR,   ONE,    ZR]
     Possible_n_vectors(4,1:3) = [  ZR,  -ONE,    ZR]
     Possible_n_vectors(5,1:3) = [  ZR,    ZR,   ONE]
     Possible_n_vectors(6,1:3) = [  ZR,    ZR,  -ONE]
     
     Possible_n_vectors(7,1:3) = [  sqrt(TWO)/TWO,    sqrt(TWO)/TWO,               ZR]     
     Possible_n_vectors(8,1:3) = [  sqrt(TWO)/TWO,   -sqrt(TWO)/TWO,               ZR]  
     Possible_n_vectors(9,1:3) = [ -sqrt(TWO)/TWO,    sqrt(TWO)/TWO,               ZR]  
     Possible_n_vectors(10,1:3)= [ -sqrt(TWO)/TWO,   -sqrt(TWO)/TWO,               ZR]  
     Possible_n_vectors(11,1:3)= [  sqrt(TWO)/TWO,               ZR,    sqrt(TWO)/TWO]  
     Possible_n_vectors(12,1:3)= [  sqrt(TWO)/TWO,               ZR,   -sqrt(TWO)/TWO]  
     Possible_n_vectors(13,1:3)= [ -sqrt(TWO)/TWO,               ZR,    sqrt(TWO)/TWO]  
     Possible_n_vectors(14,1:3)= [ -sqrt(TWO)/TWO,               ZR,   -sqrt(TWO)/TWO]  
     Possible_n_vectors(15,1:3)= [             ZR,    sqrt(TWO)/TWO,    sqrt(TWO)/TWO]  
     Possible_n_vectors(16,1:3)= [             ZR,    sqrt(TWO)/TWO,   -sqrt(TWO)/TWO]  
     Possible_n_vectors(17,1:3)= [             ZR,   -sqrt(TWO)/TWO,    sqrt(TWO)/TWO]  
     Possible_n_vectors(18,1:3)= [             ZR,   -sqrt(TWO)/TWO,   -sqrt(TWO)/TWO] 
     
     Possible_n_vectors(19,1:3)= [  ONE/sqrt(THR),    ONE/sqrt(THR),    ONE/sqrt(THR)]
     Possible_n_vectors(20,1:3)= [  ONE/sqrt(THR),    ONE/sqrt(THR),   -ONE/sqrt(THR)]
     Possible_n_vectors(21,1:3)= [  ONE/sqrt(THR),   -ONE/sqrt(THR),    ONE/sqrt(THR)]
     Possible_n_vectors(22,1:3)= [  ONE/sqrt(THR),   -ONE/sqrt(THR),   -ONE/sqrt(THR)]
     Possible_n_vectors(23,1:3)= [ -ONE/sqrt(THR),    ONE/sqrt(THR),    ONE/sqrt(THR)]
     Possible_n_vectors(24,1:3)= [ -ONE/sqrt(THR),    ONE/sqrt(THR),   -ONE/sqrt(THR)]
     Possible_n_vectors(25,1:3)= [ -ONE/sqrt(THR),   -ONE/sqrt(THR),    ONE/sqrt(THR)]
     Possible_n_vectors(26,1:3)= [ -ONE/sqrt(THR),   -ONE/sqrt(THR),   -ONE/sqrt(THR)]
endif


!############################################
! If it is a rectangular or polygonal crack.
!############################################
if (c_Crack_Type == 1) then
    !///////////
    ! Centroid.
    !///////////
    num_Cr_Poi = Each_Cr_Poi_Num(i_C)
    CR_Center_Old(1) =sum(Crack3D_Coor(i_C,1:num_Cr_Poi,1))/num_Cr_Poi
    CR_Center_Old(2) =sum(Crack3D_Coor(i_C,1:num_Cr_Poi,2))/num_Cr_Poi
    CR_Center_Old(3) =sum(Crack3D_Coor(i_C,1:num_Cr_Poi,3))/num_Cr_Poi
    
    !//////////////////////////////////////////////////////////////////
    ! Determine the normal vector of the original fracture surface. 
    ! Take the normal vector of the first discrete fracture element as 
    ! the normal vector of the original fracture.
    !//////////////////////////////////////////////////////////////////
    Nor_Vector_Old(1:3) = Crack3D_Meshed_Ele_Nor_Vector(i_C)%row(1,1:3)
    call Vector_Normalize(3,Nor_Vector_Old)
    
    allocate(Tem_Angles(num_Possible_n_vectors))
    !///////////////////////
    ! Determine each angle.
    !///////////////////////
    do i_Posi=1,num_Possible_n_vectors
        call Tool_Angle_of_Vectors_a_and_b_3D(Nor_Vector_Old,Possible_n_vectors(i_Posi,1:3),angle_1,1)
        !call Tool_Angle_of_Vectors_a_and_b_3D(Nor_Vector_Old,Possible_n_vectors(i_Posi,1:3),angle_2,2)
        !call Tool_Angle_of_Vectors_a_and_b_3D(Nor_Vector_Old,Possible_n_vectors(i_Posi,1:3),angle_3,3)
        Tem_Angles(i_Posi) = angle_1
    enddo
    
    !//////////////////////////////////////////////////////
    ! Find the vector corresponding to the smallest angle.
    !//////////////////////////////////////////////////////
    c_min_Location = MINLOC(Tem_Angles(1:num_Possible_n_vectors),1)
    Expected_n_Vector(1:3) = Possible_n_vectors(c_min_Location,1:3)
    Expected_Angle         = Tem_Angles(c_min_Location)
    
    !///////////////////////////////////////////////////////
    ! If the angle is very small, do not adjust. Just exit.
    !///////////////////////////////////////////////////////
    if(Tem_Angles(c_min_Location) <= Tol_6) then
        return
    endif
    
    !////////////////////////////////////////
    ! Adjust the coordinates of each vertex.
    !////////////////////////////////////////
    !Ref: https://math.stackexchange.com/questions/871867/rotation-matrix-of-triangle-in-3d
    !     or
    ! \theory_documents\038 Adjusting Fracture Point Coordinates Based on 
    ! Changes in Fracture Surface Normal Vectors_2022-08-02.pdf
    ! Note: The rotation algorithm here is performed around an axis passing 
    ! through the origin. For axes that do not pass through the origin
    ! (e.g., an axis passing through the centroid of the fracture surface),
    ! It is necessary to first adjust to the origin, and then adjust back to
    ! the centroid of the crack surface.
    ! Calculate rotation axis: k_vector
    call Vector_Cross_Product_3(Nor_Vector_Old(1:3),Expected_n_Vector(1:3),cross_pro_old_new(1:3))
    call Vector_Norm2(3,cross_pro_old_new(1:3),c_Norm_2_cross_pro_old_new)   
    k_vector(1:3) = cross_pro_old_new(1:3)/c_Norm_2_cross_pro_old_new
    
    ! Cycle through each crack point.
    do i_Point =1,num_Cr_Poi
        c_Point(1:3) = Crack3D_Coor(i_C,i_Point ,1:3)
        ! The following line is very important: because the rotation 
        ! algorithm here rotates around the axis k_vector that passes through the origin.
        ! So first translate to the origin.
        c_Point(1:3) = c_Point(1:3)  - CR_Center_Old(1:3)
        !c_V(1:3)     = c_Point(1:3) - CR_Center_Old(1:3)
        
        call Vector_Cross_Product_3(k_vector(1:3),c_Point(1:3),tem_Cross_Pro_Vector(1:3))
        tem_Dot_Product = DOT_PRODUCT(k_vector(1:3),c_Point(1:3))
        
        !Expected_Angle  = acos(DOT_PRODUCT(Nor_Vector_Old(1:3),Expected_n_Vector(1:3)))
         
        c_Point_New(1:3) = c_Point(1:3)*cos(Expected_Angle) + tem_Cross_Pro_Vector(1:3)*sin(Expected_Angle) +   &
                           k_vector(1:3)*tem_Dot_Product*(ONE-cos(Expected_Angle))
        
        ! The following line is very important: translate back to the centroid of the fracture surface
        c_Point_New(1:3) = c_Point_New(1:3) + CR_Center_Old(1:3)
        
        ! Update crack coordinates.
        Crack3D_Coor(i_C,i_Point ,1:3) = c_Point_New(1:3)
    enddo
endif
      
!#########################################################
! If it is a circular crack: just update the coordinates
! of the center, and leave the other variables unchanged.
!#########################################################
if (c_Crack_Type == 2) then
    !//////////////////////////////////////////
    ! Original crack surface normal direction.
    !//////////////////////////////////////////
    Nor_Vector_Old(1:3)  = Crack3D_Cir_Coor(i_C,4:6)
    call Vector_Normalize(3,Nor_Vector_Old)
    
    allocate(Tem_Angles(num_Possible_n_vectors))
    !///////////////////////
    ! Determine each angle.
    !///////////////////////
    do i_Posi=1,num_Possible_n_vectors
        call Tool_Angle_of_Vectors_a_and_b_3D(Nor_Vector_Old,Possible_n_vectors(i_Posi,1:3),angle_1,1)
        !call Tool_Angle_of_Vectors_a_and_b_3D(Nor_Vector_Old,Possible_n_vectors(i_Posi,1:3),angle_2,2)
        !call Tool_Angle_of_Vectors_a_and_b_3D(Nor_Vector_Old,Possible_n_vectors(i_Posi,1:3),angle_3,3)
        Tem_Angles(i_Posi) = angle_1
    enddo
    
    !//////////////////////////////////////////////////////
    ! Find the vector corresponding to the smallest angle.
    !//////////////////////////////////////////////////////
    c_min_Location = MINLOC(Tem_Angles(1:num_Possible_n_vectors),1)
    Expected_n_Vector(1:3) = Possible_n_vectors(c_min_Location,1:3)
    Expected_Angle         = Tem_Angles(c_min_Location)
    
    !///////////////////////////////////////////////////////
    ! If the angle is very small, do not adjust. Just exit.
    !///////////////////////////////////////////////////////
    if(Tem_Angles(c_min_Location) <= Tol_8) then
        return
    endif
    
    !////////////////////////////////////
    ! Update the crack normal direction.
    !////////////////////////////////////
    Crack3D_Cir_Coor(i_C,4:6) =  Expected_n_Vector(1:3)
    
    
endif


!############################################################
! If it is an elliptical crack, it is not supported for now.
!############################################################
if (c_Crack_Type == 3) then
    print *,'     ERROR(2022080201) :: Elliptical cracks not supported yet!'
    call Warning_Message('S',Keywords_Blank) 
endif
        

      
RETURN
END SUBROUTINE D3_Offset_Initial_Crack_to_Regular_n_Vector
