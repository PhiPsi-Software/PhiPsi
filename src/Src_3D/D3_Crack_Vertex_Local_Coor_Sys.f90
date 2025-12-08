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
 
subroutine D3_Crack_Vertex_Local_Coor_Sys(iter)
! Calculate the local coordinate system of the crack surface boundary nodes.

! real(kind=FT) Crack3D_Meshed_Vertex_x_Vector(Max_Num_Cr_3D, Max_N_Node_3D,3)! Local x-axis vector
! of 3D crack boundary points after discretization
! real(kind=FT) Crack3D_Meshed_Vertex_y_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3) ! Local y-axis
! vector of 3D crack boundary points after discretization
! real(kind=FT) Crack3D_Meshed_Vertex_z_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3) ! Local z-axis
! vector of 3D crack boundary points after discretization
!
!                   O Crack_Center
!                   | 
!                   |
!      Last O       |  V_temp  
!            \     /|\
!             \     |     O  P2
!              \    |    / 
!               \   |   / 
!          V1   /\  |  /\  V2
!                 \ | /
!                  \|/
!                   O  P1
!                 
!
!


!***************
! Public Module
!***************
use Global_Crack_Common
use Global_Crack_3D
use Global_Elem_Area_Vol
use Global_Model
use Global_Common
use Global_INTERFACE_Tool_ThetaX_ThetaY_ThetaZ_3D_rotation

!**********************
! Variable Declaration
!**********************
implicit none
integer,intent(in)::iter
integer i_C
integer num_Cr_Edges,Edge_P1,Edge_P2,i_V
integer Egde_Last_P1
real(kind=FT) Edge_P1_Point(3),Edge_P2_Point(3),Egde_Last_P1_Point(3)
real(kind=FT) Vector_V1(3),Vector_V2(3)
real(kind=FT) Vector_x(3),Vector_y(3),Vector_z(3)
real(kind=FT) ang_bisector_of_V1_V2(3)
real(kind=FT) centroid_x,centroid_y,centroid_z 
real(kind=FT) Crack_Center(3)
real(kind=FT) Vector_Crack_x(3),Vector_Crack_y(3),Vector_Crack_z(3)
real(kind=FT) ThetaX,ThetaY,ThetaZ,T_Matrx(3,3)
integer Num_Cr_Nodes
real(kind=FT) Vector_temp(3),angle_V1_Vtemp

!****************
! Variables used
!****************   

!***********************************************
! Calculate the local coordinate vectors of the 
! discrete fracture boundary points
!***********************************************
do i_C =1,num_Crack

Num_Cr_Nodes = Crack3D_Meshed_Node_num(i_C) 
centroid_x=sum(Crack3D_Meshed_Node(i_C)%row(1:Num_Cr_Nodes,1))/dble(Num_Cr_Nodes)
centroid_y=sum(Crack3D_Meshed_Node(i_C)%row(1:Num_Cr_Nodes,2))/dble(Num_Cr_Nodes)
centroid_z=sum(Crack3D_Meshed_Node(i_C)%row(1:Num_Cr_Nodes,3))/dble(Num_Cr_Nodes)

Crack_Center = [centroid_x,centroid_y,centroid_z]

num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
! Crack boundary vertex loop.
do i_V =1,num_Cr_Edges
    !##########################################
    ! STEP0: Calculate Vector_V1 and Vector_V2
    !##########################################
    Edge_P1 = Crack3D_Meshed_Outline(i_C)%row(i_V,1)
    Edge_P2 = Crack3D_Meshed_Outline(i_C)%row(i_V,2)
    
    if(i_V>=2) then
        Egde_Last_P1 =Crack3D_Meshed_Outline(i_C)%row(i_V-1,1)
    elseif(i_V==1) then
        Egde_Last_P1 =Crack3D_Meshed_Outline(i_C)%row(num_Cr_Edges,1)
    endif
    
    Edge_P1_Point = Crack3D_Meshed_Node(i_C)%row(Edge_P1,1:3)
    Edge_P2_Point = Crack3D_Meshed_Node(i_C)%row(Edge_P2,1:3)
    Egde_Last_P1_Point=Crack3D_Meshed_Node(i_C)%row(Egde_Last_P1,1:3)
    
    Vector_V1(1:3) = Edge_P2_Point - Edge_P1_Point
    call Vector_Normalize(3,Vector_V1) 
    
    Vector_V2(1:3) = Egde_Last_P1_Point - Edge_P1_Point
    call Vector_Normalize(3,Vector_V2) 

    !###########################
    ! STEP1: Calculate Vector_y
    !###########################
    Vector_y=Crack3D_Meshed_Node_Nor_Vector(i_C)%row(Edge_P1,1:3) 
    call Vector_Normalize(3,Vector_y)

    ! Note: The following is the new algorithm. 2023-08-13. First, 
    ! calculate Vector_z (along the tangent direction of the crack front), 
    ! then use the cross product to calculate Vector_x.
    ! IMPROV2023081302. Core Improvement.
    ! The new algorithm is simpler and more robust.
    !########################################
    ! STEP2: Calculate Vector_z, 2023-08-13.
    !########################################
    Vector_z = ((-Vector_V1) + Vector_V2)/TWO
    call Vector_Normalize(3,Vector_z) 
    
    !############################
    ! STEP 3: Calculate Vector_x
    !############################
    ! Vector_x and Vector_y cross to get Vector_z
    call Vector_Cross_Product_3(Vector_y,Vector_z,Vector_x)  
    call Vector_Normalize(3,Vector_x)     

    !############################
    ! Save to a global variable.
    !############################
    Crack3D_Meshed_Vertex_x_Vector(i_C)%row(i_V,1:3)=Vector_x
    Crack3D_Meshed_Vertex_y_Vector(i_C)%row(i_V,1:3)=Vector_y
    Crack3D_Meshed_Vertex_z_Vector(i_C)%row(i_V,1:3)=Vector_z
    
    !##########################################################
    ! Calculate rotation angle and rotation matrix, 2020-01-04
    !##########################################################    
    ! Determine the rotation matrix T_Matrx corresponding to the transformation matrix of the current
    ! crack coordinate system. For details, see my notes V3-P178.
    Vector_Crack_x(1:3) = [ONE,ZR,ZR]
    Vector_Crack_y(1:3) = [ZR,ONE,ZR]
    Vector_Crack_z(1:3) = [ZR,ZR,ONE]
    call Tool_ThetaX_ThetaY_ThetaZ_3D_rotation(i_V,  &
              Vector_x(1:3),Vector_y(1:3),Vector_z(1:3),  &
              Vector_Crack_x(1:3),Vector_Crack_y(1:3),Vector_Crack_z(1:3),  &
              ThetaX,ThetaY,ThetaZ,T_Matrx(1:3,1:3))
    Crack3D_Meshed_Vertex_T_Matrx(i_C)%row(i_V,1:3,1:3)=T_Matrx(1:3,1:3)
    Crack3D_Meshed_Vertex_T_Theta(i_C)%row(i_V,1:3)  =[ThetaX,ThetaY,ThetaZ]
  end do
enddo

RETURN
END SUBROUTINE D3_Crack_Vertex_Local_Coor_Sys