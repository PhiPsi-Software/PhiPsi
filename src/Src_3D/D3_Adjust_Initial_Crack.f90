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
 
SUBROUTINE D3_Adjust_Initial_Crack
! Adjust the coordinates of the initial 3D cracks.
! Objective: To prevent situations where the edges of the discrete mesh triangles on the crack
! surface
! and the edges of the elements coincide exactly, making it impossible to detect fluid nodes.
! Firstly written on 2021-08-19.



!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
      
      
      
!---------------------------     
! Variable Type Declaration
!---------------------------     
implicit none
integer i_C,i_Rand
real(kind=FT) c_rand_number(12)


      
!--------------------
! Each fracture loop
!--------------------
do i_C=1,num_Crack  
    ! If the crack has been adjusted before, proceed to the next cycle. IMPROV2022053102.
    if(Cracks_Initial_Adjusted(i_C)==1)then
        cycle
    endif
    
    !Generate random number 1 or -1   
    do i_Rand = 1,12
        call Tool_Generate_Random_1_or_Minus_1(c_rand_number(i_Rand))
    enddo
    
    !##################################
    ! Case1: Rectangular Crack Surface
    !##################################
    if(sum(abs(Crack3D_Coor(i_C,1:4,1:3)))>Tol_20)then
        Crack3D_Coor(i_C,1,1) = Crack3D_Coor(i_C,1,1) + c_rand_number(1)*Ave_Elem_L*Tol_2
        Crack3D_Coor(i_C,1,2) = Crack3D_Coor(i_C,1,2) + c_rand_number(2)*Ave_Elem_L*Tol_2
        Crack3D_Coor(i_C,1,3) = Crack3D_Coor(i_C,1,3) + c_rand_number(3)*Ave_Elem_L*Tol_2
        Crack3D_Coor(i_C,2,1) = Crack3D_Coor(i_C,2,1) + c_rand_number(4)*Ave_Elem_L*Tol_2
        Crack3D_Coor(i_C,2,2) = Crack3D_Coor(i_C,2,2) + c_rand_number(5)*Ave_Elem_L*Tol_2
        Crack3D_Coor(i_C,2,3) = Crack3D_Coor(i_C,2,3) + c_rand_number(6)*Ave_Elem_L*Tol_2
        Crack3D_Coor(i_C,3,1) = Crack3D_Coor(i_C,3,1) + c_rand_number(7)*Ave_Elem_L*Tol_2
        Crack3D_Coor(i_C,3,2) = Crack3D_Coor(i_C,3,2) + c_rand_number(8)*Ave_Elem_L*Tol_2
        Crack3D_Coor(i_C,3,3) = Crack3D_Coor(i_C,3,3) + c_rand_number(9)*Ave_Elem_L*Tol_2
        Crack3D_Coor(i_C,4,1) = Crack3D_Coor(i_C,4,1) + c_rand_number(10)*Ave_Elem_L*Tol_2
        Crack3D_Coor(i_C,4,2) = Crack3D_Coor(i_C,4,2) + c_rand_number(11)*Ave_Elem_L*Tol_2
        Crack3D_Coor(i_C,4,3) = Crack3D_Coor(i_C,4,3) + c_rand_number(12)*Ave_Elem_L*Tol_2      
    !################################
    ! Case 2: Circular crack surface
    !################################
    elseif(sum(abs(Crack3D_Cir_Coor(i_C,1:7)))>Tol_20)then     
        Crack3D_Cir_Coor(i_C,1)=Crack3D_Cir_Coor(i_C,1) + c_rand_number(1)*Ave_Elem_L*Tol_2  
        Crack3D_Cir_Coor(i_C,2)=Crack3D_Cir_Coor(i_C,2) + c_rand_number(2)*Ave_Elem_L*Tol_2
        Crack3D_Cir_Coor(i_C,3)=Crack3D_Cir_Coor(i_C,3) + c_rand_number(3)*Ave_Elem_L*Tol_2
    !#####################################
    ! Case 3: Elliptical Fracture Surface
    !#####################################
    elseif(sum(abs(Crack3D_Ellip_Coor(i_C,1:8)))>Tol_20)then 
        Crack3D_Ellip_Coor(i_C,1) = Crack3D_Ellip_Coor(i_C,1)+ c_rand_number(1)*Ave_Elem_L*Tol_2 
        Crack3D_Ellip_Coor(i_C,2) = Crack3D_Ellip_Coor(i_C,2)+ c_rand_number(2)*Ave_Elem_L*Tol_2 
        Crack3D_Ellip_Coor(i_C,3) = Crack3D_Ellip_Coor(i_C,3)+ c_rand_number(3)*Ave_Elem_L*Tol_2 
    endif
enddo


!-----------------------------------------------------------------------------------------
! Check whether the crack has been previously initialized and adjusted. IMPROV2022053102.
!-----------------------------------------------------------------------------------------
do i_C =1,num_Crack
    ! If the crack has not been previously initialized and adjusted, it is marked as initialized and
    ! adjusted.
    if(Cracks_Initial_Adjusted(i_C)==0)then
        Cracks_Initial_Adjusted(i_C)=1
    endif  
enddo   

RETURN
END SUBROUTINE D3_Adjust_Initial_Crack
