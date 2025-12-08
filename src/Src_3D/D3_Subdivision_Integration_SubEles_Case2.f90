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
 
SUBROUTINE D3_Subdivision_Integration_SubEles_Case2(c_E,c_C,num_Inters,ele_InterP_Edge,ele_InterP)
! Block-wise integral related data calculation. Case 2.
! Ref: My PhiPsi Development Notebook V1-P54-55. NEWFTU2022072701.
! 2022-07-27.
! Abandoned, 2022-07-27.

!----------------------------------
! Read the public variable module.
!----------------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_3D
      
      
!----------------------------
! Variable type declaration.
!----------------------------
implicit none
integer,intent(in):: c_E,c_C,num_Inters
integer,intent(inout):: ele_InterP_Edge(num_Inters)
real(kind=FT),intent(in):: ele_InterP(num_Inters,3)

integer Index_Vetor_4(4)
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
integer c_NN(8)  

!-----------------------------------------------
! Relevant information for the current element.
!-----------------------------------------------
c_NN(1:8)    = G_NN(1:8,c_E)
c_X_NODES = G_X_NODES(1:8,c_E)
c_Y_NODES = G_Y_NODES(1:8,c_E)  
c_Z_NODES = G_Z_NODES(1:8,c_E) 

!-----------------------------
! Sort the four edge numbers.
!-----------------------------
call Vector_Sort_Int_with_Index(num_Inters,ele_InterP_Edge(1:num_Inters),Index_Vetor_4(1:num_Inters))   

!--------------------------------------------------------------------------------
! Edge CASE 5, 6, 7, 8. At this point, it is split into two hexahedral elements.
!--------------------------------------------------------------------------------
if (ele_InterP_Edge(1)==5 .and. ele_InterP_Edge(2)==6 .and. ele_InterP_Edge(3)==7 .and. ele_InterP_Edge(4)==8) then
    ! Modify the basic data related to Gauss integration.
    num_SubEles = num_SubEles + 1
    Elems_Integration_Type(c_E,c_C) = 2
    Elems_Num_SubEles(c_E,c_C)      = 2
    Elems_Type_SubEles(c_E,c_C)     = 1
    Elems_SubEles_Index(c_E,c_C)    =  num_SubEles
    SubEles_Integ_Num(num_SubEles)  =  2*Num_Gau_P_SubInteg_6
endif

RETURN
END SUBROUTINE D3_Subdivision_Integration_SubEles_Case2
