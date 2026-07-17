!-----------------------------------------------------------
! Brief: Set up sub-elements for case-2 (4-edge intersection) integration.
!
! Parameters:
!   Input:  c_E              - element index
!           c_C              - crack index
!           num_Inters       - number of edge intersections
!           ele_InterP       - intersection point coordinates
!   In/Out: ele_InterP_Edge  - edge numbers of intersections
!
! Notes:   Splits the element into 2 hex sub-elements for the
!          specific edge configuration 5-6-7-8. Marked abandoned.
!-----------------------------------------------------------

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
