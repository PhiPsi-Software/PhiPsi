!-----------------------------------------------------------
! Brief: Build lists of new and to-be-updated 3D XFEM elements.
!
! Parameters:
!   Input: isub - current load substep index
!
! Notes:   Allocates and fills the global flags Elem_New_XFEM_Flag
!          and Elem_Update_XFEM_Flag by comparing the current
!          XFEM element state against the saved previous state.
!          No-op for the first substep.
!-----------------------------------------------------------

subroutine D3_Get_New_and_Updated_XFEM_Elements(isub)
!     After confirming the enhanced nodes, obtain the list of newly added 3D XFEM elements and the list of 3D XFEM elements whose stiffness matrices need to be updated.
!     2022-06-24. NEWFTU2022062403.


!     !---------------FEM Element List and XFEM Enriched Element List------------
!     integer num_FEM_Elem, num_XFEM_Elem ! Number of FEM elements and number of XFEM elements
!     integer, allocatable :: FEM_Elem_List(:)   ! FEM element list
!     integer, allocatable :: XFEM_Elem_List(:)  ! XFEM enhanced element list
!     integer, allocatable :: Elem_XFEM_Flag(:)  ! Used to indicate whether the element is an enhanced element
!     integer, allocatable :: Elem_Location(:,:) ! Used to store the position of the element number in the List in reverse, the first column is XFEM, the second column is FEM
!     integer, allocatable:: Elem_New_XFEM_Flag(:)     ! Used to mark whether an element is a newly added enhanced element. 2022-06-24.
!     integer, allocatable :: Elem_Update_XFEM_Flag(:)  ! Used to indicate whether an element is an enriched element that requires updating the element stiffness matrix. 2022-06-24.


!     ----------------------------
!     Read public variable module
!     ----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_3D
use Global_XFEM_Elements

!     --------------------------
!     Variable Type Declaration
!     --------------------------
implicit none
integer,intent(in)::isub
integer :: i_E

!     ---------------------------------
!     The first step is not necessary.
!     ---------------------------------
if (isub<=1)then
    return
endif

!     ------------------------
!     Variable Initialization
!     ------------------------
if (allocated(Elem_New_XFEM_Flag))then
    deallocate(Elem_New_XFEM_Flag)
endif
allocate(Elem_New_XFEM_Flag(num_Elem))
Elem_New_XFEM_Flag(1:num_Elem) = 0

if (allocated(Elem_Update_XFEM_Flag))then
    deallocate(Elem_Update_XFEM_Flag)
endif
allocate(Elem_Update_XFEM_Flag(num_Elem))
Elem_Update_XFEM_Flag(1:num_Elem) = 0

!     ----------
!     ext steps
!     ----------
! The Determine_Enriched_Nodes_3D subroutine saved Elem_Location_old and Elem_XFEM_Flag_Old.
print *,'    Getting New and to-be-Updated XFEM Elements...'

do i_E = 1,Num_Elem
    ! if the previous step was not an XFEM element and this step is an XFEM element, then it is marked
    ! as a new XFEM element.
    if(Elem_XFEM_Flag_Old(i_E)==0 .and.Elem_XFEM_Flag(i_E)==1)then
        Elem_New_XFEM_Flag(i_E) = 1
    endif
enddo



return
end subroutine D3_Get_New_and_Updated_XFEM_Elements
