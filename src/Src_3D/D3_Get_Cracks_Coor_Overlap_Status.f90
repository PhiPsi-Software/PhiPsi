!-----------------------------------------------------------
! Brief: Build a pairwise crack bounding-box overlap matrix.
!
! Parameters:
!   Output: Crack_Coor_Overlap_Status - num_Crack x num_Crack logical
!                                     matrix; true where ranges
!                                     overlap on x, y, and z
!
! Notes:   Reads Crack_Coor_Range populated by D3_Get_Cracks_Coor_Ranges
!          and tests each pair with Tool_Yes_Two_Ranges_Overlapped_Double.
!-----------------------------------------------------------

SUBROUTINE D3_Get_Cracks_Coor_Overlap_Status(Crack_Coor_Overlap_Status)
! Check whether there is an overlap in the coordinate ranges between each crack. NEWFTU2022050402.
! There is an overlap only if x, y, and z all intersect.
!2022-05-04.

!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Crack_Common
use Global_Crack_3D
use omp_lib      
      
!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
logical,intent(out):: Crack_Coor_Overlap_Status(num_Crack,num_Crack)
integer i_C,j_C
logical x_Logical_Yes,y_Logical_Yes,z_Logical_Yes


!------------
! Initialize
!------------
Crack_Coor_Overlap_Status(1:num_Crack,1:num_Crack) = .False.

!-----------------------------------------------
! If there is only one crack, exit immediately.
!-----------------------------------------------
if(num_Crack==1) then
    return
endif

!-------------
! Crack cycle
!-------------
! IMPROV2022072203. OpenMP Parallelization.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,j_C,x_Logical_Yes,y_Logical_Yes,z_Logical_Yes) SCHEDULE(static) 
do i_C =1,num_Crack
    do j_C =1,num_Crack
        if(i_C/=j_C)then
            call Tool_Yes_Two_Ranges_Overlapped_Double(Crack_Coor_Range(i_C,1,1:2),Crack_Coor_Range(j_C,1,1:2),x_Logical_Yes)  
            if (x_Logical_Yes .eqv. .False.) then
                cycle
            endif
            call Tool_Yes_Two_Ranges_Overlapped_Double(Crack_Coor_Range(i_C,2,1:2),Crack_Coor_Range(j_C,2,1:2),y_Logical_Yes)  
            if (y_Logical_Yes .eqv. .False.) then
                cycle
            endif
            call Tool_Yes_Two_Ranges_Overlapped_Double(Crack_Coor_Range(i_C,3,1:2),Crack_Coor_Range(j_C,3,1:2),z_Logical_Yes)  
            if (z_Logical_Yes .eqv. .False.) then
                cycle
            endif  
            Crack_Coor_Overlap_Status(i_C,j_C) = .True.
        endif
    enddo    
enddo
!$omp end parallel do         

RETURN
END SUBROUTINE D3_Get_Cracks_Coor_Overlap_Status
