!-----------------------------------------------------------
! Brief: Tests whether two closed double-precision intervals overlap.
!
! Parameters:
!   Input:  Range_1(2),Range_2(2) - the two intervals (low,high)
!   Output: Logical_Yes - .true. if the intervals overlap or touch
!
! Notes:   Compares max(low) with min(high); a small supplementary check
!          against Tol_11 handles the near-touching boundary case.
!-----------------------------------------------------------

subroutine Tool_Yes_Two_Ranges_Overlapped_Double(Range_1,Range_2, Logical_Yes)

!     Used to check whether the ranges of two double-precision numbers intersect (overlap). Can be used to check if coordinate ranges overlap. NEWFTU2022043001.
!     Added on 2022-04-30. 
!     Ref: https://stackoverflow.com/questions/3269434/whats-the-most-efficient-way-to-test-if-two-ranges-overlap
!     or
!     \theory_documents\032 Check Whether Two Data Ranges Intersect (Overlap) - 2022-04-30.pdf

use Global_Float_Type

implicit none
real(kind=FT),intent(in)::Range_1(2),Range_2(2)
logical,intent(out)::Logical_Yes
real(kind=FT) max_x1_y1,min_x2_y2

Logical_Yes =.False.

max_x1_y1 = max(Range_1(1),Range_2(1))
min_x2_y2 = min(Range_1(2),Range_2(2))

if(max_x1_y1 <= min_x2_y2)then
    Logical_Yes =.True.
endif

if(abs(max_x1_y1-min_x2_y2) <=Tol_11) then
    Logical_Yes =.True.
endif

return 
end SUBROUTINE Tool_Yes_Two_Ranges_Overlapped_Double          
