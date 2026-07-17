!-----------------------------------------------------------
! Brief: Sanity-check and patch the hydraulic-fracturing input
!        (forbid water in friction-type natural cracks, enforce a
!        single water-driven initial crack, etc.).
!
! Parameters:
!   Input:  (none - operates on global crack/HF state)
!   Output: (none - modifies global crack state and emits warnings)
!
! Notes:   Prints a 'Logic checking' header and aborts on fatal
!          conditions via Warning_Message.
!-----------------------------------------------------------

subroutine HF_Logic_Checking
! Some logical checks and data corrections for hydraulic fracturing analysis, the currently
! implemented inspection items include:
! (1) At the initial moment, it is obvious that friction cracks cannot contain water (only for the
! third type of natural cracks);
! (2) The current version of the program can only have one initial water-driven fracture.

!*****************************
! Read public variable module
!*****************************
use Global_Float_Type
use Global_Common
use Global_Filename
use Global_Crack
use Global_Crack_Common
use Global_HF
use Global_Model

implicit none
integer i_C,c_Na_Cr_num,Count_Fluid_Driven_Cr

print *, "    Logic checking of the HF analysis...." 

!*****************************************************
! Check item 1: Friction cracks cannot contain water.
!*****************************************************
1001 format (5X,'-- Frictioanl crack ',I4,' should not be fluid-driven!')
if(num_Na_Crack>=1 .and. Key_Na_Crack_Type==3)then
    do i_C =1,num_Crack
        ! Corresponding natural fracture number
        c_Na_Cr_num = Cracks_fric_NF_num(i_C)
        ! If the current calculated fracture number corresponds to a natural friction fracture, check
        ! whether it contains water. If it contains water, it needs to be adjusted to be water-free.
        if (c_Na_Cr_num /= 0) then
            if(Cracks_HF_State(i_C)==1)then
                Cracks_HF_State(i_C)=0
                write(*,1001) i_C
            endif
        endif
    enddo
endif

!***************************************************************************************
! Check iTerm 2: The current version of the program can only have one initial hydraulic
! fracture.
!***************************************************************************************
Count_Fluid_Driven_Cr = 0
do i_C =1,num_Crack
    if(Cracks_HF_State(i_C)==1)then
        Count_Fluid_Driven_Cr = Count_Fluid_Driven_Cr +1
    endif
enddo
if(Count_Fluid_Driven_Cr >1)then
    print *,'    -- Error:: more than one fluid-driven' // ' fractures found, which is illegal!'
    call Warning_Message('S',Keywords_Blank)
endif

!***************************************************************************************
! Check iTerm 2: The current version of the program can only have one initial hydraulic
! fracture.
!***************************************************************************************
Count_Fluid_Driven_Cr = 0
do i_C =1,num_Crack
    if(Cracks_HF_State(i_C)==1)then
        Count_Fluid_Driven_Cr = Count_Fluid_Driven_Cr +1
    endif
enddo
if(Count_Fluid_Driven_Cr >1)then
    print *,'    -- Error:: more than one fluid-driven' // ' fractures found, which is illegal!'
    call Warning_Message('S',Keywords_Blank)
endif
return
END subroutine HF_Logic_Checking
