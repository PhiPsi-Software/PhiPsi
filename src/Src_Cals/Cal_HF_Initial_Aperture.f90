!-----------------------------------------------------------
! Brief: Set the initial aperture of HF calculation points
!        for a new fracture step (zero or inherited).
!
! Parameters:
!   Input:  iter, ifra, Counter - iteration / step counters
!           Last_Cracks_CalP_Num_ifra - previous cal-point counts
!           Last_Cr_CalP_Aper_ifra    - previous apertures
!   Output: Initial_Cr_CalP_Aper     - starting apertures
!
! Notes:   Honors Key_IniPre_PassOn: 0 = zero; 1 = carry over
!          the previous step's aperture to current points.
!-----------------------------------------------------------

subroutine Cal_HF_Initial_Aperture(iter,ifra,Counter, Last_Cracks_CalP_Num_ifra, Last_Cr_CalP_Aper_ifra, &
Initial_Cr_CalP_Aper)
! Option1: If data is not inherited between each rupture step, set the opening to 0 at each rupture
! step.
! Option 2: If data is inherited between each fracture step, then read the data at the end of the
! previous iteration step and assign it to the current fracture step.
! The first iteration step

use Global_Float_Type
use Global_Common
use Global_Crack
use Global_Crack_Common
use Global_HF

implicit none

integer,intent(in)::iter,ifra,Counter
real(kind=FT),intent(in)::Last_Cr_CalP_Aper_ifra(Max_Num_Cr, &
Max_Num_Cr_CalP)
integer,intent(in)::Last_Cracks_CalP_Num_ifra(Max_Num_Cr)
real(kind=FT),intent(out)::Initial_Cr_CalP_Aper(Max_Num_Cr, Max_Num_Cr_CalP)

integer i_C,Num_Div_Points
integer :: Last_Num_Div_Points
do i_C = 1,num_Crack
    if (Cracks_HF_State(i_C) == 1) then
        if (Key_IniPre_PassOn==0)then
            Num_Div_Points = Cracks_CalP_Num(i_C)
            Initial_Cr_CalP_Aper(i_C,1:Num_Div_Points)=ZR
        elseif (Key_IniPre_PassOn==1)then
            if (ifra==1) then
                Num_Div_Points = Cracks_CalP_Num(i_C)
                Initial_Cr_CalP_Aper(i_C,1:Num_Div_Points)=ZR  
            else

                Num_Div_Points = Cracks_CalP_Num(i_C)
                Last_Num_Div_Points=Last_Cracks_CalP_Num_ifra(i_C)
                Initial_Cr_CalP_Aper(i_C,1:Last_Num_Div_Points)= Last_Cr_CalP_Aper_ifra(i_C,1:Last_Num_Div_Points)
                Initial_Cr_CalP_Aper(i_C, Last_Num_Div_Points+1:Num_Div_Points)= ZR
            endif
        end if
    end if
end do


return 
end subroutine Cal_HF_Initial_Aperture             
