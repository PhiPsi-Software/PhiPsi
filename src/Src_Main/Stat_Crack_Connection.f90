!-----------------------------------------------------------
! Brief: Determine the connectivity between discrete cracks.
!
! Parameters:
!   Input:  iter - current iteration counter
!
! Notes:   Populates the global Cracks_Cone_* arrays describing which
!          cracks connect at tips or middle junctions.
!-----------------------------------------------------------

subroutine Stat_Crack_Connection(iter)
!     Determine the connectivity of the cracks.
!
!     Variable Description (see Module_Global for details):
!     integer Cracks_Cone_Num(Max_Num_Cr)                        !Number of cracks directly connected to each crack (including tips and middle)
!     integer Cracks_Cone_Cr(Max_Num_Cr, Max_Num_Cone_Cr) ! Crack numbers directly connected to each fracture (including tips and middle parts)
!     integer Cracks_Cone_NumTipCr(Max_Num_Cr)                   !Each crack has how many crack tips connected to other cracks, either 1 or 2
!     integer Cracks_Cone_TipCrNum(Max_Num_Cr,2)                 !The connection relationship between the two crack tips of each crack and other cracks:
!     Cracks_Cone_TipType(i_C,1) -- The connected crack number of crack tip 1 of i_C
!     Cracks_Cone_TipType(i_C,2) -- The connected crack number at crack tip 2 of i_C
!     integer Cracks_Cone_TipJuEle(Max_Num_Cr,2)                 !The element number of the crack intersection point directly connecting the two tips of each crack
!     real(kind=FT) Cracks_Cone_TipJuCor(Max_Num_Cr,2,2)      !The coordinates of the crack intersection where the two tips of each crack are directly connected, which are the crack tip coordinates
!     integer Cracks_Cone_NumMidCr(Max_Num_Cr)                   !Number of cracks connected to the junction at the middle of each crack
!     integer Cracks_Cone_MidCrNum(Max_Num_Cr, Max_Num_Cone_Cr)   ! Number of cracks connected to the junction in the middle of each crack
!     integer Cracks_Cone_MidCrTip(Max_Num_Cr, Max_Num_Cone_Cr)   !The crack tip number (1 or 2) corresponding to the crack connected to the junction in the middle of each crack
!     integer Cracks_Cone_MidJuEle(Max_Num_Cr, Max_Num_Cone_Cr) ! The element number of the Junction point in the middle of each crack
!     real(kind=FT) Cracks_Cone_MidJuCor(Max_Num_Cr, Max_Num_Cone_Cr, 2) ! Coordinates of the junction point in the middle of each crack

!     ****************************
!     Read public variable module
!     ****************************
use Global_Float_Type
use Global_Crack
use Global_Crack_Common
use Global_Model
use Global_Common

!     *********************
!     Variable Declaration
!     *********************
implicit none
integer,intent(in)::iter
integer i_C,Jun_Cr,Num_Cone_Cr
integer :: Cr_num_Pr
if(Temp_Show_Message) then
    print *,"    Counting cracks' connections......"
endif

!     ************************
!     Variable Initialization
!     ************************
Cracks_Cone_Num(1:Num_Crack) =0
Cracks_Cone_Cr(1:Num_Crack,1:Max_Num_Cone_Cr)       = 0
!-----------------------------------------------------
Cracks_Cone_NumTipCr(1:Max_Num_Cr)               = 0  
Cracks_Cone_TipCrNum(1:Max_Num_Cr,1:2)           = 0
Cracks_Cone_TipJuEle(1:Max_Num_Cr,1:2)           = 0
Cracks_Cone_TipJuCor(Max_Num_Cr,1:2,1:2)         = ZR
!-----------------------------------------------------
Cracks_Cone_NumMidCr(1:Max_Num_Cr)               = 0
Cracks_Cone_MidCrNum(1:Max_Num_Cr,1:Max_Num_Cone_Cr) = 0
Cracks_Cone_MidCrTip(1:Max_Num_Cr,1:Max_Num_Cone_Cr) = 0
Cracks_Cone_MidJuEle(1:Max_Num_Cr,1:Max_Num_Cone_Cr) = 0
Cracks_Cone_MidJuCor(1:Max_Num_Cr,1:Max_Num_Cone_Cr,1:2)=ZR

!     *********************************************
!     Crack cycle
!     Algorithm: Count the junctions of each crack
!     *********************************************
do i_C=1,num_Crack
    !-----------------------------------------------------------------------------------
    ! Check whether the number of crack associations exceeds the program's scale range.
    ! If it exceeds the limit, terminate the program.
    !-----------------------------------------------------------------------------------
    if (int(maxval(Cracks_Cone_Num(1:Num_Crack))*0.9) > Max_Num_Cone_Cr) then
        write(*,1001)  
        print *, '    Error :: Fracture network is too complicated!'
        call Warning_Message('S',Keywords_Blank)           
    end if
    !------------------------------
    ! If Tip 1 is the Junction tip
    !------------------------------
    if (Crack_Tip_Type(i_C,1) ==1) then
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Get the crack numbers associated with the current crack
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Jun_Cr = Crack_Jun_CrNum(i_C,1)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! From the crack (connected at the crack tip)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Cracks_Cone_Num(i_C) = Cracks_Cone_Num(i_C) + 1
        Cracks_Cone_Cr(i_C,Cracks_Cone_Num(i_C)) =Jun_Cr
        Cracks_Cone_NumTipCr(i_C)  = Cracks_Cone_NumTipCr(i_C)+1
        Cracks_Cone_TipCrNum(i_C,1)= Jun_Cr
        ! element number and coordinates of the junction point
        Cracks_Cone_TipJuEle(i_C,1)   = Crack_Jun_Elem(i_C,1)
        Cracks_Cone_TipJuCor(i_C,1,1:2)=Crack_Tip_Coor(i_C,1,1:2)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Main Crack (Central Connection)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Cracks_Cone_Num(Jun_Cr)=Cracks_Cone_Num(Jun_Cr)+1
        Cracks_Cone_Cr(Jun_Cr,Cracks_Cone_Num(Jun_Cr))=i_C
        !Cracks_Cone_JunType(Jun_Cr,Cracks_Cone_Num(Jun_Cr))=3
        Cracks_Cone_NumMidCr(Jun_Cr) =  Cracks_Cone_NumMidCr(Jun_Cr) +1
        Cracks_Cone_MidCrNum(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr)) = i_C
        Cracks_Cone_MidCrTip(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr)) = 1
        ! element number and coordinates of the junction point
        Cracks_Cone_MidJuEle(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr))  = Crack_Jun_Elem(i_C,1)
        Cracks_Cone_MidJuCor(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr),1:2) = Crack_Tip_Coor(i_C,1,:)
    end if
    !------------------------------
    ! If Tip 2 is the Junction tip
    !------------------------------
    if (Crack_Tip_Type(i_C,2) ==1) then
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Get the crack numbers associated with the current crack
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Jun_Cr = Crack_Jun_CrNum(i_C,2)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! From the crack (connected at the crack tip)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Cracks_Cone_Num(i_C) = Cracks_Cone_Num(i_C) + 1
        Cracks_Cone_Cr(i_C,Cracks_Cone_Num(i_C)) =Jun_Cr
        Cracks_Cone_NumTipCr(i_C)  = Cracks_Cone_NumTipCr(i_C)+1
        Cracks_Cone_TipCrNum(i_C,2)   = Jun_Cr
        ! element number and coordinates of the junction point
        Cracks_Cone_TipJuEle(i_C,2)   = Crack_Jun_Elem(i_C,2)
        Cracks_Cone_TipJuCor(i_C,2,1:2)=Crack_Tip_Coor(i_C,2,1:2)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Main Crack (Central Connection)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Cracks_Cone_Num(Jun_Cr)=Cracks_Cone_Num(Jun_Cr)+1
        Cracks_Cone_Cr(Jun_Cr,Cracks_Cone_Num(Jun_Cr))=i_C
        !Cracks_Cone_JunType(Jun_Cr,Cracks_Cone_Num(Jun_Cr))=3
        Cracks_Cone_NumMidCr(Jun_Cr) = Cracks_Cone_NumMidCr(Jun_Cr) +1
        Cracks_Cone_MidCrNum(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr)) = i_C
        Cracks_Cone_MidCrTip(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr)) = 2
        ! element number and coordinates of the junction point
        Cracks_Cone_MidJuEle(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr)) = Crack_Jun_Elem(i_C,2)
        Cracks_Cone_MidJuCor(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr),1:2) = Crack_Tip_Coor(i_C,2,:)
    end if

end do
1001 format('     !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!') 

!     *************************************
!     Check if it exceeds the program size
!     *************************************
1013 format(5X,'Error :: crack ',I4,' has too many connected cracks!')
1014 format(5X,'         The avaliable maximum number is ',I4,'!')
! Each crack can connect to at most Max_Num_Cone_Cr cracks.
do i_C =1,num_Crack
    if(Cracks_Cone_NumMidCr(i_C)  > Max_Num_Cone_Cr)then
        write(*,1013) i_C
        write(*,1014) Max_Num_Cone_Cr
        call Warning_Message('S',Keywords_Blank) 
    endif
enddo

return
END subroutine Stat_Crack_Connection
