!-----------------------------------------------------------
! Brief: Save the end-of-step hydraulic-fracturing state to memory.
!
! Parameters:
!   Input:  ifra           - current fracture-step index
!           Counter_Iter   - iteration counter for the step
!           Old_total_Time - previous accumulated time per fracture
!           total_Time     - new accumulated time
!           Yes_Growth     - per-tip growth flag from the previous step
!
! Notes:   On the first call allocates the L_*_CalP snapshot arrays and
!          copies the current hydraulic-fracturing state into them for
!-----------------------------------------------------------

subroutine Write_Lfrac_Vars_to_Memory_HF(ifra,Counter_Iter, Old_total_Time,total_Time,Yes_Growth)
! Store the relevant data at the end moment of the previous fracture step into memory

!     ****************************
!     Read public variable module
!     ****************************
use Global_Float_Type
use Global_Common
use Global_Filename
use Global_Crack
use Global_HF

!     *********************
!     Variable Declaration
!     *********************
implicit none     
integer,intent(in)::ifra,Counter_Iter
real(kind=FT),intent(in)::Old_total_Time(Num_Frac),total_Time
logical,intent(in)::Yes_Growth(Max_Num_Cr,2)
logical :: c_Yes_Last_Growth

print *, '    Writting variables of last fracture step ' // 'to memory....'

!     ************************************************************************
!     If it is the end of the first rupture step, then allocate memory space.
!     ************************************************************************
if(ifra==1)then
    allocate(L_Cracks_CalP_Num(Max_Num_Cr))
    allocate(L_Cracks_CalP_Pres(Max_Num_Cr,Max_Num_Cr_CalP))
    allocate(L_Cracks_CalP_Aper(Max_Num_Cr,Max_Num_Cr_CalP))
    allocate(L_Cracks_CalP_Coors(Max_Num_Cr,Max_Num_Cr_CalP,2))
    allocate(L_Cracks_CalP_Orient(Max_Num_Cr,Max_Num_Cr_CalP))
    allocate(L_Cracks_CalP_Seg(Max_Num_Cr,Max_Num_Cr_CalP))
    allocate(L_Cracks_CalP_Elem(Max_Num_Cr,Max_Num_Cr_CalP))
    allocate(L_Cracks_CalP_Pgra(Max_Num_Cr,Max_Num_Cr_CalP))
    allocate(L_Cracks_CalP_Velo(Max_Num_Cr,Max_Num_Cr_CalP))
    allocate(L_Cracks_CalP_Quan(Max_Num_Cr,Max_Num_Cr_CalP))  
    allocate(L_Cracks_CalP_Conc(Max_Num_Cr,Max_Num_Cr_CalP))
    allocate(Map_L_Cracks_CalP_Conc(Max_Num_Cr,Max_Num_Cr_CalP))
    ! Initialized to 0
    L_Cracks_CalP_Num(1:Max_Num_Cr) = 0
    L_Cracks_CalP_Pres(1:Max_Num_Cr,1:Max_Num_Cr_CalP) = ZR
    L_Cracks_CalP_Aper(1:Max_Num_Cr,1:Max_Num_Cr_CalP) = ZR
    L_Cracks_CalP_Coors(1:Max_Num_Cr,1:Max_Num_Cr_CalP,1:2)=.0D0
    L_Cracks_CalP_Orient(1:Max_Num_Cr,1:Max_Num_Cr_CalP)=ZR
    L_Cracks_CalP_Seg(1:Max_Num_Cr,1:Max_Num_Cr_CalP)   = 0
    L_Cracks_CalP_Elem(1:Max_Num_Cr,1:Max_Num_Cr_CalP)  = 0
    L_Cracks_CalP_Pgra(1:Max_Num_Cr,1:Max_Num_Cr_CalP)  =ZR
    L_Cracks_CalP_Velo(1:Max_Num_Cr,1:Max_Num_Cr_CalP)  =ZR
    L_Cracks_CalP_Quan(1:Max_Num_Cr,1:Max_Num_Cr_CalP)  =ZR
    L_Cracks_CalP_Conc(1:Max_Num_Cr,1:Max_Num_Cr_CalP)  =ZR
    Map_L_Cracks_CalP_Conc(1:Max_Num_Cr,1:Max_Num_Cr_CalP)=ZR
endif
!     **********
!     Save Data
!     **********
if(ifra < Num_Frac)then
    L_Cracks_CalP_Num    = Cracks_CalP_Num
    L_Cracks_CalP_Pres   = Cracks_CalP_Pres
    L_Cracks_CalP_Aper   = Cracks_CalP_Aper
    L_Cracks_CalP_Coors  = Cracks_CalP_Coors
    L_Cracks_CalP_Orient = Cracks_CalP_Orient
    L_Cracks_CalP_Seg    = Cracks_CalP_Seg
    L_Cracks_CalP_Elem   = Cracks_CalP_Elem
    L_Cracks_CalP_Pgra   = Cracks_CalP_Pgra
    L_Cracks_CalP_Velo   = Cracks_CalP_Velo
    L_Cracks_CalP_Quan   = Cracks_CalP_Quan
    L_Cracks_CalP_Conc   = Cracks_CalP_Conc
endif
!     ***************************************************************
!     If it is the last breaking step, clear the relevant variables.
!     ***************************************************************
if (any(Yes_Growth).eqv..True.)then
    c_Yes_Last_Growth = .True.
    ! If no cracks have propagated, exit the program.
else
    c_Yes_Last_Growth = .False.
endif

if(ifra==Num_Frac .or. (c_Yes_Last_Growth.eqv..False..and.ifra<Num_Frac)) then
    deallocate(L_Cracks_CalP_Num)
    deallocate(L_Cracks_CalP_Pres)  
    deallocate(L_Cracks_CalP_Aper)
    deallocate(L_Cracks_CalP_Coors)
    deallocate(L_Cracks_CalP_Orient)
    deallocate(L_Cracks_CalP_Seg)  
    deallocate(L_Cracks_CalP_Elem) 
    deallocate(L_Cracks_CalP_Pgra) 
    deallocate(L_Cracks_CalP_Velo)
    deallocate(L_Cracks_CalP_Quan)  
    deallocate(L_Cracks_CalP_Conc)
    deallocate(Map_L_Cracks_CalP_Conc)
endif

return
END subroutine Write_Lfrac_Vars_to_Memory_HF
