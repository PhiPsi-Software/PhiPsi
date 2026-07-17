!-----------------------------------------------------------
! Brief: Generate the initial natural-fracture population (if the
!        user requested random generation) and process the data
!        according to the natural-fracture type (friction or
!        cementation).
!
! Parameters:
!   Input:  (none - reads global keyword and geometric state)
!   Output: (none - fills Na_Crack_Coor, num_Na_Crack, etc.)
!
! Notes:   Respects the user-defined rupture-zone bounds and
!          calls Tool_Generate_Natural_Fractures / Save_HF_Natural_
!          Cracks_Coors as needed.
!-----------------------------------------------------------

subroutine Natural_Fractures_Related_Preparation
! Preparations related to natural cracks:
! Generate initial natural fractures (if needed), and process relevant data according to the type of
! natural fractures.

!     ----------------------------
!     Read public variable module
!     ----------------------------
use Global_Float_Type
use Global_Common
use Global_Filename
use Global_Crack
use Global_Crack_Common
use Global_HF
use Global_Model

implicit none
integer i_Na_Cr,c_num_point
real(kind=FT) Ger_Min_X_Coor,Ger_Min_Y_Coor
real(kind=FT) Ger_Max_X_Coor,Ger_Max_Y_Coor

!     -----------------------------------------
!     Execute parameter initialization process
!     -----------------------------------------
print *, "    Natural fractures related operations...." 

! If the rupture zone is defined, initial inclusions and cracks will only be generated within the
! rupture zone.
if(Key_Fracture_Zone==1)then
    Ger_Min_X_Coor = Frac_Zone_MinX
    Ger_Min_Y_Coor = Frac_Zone_MinY
    Ger_Max_X_Coor = Frac_Zone_MaxX
    Ger_Max_Y_Coor = Frac_Zone_MaxY
elseif (Key_Fracture_Zone==0)then
    Ger_Min_X_Coor = Min_X_Coor
    Ger_Min_Y_Coor = Min_Y_Coor
    Ger_Max_X_Coor = Max_X_Coor
    Ger_Max_Y_Coor = Max_Y_Coor
endif

!     -------------------------------------------
!     If initial natural fractures are generated
!     -------------------------------------------
if (Key_Random_NaCr==1) then
    !*******************************************************************************************
    ! Types of natural fractures, 1: Shear fractures (extend along both sides, forming a more
    ! complex fracture network);
    ! 2: Cementation cracks (expand along one side, forming a relatively simple crack network);
    !*******************************************************************************************
    if(Key_Na_Crack_Type==1 .or. Key_Na_Crack_Type==2)then
        ! The total number of natural fractures is the number of fractures generated randomly.
        num_Na_Crack = num_Rand_Na_Crack
        call Tool_Generate_Natural_Fractures( Ger_Min_X_Coor,Ger_Max_X_Coor, Ger_Min_Y_Coor,Ger_Max_Y_Coor, &
        NaCr_Length,NaCr_Orientation, &
        NaCr_Len_Delta,NaCr_Ori_Delta, &
        num_Na_Crack,NaCr_Length*Random_NaCr_Rad_Factor, &
        Na_Crack_Coor(1:num_Na_Crack,1:2,1:2))
        ! Preserve natural cracks
        call Save_HF_Natural_Cracks_Coors(num_Na_Crack, Na_Crack_Coor(1:num_Na_Crack,1:2,1:2))
        ! Total number of natural fractures
        Each_Na_Cr_Poi_Num(1:num_Na_Crack) =2
        !********************************************************************************************
        ! Types of natural fractures, 3: True friction type, where all fractures are considered from
        ! the initial moment, taking into account the friction of all fractures.
        !********************************************************************************************
    elseif(Key_Na_Crack_Type==3)then
        ! The total number of natural fractures is the number of fractures generated randomly.
        num_Na_Crack = num_Rand_Na_Crack
        ! If a fracture zone is defined, it is obviously necessary to generate natural fractures within the
        ! fracture zone.
        if(Key_Fracture_Zone==1)then
            call Tool_Generate_Natural_Fractures( Frac_Zone_MinX,Frac_Zone_MaxX, Frac_Zone_MinY,Frac_Zone_MaxY, &
            NaCr_Length,NaCr_Orientation, &
            NaCr_Len_Delta,NaCr_Ori_Delta, &
            num_Na_Crack,NaCr_Length*Random_NaCr_Rad_Factor, &
            Na_Crack_Coor(1:num_Na_Crack,1:2,1:2))
            ! If no fracture zone is defined, natural fractures will be generated throughout the entire model.
        elseif(Key_Fracture_Zone==0)then
            call Tool_Generate_Natural_Fractures( Min_X_Coor,Max_X_Coor, Min_Y_Coor,Max_Y_Coor, &
            NaCr_Length,NaCr_Orientation, &
            NaCr_Len_Delta,NaCr_Ori_Delta, &
            num_Na_Crack,NaCr_Length*Random_NaCr_Rad_Factor, &
            Na_Crack_Coor(1:num_Na_Crack,1:2,1:2))
        endif
        ! Given the friction natural fracture number corresponding to each crack
        ! First number the non-natural fractures, then number the natural fractures.
        Cracks_fric_NF_num(1:num_Crack) = 0
        do i_Na_Cr =1,num_Na_Crack
            Cracks_fric_NF_num(num_Crack+i_Na_Cr) = i_Na_Cr
            Crack_Coor(num_Crack+i_Na_Cr,1,1:2)= Na_Crack_Coor(i_Na_Cr,1,1:2)
            Crack_Coor(num_Crack+i_Na_Cr,2,1:2)= Na_Crack_Coor(i_Na_Cr,2,1:2)

        enddo
        Each_Cr_Poi_Num(num_Crack+1:num_Crack+num_Na_Crack) =2
        num_Crack = num_Crack + num_Na_Crack

    endif 
    !     ------------------------------------------------------------------------------------------
    !     If natural fractures are not generated randomly (only the third type of natural fractures
    !     needs to
    !     be processed)
    !     ------------------------------------------------------------------------------------------
elseif (Key_Random_NaCr==0) then
    ! Then, if it is the third type of natural fracture, all natural fractures participate in the
    ! calculation from the initial moment.
    if(Key_Na_Crack_Type==3 .and. num_Na_Crack >=1)then   
        ! Given the friction natural fracture number corresponding to each crack
        ! First, assign the number 0 to non-natural fractures, and then number the natural fractures.
        Cracks_fric_NF_num(1:num_Crack) = 0
        do i_Na_Cr =1,num_Na_Crack
            c_num_point = Each_Na_Cr_Poi_Num(i_Na_Cr)
            Cracks_fric_NF_num(num_Crack+i_Na_Cr) = i_Na_Cr
            Crack_Coor(num_Crack+i_Na_Cr,1,1:c_num_point)= Na_Crack_Coor(i_Na_Cr,1,1:c_num_point)
            Crack_Coor(num_Crack+i_Na_Cr,2,1:c_num_point)= Na_Crack_Coor(i_Na_Cr,2,1:c_num_point)
            Each_Cr_Poi_Num(num_Crack+i_Na_Cr) = c_num_point
        enddo
        num_Crack = num_Crack + num_Na_Crack
    endif
endif

return
END subroutine Natural_Fractures_Related_Preparation
