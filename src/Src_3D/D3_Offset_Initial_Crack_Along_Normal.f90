!-----------------------------------------------------------
! Brief: Shift an initial 3D crack along its surface normal.
!
! Parameters:
!   Input: i_C         - crack index
!          c_Crack_Type - 1: rect/polygon, 2: circle, 3: ellipse
!          delta_L     - offset distance along the normal
!
! Notes:   For circular cracks the center is moved by delta_L*n,
!          leaving other parameters unchanged. Polygon variant is
!          a placeholder. Aborts on elliptical cracks.
!-----------------------------------------------------------

SUBROUTINE D3_Offset_Initial_Crack_Along_Normal(i_C,c_Crack_Type,delta_L)
! Adjust the initial crack along the normal direction.
! 2022-08-01.
!
! c_Crack_Type = 1  ! Rectangular crack surface or polygonal crack surface
! c_Crack_Type = 2  !Circular crack surface
! c_Crack_Type = 3  !Elliptical crack surface

!##############################        
! Read public variable module.
!##############################
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Global_Filename

!############################
! Variable Type Declaration.
!############################
implicit none
integer,intent(in)::i_C,c_Crack_Type
real(kind=FT),intent(in)::delta_L
integer num_Cr_Poi
real(kind=FT) CR_Center_Old(3),CR_Center_New(3)
real(kind=FT) n_Vector(1:3) 


!############################################
! If it is a rectangular or polygonal crack.
!############################################
if (c_Crack_Type == 1) then
    num_Cr_Poi = Each_Cr_Poi_Num(i_C)
    CR_Center_Old(1) =sum(Crack3D_Coor(i_C,1:num_Cr_Poi,1))/num_Cr_Poi
    CR_Center_Old(2) =sum(Crack3D_Coor(i_C,1:num_Cr_Poi,2))/num_Cr_Poi
    CR_Center_Old(3) =sum(Crack3D_Coor(i_C,1:num_Cr_Poi,3))/num_Cr_Poi
    
    !To be done.
    print *,'    To be done in D3_Offset_Initial_Crack_Along_Normal.f.'
endif
      
!###############################################################
! If it is a circular crack: just update the coordinates of the 
! center, and leave the other variables unchanged.
!###############################################################
if (c_Crack_Type == 2) then
    ! Original crack is circular.
    CR_Center_Old(1:3) = Crack3D_Cir_Coor(i_C,1:3)
    ! Normal direction.
    n_Vector(1:3) = Crack3D_Cir_Coor(i_C,4:6)
    ! New circular crack.
    CR_Center_New(1:3) = CR_Center_Old(1:3) + delta_L*n_Vector(1:3)
    ! Update center coordinates
    Crack3D_Cir_Coor(i_C,1:3) = CR_Center_New(1:3) 
    
    
endif

!############################################################
! If it is an elliptical crack, it is not supported for now.
!############################################################
if (c_Crack_Type == 3) then
    print *,'     ERROR(2022080202) :: Elliptical cracks not supported yet!'
    call Warning_Message('S',Keywords_Blank) 
endif
      
      
RETURN
END SUBROUTINE D3_Offset_Initial_Crack_Along_Normal
