!-----------------------------------------------------------
! Brief: Convert input data from SI to the mm-ton-MPa system.
!
! Notes:   When Key_Unit_System==2, rescales material, geometry, crack,
!          and load parameters to the internal mm-ton-MPa unit system.
!-----------------------------------------------------------

subroutine Unit_Conversion
!     Used to convert input information in the International System of Elements to the corresponding parameters in the mm-ton-MPa unit system

!     ----------------------------
!     Read public variable module
!     ----------------------------
use Global_Float_Type 
use Global_Common   
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack
use Global_Crack_Common
use Global_HF

!     --------------------------
!     Variable Type Declaration
!     --------------------------
implicit none
integer i_C,i_P

print *,'    Conversing unit system for inout parameters......'

if(Key_Unit_System==2)then
    ! Material Parameter Conversion
    Inject_Q_Val(1:20) = Inject_Q_Val(1:20)*1.0D6
    Viscosity = Viscosity/1.0D6
    Material_Para(1,1) = Material_Para(1,1)/1.0D6
    Material_Para(1,3) = Material_Para(1,3)/1.0D12
    Material_Para(1,4) = Material_Para(1,4)*1.0D3
    Material_Para(1,5) = Material_Para(1,5)/1.0D6
    Material_Para(1,6) = Material_Para(1,6)/1.0D6
    Material_Para(1,7) = Material_Para(1,7)/1.0D6
    ! Transformation of node coordinates
    Coor(1:Num_Node,1:2) = Coor(1:Num_Node,1:2)*1.0D3
    Max_X_Coor = Max_X_Coor*1.0D3
    Min_X_Coor = Min_X_Coor*1.0D3
    Max_Y_Coor = Max_Y_Coor*1.0D3
    Min_Y_Coor = Min_Y_Coor*1.0D3
    G_X_NODES(1:4,1:Num_Elem) = G_X_NODES(1:4,1:Num_Elem)*1.0D3
    G_Y_NODES(1:4,1:Num_Elem) = G_Y_NODES(1:4,1:Num_Elem)*1.0D3
    Elem_Area(1:Num_Elem) = Elem_Area(1:Num_Elem)*1.0D6
    Elem_Centroid(1:Num_Elem,1:2) = Elem_Centroid(1:Num_Elem,1:2) *1.0D3
    Max_Elem_Area = Max_Elem_Area*1.0D6
    Min_Elem_Area = Min_Elem_Area*1.0D6
    Ave_Elem_Area = Ave_Elem_Area*1.0D6
    Ave_Elem_L    = Ave_Elem_L*1.0D3
    ! Transformation of initial crack coordinates
    do i_C=1,num_Crack
        do i_P = 1,Each_Cr_Poi_Num(i_C)
            Crack_Coor(i_C,i_P,1:2)=Crack_Coor(i_C,i_P,1:2)*1.0D3
        enddo
    enddo
    ! Transformation of Initial Natural Fracture Coordinates
    do i_C=1,num_Na_Crack
        do i_P = 1,Each_Na_Cr_Poi_Num(i_C)
            Na_Crack_Coor(i_C,i_P,1:2)= Na_Crack_Coor(i_C,i_P,1:2)*1.0D3
        enddo
    enddo
    ! Conversion of external load (Note: In the mm-ton-s unit system, the element of load is still N,
    ! It itself does not need to change, but due to the model's thickness increasing from the element
    ! thickness of 1m
    ! It has been changed to a unit thickness of 1mm, so the load should be
    Foc_x(1:Num_Foc_x,2) = Foc_x(1:Num_Foc_x,2) / 1.0D3
    Foc_y(1:Num_Foc_x,2) = Foc_y(1:Num_Foc_x,2) / 1.0D3
    ! Normal and tangential stiffness of crack surface friction contact
    k_tt = k_tt /1.0D9              
    k_nn = k_nn /1.0D9               
endif

return
END subroutine Unit_Conversion
