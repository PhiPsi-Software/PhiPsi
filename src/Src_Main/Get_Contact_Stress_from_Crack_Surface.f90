!-----------------------------------------------------------
! Brief: Retrieve contact stress at a Gauss point from 3D crack data.
!
! Parameters:
!   Input:  Current_Elem    - current element index
!           GP_face_global  - global coordinates of the Gauss point
!           n_unit          - unit normal of the local face
!   Output: t_global        - retrieved contact traction vector
!
! Notes:   Finds the nearest fluid-element centroid and copies the
!          stored contact stress. Direction is reversed for explicit
!          dynamic analyses (Key_Analysis_Type==6).
!-----------------------------------------------------------

SUBROUTINE Get_Contact_Stress_from_Crack_Surface(Current_Elem,GP_face_global,n_unit, t_global)
! Get the contact stress of the crack surface. 
! First, calculate the distance between the current GP and the centroids of all contact elements.
! Then, find the minimum distance and obtain the corresponding contact stress.
!
! 2026-04-30.
!

!-----------------------------          
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Stress
use Global_Material
use Global_Crack_3D
use Global_Crack_Common
use omp_lib
use Global_XFEM_Elements
      
!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
integer, intent(in)::Current_Elem
real(kind=FT), intent(in)::GP_face_global(3),n_unit(3)
real(kind=FT), intent(out)::t_global(3)
     
integer i_C
integer i_CT_Elem
real(kind=FT) c_Center(3)
real(kind=FT),allocatable::Crack_Dis(:,:)  
real(kind=FT) Tool_Function_2Point_Dis_3D
integer best_crack,best_elem,idx(2)

!Initialized to 0
t_global(1:3) = ZR


!Get from Contact_Stress(1:num_Crack,1:Max_Max_N_FluEl_3D,1:3)

allocate(Crack_Dis(num_Crack,Max_Max_N_FluEl_3D))

Crack_Dis = TEN_15

! Circulation between each crack
do i_C=1,num_Crack
    !======================================
    ! Circulation between contact elements
    !======================================
    do i_CT_Elem = 1,Cracks_FluidEle_num_3D(i_C)
        ! Centroid of the current contact element (fluid element)
        c_Center = Cracks_FluidEle_Centroid_3D(i_C)%row(i_CT_Elem,1:3)
        Crack_Dis(i_C,i_CT_Elem) = Tool_Function_2Point_Dis_3D(GP_face_global,c_Center)
    enddo
enddo

idx = minloc(Crack_Dis, mask = Crack_Dis < TEN_14)  
if (all(idx == 0)) then
    !do nothing.
else
    best_crack = idx(1)
    best_elem = idx(2)
end if

if (allocated(Contact_Stress)) then
    t_global(1:3) = Contact_Stress(best_crack,best_elem,1:3)
!    
!    t_global(1:3) = -Contact_Stress(best_crack,best_elem,1:3)
    
!    t_global(1:3) = -Contact_Stress(best_crack,best_elem,1:3)*0.1D0

    !In explicit dynamic analysis, the direction of contact stress on the
    !contact surface needs to be reversed.
    
    if (Key_Analysis_Type == 6) then
!        t_global(1) =  ZR
!        t_global(2) =  ZR
!        t_global(3) = -t_global(3)*0.1D0
        t_global(1:3) = -Contact_Stress(best_crack,best_elem,1:3)
    end if
    
endif

deallocate(Crack_Dis)



      
RETURN
END SUBROUTINE Get_Contact_Stress_from_Crack_Surface
