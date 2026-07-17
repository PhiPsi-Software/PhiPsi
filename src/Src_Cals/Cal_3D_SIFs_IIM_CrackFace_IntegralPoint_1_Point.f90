!-----------------------------------------------------------
! Brief: Collect 1-point face Gauss data from HF fluid elements.
!
! Parameters:
!   Input:  i_Crack       - crack index
!           Elem_ID       - host solid element id
!   Output: nGP_face      - number of face Gauss points
!           GP_global     - global coordinates of each GP
!           GP_w          - integration weights (fluid-ele area)
!           GP_n_global   - outward unit normal at each GP
!
! Notes:   Uses the centroid of each HF fluid element as a
!          1-point Gauss point; skips non-HF cracks. GP positions
!          are nudged by eps_shift along the normal for stability.
!-----------------------------------------------------------

subroutine Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_1_Point(i_Crack, Elem_ID, nGP_face, GP_global, GP_w, GP_n_global)
! Get quadrature points on crack face within a given solid element, using
! pre-built HF fluid elements (centroid 1-point rule).
! 2026-01-28.
!
! OUTPUT:
!   nGP_face               : number of face Gauss points found in this Elem_ID
!   GP_global(nGP_face,3)  : global coordinates of GP (centroids)
!   GP_w(nGP_face)         : integration weights = area of corresponding fluid element
!   GP_n_global(nGP_face,3): unit normal at GP (from fluid element data)
use Global_Float_Type
use Global_Common
use Global_Crack_3D
use Global_Elem_Area_Vol

implicit none

integer, intent(in) :: i_Crack, Elem_ID
integer, intent(out) :: nGP_face
real(kind=FT), allocatable, intent(out) :: GP_global(:,:), GP_w(:), GP_n_global(:,:)

integer :: i_FluidEle, count, nFluid
real(kind=FT) :: area_i
integer, allocatable :: idx_list(:)
real(kind=FT) :: eps_shift,nvec(3)

nGP_face = 0

if (i_Crack <= 0) return
if (Crack_Type_Status_3D(i_Crack,1) /= 1) return

nFluid = Cracks_FluidEle_num_3D(i_Crack)
if (nFluid <= 0) return

count = 0
do i_FluidEle = 1, nFluid
    if (Cracks_FluidEle_EleNum_3D(i_Crack)%row(i_FluidEle) == Elem_ID) then
        count = count + 1
    endif
enddo

if (count <= 0) return

nGP_face = count
allocate(GP_global(nGP_face,3))
allocate(GP_w(nGP_face))
allocate(GP_n_global(nGP_face,3))

GP_global   = ZR
GP_w        = ZR
GP_n_global = ZR

eps_shift   = 1.0D-10*Ave_Elem_L_Enrich

count = 0
do i_FluidEle = 1, nFluid
    if (Cracks_FluidEle_EleNum_3D(i_Crack)%row(i_FluidEle) == Elem_ID) then
        area_i = Cracks_FluidEle_Area_3D(i_Crack)%row(i_FluidEle)
        nvec = Cracks_FluidEle_Vector_3D(i_Crack)%row(i_FluidEle,1:3)

        count = count + 1
        
        GP_global(count,1:3) = Cracks_FluidEle_Centroid_3D(i_Crack)%row(i_FluidEle,1:3)
        
        GP_global(count,1:3) = GP_global(count,1:3) + eps_shift *nvec

        GP_w(count) = area_i

        GP_n_global(count,1:3) = nvec
    endif
enddo

end subroutine Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_1_Point


