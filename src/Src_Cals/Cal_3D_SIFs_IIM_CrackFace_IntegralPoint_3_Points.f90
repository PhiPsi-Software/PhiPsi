!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !


subroutine Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_3_Points(i_Crack, Elem_ID, nGP_face, GP_global, GP_w, GP_n_global)
!Compute integral points of the fluid element at the crack surface.
!2026-01-28.

use Global_Float_Type
use Global_Common
use Global_Crack_3D
use Global_Elem_Area_Vol

implicit none

integer, intent(in) :: i_Crack, Elem_ID
integer, intent(out) :: nGP_face
real(kind=FT), allocatable, intent(out) :: GP_global(:,:), GP_w(:), GP_n_global(:,:)
integer :: i_FluidEle, nFluid, countTri, ig
real(kind=FT) :: x1(3), x2(3), x3(3)
real(kind=FT) :: Atri, nvec(3), nrm
real(kind=FT) :: l1, l2, l3
integer :: id1, id2, id3
integer iTri
real(kind=FT) :: area_i,eps_shift

real(kind=FT), parameter :: L(3,3) = reshape([ &
    2.0D0/3.0D0, 1.0D0/6.0D0, 1.0D0/6.0D0, &
    1.0D0/6.0D0, 2.0D0/3.0D0, 1.0D0/6.0D0, &
    1.0D0/6.0D0, 1.0D0/6.0D0, 2.0D0/3.0D0  &
], shape(L))

nGP_face = 0

if (i_Crack <= 0) return
if (Crack_Type_Status_3D(i_Crack,1) /= 1) return

nFluid = Cracks_FluidEle_num_3D(i_Crack)
if (nFluid <= 0) return

countTri = 0
do i_FluidEle = 1, nFluid
    if (Cracks_FluidEle_EleNum_3D(i_Crack)%row(i_FluidEle) == Elem_ID) then
        countTri = countTri + 1
    endif
enddo
if (countTri <= 0) return

nGP_face = 3 * countTri
allocate(GP_global(nGP_face,3))
allocate(GP_w(nGP_face))
allocate(GP_n_global(nGP_face,3))

GP_global   = ZR
GP_w        = ZR
GP_n_global = ZR
eps_shift   = 1.0D-10*Ave_Elem_L_Enrich

ig = 0
do i_FluidEle = 1, nFluid
    if (Cracks_FluidEle_EleNum_3D(i_Crack)%row(i_FluidEle) /= Elem_ID) cycle

    id1 = Cracks_FluidEle_CalP_3D(i_Crack)%row(i_FluidEle,1)     
    id2 = Cracks_FluidEle_CalP_3D(i_Crack)%row(i_FluidEle,2)  
    id3 = Cracks_FluidEle_CalP_3D(i_Crack)%row(i_FluidEle,3)
    
    x1  = Cracks_CalP_Coors_3D(i_Crack)%row(id1,1:3)
    x2  = Cracks_CalP_Coors_3D(i_Crack)%row(id2,1:3)
    x3  = Cracks_CalP_Coors_3D(i_Crack)%row(id3,1:3)
    
    nvec = Cracks_FluidEle_Vector_3D(i_Crack)%row(i_FluidEle,1:3)
    Atri = Cracks_FluidEle_Area_3D(i_Crack)%row(i_FluidEle)

    do iTri = 1, 3
        l1 = L(iTri,1)
        l2 = L(iTri,2)
        l3 = L(iTri,3)
        ig = ig + 1
        GP_global(ig,1:3) = l1*x1 + l2*x2 + l3*x3
        GP_global(ig,1:3) = GP_global(ig,1:3) + eps_shift *nvec
        GP_w(ig)            = Atri / 3.0D0
        GP_n_global(ig,1:3) = nvec
    enddo
enddo
end subroutine Cal_3D_SIFs_IIM_CrackFace_IntegralPoint_3_Points



