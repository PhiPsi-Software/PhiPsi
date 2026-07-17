!-----------------------------------------------------------
! Brief: Mesh a circular 3D crack surface using concentric rings.
!
! Parameters:
!   Input:  i_C          - crack index
!           cir_center   - circle center (x,y,z)
!           cir_vec      - circle normal vector
!           radius       - circle radius
!           elemL        - target element edge length
!   In/Out: Ele_Num_Cache - element-number cache (grows as needed)
!
! Notes:   Builds an orthonormal (a,b) basis on the circle plane,
!          then lays out a center node plus nR radial rings of
!          nTheta nodes each, with consistent CCW numbering.
!-----------------------------------------------------------

subroutine D3_Mesh_Circle_By_Rings(i_C, cir_center, cir_vec, radius, elemL, Ele_Num_Cache)
! Circular crack meshing (no crossing, consistent numbering).
!2026-02-14.
!
use Global_Float_Type
use Global_Common
use Global_Crack_3D
use Global_Cal_Ele_Num_by_Coors_3D

implicit none
integer, intent(in) :: i_C
real(kind=FT), intent(in) :: cir_center(3), cir_vec(3), radius, elemL
integer, intent(inout) :: Ele_Num_Cache

real(kind=FT) a(3), b(3), norm_a, norm_b
real(kind=FT) theta, r, dr
integer :: nR, nTheta
integer :: iRing, j
integer :: nodeCount, eleCount
integer :: in_Elem_num, i_Cr_Node
real(kind=FT) c_Kesi, c_Yita, c_Zeta
integer :: maxN
integer, allocatable :: ringStart(:)

integer :: n1, n2, n3, n4, jp
real(kind=FT) :: circumference

! Build a stable orthonormal basis (a,b) on the circle plane
call Vector_Cross_Product_3(cir_vec,[ONE,ZR,ZR],a)
if(sum(abs(a))<=Tol_20) then
    call Vector_Cross_Product_3(cir_vec,[ZR,ONE,ZR],a)
endif
call Vector_Cross_Product_3(cir_vec,a,b)
call Vector_Norm2(3,a,norm_a)
call Vector_Norm2(3,b,norm_b)
a = a/norm_a
b = b/norm_b

! Decide radial layers
nR = int(ceiling(radius/max(elemL, Tol_20)))
if(nR < 2) nR = 2
dr = radius / dble(nR)

! Decide angular divisions using outer circumference
circumference = TWO*pi*radius
nTheta = int(ceiling(circumference/max(elemL, Tol_20)))
if(nTheta < 12) nTheta = 12
if(mod(nTheta,2)/=0) nTheta = nTheta + 1

! Node layout:
! node 1 : center
! ring iRing=1..nR: nTheta nodes each, counterclockwise in (a,b) plane
allocate(ringStart(0:nR))
ringStart = 0
ringStart(0) = 1
do iRing = 1, nR
    ringStart(iRing) = 2 + (iRing-1)*nTheta
enddo

nodeCount = 1 + nR*nTheta

! Ensure memory for nodes
maxN = size(Crack3D_Meshed_Node(i_C)%row,1)
do while(nodeCount > maxN)
    call D3_Allocate_Crack_Memory(i_C,1,1)
    maxN = size(Crack3D_Meshed_Node(i_C)%row,1)
enddo

! Write nodes
Crack3D_Meshed_Node(i_C)%row(1,1:3) = cir_center(1:3)

do iRing = 1, nR
    r = dr*dble(iRing)
    do j = 1, nTheta
        theta = dble(j-1)*TWO*pi/dble(nTheta)
Crack3D_Meshed_Node(i_C)%row(ringStart(iRing)+j-1,1) = cir_center(1) + r*a(1)*cos(theta) + r*b(1)*sin(theta)
Crack3D_Meshed_Node(i_C)%row(ringStart(iRing)+j-1,2) = cir_center(2) + r*a(2)*cos(theta) + r*b(2)*sin(theta)
Crack3D_Meshed_Node(i_C)%row(ringStart(iRing)+j-1,3) = cir_center(3) + r*a(3)*cos(theta) + r*b(3)*sin(theta)
    enddo
enddo

Crack3D_Meshed_Node_num(i_C) = nodeCount

! Build triangle connectivity with a unified rule
! Rule target: triangle normal should be consistent with cir_vec.
eleCount = 0

! Center to first ring (fan), CCW: (center, j, j+1)
do j = 1, nTheta
    jp = j + 1
    if(jp > nTheta) jp = 1

    eleCount = eleCount + 1
    maxN = size(Crack3D_Meshed_Ele(i_C)%row,1)
    if(eleCount > maxN) call D3_Allocate_Crack_Memory(i_C,1,1)

    n1 = 1
    n2 = ringStart(1) + j  - 1
    n3 = ringStart(1) + jp - 1
    call D3_Ensure_Tri_Normal_Consistent(i_C, n1, n2, n3, cir_vec)
    Crack3D_Meshed_Ele(i_C)%row(eleCount,1:3) = [n1,n2,n3]
enddo

! Between rings: each quad (inner j->j+1, outer j->j+1) split into two triangles
! Consistent split:
!   T1 = (inner_j, outer_j, outer_jp)
!   T2 = (inner_j, outer_jp, inner_jp)
do iRing = 2, nR
    do j = 1, nTheta
        jp = j + 1
        if(jp > nTheta) jp = 1

        n1 = ringStart(iRing-1) + j  - 1
        n2 = ringStart(iRing)   + j  - 1
        n3 = ringStart(iRing)   + jp - 1
        n4 = ringStart(iRing-1) + jp - 1

        eleCount = eleCount + 1
        maxN = size(Crack3D_Meshed_Ele(i_C)%row,1)
        if(eleCount > maxN) call D3_Allocate_Crack_Memory(i_C,1,1)
        call D3_Ensure_Tri_Normal_Consistent(i_C, n1, n2, n3, cir_vec)
        Crack3D_Meshed_Ele(i_C)%row(eleCount,1:3) = [n1,n2,n3]

        eleCount = eleCount + 1
        maxN = size(Crack3D_Meshed_Ele(i_C)%row,1)
        if(eleCount > maxN) call D3_Allocate_Crack_Memory(i_C,1,1)
        call D3_Ensure_Tri_Normal_Consistent(i_C, n1, n3, n4, cir_vec)
        Crack3D_Meshed_Ele(i_C)%row(eleCount,1:3) = [n1,n3,n4]
    enddo
enddo

Crack3D_Meshed_Ele_num(i_C) = eleCount

! Outline detection
call D3_Find_Crack_Mesh_Outline(i_C,eleCount,nodeCount,Crack3D_Meshed_Ele(i_C)%row(1:eleCount,1:3))

! Locate nodes in background elements and compute local coordinates
do i_Cr_Node = 1, nodeCount
call Cal_Ele_Num_by_Coors_3D( Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,1), Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,2), &
Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,3), Ele_Num_Cache, in_Elem_num)

    if(in_Elem_num==0 .and. Key_Warning_Level>=3) then
        print *,'    WARN :: in_Elem_num=0 in D3_Mesh_Circle_By_Rings!'
        print *,'           Crack number:',i_C
        print *,'           Coors:',Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,1:3)
        print *,'           Maybe outside the model!'
    endif

    Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(i_Cr_Node) = in_Elem_num
    if(in_Elem_num>0) then
        call Cal_KesiYita_by_Coor_3D(Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,1:3),in_Elem_num,c_Kesi,c_Yita,c_Zeta)
        Cr3D_Meshed_Node_in_Ele_Local(i_C)%row(i_Cr_Node,1:3) = [c_Kesi,c_Yita,c_Zeta]
    endif
enddo

deallocate(ringStart)

end subroutine D3_Mesh_Circle_By_Rings