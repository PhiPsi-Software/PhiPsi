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

subroutine Tool_find_top_m_indices_of_vector(vectorA, n, m, vectorB)
use Global_Float_Type
implicit none

integer, intent(in) :: n
integer, intent(in) :: m
real(kind=FT), intent(in) :: vectorA(n)
integer, intent(out) :: vectorB(m)

real(kind=FT) :: tempA(n)
integer :: indices(n)
integer :: i, j
real(kind=FT) :: tempVal
integer :: tempIdx

if (m > n .or. m <= 0) then
    print *, "Error: m must be between 1 and n"
    return
end if

tempA = vectorA

do i = 1, n
    indices(i) = i
end do

do i = 1, n-1
    do j = 1, n-i
        if (tempA(j) < tempA(j+1)) then
            tempVal = tempA(j)
            tempA(j) = tempA(j+1)
            tempA(j+1) = tempVal
            
            tempIdx = indices(j)
            indices(j) = indices(j+1)
            indices(j+1) = tempIdx
        end if
    end do
end do

do i = 1, m
    vectorB(i) = indices(i)
end do
    
end subroutine Tool_find_top_m_indices_of_vector
