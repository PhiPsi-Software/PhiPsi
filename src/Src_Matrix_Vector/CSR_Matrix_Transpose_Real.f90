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
 
subroutine CSR_Matrix_Transpose_Real(m, n, nnz, values, col_ind, row_ptr, values_t, col_ind_t, row_ptr_t)
!2024-10-31
! Transpose of a CSR format matrix.
use Global_Float_Type  
implicit none

integer, intent(in) :: m, n, nnz
real(kind=FT), intent(in) :: values(nnz)
integer, intent(in) :: col_ind(nnz), row_ptr(m+1)

real(kind=FT), intent(out) :: values_t(nnz)
integer, intent(out) :: col_ind_t(nnz), row_ptr_t(n+1)

integer :: i, j, k, count
integer, allocatable :: row_counts(:)

row_ptr_t = 0

allocate(row_counts(n))
row_counts = 0
do i = 1, nnz
    row_counts(col_ind(i)) = row_counts(col_ind(i)) + 1
end do

row_ptr_t(1) = 1
do i = 1, n
    row_ptr_t(i+1) = row_ptr_t(i) + row_counts(i)
end do

row_counts = 0

do i = 1, m
    do j = row_ptr(i), row_ptr(i+1) - 1
        k = col_ind(j)
        count = row_ptr_t(k) + row_counts(k)
        values_t(count) = values(j)
        col_ind_t(count) = i
        row_counts(k) = row_counts(k) + 1
    end do
end do

deallocate(row_counts)
    
end subroutine CSR_Matrix_Transpose_Real