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
 
subroutine Matrix_n_x_3_Combine_Same_i_j_Lines(n, M, newM, new_n)
!
! Sort an n-row, 3-column matrix.
! The first and second columns represent positions, and the third column represents values.
! The result is sorted in ascending order by the first column, and on this basis, the second column
! is sorted in ascending order.
!
! 2024-09.
!
use Global_Float_Type 
use module_INTERFACE_Matrix_n_x_3_Check_Sort

implicit none
integer, intent(in) :: n
real(kind=FT), intent(in) :: M(n, 3)
real(kind=FT), intent(out) :: newM(n, 3)
integer, intent(out) :: new_n
integer :: i
real(kind=FT) Sorted_M(n,3)
real(kind=FT) Sorted_M_2(n,3)
integer Matrix_M1_M2(n,2)
integer idx(n)
LOGICAL check_matrix
integer count_error
integer i_Try

Sorted_M = M
do i_Try = 1,100
    Matrix_M1_M2(1:n,1) =  INT(Sorted_M(1:n,1))
    Matrix_M1_M2(1:n,2) =  INT(Sorted_M(1:n,2))
    idx = [(i, i=1,n)]
    
    call Matrix_n_x_2_Quick_Sort_Int(Matrix_M1_M2(1:n,1:2),n,idx(1:n))  

    
    do i = 1, n
        Sorted_M_2(i,1:3) = Sorted_M(idx(i),1:3) 
    enddo
        
    Sorted_M  = Sorted_M_2
        
    call Matrix_n_x_3_Check_Sort(Sorted_M, n, check_matrix,count_error)
        
    if(check_matrix)then
        exit
    endif    
enddo  

if (check_matrix .eqv. .false.) then
   write ( *, '(a)' ) 'ERROR :: Failed to sort the matrix! In Matrix_n_x_3_Combine_Same_i_j_Lines.f90!'
   call Warning_Message('S',' ')
endif


new_n = 0
newM = 0.0D0
new_n = 1
newM(new_n, :) = Sorted_M(1, :)
do i = 2, n
    if (Sorted_M(i, 1) == newM(new_n, 1) .and. Sorted_M(i, 2) == newM(new_n, 2)) then
        newM(new_n, 3) = newM(new_n, 3) + Sorted_M(i, 3)
    else
        new_n = new_n + 1
        newM(new_n, :) = Sorted_M(i, :)
    end if
end do
    
end subroutine Matrix_n_x_3_Combine_Same_i_j_Lines
