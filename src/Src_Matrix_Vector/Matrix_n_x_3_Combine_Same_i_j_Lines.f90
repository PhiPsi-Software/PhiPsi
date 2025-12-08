 
subroutine Matrix_n_x_3_Combine_Same_i_j_Lines(n, M, newM, new_n)
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
