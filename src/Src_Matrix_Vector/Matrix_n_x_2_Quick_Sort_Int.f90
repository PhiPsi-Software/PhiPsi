 
subroutine Matrix_n_x_2_Quick_Sort_Int(Matrix, n, Index)

    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: Matrix(n,2)
    integer, intent(out) :: Index(n)
    
    integer :: i
    
    do i = 1, n
        Index(i) = i
    end do
    
    call Quick_Sort(1, n)
    
contains
    recursive subroutine Quick_Sort(left, right)
        integer, intent(in) :: left, right
        integer :: i, j, pivot_idx, temp
        
        if(left >= right) return
        
        pivot_idx = (left + right)/2
        i = left
        j = right
        
        do while(i <= j)
            do while(Compare_Less(i, pivot_idx))
                i = i + 1
            end do
            
            do while(Compare_Less(pivot_idx, j))
                j = j - 1
            end do
            
            if(i <= j) then
                temp = Index(i)
                Index(i) = Index(j)
                Index(j) = temp
                
                if(pivot_idx == i) then
                    pivot_idx = j
                else if(pivot_idx == j) then
                    pivot_idx = i
                end if
                
                i = i + 1
                j = j - 1
            end if
        end do
        
        if(left < j) call Quick_Sort(left, j)
        if(i < right) call Quick_Sort(i, right)
    end subroutine Quick_Sort
    
    logical function Compare_Less(i1, i2)
        integer, intent(in) :: i1, i2
        
        if(Matrix(Index(i1),1) < Matrix(Index(i2),1)) then
            Compare_Less = .true.
            return
        else if(Matrix(Index(i1),1) > Matrix(Index(i2),1)) then
            Compare_Less = .false.
            return
        end if
        
        Compare_Less = Matrix(Index(i1),2) < Matrix(Index(i2),2)
    end function Compare_Less
    
end subroutine Matrix_n_x_2_Quick_Sort_Int