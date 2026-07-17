!-----------------------------------------------------------
! Brief: Returns the indices of the m largest entries in a real vector.
!
! Parameters:
!   Input:  vectorA(n) - input real vector of length n
!   Input:  n - length of vectorA
!   Input:  m - number of top entries to return (1 <= m <= n)
!   Output: vectorB(m) - indices into vectorA of the m largest values
!
! Notes:   Uses descending-order bubble sort on a working copy together
!          with a parallel index array so original values are preserved.
!-----------------------------------------------------------

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
