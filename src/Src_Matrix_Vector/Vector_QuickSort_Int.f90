recursive subroutine Vector_QuickSort_Int(a, first, last)
! Quick sort, firstly written on 2021-09-08.
! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
use Global_Float_Type
implicit none
integer a(*), x, t
integer first, last
integer i, j

x = a( (first+last) / 2 )
i = first
j = last
do
    do while (a(i) < x)
        i=i+1
    end do
    do while (x < a(j))
        j=j-1
    end do
    if (i >= j) exit
    t = a(i);  a(i) = a(j);  a(j) = t
    i=i+1
    j=j-1
end do
if (first < i-1) call Vector_QuickSort_Int(a, first, i-1)
if (j+1 < last)  call Vector_QuickSort_Int(a, j+1, last)
end subroutine Vector_QuickSort_Int
