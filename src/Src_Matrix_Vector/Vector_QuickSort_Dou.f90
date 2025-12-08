 
recursive subroutine Vector_QuickSort_Dou(a, first, last)
use Global_Float_Type

implicit none
real(kind=FT)  a(*), x, t
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

if (first < i-1) call Vector_QuickSort_Dou(a, first, i-1)
if (j+1 < last)  call Vector_QuickSort_Dou(a, j+1, last)
end subroutine Vector_QuickSort_Dou
