!-----------------------------------------------------------
! Brief: Apply a diagonal Jacobi-style preconditioner to a dense K.
!
! Parameters:
!   Input:     n   - matrix dimension
!   In/Out:    K   - dense n x n matrix to precondition (modified)
!   Output: Vector_Dc - scaling vector (sqrt(mean diag) / sqrt(diag))
!
! Notes:   Forms the row/column scaling from the diagonal of K so
!          that K_ij becomes Dc_i * Dc_j * K_ij. Aborts via the
!          warning message if the minimum diagonal is non-positive.
!-----------------------------------------------------------

subroutine Apply_Pre_Conditioner_to_K(n, K,Vector_Dc)

!     Apply the Pre-Conditioner to K.
!     Written first by Fang Shi on 2020-02-12.

!     ----------------------------
!     Read public variable module
!     ----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
!     --------------------------
!     Variable Type Declaration
!     --------------------------
implicit none

integer,intent(in)::n
real(kind=FT),intent(inout)::K(n,n)
real(kind=FT),intent(out)::Vector_Dc(n)
real(kind=FT) :: diagona(n)
integer i,j
real(kind=FT) max_diagona,min_diagona,tem_diagona

!Get diagonal elements of K.
call Matrix_Get_Diagonal_Elements(n,K,diagona)   
max_diagona = maxval(diagona)
min_diagona = minval(diagona)
print *,'    Max value of  diagonal element of K:',max_diagona
print *,'    Min value of  diagonal element of K:',min_diagona      

do i =1,n
    if(diagona(i)==ZR)then
        print *,i
    endif
enddo

!The min diagonal element must be positive.
if (min_diagona<=ZR)then
    Print *,'    ERROR :: in Apply_Pre_Conditioner_to_K.f!'
    Print *,'    Error code: 777'
    call Warning_Message('S',Keywords_Blank)  
endif

tem_diagona = sqrt((max_diagona+min_diagona)/TWO)
!Loop to get Vector_Dc.
do i =1,n
    Vector_Dc(i)=tem_diagona/sqrt(diagona(i))
enddo

!Apply Vector_Dc to K.
do i =1,n
    do j =1,n
        K(i,j) = K(i,j)*Vector_Dc(i)*Vector_Dc(j) 
    enddo                      
enddo

return
END subroutine Apply_Pre_Conditioner_to_K
