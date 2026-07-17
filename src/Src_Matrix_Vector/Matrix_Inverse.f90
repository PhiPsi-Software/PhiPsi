!-----------------------------------------------------------
! Brief: Invert a square matrix using Doolittle LU factorization.
!
! Parameters:
!   Input:  n - order of the matrix
!           a - input square matrix (n,n)
!   Output: c - inverse matrix (n,n)
!
! Notes:   Builds L and U from the LU factors, then solves A*x = e_k
!          for each unit vector e_k via forward and back substitution.
!-----------------------------------------------------------

subroutine Matrix_Inverse(a,c,n)   
!     Find the inverse matrix of square matrix a and output it to c.

!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
use Global_Float_Type
implicit none 
integer,intent(in)::n
real(kind=FT),intent(in):: a(n,n)
real(kind=FT),intent(out):: c(n,n)
real(kind=FT) L(n,n), U(n,n), b(n), d(n), x(n)
real(kind=FT) coeff,tem_a(n,n)
integer i, j, k

tem_a = a
L=ZR
U=ZR
b=ZR

do k=1, n-1
    do i=k+1,n
        coeff=tem_a(i,k)/tem_a(k,k)
        L(i,k) = coeff
        do j=k+1,n
            tem_a(i,j) = tem_a(i,j)-coeff*tem_a(k,j)          
        end do
    end do
end do

do i=1,n
    L(i,i) = ONE
end do
do j=1,n
    do i=1,j
        U(i,j) = tem_a(i,j)
    end do
end do

do k=1,n
    b(k)=ONE
    d(1) = b(1)
    do i=2,n
        d(i)=b(i)
        do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
        end do
    end do
    x(n)=d(n)/U(n,n)
    do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
    end do
    do i=1,n
        c(i,k) = x(i)
    end do
    b(k)=ZR
end do

return
END subroutine Matrix_Inverse



