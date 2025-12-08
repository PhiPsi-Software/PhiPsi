 
module fmin_module

use iso_fortran_env
use Global_Float_Type

implicit none

private

abstract interface
    function func(x) result(f)
    use Global_Float_Type
    implicit none
    real(kind=FT),intent(in) :: x
    real(kind=FT)            :: f
    end function func
end interface

public :: fmin

contains


function fmin(f,ax,bx,tol) result(xmin)
implicit none

procedure(func)     :: f
real(kind=FT),intent(in) :: ax
real(kind=FT),intent(in) :: bx
real(kind=FT),intent(in) :: tol
real(kind=FT)            :: xmin

real(kind=FT) :: a,b,d,e,xm,p,q,r,tol1,tol2,u,v,w
real(kind=FT) :: fu,fv,fw,fx,x
real(kind=FT) :: abs,sqrt,sign
real(kind=FT) :: c,half,sqrteps

c = (THR-sqrt(FIV))/TWO
half = ONE/TWO
sqrteps = sqrt(epsilon(ONE))

a = ax
b = bx
v = a + c*(b - a)
w = v
x = v
e = ZR
fx = f(x)
fv = fx
fw = fx

do

    xm = half*(a + b)
    tol1 = sqrteps*abs(x) + tol/THR
    tol2 = TWO*tol1


    if (abs(x - xm) <= (tol2 - half*(b - a))) then
        exit
    end if


    if (abs(e) <= tol1) then


        if (x >= xm) then
            e = a - x
        else
            e = b - x
        end if
        d = c*e

    else


        r = (x - w)*(fx - fv)
        q = (x - v)*(fx - fw)
        p = (x - v)*q - (x - w)*r
        q = TWO*(q - r)
        if (q > ZR) p = -p
        q =  abs(q)
        r = e
        e = d


        if ((abs(p) >= abs(half*q*r)) .or. (p <= q*(a - x)) .or. (p >= q*(b - x))) then


            if (x >= xm) then
                e = a - x
            else
                e = b - x
            end if
            d = c*e

        else


            d = p/q
            u = x + d


            if (((u - a) < tol2) .or. ((b - u) < tol2)) d = sign(tol1, xm - x)

        end if

    end if


    if (abs(d) >= tol1) then
        u = x + d
    else
        u = x + sign(tol1, d)
    end if
    fu = f(u)


    if (fu <= fx) then
        if (u >= x) a = x
        if (u < x) b = x
        v = w
        fv = fw
        w = x
        fw = fx
        x = u
        fx = fu
    else
        if (u < x) a = u
        if (u >= x) b = u
        if (fu <= fw .or. w == x) then
            v = w
            fv = fw
            w = u
            fw = fu
        else if (fu <= fv .or. v == x .or. v == w ) then
            v = u
            fv = fu
        end if
    end if

end do

xmin = x

end function fmin

end module fmin_module
