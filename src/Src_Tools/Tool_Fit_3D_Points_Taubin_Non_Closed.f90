subroutine Tool_Fit_3D_Points_Taubin_Non_Closed(i_C,num_Point,In_Points,Points_Flag_Inside,Out_Points)

use Global_Float_Type     
use Global_Model
use Global_Common

implicit none
integer,intent(in)::i_C,num_Point
real(kind=FT),intent(in)::In_Points(num_Point,3)
logical,intent(in)::Points_Flag_Inside(num_Point)
real(kind=FT),intent(out)::Out_Points(num_Point,3)

integer :: i_start, i_end, seg_len
integer :: i_point
integer :: num_segments

real(kind=FT) :: point_temp(num_Point,3)
integer :: n_Loops, i_Loop, i_pre, i_nex
real(kind=FT) :: lambda, mu
real(kind=FT) :: c_point(3), pre_point(3), nex_point(3)
real(kind=FT) :: pre_L, nex_L, weight_pre, weight_nex, total_weight
real(kind=FT) :: L1_vector(num_Point,3), L2_vector(num_Point,3)
logical :: found_start, in_segment

n_Loops = 20
lambda  =  0.330D0
mu      = -0.331D0

Out_Points(1:num_Point,1:3) = In_Points(1:num_Point,1:3)

num_segments = 0
in_segment = .false.

do i_point = 1, num_Point
    if (Points_Flag_Inside(i_point)) then
        if (.not. in_segment) then
            num_segments = num_segments + 1
            in_segment = .true.
        endif
    else
        if (in_segment) then
            in_segment = .false.
        endif
    endif
enddo

if (num_segments == 0) then
    return
    
elseif (num_segments > 1) then
    print *, "    Warn-2026020101: Warning in Tool_Fit_3D_Points_Taubin_Non_Closed.f90!"
    print *, '                     Number of inside segments can not > 1!'
    print *, '                     Number of inside segments detected:', num_segments
endif


i_start = 0
i_end = 0
found_start = .false.

do i_point = 1, num_Point
    if (Points_Flag_Inside(i_point)) then
        if (.not. found_start) then
            i_start = i_point
            found_start = .true.
        endif
        i_end = i_point
    else
        if (found_start) then
            exit
        endif
    endif
enddo

seg_len = i_end - i_start + 1

if (seg_len < 3) then
    return
endif


point_temp(i_start:i_end, 1:3) = Out_Points(i_start:i_end, 1:3)

L1_vector(i_start:i_end, 1:3) = ZR
L2_vector(i_start:i_end, 1:3) = ZR

do i_Loop = 1, n_Loops
    
    do i_point = i_start+1, i_end-1
        c_point = point_temp(i_point, 1:3)
        
        i_pre = i_point - 1
        i_nex = i_point + 1
        
        pre_point = point_temp(i_pre, 1:3)
        nex_point = point_temp(i_nex, 1:3)
        
        pre_L = sqrt(sum((pre_point - c_point)**2))
        nex_L = sqrt(sum((nex_point - c_point)**2))
        
        if (pre_L > Tol_15 .and. nex_L > Tol_15) then
            weight_pre = ONE / pre_L
            weight_nex = ONE / nex_L
            total_weight = weight_pre + weight_nex
            
            L1_vector(i_point, 1:3) = &
                (weight_pre*pre_point + weight_nex*nex_point) / total_weight - c_point
        else
            L1_vector(i_point, 1:3) = ZR
        endif
    enddo
    
    do i_point = i_start+1, i_end-1
        point_temp(i_point, 1:3) = point_temp(i_point, 1:3) + lambda * L1_vector(i_point, 1:3)
    enddo
    
    do i_point = i_start+1, i_end-1
        c_point = point_temp(i_point, 1:3)
        
        i_pre = i_point - 1
        i_nex = i_point + 1
        
        pre_point = point_temp(i_pre, 1:3)
        nex_point = point_temp(i_nex, 1:3)
        
        pre_L = sqrt(sum((pre_point - c_point)**2))
        nex_L = sqrt(sum((nex_point - c_point)**2))
        
        if (pre_L > Tol_15 .and. nex_L > Tol_15) then
            weight_pre = ONE / pre_L
            weight_nex = ONE / nex_L
            total_weight = weight_pre + weight_nex
            
            L2_vector(i_point, 1:3) = &
                (weight_pre*pre_point + weight_nex*nex_point) / total_weight - c_point
        else
            L2_vector(i_point, 1:3) = ZR
        endif
    enddo
    
    do i_point = i_start+1, i_end-1
        point_temp(i_point, 1:3) = point_temp(i_point, 1:3) + mu * L2_vector(i_point, 1:3)
    enddo
    
enddo

Out_Points(i_start+1:i_end-1, 1:3) = point_temp(i_start+1:i_end-1, 1:3)


return 
end subroutine Tool_Fit_3D_Points_Taubin_Non_Closed