!-----------------------------------------------------------
! Brief: Removes outliers from a 1D data series using an n-sigma filter.
!
! Parameters:
!   Input:  Data_Points, num_points - data to denoise (modified in place)
!   Input:  Denoise_method - currently only method 1 (n-sigma) is supported
!   Input:  n_Sigma        - sigma multiplier for the rejection band
!   Input:  logical_Close  - closed-curve flag (unused by method 1)
!
! Notes:   Outliers above the mean are clipped to mean+SD and below to mean-SD
!          while non-outliers are left unchanged.
!-----------------------------------------------------------

subroutine Tool_Denoise_Data(Data_Points,num_points, Denoise_method, n_Sigma, logical_Close)

!     Used for denoising data. One-dimensional curve data.
!     Added on 2022-04-26.
!     Denoise = 1, standard deviation method.
!                  Ref:https://medium.com/geekculture/curve-smoothing-and-outliers-removal-using-c-d3bf6e2fbc78 
!                  or
!     \theory_documents\030 Curve Denoising (Deburring) and Smoothing Treatment-2022-04-26.pdf

!......................
! Variable Declaration
!......................
use Global_Float_Type
use Global_Common
implicit none
integer,intent(in)::num_points,Denoise_method
integer,intent(in)::n_Sigma
logical,intent(in)::logical_Close
real(kind=FT),intent(inout)::Data_Points(num_points)
real(kind=FT) :: in_Data(num_points)
integer :: i
real(kind=FT) c_average,SD,c_x,low_limit,up_limit
integer :: c_count

in_Data = Data_Points


if (Denoise_method==1) then
    c_average = sum(in_Data)/dble(num_points)
    call Tool_Standard_Deviation(in_Data,num_points,SD)


    c_count = 0 
    do i=1,num_points
        c_x = in_Data(i)
        low_limit = c_average - n_Sigma*SD
        up_limit  = c_average + n_Sigma*SD 
        if((c_x < low_limit) .or. (c_x > up_limit)) then
            if (c_x > c_average ) then 
                Data_Points(i) = c_average + SD
            else
                Data_Points(i) = c_average - SD
            endif
            c_count =  c_count +1
        endif
    enddo  

endif

return 
end subroutine Tool_Denoise_Data             
