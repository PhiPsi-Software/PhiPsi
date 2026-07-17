!-----------------------------------------------------------
! Brief: Smooth a 1D data series using a moving average filter.
!
! Parameters:
!   In/Out: Data_Points             - Series to be smoothed in place
!   Input:  num_points              - Number of samples
!   Input:  Smooth_method           - Smoothing method code (1 = MA)
!   Input:  Moving_Average_Method_n - Half-width of the MA window
!   Input:  logical_Close           - True for periodic wrap-around
!
! Notes:   Averages 2n+1 neighbouring values per sample; supports
!   closed curves via replicated end padding.
!-----------------------------------------------------------

subroutine Tool_Smooth_Data(Data_Points,num_points,Smooth_method, Moving_Average_Method_n,logical_Close)
!     Used for smoothing data. One-dimensional curve data.
!     Added on 2022-04-25.
!     Smooth_method = 1, moving average method, using the values from previous and subsequent moments
!     A total of 2n + 1 values are averaged (\theory_documents\029 Moving Average Smoothing Method-2022-04-25.docx).
!     

!......................
! Variable Declaration
!......................
use Global_Float_Type
use Global_Common
implicit none
integer,intent(in)::num_points,Smooth_method
integer,intent(in)::Moving_Average_Method_n
logical,intent(in)::logical_Close
real(kind=FT),intent(inout)::Data_Points(num_points)
real(kind=FT) in_Data(num_points),out_Data(num_points)
integer :: n
integer i,j,k
integer :: twice_n1
real(kind=FT),allocatable::tem_Data(:)
real(kind=FT) :: c_average

in_Data = Data_Points
n       = Moving_Average_Method_n


twice_n1 = 2*n+1

if(twice_n1 >= num_points) then
    print *,'    Error :: n is too big or num_points is too small!'
    print *,'             in Tool_Smooth_Data.f!'
    print *,'             n:',n
    print *,'             2*n+1:',2*n+1
    print *,'             num_points:',num_points
    call Warning_Message('S',Keywords_Blank)   
endif

if (Smooth_method==1) then
    allocate(tem_Data(-num_points:2*num_points))
    if (logical_Close .eqv. .False.) then
        tem_Data(-num_points:0) = ZR
        tem_Data(1:num_points)  = in_Data(1:num_points)
        tem_Data((num_points+1):2*num_points) = ZR
    elseif (logical_Close .eqv. .True.) then
        tem_Data(-num_points+1:0)  = in_Data(1:num_points)
        tem_Data(1:num_points)     = in_Data(1:num_points)
        tem_Data((num_points+1):2*num_points) = in_Data(1:num_points)
    endif
    do i=1,num_points
        c_average   = sum(tem_Data((i-n):(i+n)))/dble(twice_n1)
        out_Data(i) = c_average 
    enddo
    deallocate(tem_Data)
endif

Data_Points = out_Data

return 
end subroutine Tool_Smooth_Data             
