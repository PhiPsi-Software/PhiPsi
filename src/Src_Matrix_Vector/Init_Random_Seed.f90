!-----------------------------------------------------------
! Brief: Seed the Fortran intrinsic RANDOM_NUMBER from system clock.
!
! Notes:   Uses SYSTEM_CLOCK to derive a seed vector, ensuring a different
!          random sequence on every run.
!-----------------------------------------------------------

subroutine Init_Random_Seed()
! Subroutine init_random_seed() that seeds based on the system time
use Global_Float_Type
implicit none    
integer :: i, n, clock
integer, dimension(:), allocatable :: seed

call random_seed(size = n)
allocate(seed(n))

call system_clock(count=clock)

seed = clock + 37 * (/ (i - 1, i = 1, n) /)
call random_seed(put = seed)

deallocate(seed)
END subroutine Init_Random_Seed
