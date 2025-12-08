 
      SUBROUTINE Init_Random_Seed()
      use Global_Float_Type
      implicit none    
      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed
   
      call random_seed(size = n)
      ALLOCATE(seed(n))
   
      call system_clock(count=clock)
   
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)
   
      DEALLOCATE(seed)
      END SUBROUTINE Init_Random_Seed