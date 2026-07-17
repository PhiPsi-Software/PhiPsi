!-----------------------------------------------------------
! Brief: Count duplicate rows of a 2-column integer matrix in parallel.
!
! Parameters:
!   Input:  m           - number of rows in Matrix
!           Matrix      - integer matrix with two columns (m,2)
!           key_process - whether to print progress (unused)
!   Output: Count_Row   - per-row duplicate count
!           Num_Once    - number of rows that appear only once
!
! Notes:   Parallelized with OpenMP over the rows; each row is
!          compared to the entire matrix via masked count.
!-----------------------------------------------------------

subroutine Matrix_Count_Row_Int_m_x_2(m,Matrix, Count_Row,Num_Once, key_process)
! Count the occurrences of each row in the matrix (m,2) and store them in Count_Row.
! 2020-03-27: Add openMP, disable key_process.    

use Global_Float_Type    
use CMD_Progress
use Global_Common
use omp_lib
implicit none

integer,intent(in)::m
integer,intent(in)::key_process
integer,intent(in)::Matrix(m,2)
integer,intent(out)::Count_Row(m)
integer,intent(out)::Num_Once

integer :: i

type(CLS_CMD_Progress)::Progress

Count_Row(1:m) = 0
!$OMP PARALLEL do DEFAULT(SHARED) private(i) SCHEDULE(static)  
do i=1,m    
Count_Row(i) =count((Matrix(1:m,1) .eq. Matrix(i,1)) .and. (Matrix(1:m,2) .EQ. Matrix(i,2)))
end do      
!$omp end parallel do 

Num_Once= 0
Num_Once=count(Count_Row .eq. 1)

return
END subroutine Matrix_Count_Row_Int_m_x_2
    


