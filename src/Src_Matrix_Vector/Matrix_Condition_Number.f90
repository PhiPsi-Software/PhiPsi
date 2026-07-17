!-----------------------------------------------------------
! Brief: Compute the 1-norm condition number of a sparse matrix.
!
! Parameters:
!   Input:  n          - order of the square matrix
!           Con_NNZ    - number of nonzeros in the matrix
!           Cond_A     - values of nonzero entries (size Con_NNZ)
!           Cond_IRN   - row indices of nonzero entries
!           Cond_JCN   - column indices of nonzero entries
!   Output: Condition_Num_Norm2 - 1-norm condition number
!
! Notes:   Wraps the HSL MC75AD sparse direct condition-number
!          estimator; on negative INFO it prints an error and leaves
!          the output unchanged.
!-----------------------------------------------------------

subroutine Matrix_Condition_Number(n,Con_NNZ, Cond_A, Cond_IRN, Cond_JCN,Condition_Num_Norm2)
!     Calculate the condition number of matrix.
!     Fist written on 2020-02-10 by Fang Shi.
!     http://www.hsl.rl.ac.uk/catalogue/mc75.html

use Global_Float_Type
use Global_Common
use Global_Model

implicit none
integer,intent(in)::n,Con_NNZ
integer,intent(in)::Cond_IRN(1:Con_NNZ)
integer,intent(in)::Cond_JCN(1:Con_NNZ)
real(kind=FT),intent(in)::Cond_A(1:Con_NNZ)
real(kind=FT),intent(out)::Condition_Num_Norm2
integer INFO(5),ICNTL(5)
real(kind=FT) :: COND(2)
integer LA,LIW,LW
integer,allocatable::IRN(:), JCN(:), IW(:)
integer,allocatable::K_Ptr(:)
real(kind=FT),allocatable::A(:),W(:)

LA=1000*Con_NNZ
LIW=14*(Con_NNZ-1)+7
LW=8*(Con_NNZ-1)

allocate(A(LA))
allocate(W(LW))
allocate(IRN(LA))
allocate(JCN(LA))
allocate(IW(LIW))

A(1:Con_NNZ)   = Cond_A(1:Con_NNZ)
IRN(1:Con_NNZ) = Cond_IRN(1:Con_NNZ)
JCN(1:Con_NNZ) = Cond_JCN(1:Con_NNZ)

call MC75ID(ICNTL)

call MC75AD(n,Con_NNZ,LA,A,IRN,JCN,COND,LIW,IW,LW,W,ICNTL,INFO)
if (INFO(1) .ge. 0) then
    Condition_Num_Norm2 =COND(1)
else
    print *, '    ERROR :: in Matrix_Condition_Number.f! ', INFO(1)
endif

END subroutine Matrix_Condition_Number

