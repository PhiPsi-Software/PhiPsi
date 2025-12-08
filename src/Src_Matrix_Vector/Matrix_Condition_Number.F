 
      SUBROUTINE Matrix_Condition_Number(n,Con_NNZ,
     &                        Cond_A,
     &                        Cond_IRN,
     &                        Cond_JCN,Condition_Num_Norm2)

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
      real(kind=FT) COND(2)
      integer LA,LIW,LW
      integer,allocatable::IRN(:), JCN(:), IW(:)
      integer,allocatable::K_Ptr(:)
      real(kind=FT),allocatable::A(:),W(:)
      
      LA=1000*Con_NNZ
      LIW=14*(Con_NNZ-1)+7
      LW=8*(Con_NNZ-1)
      
      ALLOCATE(A(LA))
      ALLOCATE(W(LW))
      ALLOCATE(IRN(LA))
      ALLOCATE(JCN(LA))
      ALLOCATE(IW(LIW))
     
      A(1:Con_NNZ)   = Cond_A(1:Con_NNZ)
      IRN(1:Con_NNZ) = Cond_IRN(1:Con_NNZ)
      JCN(1:Con_NNZ) = Cond_JCN(1:Con_NNZ)
      
      CALL MC75ID(ICNTL)
      
      CALL MC75AD(n,Con_NNZ,LA,A,IRN,JCN,COND,LIW,IW,LW,W,ICNTL,INFO)
      IF (INFO(1) .GE. 0) THEN
        Condition_Num_Norm2 =COND(1)
      ELSE
        print *, '    ERROR :: in Matrix_Condition_Number.f! ', INFO(1)
      ENDIF
      
      END SUBROUTINE Matrix_Condition_Number
      
