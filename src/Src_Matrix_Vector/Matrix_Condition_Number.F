!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
      SUBROUTINE Matrix_Condition_Number(n,Con_NNZ,
     &                        Cond_A,
     &                        Cond_IRN,
     &                        Cond_JCN,Condition_Num_Norm2)
c     Calculate the condition number of matrix.
c     Fist written on 2020-02-10 by Fang Shi.
c     http://www.hsl.rl.ac.uk/catalogue/mc75.html
c     Use Condition_Number_HSL_ddeps.f + Condition_Number_HSL_mc75d.f

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
      
