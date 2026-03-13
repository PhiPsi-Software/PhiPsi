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
 
      SUBROUTINE Matrix_Unique_Row_Int(m,n,m_Finish,Matrix,
     &                             Uniqued_Matrix,Uniqued_m,
     &                             Uni_Mat_Count)   
c     Operation target: from the 1st row to the m_Finish-th row of the matrix
C     Remove duplicate rows from the original matrix and store them in Uniqued_Matrix as integers.
c     m, n: Original matrix dimensions
c     Uniqued_Matrix: New Matrix
c     Uniqued_m: Number of usable data rows in the new matrix
c     Uni_Mat_Count: Number of times each row in the new matrix is repeated
      
      use Global_Float_Type
      implicit none
      integer,intent(in):: m,n,m_Finish
      integer,intent(in):: Matrix(m,n)
      integer,intent(out):: Uniqued_m
      integer,intent(out):: Uniqued_Matrix(m,n)
      integer,intent(out):: Uni_Mat_Count(m)
      
      integer tem(n)
      integer i,j,k
      
      Uniqued_Matrix(1:m,1:n) = 0

      Uni_Mat_Count(1:m)=1
      
      k = 1
      Uniqued_Matrix(1,:) = Matrix(1,:)
      
      outer: do i=2,m_Finish
          do j=1,k
              tem = Uniqued_Matrix(j,:) - Matrix(i,:)
              if ((maxval(tem).eq. 0).and.
     &            (minval(tem).eq. 0)) then
                  Uni_Mat_Count(j) = Uni_Mat_Count(j)+1
                  cycle outer
              end if
          end do
          k = k + 1
          Uniqued_Matrix(k,:) = Matrix(i,:)
      end do outer
                  
      Uniqued_m = k            

      return
      END SUBROUTINE Matrix_Unique_Row_Int
    


