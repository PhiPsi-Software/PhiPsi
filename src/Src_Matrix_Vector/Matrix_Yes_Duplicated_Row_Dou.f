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
 
      SUBROUTINE Matrix_Yes_Duplicated_Row_Dou(m,n,Matrix,Yes_Dup,
     &                                         num_Dup)   
c     Determine whether the matrix has duplicate rows.
c     Written based on Matrix_Unique_Row_Dou.f.
c     Added on 2022-04-23.

      use Global_Float_Type
      implicit none
      integer,intent(in):: m,n
      real(kind=FT),intent(in):: Matrix(m,n)
      logical,intent(out):: Yes_Dup
      integer,intent(out):: num_Dup
      integer Uniqued_m,Uni_Mat_Count(m)
      real(kind=FT) Uniqued_Matrix(m,n),tem(n)
      integer i,j,k
      
      Uniqued_Matrix(1:m,1:n) = ZR
      Yes_Dup = .False.
      Uni_Mat_Count(1:m)=1
      
      k = 1
      Uniqued_Matrix(1,:) = Matrix(1,:)
      
      outer: do i=2,m
          do j=1,k
              tem = Uniqued_Matrix(j,:) - Matrix(i,:)
              if ((abs(maxval(tem))<=Tol_11) .and.
     &            (abs(minval(tem))<=Tol_11)) then
                  Uni_Mat_Count(j) = Uni_Mat_Count(j)+1
                  Yes_Dup = .True.
                  cycle outer
              end if
          end do
          k = k + 1
          Uniqued_Matrix(k,:) = Matrix(i,:)
      end do outer
                  
      Uniqued_m = k    
      
      num_Dup = m-Uniqued_m

      return
      END SUBROUTINE Matrix_Yes_Duplicated_Row_Dou
    


