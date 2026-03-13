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
 
      SUBROUTINE Matrix_Inverse(a,c,n)   
C     Find the inverse matrix of square matrix a and output it to c.
       
      !============================================================
      ! Inverse matrix
      ! Method: Based on Doolittle LU factorization for Ax=b
      ! Alex G. December 2009
      !-----------------------------------------------------------
      ! input ...
      ! a(n,n) - array of coefficients for matrix A
      ! n      - dimension
      ! output ...
      ! c(n,n) - inverse matrix of A
      ! comments ...
      ! the original matrix a(n,n) will be destroyed 
      ! during the calculation
      !===========================================================
      use Global_Float_Type
      implicit none 
      integer,intent(in)::n
      real(kind=FT),intent(in):: a(n,n)
      real(kind=FT),intent(out):: c(n,n)
      real(kind=FT) L(n,n), U(n,n), b(n), d(n), x(n)
      real(kind=FT) coeff,tem_a(n,n)
      integer i, j, k
      
      tem_a = a
      L=ZR
      U=ZR
      b=ZR

      do k=1, n-1
          do i=k+1,n
              coeff=tem_a(i,k)/tem_a(k,k)
              L(i,k) = coeff
              do j=k+1,n
                  tem_a(i,j) = tem_a(i,j)-coeff*tem_a(k,j)          
              end do
          end do
      end do

      do i=1,n
        L(i,i) = ONE
      end do
      do j=1,n
          do i=1,j
              U(i,j) = tem_a(i,j)
          end do
      end do

      do k=1,n
          b(k)=ONE
          d(1) = b(1)
          do i=2,n
              d(i)=b(i)
              do j=1,i-1
                 d(i) = d(i) - L(i,j)*d(j)
              end do
          end do
        x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
              x(i) = d(i)
              do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
              end do
          x(i) = x(i)/u(i,i)
          end do
          do i=1,n
              c(i,k) = x(i)
          end do
        b(k)=ZR
      end do
        
      return
      END SUBROUTINE Matrix_Inverse
    


