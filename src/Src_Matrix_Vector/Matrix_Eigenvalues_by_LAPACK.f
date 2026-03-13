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
 
      SUBROUTINE Matrix_Eigenvalues_by_LAPACK(n,Matrix,WR,WI)
C     The eigenvalues of a square matrix can be calculated, and the eigenvectors can also be calculated (can be calculated but not calculated).
c     Call the LAPACK library DGEEV subroutine to calculate eigenvalues and eigenvectors
c
c      (1) Official online documentation for the DGEEV subroutine in the LAPACK library:
C         http://www.netlib.org/lapack/explore-html/d9/d28/dgeev_8f.html
C         DGEEV computes for an N-by-N real nonsymmetric matrix A, the
C         eigenvalues and, optionally, the left and/or right eigenvectors.
C         The right eigenvector v(j) of A satisfies:
C                       A * v(j) = lambda(j) * v(j),where lambda(j) is its eigenvalue.
c     (2) Example of using the DGEEV subroutine from the LAPACK library:
C         https://www.nag.com/lapack-ex/node87.html
C         CALL DGEEV('No left vectors','Vectors (right)',N,A,LDA,WR,WI,
C                     DUMMY,1,VR,LDVR,WORK,LWORK,INFO)
      use Global_Float_Type
      implicit none

      integer,intent(in)::n
      real(kind=FT),intent(in)::Matrix(n,n)
      real(kind=FT),intent(out)::WR(n),WI(n)

      integer out_INFO,LWORK
      real(kind=FT) VR(n,n)
      integer LDVR
      real(kind=FT) DUMMY(1,N)
      integer lda
      real(kind=FT),ALLOCATABLE::Lapack_work(:)
      
      
      lda     = n
      LDVR    = n
      WR(1:n) = ZR
      WI(1:n) = ZR
      LWORK = 100*n
      ALLOCATE(Lapack_work(LWORK))
      
      call dgeev ('N','N', n,     Matrix, lda,
     &             WR, WI, DUMMY, 1, VR,
     &             LDVR,Lapack_work, LWORK, out_INFO)
     
      DEALLOCATE(Lapack_work)
      
      
      return
      END SUBROUTINE Matrix_Eigenvalues_by_LAPACK



