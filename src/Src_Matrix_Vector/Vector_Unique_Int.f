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
 
      SUBROUTINE Vector_Unique_Int(n,n_Op_F,Vector,
     &                             Uniqued_Vec,Uniqued_n)   
C     Remove duplicate elements from the array, generate a new array, and store it in Uniqued_Vec, integer type
c     n: length of the original array
c     n_Op_F: The end operation position of the original array, that is, the Vector(1:n_Op_F) of the processed array
c     Uniqued_Vec: New Array
c     Uniqued_n: Length of usable data in the new array
      use Global_Float_Type
      implicit none
      integer,intent(in):: n,n_Op_F
      integer,intent(in):: Vector(n)
      integer,intent(out)::Uniqued_Vec(n)
      integer,intent(out)::Uniqued_n
      
      integer i,j,k
      Uniqued_Vec(1:n) = 0
      k = 1
      Uniqued_Vec(1) = Vector(1)
          
      outer: do i=2,n_Op_F
          do j=1,k
              if (Uniqued_Vec(j) .eq. Vector(i)) then
                  cycle outer
              end if
          end do
          k = k + 1
          Uniqued_Vec(k) = Vector(i)
      end do outer
             
      Uniqued_n   = k

      return
      END SUBROUTINE Vector_Unique_Int
    


