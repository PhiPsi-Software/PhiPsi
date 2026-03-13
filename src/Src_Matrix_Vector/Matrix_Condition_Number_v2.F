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
 
      SUBROUTINE Matrix_Condition_Number_v2(norm_Type,n,Matrix,
     &                                      Condition_Num_Norm)
c     Calculate the condition number of a matrix. Suitable for calculating the condition number of smaller matrices.
c     Fist written on 2022-07-10 by Fang Shi.
c     Test passed.
c     Note: If the matrix does not have an inverse, the output may be NAN.

      use Global_Float_Type
      
      implicit none
      integer,intent(in)::n,norm_Type
      real(kind=FT),intent(in)::Matrix(n,n)
      real(kind=FT),intent(out)::Condition_Num_Norm
      
      real(kind=FT) Norm_2_K,Norm_2_K_inv
      real(kind=FT) Norm_1_K,Norm_1_K_inv
      real(kind=FT) Matrix_inv(n,n)
      
      
      call Matrix_Inverse(Matrix,Matrix_inv,n)   
      
      if (norm_Type==2) then
          call Matrix_Norm2(n,Matrix,Norm_2_K) 

          call Matrix_Norm2(n,Matrix_inv,Norm_2_K_inv) 
          
          Condition_Num_Norm =Norm_2_K* Norm_2_K_inv
      endif

      if (norm_Type==1) then
          call Matrix_Norm1(n,Matrix,Norm_1_K) 

          call Matrix_Norm1(n,Matrix_inv,Norm_1_K_inv) 
          
          Condition_Num_Norm =Norm_1_K* Norm_1_K_inv
      endif      
      
      return
      END SUBROUTINE Matrix_Condition_Number_v2
      
