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
 
      subroutine Tool_Get_Normal_Vector_of_Points_3D(
     &                 In_Points,num_Point,n_Vector)
      ! Calculate the outward normal vector of the plane on which the 3D point In_Points lies.
      ! For information on obtaining an out-of-plane normal vector, refer to the following materials
      !Theory Ref-https://www.freesion.com/article/7157945105/
      !            or
      ! \theory_documents\026 Optimal Spatial Circle Fitting for 3D Discrete Points and
      ! Implementation-2021-10-30.pdf
      !2021-10-30.
      
      use Global_Float_Type
      
      implicit none
      integer,intent(in)::num_Point
      real(kind=FT),intent(in)::In_Points(num_Point,3)
      real(kind=FT),intent(out)::n_Vector(3)
      real(kind=FT) M(num_Point,3),MTM_Inv(3,3),MTM(3,3),L(num_Point)
      real(kind=FT) norm_n_Vector
      
      M = In_Points
      MTM = matmul(TRANSPOSE(M),M)
      call Matrix_Inverse_3x3(MTM,MTM_Inv)  
      L(1:num_Point) = ONE
      n_Vector =  matmul(matmul(MTM_Inv,TRANSPOSE(M)),L)
      call Vector_Norm2(3,n_Vector,norm_n_Vector) 
      n_Vector =  n_Vector /norm_n_Vector
      
      
      return 
      end SUBROUTINE Tool_Get_Normal_Vector_of_Points_3D                        
