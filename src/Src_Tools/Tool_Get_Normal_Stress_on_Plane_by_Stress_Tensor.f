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
 
      subroutine Tool_Get_Normal_Stress_on_Plane_by_Stress_Tensor(
     &  Stress,n,Normal_Stress)
c     Calculate the normal stress on any oblique section based on the stress.
c     Ref: \theory_documents\036 Calculation of Normal Stress on Any Oblique Section of 3D Elements_P5_2022-06-04.pdf
c     2022-06-04.
      use Global_Float_Type
      use Global_Common
      
      implicit none
      real(kind=FT),intent(in)::Stress(6),n(3)
      real(kind=FT),intent(out)::Normal_Stress
      real(kind=FT) Stress_Matrx(3,3),t_Vector(3)
      
      real(kind=FT) n_vector(3)
      
      
      n_vector = n
      
      call Vector_Normalize(3,n_vector)
      
      Stress_Matrx(1,1:3) = [Stress(1),Stress(4),Stress(6)]
      Stress_Matrx(2,1:3) = [Stress(4),Stress(2),Stress(5)]
      Stress_Matrx(3,1:3) = [Stress(6),Stress(5),Stress(3)]
      
      t_Vector = MATMUL(Stress_Matrx,n_vector)
      
      Normal_Stress = DOT_PRODUCT(t_Vector,n_vector)
      
      return 
      end SUBROUTINE Tool_Get_Normal_Stress_on_Plane_by_Stress_Tensor          
