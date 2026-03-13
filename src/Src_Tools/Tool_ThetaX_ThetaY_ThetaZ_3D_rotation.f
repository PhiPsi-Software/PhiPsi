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
 
      subroutine Tool_ThetaX_ThetaY_ThetaZ_3D_rotation(i_V,     
     &                    Vector_Origi_1,
     &                    Vector_Origi_2,
     &                    Vector_Origi_3,
     &                    Vector_Crack_1,
     &                    Vector_Crack_2,
     &                    Vector_Crack_3,
     &                    ThetaX,ThetaY,ThetaZ,T_Matrix)
C     After the 3D spatial coordinate system is rotated, its vectors also transform accordingly. 
c     For example: from Vector_Origi_1(1:3) to Vector_Crack_1(1:3)
c     Vector_Origi_2(1:3) to Vector_Crack_2(1:3)
c     From Vector_Origi_3(1:3) to Vector_Crack_3(1:3)
c     The purpose of this subroutine is to determine the corresponding rotation matrix 
c     T_Matrix based on known transformations.

      !......................
      ! Variable Declaration
      !......................
      use Global_Float_Type   
      use Global_Crack_3D
      
      implicit none
      integer,intent(in)::i_V
      real(kind=FT),intent(in)::     
     &         Vector_Origi_1(3),Vector_Origi_2(3),Vector_Origi_3(3),
     &         Vector_Crack_1(3),Vector_Crack_2(3),Vector_Crack_3(3)
      real(kind=FT),intent(out):: ThetaX,ThetaY,ThetaZ,T_Matrix(3,3)
      real(kind=FT) c_T(3,3)      
      T_Matrix(1:3,1:3) = ZR
      ThetaX = ZR
      ThetaY = ZR
      ThetaZ = ZR

      
      call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_1,
     &                                          Vector_Crack_1,c_T(1,1))
      call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_1,
     &                                          Vector_Crack_2,c_T(1,2)) 
      call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_1,
     &                                          Vector_Crack_3,c_T(1,3))     
      call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_2,
     &                                          Vector_Crack_1,c_T(2,1))
      call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_2,
     &                                          Vector_Crack_2,c_T(2,2)) 
      call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_2,
     &                                          Vector_Crack_3,c_T(2,3))     
      call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_3,
     &                                          Vector_Crack_1,c_T(3,1))
      call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_3,
     &                                          Vector_Crack_2,c_T(3,2)) 
      call Tool_Cos_Value_of_Vectors_a_and_b_3D(Vector_Origi_3,
     &                                          Vector_Crack_3,c_T(3,3))  
      T_Matrix = c_T
             
      
      return 
      end SUBROUTINE Tool_ThetaX_ThetaY_ThetaZ_3D_rotation                   
