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
 
      SUBROUTINE Convert_DISP_to_Cylinder_3D(n_FD,in_Disp,out_Disp)
c     Transform node displacements from the Cartesian coordinate system to the cylindrical coordinate system. 3D.
c     Ref: \theory_documents\024.2 Stress Tensor and Strain Tensor Transformation from Cartesian Coordinate System to Cylindrical Coordinate System --- 2.4-2021-09-10.html
c     2021-09-11.

c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Model
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      integer,intent(in)::n_FD
      real(kind=FT),intent(in)::  in_Disp(n_FD)
      real(kind=FT),intent(out):: out_Disp(n_FD)
      integer i_E,i_N
      real(kind=FT) c_Theta,c_Disp(1:3),sTheta,cTheta
      real(kind=FT) tem_Disp(1:3)
      
      print *,'    Converting disp from Cartesian CS cylindical CS...'
      do i_N = 1,Num_Node
            c_Disp(1)  = in_Disp(3*i_N-2)
            c_Disp(2)  = in_Disp(3*i_N-1)
            c_Disp(3)  = in_Disp(3*i_N-0)
            
            c_Theta = Theta_Cartesian_to_Cylinder_Node(i_N)
            
            sTheta      =  sin(c_Theta)
            cTheta      =  cos(c_Theta)
            tem_Disp(1) =  cTheta*c_Disp(1) + sTheta*c_Disp(2)
            tem_Disp(2) = -sTheta*c_Disp(1) + cTheta*c_Disp(2)
            tem_Disp(3) = c_Disp(3)
            
            out_Disp(3*i_N-2) = tem_Disp(1)
            out_Disp(3*i_N-1) = tem_Disp(2)
            out_Disp(3*i_N-0) = tem_Disp(3)

      end do
      RETURN
      END SUBROUTINE Convert_DISP_to_Cylinder_3D
