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
 
      SUBROUTINE D3_Flip_Crack_Mesh_Outline(i_C)
c     Outline of the ordering and other information of discrete nodes at the front edge vertices of flipped 3D crack surfaces.
c     2022-07-13. Ref: My PhiPsi Development Notebook V1-P41.

c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Crack_3D
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      integer,intent(in)::i_C
      integer i_out,Outline_num,c_Max_N_Node_3D
      integer,ALLOCATABLE::Temp(:,:) 
      
      c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
      
c     ---------------------
C     Main program section
c     ---------------------
      ALLOCATE(Temp(c_Max_N_Node_3D,4))
      Temp(1:c_Max_N_Node_3D,1:4)= 
     &            Crack3D_Meshed_Outline(i_C)%row(1:c_Max_N_Node_3D,1:4)
      
      Outline_num = Crack3D_Meshed_Outline_num(i_C)

      do i_out =1,Outline_num
        Crack3D_Meshed_Outline(i_C)%row(i_out,1)=
     &             Temp(Outline_num-i_out+1,2) 
        Crack3D_Meshed_Outline(i_C)%row(i_out,2)=
     &             Temp(Outline_num-i_out+1,1) 
        Crack3D_Meshed_Outline(i_C)%row(i_out,3)=
     &             Temp(Outline_num-i_out+1,3)      
      enddo

      deALLOCATE(Temp)

      RETURN
      END SUBROUTINE D3_Flip_Crack_Mesh_Outline
