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
 
      SUBROUTINE D3_Find_Crack_Mesh_Outline(i_C,num_of_ele,num_of_nodes,
     &                               c_Crack3D_Meshed_Ele)
c     Find the outer boundary of the crack surface mesh and store it in the global variable Crack3D_Meshed_Outline(Max_Num_Cr,1000,3) for each crack.
c     The number of boundary lines is stored in the global variable Crack3D_Meshed_Outline_num(Max_Num_Cr).
c     
c     Firstly written by Fang Shi on 2019-03-28.
C     integer Crack3D_Meshed_Outline(Max_Num_Cr_3D, Max_N_Node_3D, 4)
                                                                        ! Data 1 is the first point on the boundary line of the crack front edge
                                                                        ! Data 2 is the second point on the leading edge boundary line of the crack
                                                                        ! Data 3 corresponds to the discrete fracture element number
                                                                        ! Data 4 is used to mark whether the two points of the boundary line are allowed to extend,
                                                                        ! extending in very small steps (2021-08-20)
c     integer Crack3D_Meshed_Outline_num(Max_Num_Cr)


      ! The following variables Phi and Psi were not used
c     real(kind=FT) Crack3D_Meshed_Outline_Vertex(Max_Num_Cr_3D, Max_N_Node_3D, 3)
c     integer Crack3D_Meshed_Outline_Vertex_Ele_num(Max_Num_Cr_3D, Max_N_Node_3D)

c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Crack_3D
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      integer,intent(in)::i_C,num_of_ele,num_of_nodes
      integer,intent(in)::c_Crack3D_Meshed_Ele(num_of_ele,3)
      integer tem1(3*num_of_ele,2)
      integer tem2(3*num_of_ele)
      integer i,j,all_num_Outline,Out_num_Outline
      integer,ALLOCATABLE::All_Outline(:,:)
      integer,ALLOCATABLE::Temp_Outline(:,:) 
      integer location
      logical Yes_In_1,Yes_In_2


c     ---------------------
C     Main program section
c     ---------------------
      ! The three sides of the element are stored in the temporary variable tem1
      tem1(1:num_of_ele,1)   = c_Crack3D_Meshed_Ele(1:num_of_ele,1)
      tem1(1:num_of_ele,2)   = c_Crack3D_Meshed_Ele(1:num_of_ele,2)
      tem1(num_of_ele+1:2*num_of_ele,1)   = 
     &                              c_Crack3D_Meshed_Ele(1:num_of_ele,2)
      tem1(num_of_ele+1:2*num_of_ele,2)   = 
     &                              c_Crack3D_Meshed_Ele(1:num_of_ele,3)
      tem1(2*num_of_ele+1:3*num_of_ele,1) = 
     &                              c_Crack3D_Meshed_Ele(1:num_of_ele,3)
      tem1(2*num_of_ele+1:3*num_of_ele,2) = 
     &                              c_Crack3D_Meshed_Ele(1:num_of_ele,1)
      
      ! Sort
      call Matrix_Sort_Int(3*num_of_ele,2,tem1)
      
      ! Count the occurrences of each row in the matrix and store them in tem2
      call Matrix_Count_Row_Int(3*num_of_ele,2,
     &                          tem1,tem2,all_num_Outline)
      
      ! Extract edges that appear only once
      ALLOCATE(All_Outline(all_num_Outline,2))
      ALLOCATE(Temp_Outline(all_num_Outline,2))
      j=0
      do i=1,3*num_of_ele
          if (tem2(i).eq.1)then
              j=j+1
              All_Outline(j,:) = tem1(i,:)
          end if
      end do
      
      ! Connected end to end, containing only the outer boundary
      call Tool_Sort_by_End_to_End(all_num_Outline,all_num_Outline,
     &                             All_Outline,Temp_Outline,
     &                             Out_num_Outline)
      
      Crack3D_Meshed_Outline_num(i_C) = Out_num_Outline
      
      ! Find the elements corresponding to the outer boundary of the crack
      do i=1,Out_num_Outline
          Crack3D_Meshed_Outline(i_C)%row(i,1:2)= Temp_Outline(i,1:2)
          do j=1,num_of_ele
              call Vector_Location_Int(3,c_Crack3D_Meshed_Ele(j,1:3),
     &                            Temp_Outline(i,1),location,Yes_In_1) 
              call Vector_Location_Int(3,c_Crack3D_Meshed_Ele(j,1:3),
     &                            Temp_Outline(i,2),location,Yes_In_2)      
              ! If the two nodes on the current outer boundary belong to a certain element, then that element is
              ! the one we are looking for.
              if (Yes_In_1 .and. Yes_In_2)then
                  Crack3D_Meshed_Outline(i_C)%row(i,3)  = j
                  exit
              end if
          end do
      end do      
      if (allocated(All_Outline))deallocate(All_Outline)
      if (allocated(Temp_Outline))deallocate(Temp_Outline)
 
      
      
      RETURN
      END SUBROUTINE D3_Find_Crack_Mesh_Outline
