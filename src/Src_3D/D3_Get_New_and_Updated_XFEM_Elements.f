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
 
      SUBROUTINE D3_Get_New_and_Updated_XFEM_Elements(isub)
c     After confirming the enhanced nodes, obtain the list of newly added 3D XFEM elements and the list of 3D XFEM elements whose stiffness matrices need to be updated.
c     2022-06-24. NEWFTU2022062403.


C
C     integer num_FEM_Elem, num_XFEM_Elem
C     integer, ALLOCATABLE :: FEM_Elem_List(:)
C     integer, ALLOCATABLE :: XFEM_Elem_List(:)
C     integer, ALLOCATABLE :: Elem_XFEM_Flag(:)
C     integer, ALLOCATABLE :: Elem_Location(:,:)
C     integer, ALLOCATABLE:: Elem_New_XFEM_Flag(:)
C     integer, ALLOCATABLE :: Elem_Update_XFEM_Flag(:)


c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Crack_3D
      use Global_XFEM_Elements
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      integer,intent(in)::isub
      integer i_E
      
c     ---------------------------------
C     The first step is not necessary.
c     ---------------------------------
      if (isub<=1)then
          return
      endif

c     ------------------------
C     Variable Initialization
c     ------------------------
      if (allocated(Elem_New_XFEM_Flag))then
          DEALLOCATE(Elem_New_XFEM_Flag)
      endif
      allocate(Elem_New_XFEM_Flag(num_Elem))
      Elem_New_XFEM_Flag(1:num_Elem) = 0
      
      if (allocated(Elem_Update_XFEM_Flag))then
          DEALLOCATE(Elem_Update_XFEM_Flag)
      endif
      allocate(Elem_Update_XFEM_Flag(num_Elem))
      Elem_Update_XFEM_Flag(1:num_Elem) = 0

c     ----------
C     ext steps
c     ----------
      ! The Determine_Enriched_Nodes_3D subroutine saved Elem_Location_old and Elem_XFEM_Flag_Old.
      print *,'    Getting New and to-be-Updated XFEM Elements...'  
      
      do i_E = 1,Num_Elem
          ! If the previous step was not an XFEM element and this step is an XFEM element, then it is marked
          ! as a new XFEM element.
          if(Elem_XFEM_Flag_Old(i_E)==0 .and.Elem_XFEM_Flag(i_E)==1)then
              Elem_New_XFEM_Flag(i_E) = 1
          endif
      enddo
      

      
      RETURN
      END SUBROUTINE D3_Get_New_and_Updated_XFEM_Elements
