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
 
      SUBROUTINE Boundary_Cond_HF(iter,
     &                            freeDOF_HF,num_free_CalP,
     &                            Local_freeDOF_HF) 
c     Boundary condition (the water pressure at the crack tip is 0, so the water pressure degree of freedom at the crack tip is constrained)

      ! num_Tol_CalP: Total number of calculation points
      ! num_free_CalP: Number of calculation points after removing constraints
      
c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Crack
      use Global_Crack_Common
      use Global_Material
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      
      integer,intent(in)::iter
      integer,intent(out)::freeDOF_HF(num_Tol_CalP_Water),num_free_CalP
      integer,intent(out)::Local_freeDOF_HF(num_Tol_CalP_Water)
      
      integer i_CalP,i_C,i_fd
      integer fixedDOF(num_Tol_CalP_Water),num_Count
      
      print *,'    Dealing with boundary condition of HF...'   
      
      !Initialize vector of fixed DOFs.
      freeDOF_HF(1:num_Tol_CalP_Water) = 0
      fixedDOF(1:num_Tol_CalP_Water)   = 0
      
      num_Count = 0
      do i_C = 1,num_Crack
          ! If the current crack is driven by fracturing fluid, then:
          if (Cracks_HF_State(i_C) == 1) then  
              ! Current Crack Calculation Point Cycle
              do i_CalP=1,Cracks_CalP_Num(i_C)
                  num_Count = num_Count + 1
                  if (Cracks_CalP_Type(i_C,i_CalP,1) == 5) then
                                                                 ! See Cal_HF_Crack_Points_Info_Linear(isub) for details
                      ! If viscosity dominates and the keyword for setting crack tip water pressure to 0 is enabled, then
                      ! set the crack tip water pressure to 0.
                      if(Global_K_m < ONE .and.
     &                   Key_Tip_Pres_Zero==1) then
                          fixedDOF(num_Count) =  num_Count  
                      endif
                  end if

              end do
          end if
      end do
      
      num_free_CalP = 0              
      do i_fd =1,num_Tol_CalP_Water
          if ( .not.(ANY( fixedDOF(1:num_Count) .eq. i_fd)) ) then
              num_free_CalP = num_free_CalP +1
              freeDOF_HF(num_free_CalP) = Total_FD+i_fd
              Local_freeDOF_HF(num_free_CalP) = i_fd
          end if
      end do
      
      RETURN
      END SUBROUTINE Boundary_Cond_HF
