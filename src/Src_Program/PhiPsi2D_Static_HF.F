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
 
      SUBROUTINE PhiPsi2D_Static_HF
      ! Hydraulic Fracturing Static Analysis Main Program
      ! Supports Newton-Raphson iteration with a linear search backtracking algorithm to improve
      ! convergence stability (see the description in my notes V3-P146)
      ! Note: (1) Supports crack closure detection
      !       (2) Support T-shaped cracks
      !       (3) Support multiple break steps
      !       (4) Proppants are not supported for now, but will be added in the future
      !       (5) The more advanced iterator based on this iterator is Iterator No. 9.
      !       (6) Support mm-ton-s unit system
      !       (7) Support stress (treatment strategy, linear superposition)
      !       (8) The most powerful solver
      ! Some important variables related to degrees of freedom:
                  ! Total_FD:           Total number of degrees of freedom of the solid
                  ! num_FreeD:          Total number of active solid degrees of freedom
                  ! num_Tol_CalP_Water: Number of fluid degrees of freedom
                  ! num_free_CalP:      Number of degrees of freedom of the active fluid
                  ! num_ALLDOF = Total_FD   num_Tol_CalP_Water: Total number of degrees of freedom for fluids and solids
                  ! num_Total_FD = num_FreeD num_free_CalP:     Total number of active fluid and solid degrees of freedom

c     ----------------------------
c     Read Public Variable Module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Filename
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Crack
      use Global_Crack_Common
      use Global_DISP
      use Global_HF
      use Global_Stress
      use Global_POST
      use Global_Contact
      use Global_Material

c     ---------------------------------------------------------------------------------
c     Read subroutine interface module (activate compiler parameter consistency check)
c     ---------------------------------------------------------------------------------
      use Global_Inter_Cal_HF_Initial_Pressure
      use Global_Inter_Cal_HF_Resid
      use Global_Inter_Cal_HF_Jacobian_NR
      use Global_Inter_Matrix_Solve_LSOE
      use Global_Inter_Cal_HF_Line_Searching
      use Global_Inter_Determine_Contact_State_by_Iteration
      use Global_Inter_Assemble_Stiffness_Matrix_SPARS_XFEM
      use Global_Inter_Matrix_Solve_LSOE_Sparse
      use Global_Inter_Cal_HF_Resid_SPARS
      use Global_Inter_Cal_HF_Jacobian_NR_SPARS

c     --------------------------
c     Variable Type Declaration
c     --------------------------
      implicit none

      integer ifra,iter
      real(kind=FT) Lambda
      integer,ALLOCATABLE::freeDOF(:),freeDOF_HF(:),Total_freeDOF(:)
      integer,ALLOCATABLE::Local_freeDOF_HF(:)
      real(kind=FT),ALLOCATABLE:: globalK(:,:),F_U(:),F(:),F_U_Ini(:)
      real(kind=FT),ALLOCATABLE:: ori_globalK(:,:)
      real(kind=FT),ALLOCATABLE:: F_CT(:)
      real(kind=FT),ALLOCATABLE:: Coupled_Q(:,:)
      real(kind=FT),ALLOCATABLE:: H(:,:)
      real(kind=FT),ALLOCATABLE:: c_S(:)
      real(kind=FT),ALLOCATABLE:: Last_DISP(:)
      real(kind=FT),ALLOCATABLE:: Last_Last_DISP(:)
      real(kind=FT),ALLOCATABLE:: cLast_Last_DISP(:)
      real(kind=FT),ALLOCATABLE:: Real_DISP(:),POST_DISP(:)
      real(kind=FT),ALLOCATABLE:: R(:)
      real(kind=FT),ALLOCATABLE:: NR_Deri(:,:)
      real(kind=FT),ALLOCATABLE:: Last_CalP_Pres(:)
      real(kind=FT),ALLOCATABLE:: CalP_Pres(:)
      real(kind=FT),ALLOCATABLE:: CalP_Pres_Temp(:)
      real(kind=FT),ALLOCATABLE:: DISP_Temp(:)
      real(kind=FT),ALLOCATABLE:: Added_DISP(:)
      real(kind=FT),ALLOCATABLE:: delta_x(:)
      real(kind=FT),ALLOCATABLE:: Contact_x(:)
      real(kind=FT),ALLOCATABLE:: Initial_DISP(:)
      real(kind=FT),ALLOCATABLE:: Contact_DISP(:)
      real(kind=FT),ALLOCATABLE:: x(:)
      integer num_FreeD
      integer i_C,Counter_Iter
      real(kind=FT) delta_Time,total_Time
      real(kind=FT) Last_Cr_CalP_Aper(Max_Num_Cr,Max_Num_Cr_CalP),
     &    Last_Cr_CalP_Pres(Max_Num_Cr,Max_Num_Cr_CalP),
     &    Last_Cr_CalP_Aper_ifra(Max_Num_Cr,Max_Num_Cr_CalP),
     &    Last_Cr_CalP_Pres_ifra(Max_Num_Cr, Max_Num_Cr_CalP),
     &    Temp_Cr_CalP_Aper(Max_Num_Cr,Max_Num_Cr_CalP)
      integer Last_Cracks_CalP_Num_ifra(Max_Num_Cr)
      logical Yes_Conv
      integer num_free_CalP,
     &        num_Count,i_Calp,
     &        num_ALlDOF
      integer num_Total_FD,Total_Num_G_P
      logical Yes_Growth(Max_Num_Cr,2)
      logical Yes_Last_Growth
      integer num_CalP
      integer Last_Total_Freedom_ifra
      real(kind=FT) Max_vm_Node,Max_vm_Gauss
      real(kind=FT) Old_total_Time(Num_Frac)
      real(kind=FT) Temp_total_Time
      integer date_time(8)
      integer(LIT) curent_time
      character*10  current_data
      real(kind=FT) ifra_CPU_time(Num_Frac)
      real(kind=FT) ifra_LS_CPU_time(Num_Frac)
      real(kind=FT) ifra_CT_CPU_time(Num_Frac)
      integer(LIT) c_LS_time_s,c_LS_time_f
      integer(LIT) c_CT_time_s,c_CT_time_f
      real(kind=FT) Total_LS_CPU_time,Total_CT_CPU_time
      real(kind=FT) ifra_HF_time(Num_Frac)
      integer ifra_num_iter(Num_Frac)
      real(kind=FT) ifra_Max_Press(Num_Frac)
      real(kind=FT) ifra_Length(Num_Frac)
      ! Used to determine the contact condition of the fracture surface
      real(kind=FT) trans_factor
      integer Num_Line_Search
      integer t_num_CalP
      real(kind=FT),ALLOCATABLE:: F_U_Backup(:)
      real(kind=FT),ALLOCATABLE::Cracks_CalP_Tan_Aper(:,:)

c     ---------------------------
c     Formatted output statement
c     ---------------------------
 1001 FORMAT(' >>  MNR iteration ',I3,' of fracturing step'
     &                     ,I3,'(',I3,') started:')
 1023 FORMAT('     Sum of globalK: ',E15.8)
 1041 FORMAT(5X,'Max Von-mises stress of all nodes:   ',F18.6,' MPa')
 1042 FORMAT(5X,'Max Von-mises stress of Gauss points:   ',F18.6,' MPa')
 1131 FORMAT(5X,'KI and KII of crack ',I5,' tip 1 are ', F18.6,
     &                 ' and ',F12.6,' MPa.m^(1/2)')
 1132 FORMAT(5X,'KI and KII of crack ',I5,' tip 2 are ', F18.6,
     &                 ' and ',F12.6,' MPa.m^(1/2)')
 1201 FORMAT(5X,'Computing counter: ',I6)
 5001 FORMAT('     Elapsed CPU time：',F7.3,' mins')



c     ----------------------------------
c     Variable Initialization and Setup
c     ----------------------------------
      ifra_CPU_time(1:Num_Frac) =ZR
      ifra_LS_CPU_time(1:Num_Frac) = ZR
      ifra_CT_CPU_time(1:Num_Frac) = ZR
      ifra_HF_time(1:Num_Frac)     = ZR
      ifra_num_iter(1:Num_Frac) =0
      ifra_Max_Press(1:Num_Frac)=0
      ifra_Length(1:Num_Frac)   =0
      Num_Line_Search = 0
      Total_LS_CPU_time            = ZR
      Total_CT_CPU_time            = ZR
      ! Get Start Time
      call Tool_Get_Current_Time(current_data,date_time,S_time)

c     --------------------------------------------------
c     Calculate the material matrix D for each material
c     --------------------------------------------------
      call Get_Material_Matrix

c     -------------------------------------------------------------
c     Handling natural fractures related: Generate initial natural 
c     fractures (if needed), and based on the natural fractures
c     Perform relevant data processing based on the type
c     -------------------------------------------------------------
      call Natural_Fractures_Related_Preparation
      
c     --------------------------------------------------------------------------------------------
c     Some logical checks and data corrections:
c     (1) Friction cracks cannot contain water (only applicable to third-type natural fractures).
c     (2) The current version of the program can only have one initial water-driven fracture.
c     --------------------------------------------------------------------------------------------
      call HF_Logic_Checking

c     -----------------------------------------------------------------------------------------
c     Calculate the traditional degree-of-freedom displacement field under 
c     the action of in-situ stress (far-field stress), and calculate the in-situ stress level.
c     -----------------------------------------------------------------------------------------
      if(Key_InSitu_Strategy /= 0)then
          call Cal_InSitu_Stress
      endif

c     --------------------
c     Fracture Step Cycle
c     --------------------
      Counter_Iter   = 0
      ! If each step of the break inherits the calculation results of the previous step, then total_Time
      ! is only set to zero at the initial moment.
      ! For more information, see my notes V3-P78
      if (Key_IniPre_PassOn == 1) then
          total_Time     = ZR
      end if
      ! Old_total_Time is used to indicate the time corresponding to the end of the previous rupture step.
      Old_total_Time(1:Num_Frac) = ZR

      ! Marked for XFEM analysis
      Yes_XFEM = .True.

      do ifra = 1,Num_Frac   
          ! If each rupture step does not inherit the calculation results from the previous step, then at the
          ! start of each rupture step, total_Time is set to zero.
          ! For more information, see my notes V3-P78
          if (Key_IniPre_PassOn == 0) then
              total_Time     = ZR
          end if
          ! Load Factor
          Lambda = ONE
          !-------------------------------
          ! Newton-Raphson Iteration Loop
          !-------------------------------
          do iter = 1,Max_MNR_Iter
              print *, "  "
              WRITE(*,1001) iter,ifra,Num_Frac
              Counter_Iter   =   Counter_Iter +1
              write(*,1201) Counter_Iter
              ifra_Data_iter_num(ifra) = Counter_Iter
              !#########################################################
              ! In the initial iteration step, the enhancement nodes, 
              ! K matrix, and Q matrix only need to be calculated once.
              !#########################################################
              if (iter == 1) then
                  !**************************************************
                  ! Initialize variables (set water pressure and 
                  ! valve opening at each calculation point to zero)
                  !**************************************************
                  Cracks_CalP_Pres(1:Max_Num_Cr,1:Max_Num_Cr_CalP)=ZR
                  Cracks_CalP_Aper(1:Max_Num_Cr,1:Max_Num_Cr_CalP)=ZR


                  !**************************
                  ! Confirm enhancement node
                  !**************************
                  call Determine_Enriched_Nodes(ifra,Counter_Iter)

                  !*********************************************************************************
                  ! Assign numbers to the enhanced nodes and calculate the total degrees of freedom
                  !*********************************************************************************
                  call Number_Enriched_Nodes(Counter_Iter)
                  
                  !********************************************************************
                  ! Save crack-related files (including enhanced node numbering c_POS)
                  !********************************************************************
                  call Save_Files_Crack(Counter_Iter)
                  
                  !**********************
                  ! Save VTK CRACK file.
                  !**********************
                  call Save_vtk_file_for_Crack(Counter_Iter)
      
                  !**********************************************************************************
                  ! Statistically analyze the connectivity between each crack and other cracks, only
                  ! in
                  ! Necessary for hydraulic fracturing analysis
                  !**********************************************************************************
                  call Stat_Crack_Connection(Counter_Iter)
                  
                  !******************************
                  ! Consider boundary conditions
                  !******************************
                  ALLOCATE(freeDOF(Total_FD))
                  call Boundary_Cond(Total_FD,Counter_Iter,
     &                                 freeDOF,num_FreeD)
                  
                  !****************************************************
                  ! Calculate data related to crack calculation points
                  !****************************************************
                  call Cal_HF_Crack_Points_Info_Linear(Counter_Iter)

                  !**********************************************************************************
                  ! Save the coordinates, direction, aperture, and pressure of the crack calculation
                  ! points
                  !**********************************************************************************
                  call Save_Files_Cr_CalP(Counter_Iter)
                  
                  !***********************************************************************
                  ! Hydraulic boundary condition: Constrain the hydraulic degree 
                  ! of freedom at the crack tip (hydraulic pressure atthe crack tip is 0)
                  !***********************************************************************
                  ALLOCATE(freeDOF_HF(num_Tol_CalP_Water))
                  ALLOCATE(Local_freeDOF_HF(num_Tol_CalP_Water))
                  call Boundary_Cond_HF(Counter_Iter,freeDOF_HF,
     &                                  num_free_CalP,Local_freeDOF_HF)
                  
                  !********************************************************************************
                  ! Number vector Total_freeDOF corresponding to all solid and hydrostatic degrees
                  ! of freedom
                  !********************************************************************************
                  ALLOCATE(Total_freeDOF(Total_FD+num_Tol_CalP_Water))
                  num_Total_FD = num_FreeD + num_free_CalP
                  Total_freeDOF(1:num_Total_FD) = 0
                  Total_freeDOF(1:num_FreeD) = freeDOF(1:num_FreeD)
                  Total_freeDOF(num_FreeD+1:num_Total_FD) =
     &                                     freeDOF_HF(1:num_free_CalP)

                  !*************************************
                  ! Show degrees of freedom information
                  !*************************************
 1101             FORMAT(5X,'Total DOF of solid:       ',I8)
 1102             FORMAT(5X,'Total active DOF of solid:',I8)
 1103             FORMAT(5X,'Total CPC:                ',I8)
 1104             FORMAT(5X,'Total active CPC:         ',I8)
                  WRITE(*,1101) Total_FD
                  WRITE(*,1102) num_FreeD
                  WRITE(*,1103) num_Tol_CalP_Water
                  WRITE(*,1104) num_free_CalP

                  !********************************************
                  ! Allocate other relevant dynamic data space
                  !********************************************
                  ALLOCATE(F_U(Total_FD))
                  ALLOCATE(F_U_Backup(Total_FD))
#ifndef github
                  if(Key_Proppant==1) then
                      ALLOCATE(F_U_Ini(Total_FD))
                  endif
#endif
                  ALLOCATE(Last_DISP(Total_FD))
                  ALLOCATE(Last_Last_DISP(Total_FD))
                  ALLOCATE(cLast_Last_DISP(Total_FD))
                  ALLOCATE(DISP(Total_FD))
                  ALLOCATE(DISP_Temp(Total_FD))
                  ALLOCATE(Added_DISP(Total_FD))
                  ALLOCATE(Real_DISP(Total_FD),POST_DISP(Total_FD))
                  ALLOCATE(Coupled_Q(Total_FD,num_Tol_CalP_Water))
                  ALLOCATE(F(Total_FD))
                  ALLOCATE(F_CT(Total_FD))
                  num_ALlDOF = Total_FD+num_Tol_CalP_Water
                  ALLOCATE(R(num_ALlDOF))
                  ALLOCATE(Last_CalP_Pres(num_Tol_CalP_Water))
                  ALLOCATE(CalP_Pres(num_Tol_CalP_Water))
                  ALLOCATE(CalP_Pres_Temp(num_Tol_CalP_Water))
                  ALLOCATE(H(num_Tol_CalP_Water,num_Tol_CalP_Water))
                  ALLOCATE(c_S(num_Tol_CalP_Water))
                  ALLOCATE(NR_Deri(num_ALlDOF,num_ALlDOF))
                  ALLOCATE(delta_x(num_FreeD+num_free_CalP))
                  ALLOCATE(Contact_x(num_FreeD+num_free_CalP))
                  ALLOCATE(Initial_DISP(Total_FD))
                  ALLOCATE(Contact_DISP(Total_FD))
                  ! Related to MNR iterative line search and backtracking algorithm
                  ALLOCATE(x(num_FreeD+num_free_CalP))

                  !***************************
                  ! Initialize some variables
                  !***************************
                  CalP_Pres(1:num_Tol_CalP_Water)      = ZR
                  DISP_Temp(1:Total_FD)                = ZR
                  Added_DISP(1:Total_FD)               = ZR
                  CalP_Pres_Temp(1:num_Tol_CalP_Water) = ZR
                  F(1:Total_FD)                        = ZR
                  F_U(1:Total_FD)                      = ZR

                  !*********************************************************************************
                  ! If the transport of proppant is considered, then the proppant from the previous
                  ! fracturing step
                  ! Concentration mapped to the calculation points corresponding to the new
                  ! fractures
                  !*********************************************************************************
#ifndef github
                  if(Key_Propp_Trans==1 .and. ifra >1)then
                      call Cal_HF_Mapping_ifrac_old2new(ifra,
     &                     Counter_Iter,
     &                     Old_total_Time,total_Time,Yes_Growth,
     &                     L_Cracks_CalP_Conc,
     &                     Map_L_Cracks_CalP_Conc)
                  endif
#endif                  
                  !************************************************************************
                  ! Assembled general load F_U (only executed in the first iteration step)
                  !************************************************************************
                  call Force_Vector(Total_FD,Counter_Iter,.False.,
     &                              Lambda,F_U)
                  F_U_Backup = F_U
                  ! If the effect of the in-situ stress field is taken into account, the external load is zero.
                  if(Key_InSitu_Strategy==1)then
                      F_U(1:Total_FD) = ZR
                  endif
#ifndef github
                  if(Key_Proppant==1)then
                      F_U_Ini = F_U
                  endif
#endif
                  ! If the method of dealing with ground stress is: shielding ground stress, then the ground stress is
                  ! 0.
                  if(Key_InSitu_Strategy==3)then
                      F_U(1:Total_FD) = ZR
                  endif
                  
                  !******************************************************************************
                  ! Calculate the coupling matrix Q (performed only in the first iteration step)
                  !******************************************************************************
                  call Cal_HF_Matrix_Q_Linear(Counter_Iter,
     &                                         Coupled_Q,Total_FD)

                  !******************************************************************************
                  ! Initial iteration step calculation, calculating node displacements under the
                  ! applied test pressure
                  ! Initial test pressure, depending on whether the calculation results of the
                  ! previous fracture step are inherited
                  ! It is also different
                  !******************************************************************************
                  ! Given the initial pressure of each crack
                  call Cal_HF_Initial_Pressure(iter,ifra,Counter_Iter,
     &                     Yes_Growth,Last_Cracks_CalP_Num_ifra,
     &                     Last_Cr_CalP_Pres_ifra,Last_Cr_CalP_Pres)
                  
                  ! Transform into the pressure vector for all calculation points (including only hydraulic fractures)
                  num_Count = 0
                  do i_C = 1,num_Crack
                      ! If the current crack is driven by fracturing fluid, then:
                      if (Cracks_HF_State(i_C) == 1) then
                          do i_CalP=1,Cracks_CalP_Num(i_C)
                              num_Count = num_Count + 1
                              Last_CalP_Pres(num_Count) =
     &                                    Last_Cr_CalP_Pres(i_C,i_CalP)
                          end do
                      end if
                  end do

                  ! Assemble the load vector for calculating the initial displacement
                  F(freeDOF(1:num_FreeD)) =
     &                  MATMUL(Coupled_Q(freeDOF(1:num_FreeD),
     &                         Local_freeDOF_HF(1:num_free_CalP)),
     &               thick(1)*
     &               Last_CalP_Pres(Local_freeDOF_HF(1:num_free_CalP)))
     &                  + F_U(freeDOF(1:num_FreeD))
     
                  !***************************
                  ! Assembly stiffness matrix
                  !***************************
                  ALLOCATE(globalK(Total_FD,Total_FD))
                  print *,'    Assembling the K...'
                  call Assemble_Stiffness_Matrix_XFEM(Counter_Iter,
     &                             globalK,Total_FD,Total_Num_G_P)
                  
                  if(Key_Contact/=0)then
                       ALLOCATE(ori_globalK(Total_FD,Total_FD))
                       ori_globalK = globalK
                  endif
                  write(*,1023) sum(globalK)
                  
                  
 1301             FORMAT(5X,'Total num of Gauss Points:     ',I8)
                  WRITE(*,1301) Total_Num_G_P

                  !************************
                  ! Solve for displacement
                  !************************
                  call Matrix_Solve_LSOE(0,1,Key_SLOE,
     &                         globalK(freeDOF(1:num_FreeD),
     &                                 freeDOF(1:num_FreeD)),
     &                         F(freeDOF(1:num_FreeD)),
     &                         delta_x(1:num_FreeD),num_FreeD)
                  ! Initialize displacement vector
                  DISP(1:Total_FD) = ZR
                  ! Get displacement vector
                  Initial_DISP = ZR
                  Initial_DISP(freeDOF(1:num_FreeD)) =
     &                                   delta_x(1:num_FreeD)
                  DISP(freeDOF(1:num_FreeD))         =
     &                                   delta_x(1:num_FreeD)

                  !*****************************************************************************
                  ! Determine the contact state of the fracture surface 
                  ! (through iteration), and return the new stiffness matrix and displacement
                  !  ----------------
                  ! Note: Sparse matrix section not yet completed_2016-04-04
                  !*****************************************************************************
                  ! Contact time tracking timestamp start
                  call Tool_Get_Current_Time(current_data,
     &                                       date_time,c_CT_time_s)
                  if (Key_Contact /= 0) then
                      Contact_DISP(1:Total_FD) =ZR
                      ! Displacement for contact iteration
                      Contact_DISP(freeDOF(1:num_FreeD))=
     &                                      delta_x(1:num_FreeD)
                      call Determine_Contact_State_by_Iteration(
     &                     iter,ifra,Counter_Iter,Contact_DISP,Total_FD,
     &                     Usual_Freedom,Enrich_Freedom,
     &                     freeDOF,num_FreeD,F,ori_globalK,globalK)
                      ! Return the displacement calculated after the new contact consideration
                      DISP = Contact_DISP
                  endif
                  
                  ! Contact time statistics timestamp ended
                  call Tool_Get_Current_Time(current_data,
     &                                       date_time,c_CT_time_f)
                  Total_CT_CPU_time =  Total_CT_CPU_time +
     &                (dble(c_CT_time_f)-dble(c_CT_time_s))
                  ifra_CT_CPU_time(ifra) = Total_CT_CPU_time

                  !***********************************
                  ! Calculate the initial crack width
                  !***********************************
                  call Cal_Crack_Aperture(Counter_Iter,DISP)
                  call Print_to_Screen_Value_Crack1(4,'WIDTH')
                  call Cal_HF_Initial_Aperture(iter,ifra,Counter_Iter,
     &                  Last_Cracks_CalP_Num_ifra,
     &                  Last_Cr_CalP_Aper_ifra,Last_Cr_CalP_Aper)
                  
                  !***********************************************************
                  ! Calculate delta_Time (based on overall mass conservation)
                  !***********************************************************
                  call Cal_HF_delta_Time_Linear(Counter_Iter,
     &                  Last_Cr_CalP_Aper,total_Time,delta_Time)
 

                  total_Time = total_Time+delta_Time
                  print*,'    Delta time: ',delta_Time,delta_Time/Con_60
                  print*,'    Current time: ',
     &                               total_Time,total_Time/Con_60
                  ! Save the fracturing time corresponding to the current iteration step
                  call Save_HF_time(1,ifra,Counter_Iter,total_Time)
                  
                  !************************************************************************
                  ! Calculate matrix H (requires the crack opening from the previous step)
                  !************************************************************************
                  call Cal_HF_Matrix_H_Linear(ifra,Counter_Iter,H,
     &                                 Cracks_CalP_Aper,total_Time)

                  !********************
                  ! Calculate vector S
                  !********************
                  call Cal_HF_S(Counter_Iter,c_S,total_Time)

                  !********************************************************************************
                  ! Calculate the residual (since it involves K, we need to consider both the full 
                  ! matrix and the sparse matrix cases)
                  !********************************************************************************
                  Last_DISP = DISP
                  print *,'    Calculating residual vector R......'
                  call Cal_HF_Resid(ifra,iter,Counter_Iter,R,
     &                 Total_FD,num_FreeD,num_free_CalP,globalK,
     &                 Coupled_Q,F_U,Last_DISP,Last_Last_DISP,
     &                 freeDOF,freeDOF_HF,Local_freeDOF_HF,
     &                 Last_CalP_Pres,delta_Time,H,c_S)
     
                  
                  !*****************************************************
                  ! Assemble Newton-Raphson Iterative Derivative Matrix
                  !*****************************************************
                  call Cal_HF_Jacobian_NR(Counter_Iter,
     &                   NR_Deri,Total_FD,num_FreeD,
     &                   num_free_CalP,globalK,Coupled_Q,freeDOF,
     &                   freeDOF_HF,Local_freeDOF_HF,delta_Time,H)
                  
                  !******************************************************************************
                  ! Using the penalty function method to correct the Jacobian derivative matrix,
                  ! In order to take into account the condition that the water pressure at the
                  ! junction is the same
                  !******************************************************************************
                  if (Num_JunPair > 0) then
                      call Cal_HF_Jacobian_PenaltyFunc(Counter_Iter,
     &                        NR_Deri,Total_FD,num_FreeD,
     &                        num_free_CalP,freeDOF,freeDOF_HF)
                  end if
                  
                  !******************************************************
                  ! Calculate the iteration increment for the first step
                  !******************************************************
                  ! Timestamp Start
                  call Tool_Get_Current_Time(current_data,
     &                                       date_time,c_LS_time_s)
   
                  call Matrix_Solve_LSOE(0,1,Key_SLOE,
     &                 NR_Deri(Total_freeDOF(1:num_Total_FD),
     &                         Total_freeDOF(1:num_Total_FD)),
     &                 R(Total_freeDOF(1:num_Total_FD)),delta_x,
     &                 num_FreeD+num_free_CalP)
                  !********************************************
                  ! Update the Last variable for linear search
                  !********************************************
                  Temp_Cr_CalP_Aper = Cracks_CalP_Aper
                  Temp_total_Time   = total_Time
                  
                  !*********************************************
                  ! MNR Iterative Line Search and Backtracking!
                  !*********************************************
                  if(Key_HF_LineSearch==1)then
                    call Cal_HF_Line_Searching(ifra,iter,Counter_Iter,
     &                 num_ALlDOF,
     &                 Total_FD,num_FreeD,num_free_CalP,Total_freeDOF,
     &                 num_Tol_CalP_Water,num_Total_FD,freeDOF_HF,
     &                 Local_freeDOF_HF,NR_Deri,delta_x,
     &                 R,Initial_DISP,DISP,freeDOF,
     &                 Last_CalP_Pres,globalK,Coupled_Q,F_U,
     &                 Temp_Cr_CalP_Aper,Temp_total_Time,
     &                 Num_Line_Search,total_Time,x)
                  endif
                  
                  ! Timestamp ended
                  call Tool_Get_Current_Time(current_data,
     &                                       date_time,c_LS_time_f)
                  Total_LS_CPU_time =  Total_LS_CPU_time +
     &                (dble(c_LS_time_f)-dble(c_LS_time_s))
                  ifra_LS_CPU_time(ifra) = Total_LS_CPU_time

                  !******************************************
                  ! Obtain the new displacement and pressure
                  !******************************************
                  if(Key_HF_LineSearch==1)then
                      DISP(freeDOF(1:num_FreeD)) =  x(1:num_FreeD)
                      CalP_Pres(Local_freeDOF_HF(1:num_free_CalP)) =
     &                              x(num_FreeD+1:num_Total_FD)
                  elseif(Key_HF_LineSearch==0)then
                      DISP(freeDOF(1:num_FreeD)) =
     &                          DISP(freeDOF(1:num_FreeD)) -
     &                          delta_x(1:num_FreeD)
                      CalP_Pres(Local_freeDOF_HF(1:num_free_CalP)) =
     &                Last_CalP_Pres(Local_freeDOF_HF(1:num_free_CalP))-
     &                delta_x(num_FreeD+1:num_FreeD+num_free_CalP)
                  endif
                  ! Water pressure cannot be negative
                  if(Key_HF_Del_Ne_Pres==1)then
                      call Vector_ZeroOut_Neg_value(CalP_Pres,
     &                                              num_Tol_CalP_Water)
                  endif
                  !************************************************************************
                  ! Transform the pressure vectors of all calculation points to each crack
                  !************************************************************************
                  call Tool_HF_Pres_Transform(Counter_Iter,CalP_Pres,
     &                                        Cracks_CalP_Pres)
                  call Print_to_Screen_Value_Crack1(4,'PRESS')
                  
                  !************************************
                  ! Calculate the flow rate at a point
                  !************************************
                  call Cal_HF_Flow_Quan(ifra,iter,
     &                                      Counter_Iter,total_time)
                  call Print_to_Screen_Value_Crack1(4,'FLUXX')
                  
                  !*******************
                  ! Save displacement
                  !*******************
                  call Save_Disp(Counter_Iter,1)
                  
                  !**********************************************************************************
                  ! Save the coordinates, direction, aperture, and pressure of the crack calculation
                  ! points
                  !**********************************************************************************
                  call Save_Files_Cr_CalP(Counter_Iter)
                  
                  !*************************************************
                  ! Update the Last variable for the next iteration
                  !*************************************************
                  Last_Cr_CalP_Aper = Temp_Cr_CalP_Aper
                  Last_Cr_CalP_Pres = Cracks_CalP_Pres
                  Last_DISP         = DISP
                  total_Time        = Temp_total_Time
              !#########################################################!######################################
              ! Next iteration step: calculate the H matrix and S vector (considering leakage), then compute
              ! time increment, group set derivative matrix (only need the group set the delta_t*H part in
              ! the lower right corner), calculating residuals, etc.
              !#########################################################!######################################
              elseif (iter >= 2) then
                  !******************************************************************************************
                  ! Determine the contact state of the fracture surface (through iteration), 
                  ! and return the new stiffness matrix and displacement
                  !  --------------------------
                  ! Note, for HF analysis, there are three contact detection schemes: Key_HF_Cont_Scheme=1,
                  ! Then contact detection is performed only during the first fluid-structure 
                  ! iteration; =2, detection in the first two steps; =3,detection at every step
                  !******************************************************************************************
                  ! Contact time tracking timestamp start
                  call Tool_Get_Current_Time(current_data,
     &                                       date_time,c_CT_time_s)
                  if (Key_Contact /= 0) then
                      if (Key_HF_Cont_Scheme==1) then
                          goto 50
                      endif
                      if (Key_HF_Cont_Scheme==2) then
                          if(iter > 2)then
                              goto 50
                          endif
                      endif
                      ! Assembled load vector (including hydrostatic load)
                      F_CT(1:Total_FD) = ZR
                      F_CT(freeDOF(1:num_FreeD)) =
     &                  MATMUL(Coupled_Q(freeDOF(1:num_FreeD),
     &                         Local_freeDOF_HF(1:num_free_CalP)),
     &                   thick(1)*
     &                    CalP_Pres(Local_freeDOF_HF(1:num_free_CalP)))
     &                    + F_U(freeDOF(1:num_FreeD))
                      Contact_DISP(1:Total_FD) =ZR
                      ! Displacement for contact iteration
                      Contact_DISP = DISP
                      call Determine_Contact_State_by_Iteration(
     &                     iter,ifra,Counter_Iter,Contact_DISP,Total_FD,
     &                     Usual_Freedom,Enrich_Freedom,
     &                     freeDOF,num_FreeD,F_CT,ori_globalK,globalK)
                      
                      !write(*,1023) sum(globalK)
                      ! Return the displacement calculated after the new contact consideration
                      DISP = Contact_DISP
                  endif
                  ! Contact time statistics timestamp ended
                  call Tool_Get_Current_Time(current_data,
     &                                       date_time,c_CT_time_f)
                  Total_CT_CPU_time =  Total_CT_CPU_time +
     &                (dble(c_CT_time_f)-dble(c_CT_time_s))
                  ifra_CT_CPU_time(ifra) = Total_CT_CPU_time

   50             continue

                  !********************************************************************************
                  ! Read the displacement from two steps ago, Last_Last_DISP
                  ! Description: The water pressure from the previous step Last_Cr_CalP_Pres, from
                  ! the previous step
                  ! Opening Last_Cr_CalP_Pres, and the previous step's displacement Last_DISP
                  ! It has already been given at the end of the previous calculation!
                  !******************************************************************************** 
                  if(iter==2) then
                      Last_Last_DISP(1:Total_FD)=Initial_DISP
                  elseif (iter >= 3) then
                      
                      call Read_Disp_to(Counter_Iter-2,Last_Last_DISP,
     &                                  Total_FD)
                  end if
                  ! Calculate the new crack width
                  call Cal_Crack_Aperture(Counter_Iter,DISP)
                  call Print_to_Screen_Value_Crack1(4,'WIDTH')
                  
                  !***********************************************************
                  ! Calculate delta_Time (based on overall mass conservation)
                  !***********************************************************
                  call Cal_HF_delta_Time_Linear(Counter_Iter,
     &                     Last_Cr_CalP_Aper,total_Time,delta_Time)
                  

                  total_Time = total_Time + delta_Time
                  print*,'    Delta time: ',delta_Time,delta_Time/Con_60
                  print*,'    Current time: ',total_Time,
     &                                        total_Time/Con_60
                  ! Save the fracturing time corresponding to the current iteration step
                  call Save_HF_time(1,ifra,Counter_Iter,total_Time)
                  
                  !************************************************************************
                  ! Calculate matrix H (matrix H is directly related to fracture aperture)
                  !************************************************************************
                  call Cal_HF_Matrix_H_Linear(ifra,Counter_Iter,H,
     &                                Cracks_CalP_Aper,total_time)
                  
                  !*********************************************************************
                  ! Calculate vector S (S needs to be updated when considering leakage)
                  !*********************************************************************
                  call Cal_HF_S(Counter_Iter,c_S,total_Time)
                  
                  !******************************************************
                  ! Calculate the residual of the current iteration step
                  !******************************************************
                  print *,'    Calculating residual vector R......'
                  call Cal_HF_Resid(ifra,iter,Counter_Iter,R,
     &                              Total_FD,num_FreeD,
     &                              num_free_CalP,globalK,Coupled_Q,F_U,
     &                              Last_DISP,Last_Last_DISP,
     &                              freeDOF,freeDOF_HF,
     &                              Local_freeDOF_HF,CalP_Pres,
     &                              delta_Time,H,c_S)
                  
                  !***************************************************************************
                  ! Modify the Newton-Raphson iteration derivative matrix (bottom-right part)
                  !***************************************************************************
                  call Cal_HF_Modify_Jacobian(Counter_Iter,NR_Deri,
     &                        Total_FD,num_free_CalP,freeDOF_HF,
     &                        Local_freeDOF_HF,delta_Time,H)

                  !******************************************************
                  ! in order to account for the condition that the water
                  ! pressure at the junction points is the same.
                  !******************************************************
                  if (Num_JunPair > 0) then
                      call Cal_HF_Jacobian_PenaltyFunc(Counter_Iter,
     &                      NR_Deri,Total_FD,num_FreeD,num_free_CalP,
     &                      freeDOF,freeDOF_HF)
                  end if
                  
                  !***************************************************
                  ! Calculate the displacement and water 
                  ! pressure increment for the current iteration step
                  !***************************************************
                  ! Timestamp Start
                  call Tool_Get_Current_Time(current_data,
     &                                       date_time,c_LS_time_s)
                  call Matrix_Solve_LSOE(0,1,Key_SLOE,
     &                   NR_Deri(Total_freeDOF(1:num_Total_FD),
     &                           Total_freeDOF(1:num_Total_FD)),
     &                   R(Total_freeDOF(1:num_Total_FD)),delta_x,
     &                   num_FreeD+num_free_CalP)
                  ! Timestamp ended
                  call Tool_Get_Current_Time(current_data,
     &                                       date_time,c_LS_time_f)
                  Total_LS_CPU_time =  Total_LS_CPU_time +
     &                (dble(c_LS_time_f)-dble(c_LS_time_s))
                  ifra_LS_CPU_time(ifra) = Total_LS_CPU_time


                  !********************************************
                  ! Update the Last variable for linear search
                  !********************************************
                  Temp_Cr_CalP_Aper = Cracks_CalP_Aper
                  Temp_total_Time   = total_Time

                  !*********************************************
                  ! MNR iterative line search and backtracking!
                  !*********************************************
                  ! Timestamp Start
                  call Tool_Get_Current_Time(current_data,
     &                                       date_time,c_LS_time_s)
                  if(Key_HF_LineSearch==1)then
                    call Cal_HF_Line_Searching(ifra,iter,Counter_Iter,
     &                 num_ALlDOF,
     &                 Total_FD,num_FreeD,num_free_CalP,Total_freeDOF,
     &                 num_Tol_CalP_Water,num_Total_FD,freeDOF_HF,
     &                 Local_freeDOF_HF,NR_Deri,delta_x,
     &                 R,Initial_DISP,DISP,freeDOF,
     &                 CalP_Pres,globalK,Coupled_Q,F_U,
     &                 Temp_Cr_CalP_Aper,Temp_total_Time,
     &                 Num_Line_Search,total_Time,x)
                  endif
                  
                  call Tool_Get_Current_Time(current_data,
     &                                       date_time,c_LS_time_f)
                  ifra_LS_CPU_time(ifra) = ifra_LS_CPU_time(ifra)+
     &                (dble(c_LS_time_f)-dble(c_LS_time_s))

                  !*************************************************
                  ! Obtain new displacement and pressure (SOR, 
                  ! successive over-relaxation iterative algorithm)
                  !*************************************************
                  ! No slack
                  if(Key_HF_LineSearch==1)then
                      DISP(freeDOF(1:num_FreeD)) =  x(1:num_FreeD)
                      CalP_Pres(Local_freeDOF_HF(1:num_free_CalP)) =
     &                              x(num_FreeD+1:num_Total_FD)
                  elseif(Key_HF_LineSearch==0)then
                      DISP(freeDOF(1:num_FreeD)) =
     &                          DISP(freeDOF(1:num_FreeD)) -
     &                          delta_x(1:num_FreeD)
                      CalP_Pres(Local_freeDOF_HF(1:num_free_CalP)) =
     &                    CalP_Pres(Local_freeDOF_HF(1:num_free_CalP))-
     &                delta_x(num_FreeD+1:num_FreeD+num_free_CalP)
                  endif
                  ! Water pressure cannot be negative
                  if(Key_HF_Del_Ne_Pres==1)then
                      call Vector_ZeroOut_Neg_value(CalP_Pres,
     &                                              num_Tol_CalP_Water)
                  endif

                  !************************************************************************
                  ! Transform the pressure vectors of all calculation points to each crack
                  !************************************************************************
                  call Tool_HF_Pres_Transform(Counter_Iter,CalP_Pres,
     &                                        Cracks_CalP_Pres)
#ifndef Silverfrost
                  if(Key_Unit_System==1)then
                      print *,'    Max pressure of crack 1(MPa):',
     &                 maxval(Cracks_CalP_Pres(1,1:Cracks_CalP_Num(1)))
     &                /Cn_M
                  elseif(Key_Unit_System==2)then
                      print *,'    Max pressure of crack 1(MPa):',
     &                 maxval(Cracks_CalP_Pres(1,1:Cracks_CalP_Num(1)))
                  endif
#endif                  

                  !************************************
                  ! Calculate the flow rate at a point
                  !************************************
                  call Cal_HF_Flow_Quan(ifra,iter,
     &                                      Counter_Iter,total_time)
#ifndef Silverfrost
                  if(Key_Unit_System==1)then
                      print *,'    Max flux of crack 1(m*2/s):',
     &                 maxval(Cracks_CalP_Quan(1,1:Cracks_CalP_Num(1)))
                  elseif(Key_Unit_System==2)then
                      print *,'    Max pressure of crack 1(m*2/s):',
     &                 maxval(Cracks_CalP_Quan(1,1:Cracks_CalP_Num(1)))
     &                 /Cn_M
                  endif
#endif
                  !*******************
                  ! Save displacement
                  !*******************
                  call Save_Disp(Counter_Iter,1)
                  
                  !**********************************************
                  ! Save the coordinates, direction, aperture, 
                  ! and pressure of the crack calculation points
                  !**********************************************
                  call Save_Files_Cr_CalP(Counter_Iter)
                  
                  !****************************************************
                  ! Newton-Raphson Iteration Convergence Determination
                  !****************************************************
                  call Cal_HF_CheckConv(Counter_Iter,
     &                  Last_Cr_CalP_Aper,Last_Cr_CalP_Pres,Yes_Conv)
                  
                  !*************************************************
                  ! Update the Last variable for the next iteration
                  !*************************************************
                  Last_Cr_CalP_Aper = Temp_Cr_CalP_Aper
                  Last_Cr_CalP_Pres = Cracks_CalP_Pres
                  Last_DISP         = DISP
                  total_Time        = Temp_total_Time

                  !*************************************************************
                  ! Determine the new number of cracks and update the next step
                  !num_Crack_Log(Counter_Iter+1)
                  !*************************************************************
                  num_Crack              = num_Crack
                  num_Crack_Log(Counter_Iter+1)  = num_Crack
                  
                  !*********************************************************
                  ! If it converges, exit the Newton-Raphson iteration loop
                  !*********************************************************
                  if (Yes_Conv) then
                      Old_total_Time(ifra) = total_Time
                      goto 100
                  end if
                  !******************************************************
                  ! If it does not converge, then terminate the program.
                  !******************************************************
                  if ((iter ==Max_MNR_Iter)
     &                .and. (Yes_Conv.eqv..False.)) then
                      print *,'    PhiPsi failed to converge!'
                      !call Warning_Message('S',Keywords_Blank)
                  end if
              end if
          end do

  100     continue

          !-----------------------------------
          ! Save or calculate some other data
          !-----------------------------------
     
          ifra_num_iter(ifra)= iter
          ifra_HF_time(ifra) = total_Time
          ! Calculate 1 fracture length
          call Tool_Length_of_Crack(1,ifra_Length(ifra))

          !-----------------------------------------------------------------------
          ! Calculate proppant concentration (upwind method), added on 2015-11-06
          !-----------------------------------------------------------------------
#ifndef github
          if(Key_Propp_Trans==1)then
              ! Calculate proppant concentration
              call Cal_HF_Proppant_c(ifra,Counter_Iter,
     &                        Old_total_Time,total_Time,Yes_Growth)
              ! Maintain proppant concentration
              call Save_Files_Cr_CalP_Propp_c(ifra,Counter_Iter)
              ! Calculate and save the wpnp file for use in static support crack width calculations
              call Cal_HF_Proppant_wpnp(ifra,Counter_Iter)
              call Save_Files_Cr_CalP_Propp_wpnp(ifra,Counter_Iter)
          end if
#endif
          !--------------------------------------------------------------
          ! Save the crack-related files corresponding to the 
          ! convergence steps (including the enhanced node number c_POS)
          ! Purpose: Needed for post-processing
          !--------------------------------------------------------------
          call Save_Files_Crack(Counter_Iter)
          
          !**********************
          ! Save VTK CRACK file.
          !**********************
          call Save_vtk_file_for_Crack(Counter_Iter)

          !----------------------------------------------------------------------------------------
          ! The water pressure and opening at the end of the current fracture step are saved to be
          ! used as the initial values for the next fracture step.
          !----------------------------------------------------------------------------------------
          Last_Total_Freedom_ifra = Total_FD
          do i_C = 1,num_Crack
              num_CalP = Cracks_CalP_Num(i_C)
              Last_Cracks_CalP_Num_ifra(i_C) = num_CalP
              Last_Cr_CalP_Aper_ifra(i_C,1:num_CalP) =
     &               Cracks_CalP_Aper(i_C,1:num_CalP)
              Last_Cr_CalP_Pres_ifra(i_C,1:num_CalP) =
     &               Cracks_CalP_Pres(i_C,1:num_CalP)
          end do

          !-----------------------------------------------------------------------------------
          ! Save the loads corresponding to each solid degree of freedom 
          ! at the convergence step, so that load arrows can be drawn during post-processing.
          !-----------------------------------------------------------------------------------
          if (Key_Post_S_Dof_F == 1) then
              F(1:Total_FD) = ZR
              F(freeDOF(1:num_FreeD)) =
     &           MATMUL(Coupled_Q(freeDOF(1:num_FreeD),
     &                Local_freeDOF_HF(1:num_free_CalP)),
     &           thick(1)*CalP_Pres(Local_freeDOF_HF(1:num_free_CalP)))
     &                    + F_U(freeDOF(1:num_FreeD))
              call Save_Dof_Force(Counter_Iter,F,Total_FD)
          end if
     
          !----------------------------------------------------------
          ! If necessary, calculate and save the tangential relative 
          ! displacement at each crack calculation point.
          !----------------------------------------------------------
          if (Key_Post_S_TanDisp==1)then
             print *,'    Cal and save tangent relative disp of cracks.'
             allocate(Cracks_CalP_Tan_Aper(num_Crack,Max_Num_Cr_CalP))
             call Cal_Crack_Tan_Relative_Disp(Counter_Iter,DISP,
     &                                        Cracks_CalP_Tan_Aper)
             call Save_Crack_Tan_Relative_Disp(Counter_Iter,
     &                                         Cracks_CalP_Tan_Aper)
             deallocate(Cracks_CalP_Tan_Aper)
          endif

          !--------------------------------------------------------------------------------
          ! If needed, calculate and save the Gauss point coordinates for post-processing.
          !--------------------------------------------------------------------------------
          if (Key_Post_CS_G_Coor==1) then
              ALLOCATE(Gauss_CoorX(Total_Num_G_P))
              ALLOCATE(Gauss_CoorY(Total_Num_G_P))
              ! Obtain Gauss point coordinates (Note: Actually, during the process of assembling the stiffness
              ! matrix
              ! The Gaussian point coordinates have already been indirectly calculated, recalculated here)
              call Cal_Gauss_Coors(Counter_Iter,Total_Num_G_P,
     &                             Gauss_CoorX,Gauss_CoorY)
              ! Save Gauss point coordinates
              call Save_Gauss_Coors(Counter_Iter,Total_Num_G_P,
     &                              Gauss_CoorX,Gauss_CoorY)
          end if

          !-------------------------------------------------------------------------------
          ! Calculate the displacement field under the actual load and save it for MATLAB
          ! post-processing.
          !-------------------------------------------------------------------------------
          if(Key_InSitu_Strategy == 2)then
              POST_DISP = DISP
              goto 300
          endif
          ! If there is no ground stress, then directly jump to statement 300.
          if(abs(InSitu_x) < 0.001D6 .and. abs(InSitu_y) < 0.001D6)then
              POST_DISP = DISP
              goto 300
          endif
          ! Water pressure load combined with far-field stress
          t_num_CalP = 0
          do i_C=1,num_Crack
              if(Cracks_HF_State(i_C)==1)then
                  do i_CalP = 1,Cracks_CalP_Num(i_C)
                      t_num_CalP=t_num_CalP+1
                      CalP_Pres(t_num_CalP)= CalP_Pres(t_num_CalP)+
     &                          Cracks_CalP_Remo_Strs(i_C,i_CalP)*ONE
                  end do
              endif
          end do
          ! Assemble the actual water pressure load (Note: it is necessary to include
          ! the water pressure at the calculation point at the crack tip, so it cannot be.
          ! Local_freeDOF_HF(1:num_free_CalP), but 1:num_Tol_CalP_Water)
          F(freeDOF(1:num_FreeD))=MATMUL(Coupled_Q(freeDOF(1:num_FreeD),
     &                         1:num_Tol_CalP_Water),
     &       thick(1)*CalP_Pres(1:num_Tol_CalP_Water))
     &                  + F_U_Backup(freeDOF(1:num_FreeD))

          ! If contact is not considered, the actual displacement field can be used directly for
          ! post-processing.
          if (Key_Contact == 0) then
              POST_DISP = Real_DISP
          ! If contact is considered, determine the contact status of the crack surface (through iteration)
          ! and return the new stiffness matrix and displacement.
          ! Because the contact state of crack surfaces under actual loading may be completely different from
          ! the situation without in-situ stress.
          elseif (Key_Contact /= 0) then
              Contact_DISP(1:Total_FD) =ZR
              ! Displacement for contact iteration
              Contact_DISP = Real_DISP
              call Determine_Contact_State_by_Iteration(
     &             iter,ifra,Counter_Iter,Contact_DISP,Total_FD,
     &             Usual_Freedom,Enrich_Freedom,
     &             freeDOF,num_FreeD,F,ori_globalK,globalK)
              ! Return the displacement calculated after the new contact consideration
              POST_DISP = Contact_DISP
          endif
          
  300     continue

          ! Save the real displacement field for post-processing in Matlab. Note: the DISP involved in
          ! iterative calculations is not the real displacement field, but only
          ! displacement field under net pressure
          call Save_Disp_Real_to_POST(Counter_Iter,POST_DISP)

          !-------------------------------------------------------------------------------------
          ! If necessary, calculate and save the Gauss point displacements for post-processing.
          !-------------------------------------------------------------------------------------
          if (Key_Post_CS_G_Disp==1) then
              ALLOCATE(DISP_x_Gauss(Total_Num_G_P))
              ALLOCATE(DISP_y_Gauss(Total_Num_G_P))
              ! Obtain Gauss point displacement
              call Get_Gauss_Disps(Counter_Iter,POST_DISP,Total_Num_G_P)
              ! Save Gauss point displacement
              call Save_Gauss_Disps(Counter_Iter,Total_Num_G_P)
          end if
          !-----------------------------------------------------
          ! If necessary, calculate and save the node stresses.
          !-----------------------------------------------------
          if (Key_Post_CS_N_Strs ==1) then
              ! Allocate memory space for node stress
              ALLOCATE(Stress_xx_Node(num_Node))
              ALLOCATE(Stress_yy_Node(num_Node))
              ALLOCATE(Stress_xy_Node(num_Node))
              ALLOCATE(Stress_vm_Node(num_Node))
              ! Calculate the nodal stress and store it in the global variables Stress_xx_Node, Stress_yy_Node,
              ! Stress_xy_Node, Stress_vm_Node
              call Get_Node_Stress_XFEM(Counter_Iter,POST_DISP)
              ! Maximum Mises stress at the screen output node
              Max_vm_Node=maxval(Stress_vm_Node(1:num_Node))/Cn_M
              if(Key_Unit_System==1)then
                  WRITE(*,1041) Max_vm_Node
              elseif(Key_Unit_System==2)then
                  WRITE(*,1041) Max_vm_Node*Cn_M
              endif
              ! Save node stress
              call Save_Stress_Node(Counter_Iter,1)
          end if
          !-----------------------------------------------------------------------------------
          ! If needed, calculate and save the Gauss point stress for post-processing display.
          !-----------------------------------------------------------------------------------
          if (Key_Post_CS_G_Strs==1) then
              ALLOCATE(Stress_xx_Gauss(Total_Num_G_P))
              ALLOCATE(Stress_yy_Gauss(Total_Num_G_P))
              ALLOCATE(Stress_xy_Gauss(Total_Num_G_P))
              ALLOCATE(Stress_vm_Gauss(Total_Num_G_P))
              ! Calculate the nodal stress and store it in the global variables Stress_xx_Node, Stress_yy_Node,
              ! Stress_xy_Node, Stress_vm_Node
              call Get_Gauss_Stress_XFEM(Counter_Iter,POST_DISP)
              ! Display the maximum Mises stress at Gauss points on the screen
              Max_vm_Gauss=maxval(Stress_vm_Gauss(1:Total_Num_G_P))/Cn_M
              if(Key_Unit_System==1)then
                  WRITE(*,1042) Max_vm_Gauss
              elseif(Key_Unit_System==2)then
                  WRITE(*,1042) Max_vm_Gauss*Cn_M
              endif
              ! Save Gauss point stress
              call Save_Gauss_Stress(Counter_Iter,Total_Num_G_P)
          end if

          !---------------------------------------------------------------------------------------
          ! If crack propagation is allowed, calculate the stress intensity factor and record the
          ! stress intensity factor.
          !---------------------------------------------------------------------------------------
          if (Key_Propagation ==1) then
              if (Key_SIFs_Method==1) then
                  call Cal_SIFs_DIM(Counter_Iter,POST_DISP)
              elseif (Key_SIFs_Method==2) then
                  call Cal_SIFs_IIM(Counter_Iter,.True.,POST_DISP)
              end if
              do i_C = 1,num_Crack
                  if(Key_Unit_System==1)then
                      WRITE(*,1131)i_C,KI(i_C,1)/Cn_M,KII(i_C,1)/Cn_M
                      WRITE(*,1132)i_C,KI(i_C,2)/Cn_M,KII(i_C,2)/Cn_M
                  elseif(Key_Unit_System==2)then
                      trans_factor = sqrt(1000.0D0)
                      WRITE(*,1131)i_C,KI(i_C,1)/trans_factor,
     &                                 KII(i_C,1)/trans_factor
                      WRITE(*,1132)i_C,KI(i_C,2)/trans_factor,
     &                                 KII(i_C,2)/trans_factor
                  endif
              end do
              call Save_SIFs_KI_and_KII(Counter_Iter)
          end if

          !---------------------------
          ! Save VTK file, 2021-07-16
          !---------------------------
          call Save_vtk_file(Counter_Iter)    
          
          !---------------------------------------------------------------------
          ! Displays the time consumed, and saves or calculates some other data
          !---------------------------------------------------------------------
          call Tool_Get_Current_Time(current_data,date_time,
     &                               curent_time)
          WRITE(*,5001) (dble(curent_time)-dble(S_time))/Con_60
          ifra_CPU_time(ifra) = (dble(curent_time)-dble(S_time))

          ! Solving the real displacement field by solving a system of linear equations
          ! Timestamp Start
          call Tool_Get_Current_Time(current_data,
     &                                       date_time,c_LS_time_s)
          call Matrix_Solve_LSOE(0,1,Key_SLOE,
     &          globalK(freeDOF(1:num_FreeD),freeDOF(1:num_FreeD)),
     &          F(freeDOF(1:num_FreeD)),delta_x(1:num_FreeD),num_FreeD)
          ! Timestamp ended
          call Tool_Get_Current_Time(current_data,
     &                                       date_time,c_LS_time_f)
          Total_LS_CPU_time =  Total_LS_CPU_time +
     &                (dble(c_LS_time_f)-dble(c_LS_time_s))
          ifra_LS_CPU_time(ifra) = Total_LS_CPU_time

          Real_DISP(1:Total_FD) = ZR
          Real_DISP(freeDOF(1:num_FreeD))=delta_x(1:num_FreeD)
          
          !------------------------------------------------------------------------------------------
          ! If the crack is allowed to extend, then determine whether it extends.
          ! If an extension occurs, calculate the new crack tip coordinates, add new crack segments,
          ! and update.
          ! Crack coordinates
          !------------------------------------------------------------------------------------------
          if (Key_Propagation ==1) then
              call Check_Crack_Grows(ifra,Counter_Iter,Yes_Growth)
              ! If a crack has propagated, update the crack propagation flag: Yes_Last_Growth
              if (any(Yes_Growth).eqv..True.)then
                  Yes_Last_Growth = .True.
              ! If no cracks have propagated, exit the program.
              else
                  Yes_Last_Growth = .False.
                  if(ifra < Num_Frac)then
                      print *,'    -----------------------------'
     &                                //'-----------------------'
                      print *,'    Warning :: No crack propapated, '
     &                           //  'program was ended!'
                      print *,'    -----------------------------'
     &                                //'-----------------------'
                  elseif(ifra == Num_Frac)then
                      print *,'    |<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|'
                      print *,'    |    All propagation steps done!  |'
                      print *,'    |<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|'
                  endif
                  goto 199
              endif
          end if

  199     continue
          ! Total number of linear searches (for testing)
          print *,'    Totoal number of line search:',Num_Line_Search

          !--------------------------------------------------------------
          ! If the contact condition of the crack surface is considered,
          ! save the element contact state elcs file.
          !--------------------------------------------------------------
          if(Key_Contact /=0) then
              call Save_Ele_Contac_State(Counter_Iter)
          end if

          !------------------------------------------------------------
          ! If considering the proppant, save the proppant coordinates 
          ! for easier post-processing display.
          !------------------------------------------------------------
#ifndef github
          if(Key_Proppant==1) then
              call Save_Ele_Proppant_Coor(Counter_Iter)
          end if
#endif
          !---------------------------------------------------------------
          ! Update or save the relevant variables of the current fracture 
          ! step into memory for use in the next fracture step
          !---------------------------------------------------------------
          Counter_Num_iFrac(ifra) = Counter_Iter
          if(Num_Frac>=2) then
              call Write_Lfrac_Vars_to_Memory_HF(ifra,Counter_Iter,
     &                        Old_total_Time,total_Time,Yes_Growth)
          end if

          !---------------------------------------------------------------------
          ! Save the coordinates of the natural cracks (only need to save once)
          !---------------------------------------------------------------------
          if(ifra==1 .and. num_Na_Crack > 0)then
              call Save_HF_Natural_Cracks_Coors(num_Na_Crack,
     &                  Na_Crack_Coor(1:num_Na_Crack,1:2,1:2))
          endif

          !--------------------------------------------------------------------------------------
          ! Save the water pressure and corresponding time at each rupture injection point (injp
          ! file)
          !--------------------------------------------------------------------------------------
          call Save_HF_Injection_Pressure(ifra,total_Time)

          !-------------------------------------------------------------------------------------
          ! Save the CPU time corresponding to each fracturing step, fracturing time, length of 
          ! fracture No. 1, and the number of iterations for each fracturing step,
          ! The maximum water pressure of crack No. 1 corresponding to each rupture step. 
          ! (Note: Generally, it only needs to be saved at the end of the program,
          ! But the program might terminate, so the safe approach is to save after each break).
          !-------------------------------------------------------------------------------------
  401     FORMAT(5X,'Pecentage of LS cpu time:   ',F8.4,'%')
  402     FORMAT(5X,'Pecentage of CT cpu time:   ',F8.4,'%')
          write(*,401)(ifra_LS_CPU_time(ifra)/ifra_CPU_time(ifra))*Cn_H
          write(*,402)(ifra_CT_CPU_time(ifra)/ifra_CPU_time(ifra))*Cn_H
          call Save_Ifrac_CPU_time(ifra_CPU_time,Num_Frac)
          call Save_Ifrac_LS_CPU_time(ifra_LS_CPU_time,Num_Frac)
          call Save_Ifrac_CT_CPU_time(ifra_CT_CPU_time,Num_Frac)
          call Save_Ifrac_HF_time(ifra_HF_time,Num_Frac)
          call Save_Ifrac_Num_Iter(ifra_num_iter,Num_Frac)
          call Save_Ifrac_Max_Pres_Crack1(ifra_Max_Press,Num_Frac)
          call Save_Ifrac_Length_Crack1(ifra_Length,Num_Frac)

          !------------------------------------------
          ! Clear the dynamic array and memory space
          !------------------------------------------
          DEALLOCATE(globalK)
          if(Key_Contact/=0)then
              DEALLOCATE(ori_globalK)
          endif
          DEALLOCATE(freeDOF,freeDOF_HF)
          DEALLOCATE(Local_freeDOF_HF);   DEALLOCATE(Total_freeDOF)
          DEALLOCATE(F_U,F_U_Backup,F_CT)
          DEALLOCATE(DISP,DISP_Temp,Last_DISP,Real_DISP,POST_DISP)
          DEALLOCATE(Last_Last_DISP,cLast_Last_DISP)
          DEALLOCATE(Coupled_Q,R);        DEALLOCATE(Last_CalP_Pres)
          DEALLOCATE(CalP_Pres,CalP_Pres_Temp)
          DEALLOCATE(H,c_S);              DEALLOCATE(NR_Deri)
          DEALLOCATE(delta_x);            DEALLOCATE(Contact_x)
          DEALLOCATE(Initial_DISP);       DEALLOCATE(F)
          DEALLOCATE(x)
          DEALLOCATE(Contact_DISP)
          DEALLOCATE(Added_DISP)
          ! Clear the memory space for node stress
          if (Key_Post_CS_N_Strs ==1) then
              DEALLOCATE(Stress_xx_Node,Stress_yy_Node)
              DEALLOCATE(Stress_xy_Node,Stress_vm_Node)
          end if
          ! Clear Gauss point coordinate memory
          if (Key_Post_CS_G_Coor==1) then
              DEALLOCATE(Gauss_CoorX,Gauss_CoorY)
          endif
          ! Clear Gauss point displacement memory
          if (Key_Post_CS_G_Disp==1) then
              DEALLOCATE(DISP_x_Gauss,DISP_y_Gauss)
          endif
          ! Clear Gauss point stress memory space
          if (Key_Post_CS_G_Strs==1) then
              DEALLOCATE(Stress_xx_Gauss,Stress_yy_Gauss)
              DEALLOCATE(Stress_xy_Gauss,Stress_vm_Gauss)
          endif
          ! Clear memory space related to crack surface contact iterations and proppant support iterations
#ifndef github
          if(Key_Proppant==1)then
              DEALLOCATE(F_U_Ini)
          endif
#endif
          ! If all cracks stop propagating, terminate the program.
          if(Yes_Last_Growth .eqv. .False.)then
              goto 200
          endif
      end do

  200 continue

      ! Clear the memory public variables related to materials
      DEALLOCATE(D,S,St,Sc,T_Alpha,KIc,E,v,thick,density)


      RETURN
      END SUBROUTINE PhiPsi2D_Static_HF
