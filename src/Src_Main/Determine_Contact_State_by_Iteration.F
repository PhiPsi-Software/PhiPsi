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
 
      SUBROUTINE Determine_Contact_State_by_Iteration(
     &             iter,ifra,Counter_Iter,
     &             Contact_DISP,c_Total_FD,
     &             usual_FD,enrich_FD,
     &             c_freeDOF,c_num_freeDOF,c_F,
     &             ori_globalK,CT_Jacobian)
                                                     ! CT_Jacobian is the stiffness matrix corresponding to the contact state
     
c     This subroutine determines the contact status of the fracture surface through iterative calculation.
c     Contact detection algorithm for crack surfaces Keywords: Key_Contact
c                   1: Penalty method.
c                   2: Reduced Penalty method (under testing).
c                   3: Lagrange multipliers technique (temporarily unavailable).
c                   4: Based on linear complementarity theory (temporarily unavailable).
                                  
c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Filename
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      use Global_Contact
      use Global_Inter_Cal_Contact_Red_Resid
      
c     ---------------------------------------------------------------------------------
c     Read subroutine interface module (activate compiler parameter consistency check)
c     ---------------------------------------------------------------------------------
      use Global_Inter_Get_Contact_State_Eles_IN
      
c     --------------------------
c     Variable Type Declaration
c     --------------------------
      implicit none
      integer,intent(in)::iter,ifra,Counter_Iter
      integer,intent(in)::c_Total_FD,c_num_freeDOF
      integer,intent(in)::usual_FD,enrich_FD
      real(kind=FT),intent(inout)::Contact_DISP(c_Total_FD)
      integer,intent(in)::c_freeDOF(c_Total_FD)
      real(kind=FT),intent(in)::c_F(c_Total_FD)
      real(kind=FT),intent(in)::ori_globalK(c_Total_FD,
     &                                         c_Total_FD)
      real(kind=FT),intent(out)::CT_Jacobian(c_Total_FD,
     &                                          c_Total_FD)
      integer i_C,i_E,i_N,c_Edge_Elem
      logical Yes_Contact
      ! ------------Variables Related to Penalty Function Method-------------
      integer i_NR_P
      real(kind=FT) R_PSI(c_Total_FD)
      real(kind=FT) Last_R_PSI(c_Total_FD)
      real(kind=FT) NR_DISP(c_Total_FD)
      real(kind=FT) PC_Gauss_x(num_Crack,Max_Num_Cr_CalP-1,2)
      real(kind=FT) PC_Gauss_y(num_Crack,Max_Num_Cr_CalP-1,2)
      integer CT_State_Gauss(num_Crack,Max_Num_Cr_CalP-1,2)
                                                                    ! = 0, Separation; = 1, Bonding; = 2, Sliding
      real(kind=FT) Kn,Kt
      real(kind=FT) Kn_Gauss(num_Crack,Max_Num_Cr_CalP-1,2)
      real(kind=FT) Kt_Gauss(num_Crack,Max_Num_Cr_CalP-1,2)
      real(kind=FT) delta_u_a(c_Total_FD)
      real(kind=FT) Conve_Tolerance,Conve_Factor
      real(kind=FT) tem_DISP(c_num_freeDOF),fric_mu
      logical Yes_Conve
      real(kind=FT) Contact_DISP_0(c_Total_FD)
      real(kind=FT),ALLOCATABLE:: Part1(:,:), Part2(:) 
      real(kind=FT),ALLOCATABLE:: K_uu_Inv(:,:),K_uu(:,:)
      real(kind=FT),ALLOCATABLE:: K_ee(:,:),K_ue(:,:),J_red(:,:)
      real(kind=FT),ALLOCATABLE:: ori_K_ee(:,:)
      real(kind=FT),ALLOCATABLE:: U_e(:),delta_U_e(:),R_red(:)
      real(kind=FT),ALLOCATABLE:: Last_R_red(:)
      real(kind=FT),ALLOCATABLE:: tem_DISP_red(:)
      real(kind=FT),ALLOCATABLE:: U_e_0(:),U_u(:),red_tem1(:)
      real(kind=FT),ALLOCATABLE:: U_all_0(:),U_all(:),U_all_last(:)
      real(kind=FT),ALLOCATABLE:: delta_U_all(:)
      integer n_freeDOF_FEM
      integer,ALLOCATABLE:: freeDOF_FEM(:)
      real(kind=FT) Saved_Conv_Factor(Max_Contact_Iter)
      logical Yes_Oscill
      
c     ---------------------------
c     Formatted output statement
c     ---------------------------
 1997 FORMAT(5X,'########################################!')  
 1998 FORMAT(5X,'########################################!')  
 2002 FORMAT(5X,'  PhiPsi failed to converge for contact iteration!')   
 2003 FORMAT(5X,'  WARNING :: Oscillation detected!')   
 4001 FORMAT(5X,'+++  Contact NR iteration ',I3,' of ',I3,' started:')
 4012 FORMAT(5X,'No penetration detected, leaving contact iteration!')
 4022 FORMAT(10X,'Convergence factor is', E16.7 ,'.')  
 4032 FORMAT(5X,'Contact iteration done after ',I3,' tries!') 
 4033 FORMAT(5X,'Number of contact elements is ',I5) 
 5001 FORMAT(10X,'Sticking elements:',I5,'| Sliding elements:',I5) 
      print *,'    Determine contact states by iteration......' 
      
      
      
c     ------------------------------------
c     Initialization and Data Preparation
c     ------------------------------------
      Contact_DISP_0(1:c_Total_FD) = ZR
      Contact_DISP_0 = Contact_DISP
      Saved_Conv_Factor(1:Max_Contact_Iter)  = ZR
      
c     -----------------------
c     Key control parameters
c     -----------------------
      ! Penalty function method related (applicable to contact iterations 1 and 2)
      ! Penalty stiffness (reducing it to 1.0D11 and 1.0D10 results in better convergence stability, but
      ! is inaccurate)
      kn = kn_Cont_Penalty
      kt = kt_Cont_Penalty
      ! Coefficient of friction
      fric_mu = fric_mu_Cont
      ! Convergence tolerance
      Conve_Tolerance = Conve_Tol_Penalty
      
c     ---------------------------------------------------------------------------------
c     Select different contact algorithms based on the keyword Key_Contact
c     These are further divided into conventional algorithms and reduction algorithms.
c     ---------------------------------------------------------------------------------
      select case (Key_Contact)
      !!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !case(1): Penalty function method and Newton-Raphson iteration
      !---------------------------------------------------------------
      ! Note: (1) The crack-tip enhancement element corresponds to the part of the crack segment that
      ! participates in the iterative calculation (previous older version
      ! Originally not included in the calculation);
      !!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      case(1)
          ! Newton-Raphson Contact Iteration
          do i_NR_P = 1,Max_Contact_Iter
              write(*,4001) i_NR_P,Max_Contact_Iter
              !#####################################################################################
              ! Special handling for the first iteration step: If it is the first NR iteration, the
              ! crack needs to be determined first
              ! If there is no contact on the crack surface in the joint, the iteration ends
              ! immediately.
              !#####################################################################################
              if(i_NR_P==1)then 
                  NR_DISP(1:c_Total_FD)    = ZR
                  delta_u_a(1:c_Total_FD)  = ZR    
                  Last_R_PSI(1:c_Total_FD) = ZR
                  Elem_Conta_Sta(1:Num_Elem,1:Max_Num_Cr) = 0 
                  Kt_Gauss(1:num_Crack,1:Max_Num_Cr_CalP-1,1:2) = kt
                  Kn_Gauss(1:num_Crack,1:Max_Num_Cr_CalP-1,1:2) = kn
                  call Cal_Contact_Contact_State_Gauss(
     &                    iter,ifra,Counter_Iter,i_NR_P,
     &                    Contact_DISP,Yes_Contact,Elem_Conta_Sta,
     &                    CT_State_Gauss)    
                  ! If no embedding occurs on any of the crack surfaces, exit the iteration.
                  if(Yes_Contact.eqv..False.)then
                      write(*,4012)
                      CT_Jacobian = ori_globalK
                      goto 9999
                  endif
              endif
              
              !##################################################################################
              ! Calculate the components of the Gauss point contact force in the x and y directions, and update
              ! and save the Gauss point contact status in the process
              !------------------------------------------------
              ! First, obtain the contact status of the Gauss points on each fracture surface and store it in
              ! CT_State_Gauss and Elem_Conta_Sta (for post-processing).
              !##################################################################################
              print *,'         Get contact force and update Kt.'
              call Cal_Contact_PN_and_PT(iter,ifra,Counter_Iter,
     &                  i_NR_P,c_Total_FD,c_num_freeDOF, 
     &                  Kn,Kn_Gauss,Kt_Gauss,fric_mu,NR_DISP,delta_u_a,
     &                  CT_State_Gauss,
     &                  Elem_Conta_Sta,
     &                  PC_Gauss_x,PC_Gauss_y)
              write(*,5001) count(Elem_Conta_Sta==1),
     &                      count(Elem_Conta_Sta==2)
     
              !######################################################################################
              ! Assemble the contact iteration Jacobian matrix, calculated based on ori_globalK (can
              ! reduce computation)
              !######################################################################################
              print *,'         Assemble Jacobian matrix.'
              call Cal_Contact_Jacobian(     
     &           iter,ifra,Counter_Iter,i_NR_P,
     &           c_Total_FD,c_num_freeDOF,ori_globalK,
     &           c_freeDOF,Kn,Kn_Gauss,Kt_Gauss,
     &           CT_State_Gauss,CT_Jacobian)
              
              
              !##########################
              ! Calculate residuals, PSI
              !##########################
              print *,'         Get residual.'
              !------------------------------------------------------------
              ! option 1: Calculate internal forces through Bσ integration
              !------------------------------------------------------------
              call Cal_Contact_Resid(
     &           iter,ifra,Counter_Iter,i_NR_P,
     &           c_Total_FD,c_num_freeDOF,c_F,
     &           NR_DISP,c_freeDOF,PC_Gauss_x,PC_Gauss_y,R_PSI)
              
              !-------------------------------------------------------------------------------
              ! Option 2: Directly calculate internal forces through Ku (this option does not
              ! support crack reinforcement)
              !-------------------------------------------------------------------------------
C             call Cal_Contact_Resid_option2(
C    &           iter,ifra,Counter_Iter,i_NR_P,
C    &           c_Total_FD,c_num_freeDOF,c_F,
C    &           CT_Jacobian,ori_globalK,      
C    &           NR_DISP,c_freeDOF,PC_Gauss_x,PC_Gauss_y,R_PSI)
              
              !############################################
              ! Solve for displacement increments (△u, △a)
              !############################################
              print *,'         Get displacemnt increament.'
              delta_u_a(1:c_Total_FD)  =ZR
              call Matrix_Solve_LSOE(5,1,Key_SLOE,    
     &                    CT_Jacobian(c_freeDOF(1:c_num_freeDOF),
     &                                c_freeDOF(1:c_num_freeDOF)),
     &                   -R_PSI(c_freeDOF(1:c_num_freeDOF)),
     &                    tem_DISP,c_num_freeDOF)
              delta_u_a(c_freeDOF(1:c_num_freeDOF)) = tem_DISP
              
              !######################
              ! Update displacement.
              !######################
              print *,'         Update displacemnt.'
              
              NR_DISP = NR_DISP + delta_u_a
              
              
              !########################
              ! Check for convergence.
              !########################
              print *,'         Check convergence.'
              call Cal_Contact_Conve_Factor(
     &                iter,ifra,Counter_Iter,i_NR_P,Conve_Tolerance,
     &                c_Total_FD,c_freeDOF,c_num_freeDOF,
     &                c_F,R_PSI,Last_R_PSI,
     &                delta_u_a,NR_DISP,Contact_DISP_0,
     &                Yes_Conve,Conve_Factor)
              write(*,4022) Conve_Factor
              Saved_Conv_Factor(i_NR_P) =Conve_Factor
              ! Exit the iteration loop if convergence is achieved.
              if(Yes_Conve)then
                  write(*,1997)
                  write(*,4032) i_NR_P
                  write(*,4033) count(Elem_Conta_Sta/=0)
                  write(*,1997)
                  exit
              endif
              ! Check whether the convergence factor has oscillated (when the number of iterations is greater than
              ! or equal to 6)
              if (i_NR_P>=6)then
                  ! Determine whether oscillation has occurred based on the last six convergence values.
                  print *,'         Check oscillation.'
                  call Tool_Check_Oscillation_by_6_Variables(
     &                 Saved_Conv_Factor(i_NR_P-5:i_NR_P),Yes_Oscill)
                  if (Yes_Oscill) then
                      write(*,1998)
                      write(*,2003)
                      write(*,1998)
                      exit
                  endif
              endif
              ! Update the residual from the previous step
              Last_R_PSI = R_PSI
          enddo
          
          ! If all iterations are completed and it still hasn't converged, the program will terminate.
          if(.not.Yes_Conve)then
              write(*,2002)
              !call Warning_Message('S',Keywords_Blank)
          endif
          
          ! Save output variable
          Contact_DISP = NR_DISP
          
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !case(2): Reduced Penalty function method (under testing).
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      case(2)
#ifndef github
          !*****************************************************************
          !First Step: data preparation.
          !usual_FD and enrich_FD are num of FEM DOFs and num of XFEM DOFs.
          !*****************************************************************
          ! Number of degrees of freedom for traditional degrees of freedom
          n_freeDOF_FEM = c_num_freeDOF - enrich_FD
          ! Conventional degrees of freedom numbering
          ALLOCATE(freeDOF_FEM(n_freeDOF_FEM))  
          freeDOF_FEM(1:n_freeDOF_FEM) = c_freeDOF(1:n_freeDOF_FEM)
          ! Allocate memory space for other variables
          ALLOCATE(Part1(enrich_FD,enrich_FD))
          ALLOCATE(Part2(enrich_FD))
          ALLOCATE(K_uu_Inv(n_freeDOF_FEM,n_freeDOF_FEM))
          ALLOCATE(K_ee(enrich_FD,enrich_FD))         
          ALLOCATE(ori_K_ee(enrich_FD,enrich_FD))
          ALLOCATE(J_red(enrich_FD,enrich_FD))
          ALLOCATE(K_ue(n_freeDOF_FEM,enrich_FD))
          ALLOCATE(U_e_0(enrich_FD),U_e(enrich_FD))
          ALLOCATE(delta_U_e(enrich_FD))
          ALLOCATE(R_red(enrich_FD))
          ALLOCATE(Last_R_red(enrich_FD))       
          ALLOCATE(tem_DISP_red(enrich_FD))  
          ALLOCATE(U_u(c_Total_FD))  
          ALLOCATE(red_tem1(n_freeDOF_FEM))      
          ALLOCATE(U_all_0(c_Total_FD),U_all(c_Total_FD),
     &             U_all_last(c_Total_FD),delta_U_all(c_Total_FD))
          K_uu_Inv(1:n_freeDOF_FEM,1:n_freeDOF_FEM) = 
     &               ori_globalK(freeDOF_FEM(1:n_freeDOF_FEM),
     &                           freeDOF_FEM(1:n_freeDOF_FEM))
          print *,'      Get the inverse matrix of Kuu......'
#ifndef Silverfrost
          call Matrix_Inverse_Lapack(K_uu_Inv(1:n_freeDOF_FEM,
     &                                        1:n_freeDOF_FEM),
     &                               n_freeDOF_FEM)  
#endif
#ifdef Silverfrost
          ALLOCATE(K_uu(n_freeDOF_FEM,n_freeDOF_FEM))
          call Matrix_Inverse(K_uu(1:n_freeDOF_FEM,
     &                                 1:n_freeDOF_FEM),
     &                        K_uu_Inv(1:n_freeDOF_FEM,
     &                                 1:n_freeDOF_FEM),
     &                        n_freeDOF_FEM)
          deallocate(K_uu) 
#endif
          ! Check whether the calculation of the inverse matrix is correct
          print *,'      Get Kue and Kee......'
          K_ue(1:n_freeDOF_FEM,1:enrich_FD) = 
     &               ori_globalK(freeDOF_FEM(1:n_freeDOF_FEM),
     &                           usual_FD+1:c_Total_FD)
          K_ee(1:enrich_FD,1:enrich_FD) = 
     &          ori_globalK(usual_FD+1:c_Total_FD,usual_FD+1:c_Total_FD)
          !Invariant Part1 and Part2.
          Part1 = MATMUL(transpose(K_ue),MATMUL(K_uu_Inv,K_ue))
          Part2 = MATMUL(MATMUL(transpose(K_ue),K_uu_Inv),
     &                   c_F(freeDOF_FEM(1:n_freeDOF_FEM)))
          
          !**************************
          !Newton-Raphson iteration.
          !**************************
          do i_NR_P = 1,Max_Contact_Iter
              write(*,4001) i_NR_P,Max_Contact_Iter
              !#####################################################################################
              ! Special handling for the first iteration step: If it is the first NR iteration, the
              ! crack needs to be determined first
              ! If there is no contact on the crack surface in the joint, the iteration ends
              ! immediately.
              !#####################################################################################
              if(i_NR_P==1)then 
                  ori_K_ee = K_ee
                  U_e(1:enrich_FD) = ZR
                  U_e_0(1:enrich_FD) =  Contact_DISP(Usual_Freedom+1:
     &                                               c_Total_FD)
                  U_all_0 = Contact_DISP
                  !U_all_last = Contact_DISP 
                  U_all_last(1:c_Total_FD) = ZR
                  !U_all = Contact_DISP  
                  U_all(1:c_Total_FD) = ZR
                  delta_U_all(1:c_Total_FD) = ZR
                  delta_U_e(1:enrich_FD)  = ZR
                  Last_R_red(1:enrich_FD) = ZR
                  
                  Elem_Conta_Sta(1:Num_Elem,1:Max_Num_Cr) = 0 
                  Kt_Gauss(1:num_Crack,1:Max_Num_Cr_CalP-1,1:2) = kt
                  call Cal_Contact_Red_Contact_State_Gauss(
     &                    iter,ifra,Counter_Iter,i_NR_P,
     &                    U_e_0,Yes_Contact,Elem_Conta_Sta,
     &                    CT_State_Gauss)    
                  ! If no embedding occurs on any of the crack surfaces, exit the iteration (go to the end of the
                  ! subroutine)
                  if(Yes_Contact.eqv..False.)then
                      write(*,4012)
                      CT_Jacobian = ori_globalK
                      goto 9999
                  endif
              endif
              
              !##################################################################################
              ! Calculate the components of the contact force at the Gauss points of the contact element in the x
              ! and y directions, and determine Kn and Kt (normal and tangential penalty stiffness) according to
              ! the contact criteria.
              ! And conveniently update and save the Gauss point contact status
              !------------------------------------------------
              ! First, obtain the contact status of the Gauss points on each crack surface and store it in
              ! CT_State_Gauss and Elem_Conta_Sta (for post-processing).
              !##################################################################################
              print *,'         Get contact force and update Kt.'
              call Cal_Contact_PN_and_PT(iter,ifra,Counter_Iter,
     &                  i_NR_P,c_Total_FD,c_num_freeDOF, 
     &                  Kn,Kn_Gauss,Kt_Gauss,fric_mu,U_all,delta_U_all,
     &                  CT_State_Gauss,
     &                  Elem_Conta_Sta,
     &                  PC_Gauss_x,PC_Gauss_y)
     
              write(*,5001) count(Elem_Conta_Sta==1),
     &                      count(Elem_Conta_Sta==2)

              !###########################
              ! Calculate residual, R_red
              !###########################
              print *,'         Get residual for the reduced scheme.'
              call Cal_Contact_Red_Resid(
     &                 iter,ifra,Counter_Iter,i_NR_P,
     &                 c_Total_FD,c_num_freeDOF,n_freeDOF_FEM,c_F,
     &                 freeDOF_FEM,
     &                 U_e,U_all,c_freeDOF,
     &                 Part1,Part2,K_ee,K_uu_Inv,K_ue,
     &                 PC_Gauss_x,PC_Gauss_y,R_red)
     
      
              !######################################################################################
              ! Update K_ee based on the contact state and Kn, Kt, and calculate the Jacobian matrix
              ! J_red
              !######################################################################################
              print *,'         Update Kee and Jacobian matrix.'
              call Cal_Contact_Red_Update_Kee_and_Jacobian(     
     &                 iter,ifra,Counter_Iter,i_NR_P,
     &                 c_Total_FD,c_num_freeDOF,ori_globalK,
     &                 c_freeDOF,Kn,Kt_Gauss,CT_State_Gauss,
     &                 Part1,Part2,
     &                 ori_K_ee,K_ee,J_red)
              print *,'         sum(K_ee)/1.0D9:',sum(K_ee)/1.0D9             
     
              !###################################
              ! Solve for displacement increment.
              !###################################
              print *,'         Get displacemnt increament.'
              delta_U_e(1:enrich_FD)  =ZR
              call Matrix_Solve_LSOE(7,1,Key_SLOE,    
     &                    J_red,
     &                    R_red,
     &                    tem_DISP_red,enrich_FD)
              delta_U_e = tem_DISP_red
              
              !######################
              ! Update displacement.
              !######################
              print *,'         Update displacemnt.'
              U_e = U_e + delta_U_e
              
              !###########################
              ! Solve for U_u.
              !--------------------------
              !Kuu*Uu + Kue*Ue - Fu =0
              !###########################
              U_u(1:c_Total_FD) = ZR
              ! Output displacement
              red_tem1 = c_F(freeDOF_FEM(1:n_freeDOF_FEM)) - 
     &                   MATMUL(K_ue,U_e)
              U_u(freeDOF_FEM(1:n_freeDOF_FEM)) = 
     &                               MATMUL(K_uu_Inv,red_tem1)
              U_all(1:usual_FD)  = U_u
              U_all(usual_FD+1:c_Total_FD) = U_e
              delta_U_all = U_all - U_all_last
              
              !#####################################################################################
              ! Check convergence (the convergence displacement algorithm must be determined by the
              ! full displacement, not just
              ! By enhancing the displacement, because the displacements of the FEM degrees of
              ! freedom are also changing.
              !#####################################################################################
              print *,'         Check convergence.'
              call Cal_Contact_Red_Conve_Factor(
     &                iter,ifra,Counter_Iter,i_NR_P,Conve_Tolerance,
     &                c_Total_FD,c_F,R_red,Last_R_red,
     &                U_all_last,U_all,U_all_0,
     &                Yes_Conve,Conve_Factor)
              write(*,4022) Conve_Factor   
              ! Exit the iteration loop if convergence is achieved.
              if(Yes_Conve)then
              !if(i_NR_P == Max_Contact_Iter)then
                  write(*,1997)
                  write(*,4032) i_NR_P
                  write(*,4033) count(Elem_Conta_Sta/=0)
                  write(*,1997)
                  exit
              endif
              ! Update the residual from the previous step
              Last_R_red = R_red         
              ! Update the full displacement from the previous step
              U_all_last = U_all
      
          enddo
          
          !****************************************************
          ! Obtain the full stiffness matrix and displacement.
          !****************************************************
          CT_Jacobian = ori_globalK
          CT_Jacobian(usual_FD+1:c_Total_FD,
     &                usual_FD+1:c_Total_FD) = K_ee
          Contact_DISP = U_all
          
          !**********************
          !Clear temperary data.
          !**********************
          DEALLOCATE(Part1,Part2)
          DEALLOCATE(K_uu_Inv,K_ue,K_ee,ori_K_ee,J_red)
          DEALLOCATE(U_e_0,U_e,delta_U_e,R_red,Last_R_red)
          DEALLOCATE(tem_DISP_red)
          DEALLOCATE(delta_U_all)
          DEALLOCATE(U_u)
#endif
      end select

 9999 continue
 
      
      RETURN
      END SUBROUTINE Determine_Contact_State_by_Iteration
