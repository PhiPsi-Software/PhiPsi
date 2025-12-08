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
 
SUBROUTINE Determine_Contact_State_by_Iteration_3D( &
           iter,ifra,Counter_Iter, &
           Contact_DISP,c_Total_FD,   &
           usual_FD,enrich_FD,&
           c_freeDOF,c_num_freeDOF,c_F,&
           ori_globalK,CT_Jacobian)
                                             ! CT_Jacobian is the stiffness matrix corresponding to the contact state

! This subroutine determines the contact status of the fracture surface through iterative
! calculation.
! Contact detection algorithm for crack surfaces Keywords: Key_Contact
!                   1: Penalty method.
!                   2: Reduced Penalty method (under testing).
!                   3: Lagrange multipliers technique (temporarily unavailable).
!                   4: Based on linear complementarity theory (temporarily unavailable).
!
!   
   
!/////////////////////////////
! Read public variable module
!/////////////////////////////
use Global_Float_Type
use Global_Common
use Global_Filename
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Global_Contact
use Global_POST
use Global_Inter_Cal_Contact_Red_Resid

!////////////////////////////////////////////////////////
! Read subroutine interface module (activate compilation
! Device parameter consistency check
!////////////////////////////////////////////////////////
use Global_Cal_Contact_Resid_3D
use Global_Cal_Contact_Contact_State_Gauss_3D
use Global_Cal_Contact_PN_and_PT_3D
use Global_Cal_Contact_Conve_Factor

!///////////////////////////
! Variable Type Declaration
!///////////////////////////
implicit none
integer,intent(in)::iter,ifra,Counter_Iter
integer,intent(in)::c_Total_FD,c_num_freeDOF
integer,intent(in)::usual_FD,enrich_FD
real(kind=FT),intent(inout)::Contact_DISP(c_Total_FD)
integer,intent(in)::c_freeDOF(c_Total_FD)
real(kind=FT),intent(in)::c_F(c_Total_FD)
real(kind=FT),intent(in)::ori_globalK(c_Total_FD,c_Total_FD)
real(kind=FT),intent(out)::CT_Jacobian(c_Total_FD,c_Total_FD)
integer i_C,i_E,i_N,c_Edge_Elem
logical Yes_Contact
! ------------Variables Related to Penalty Function Method-------------
integer i_NR_P
real(kind=FT) R_PSI(c_Total_FD)
real(kind=FT) Last_R_PSI(c_Total_FD)
real(kind=FT) NR_DISP(c_Total_FD)
real(kind=FT) PC_Gauss_x(num_Crack,Max_Max_N_FluEl_3D)
real(kind=FT) PC_Gauss_y(num_Crack,Max_Max_N_FluEl_3D)
real(kind=FT) PC_Gauss_z(num_Crack,Max_Max_N_FluEl_3D)
integer CT_State_Gauss(num_Crack,Max_Max_N_FluEl_3D)
                                                            ! = 0, Separation; = 1, Bonding; = 2, Sliding
real(kind=FT) Kn,Kt
real(kind=FT) Kn_Gauss(num_Crack,Max_Max_N_FluEl_3D)
real(kind=FT) Kt1_Gauss(num_Crack,Max_Max_N_FluEl_3D)
real(kind=FT) Kt2_Gauss(num_Crack,Max_Max_N_FluEl_3D)
real(kind=FT) delta_u_a(c_Total_FD)
real(kind=FT) Conve_Tolerance,Conve_Factor
real(kind=FT) tem_DISP(c_num_freeDOF),fric_mu
logical Yes_Conve
real(kind=FT) Contact_DISP_0(c_Total_FD)

real(kind=FT) Saved_Conv_Factor(Max_Contact_Iter)
logical Yes_Oscill
character(200) c_File_name_1,c_File_name_2,c_File_name_3
character(5) temp
real(kind=FT) tem_1 ,tem_2
integer i,j

!////////////////////////////    
! Formatted output statement
!////////////////////////////  
1997 FORMAT(7X,'\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!')  
1998 FORMAT(7X,'\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!')  
2002 FORMAT(7X,'  PhiPsi failed to converge for contact iteration!')   
2003 FORMAT(7X,'  WARNING :: Oscillation detected!')   
4001 FORMAT(7X,'+++  Contact NR iteration ',I3,' of ',I3,' started:') 
4012 FORMAT(7X,'No penetration detected, leaving contact iteration!') 
4022 FORMAT(12X,'Convergence factor is', E12.4 ,' | ',E12.4,'.')  
4032 FORMAT(7X,'Contact iteration done after ',I3,' tries!') 
4033 FORMAT(7X,'Number of contact elements is ',I5) 
5001 FORMAT(12X,'Sticking solid elements:  ',I5,'| Sliding solid elements:  ',I5) 
5002 FORMAT(12X,'Sticking contact elements:',I5,'| Sliding contact elements:',I5) 
6001 FORMAT(12X,'SUM(ABS(R_PSI)):', E12.4)  
print *,'    Determine contact states by iteration......' 


!/////////////////////////////////////
! Initialization and Data Preparation
!/////////////////////////////////////
Contact_DISP_0(1:c_Total_FD) = ZR
Contact_DISP_0 = Contact_DISP
Saved_Conv_Factor(1:Max_Contact_Iter)  = ZR

!IMPROV2022081501.
IF(ALLOCATED(Elem_Conta_Sta)) DEALLOCATE(Elem_Conta_Sta)  
ALLOCATE(Elem_Conta_Sta(Num_Elem,num_Crack))
Elem_Conta_Sta(1:Num_Elem,1:num_Crack)=0

!////////////////////////
! Key control parameters
!////////////////////////
! Penalty function method related (applicable to contact iterations 1 and 2)
!-------
! Penalty stiffness (reducing it to 1.0D11 and 1.0D10 results in better convergence stability, but
! it is inaccurate)
kn = kn_Cont_Penalty
kt = kt_Cont_Penalty
! Coefficient of friction
fric_mu = fric_mu_Cont
! Convergence tolerance
Conve_Tolerance = Conve_Tol_Penalty

!//////////////////////////////////////////////////////////////////////////////////
!  
! Select different contact algorithms based on the keyword Key_Contact
! These are further divided into conventional algorithms and reduction algorithms.
!
!//////////////////////////////////////////////////////////////////////////////////
select case (Key_Contact)
!
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!case(1): Penalty function method and Newton-Raphson iteration
!
! Note: (1) The crack tip enrichment element corresponds to the part of the crack segment that
! participates in the iterative calculation (in previous older versions, it did not participate in
! the calculation);
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
case(1)
  ! Newton-Raphson Contact Iteration
  do i_NR_P = 1,Max_Contact_Iter
      write(*,4001) i_NR_P,Max_Contact_Iter
      !###########################################################################################
      ! Special handling for the first iteration step: If it is the first NR iteration, the crack
      ! needs to be determined first
      ! If there is no contact on the crack surface in the joint, the iteration ends immediately.
      !###########################################################################################
      if(i_NR_P==1)then 
          ! The initial value for the Newton-Raphson iteration is assumed to be a crack width of 0.
          NR_DISP(1:c_Total_FD)    = ZR
          delta_u_a(1:c_Total_FD)  = ZR
          ! NR_DISP(1:c_Total_FD) = Contact_DISP !Note: NR_DISP is the displacement involved in the
          ! Newton-Raphson contact iteration calculation
          !delta_u_a(1:c_Total_FD)  = ZR                  
          Last_R_PSI(1:c_Total_FD) = ZR
          
          !Elem_Conta_Sta(1:Num_Elem,1:Max_Num_Cr_3D) = 0  
          Elem_Conta_Sta(1:Num_Elem,1:num_Crack) = 0
          
          Kn_Gauss(1:num_Crack,1:Max_Max_N_FluEl_3D) = kn
          Kt1_Gauss(1:num_Crack,1:Max_Max_N_FluEl_3D) = kt
          Kt2_Gauss(1:num_Crack,1:Max_Max_N_FluEl_3D) = kt
          call Cal_Contact_Contact_State_Gauss_3D(  &
                 iter,ifra,Counter_Iter,i_NR_P, &
                  Contact_DISP,Yes_Contact,Elem_Conta_Sta, &
                  CT_State_Gauss)    
          ! If no embedding occurs on any of the crack surfaces, exit the iteration.
          if(Yes_Contact .eqv. .False.)then
              write(*,4012)
              CT_Jacobian = ori_globalK
              goto 9999
          endif
      endif

      !#############################################################################################
      ! Calculate the components of the contact force at the Gauss points of the contact element in
      ! the x, y, and z directions, and update and save the contact status at the Gauss points.
      !#############################################################################################
      print *,'           Get contact force and update Kt.'
      call Cal_Contact_PN_and_PT_3D(iter,ifra,Counter_Iter, &
                i_NR_P,c_Total_FD,c_num_freeDOF,  &
                Kn,Kn_Gauss,Kt1_Gauss,Kt2_Gauss, &
                fric_mu,NR_DISP,delta_u_a, &
                CT_State_Gauss,Elem_Conta_Sta, &
                PC_Gauss_x,PC_Gauss_y,PC_Gauss_z,8)
      write(*,5001) count(Elem_Conta_Sta==1),count(Elem_Conta_Sta==2)
      write(*,5002) count(CT_State_Gauss==1),count(CT_State_Gauss==2)
      
      !#############################################################################################
      ! Assemble the contact iteration Jacobian matrix, calculated based on ori_globalK (can reduce
      ! computation)
      !#############################################################################################
      print *,'           Assemble Jacobian matrix.'
      call Cal_Contact_Jacobian_3D(iter,ifra,Counter_Iter,i_NR_P, &
         c_Total_FD,c_num_freeDOF,ori_globalK, &
         c_freeDOF,Kn,Kn_Gauss,Kt1_Gauss,Kt2_Gauss, &
         CT_State_Gauss,CT_Jacobian)
      
      
      !##########################
      ! Calculate residuals, PSI
      !##########################
      print *,'           Get residual.'
      !------------------------------------------------------------
      ! option 1: Calculate internal forces through Bσ integration
      !------------------------------------------------------------
      call Cal_Contact_Resid_3D(iter,ifra,Counter_Iter,i_NR_P, &
         c_Total_FD,c_num_freeDOF,c_F,NR_DISP,c_freeDOF, &
         PC_Gauss_x,PC_Gauss_y,PC_Gauss_z,R_PSI) 
      write(*,6001) sum(abs(R_PSI))
      
      !############################################
      ! Solve for displacement increments (△u, △a)
      !############################################
      print *,'           Get displacemnt increament.'
      delta_u_a(1:c_Total_FD)  =ZR
      call Matrix_Solve_LSOE(7,1,Key_SLOE,     &
                  CT_Jacobian(c_freeDOF(1:c_num_freeDOF),c_freeDOF(1:c_num_freeDOF)), &
                 -R_PSI(c_freeDOF(1:c_num_freeDOF)), &
                  tem_DISP,c_num_freeDOF)
      delta_u_a(c_freeDOF(1:c_num_freeDOF)) = tem_DISP
      
      
      !######################
      ! Update displacement.
      !######################
      print *,'           Update displacemnt.'
      NR_DISP = NR_DISP + delta_u_a
      
      
      !########################
      ! Check for convergence.
      !########################
      print *,'           Check convergence.'
      call Cal_Contact_Conve_Factor(iter,ifra,Counter_Iter,i_NR_P,Conve_Tolerance, &
              c_Total_FD,c_freeDOF,c_num_freeDOF,c_F,R_PSI,Last_R_PSI, &
              delta_u_a,NR_DISP,Contact_DISP_0,Yes_Conve,Conve_Factor)
      write(*,4022) Conve_Factor,Conve_Tolerance
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
          ! Determine whether oscillation has occurred based on the last six convergence values
          print *,'           Check oscillation.'
          call Tool_Check_Oscillation_by_6_Variables(Saved_Conv_Factor(i_NR_P-5:i_NR_P),Yes_Oscill)
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
      !all Warning_Message('S',Keywords_Blank)
  endif
  
  ! Save output variable
  Contact_DISP = NR_DISP
  
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! case(2): Reduced Penalty function method (not implemented).
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
case(2)
  
end select


9999 continue


!////////////////////////////////
! Save data for post-processing.
!////////////////////////////////
if(Key_Save_Nothing ==0 .and. Key_Simple_Post==0 .and. Key_Contact==1) then
    !
    ! Save the fluid element (contact element) Gauss point contact force file for post-processing.
    !
    write(temp,'(I5)') iter
    c_File_name_1=trim(Full_Pathname)//'.cfrx_'//ADJUSTL(temp)     
    c_File_name_2=trim(Full_Pathname)//'.cfry_'//ADJUSTL(temp)    
    c_File_name_3=trim(Full_Pathname)//'.cfrz_'//ADJUSTL(temp)  
    open(801,file=c_File_name_1,status='unknown') 
    open(802,file=c_File_name_2,status='unknown') 
    open(803,file=c_File_name_3,status='unknown') 
    do i_C=1,num_Crack
        write(801, '(5000E20.12)') (PC_Gauss_x(i_C,j),j=1,Cracks_FluidEle_num_3D(i_C))       
        write(802, '(5000E20.12)') (PC_Gauss_y(i_C,j),j=1,Cracks_FluidEle_num_3D(i_C))      
        write(803, '(5000E20.12)') (PC_Gauss_z(i_C,j),j=1,Cracks_FluidEle_num_3D(i_C))         
    end do
    close(801)          
    close(802)    
    close(803)   
    !
    ! Save the fluid element (contact element) Gauss point contact state file for post-processing.
    ! 2023-09-23. NEWFTU2023092301.
    !
    write(temp,'(I5)') iter
    c_File_name_1=trim(Full_Pathname)//'.csce_'//ADJUSTL(temp)     
    open(901,file=c_File_name_1,status='unknown') 
    do i_C=1,num_Crack
        write(901, '(5000I3)') (CT_State_Gauss(i_C,j),j=1,Cracks_FluidEle_num_3D(i_C))    
    end do
    close(901)          
endif

RETURN
END SUBROUTINE Determine_Contact_State_by_Iteration_3D
