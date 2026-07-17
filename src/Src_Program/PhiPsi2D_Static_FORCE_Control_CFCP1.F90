!-----------------------------------------------------------
! Brief: 2D quasi-static driver using cohesive force-control algorithm CFCP1.
!
! Parameters:
!   Input:  none
!   Output: none
!   In/Out: none
!
! Notes:   Key_Force_Control==5; supports cross cracks and assumes at
!          least one initial crack. Iteratively adjusts the equivalent
!          K to drive cohesive crack growth; contact is temporarily off.
!-----------------------------------------------------------

SUBROUTINE PhiPsi2D_Static_FORCE_Control_CFCP1
! 2D quasi-static program and Key_Force_Control == 5 (Cohesive crack load control algorithm 1 
! (2010_MS_Numerical analysis of cohesive crack growth using XFEM.pdf_P31))
!CFCP=1
! Supports cross cracks
! Contains at least one crack
! Temporarily not available for contact


!-----------------------------          
! Read public variable module
!----------------------------- 
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
use Global_Inclusion
use Global_Plasticity
use Global_Cohesive
use Global_Cross

!------------------------------------------------          
! Read subroutine interface module (activate compiler parameter consistency check)
!------------------------------------------------ 
use Global_Inter_Matrix_Solve_LSOE
use Global_Inter_Determine_Contact_State_by_Iteration


!-----------------------------       
! Variable Type Declaration
!----------------------------- 
implicit none
integer igrow
logical Yes_Last_Growth
real(kind=FT) Lambda
integer,ALLOCATABLE::freeDOF(:)
real(kind=FT),ALLOCATABLE:: globalK(:,:),globalF(:)
real(kind=FT),ALLOCATABLE:: ori_globalK(:,:)
real(kind=FT),ALLOCATABLE:: Coh_Jacobian(:,:)
real(kind=FT),ALLOCATABLE:: Coh_DISP(:)
real(kind=FT),ALLOCATABLE:: Coh_R(:),Last_Coh_R(:)
real(kind=FT),ALLOCATABLE:: tem_DISP(:)
real(kind=FT),ALLOCATABLE:: Delta_U(:)
! real(kind=FT), ALLOCATABLE :: dsdePl(:,:,:,:)   ! D matrix for each element at each Gauss point
real(kind=FT),ALLOCATABLE::Cracks_CalP_Tan_Aper(:,:)
!------
integer num_FreeD
real(kind=FT) Max_Disp_x,Max_Disp_y,Min_Disp_x,Min_Disp_y
integer i_C,i_E,i_G,i_Tip
logical Yes_Growth(Max_Num_Cr,2)
integer ifra,Total_Num_G_P
real(kind=FT) Max_vm_Node,Max_vm_Gauss
integer i_Contact
logical Yes_Contact_Con
integer i,Result_File_Num
real(kind=FT) det_globalK,Max_Eq_Pla_Strain
logical Yes_Bou

integer iter,Counter_Iter
integer Max_Equ_Iteration
logical Yes_Conve,Yes_Plastic
real(kind=FT) Conve_Tolerance,Conve_Factor
real(kind=FT) Kequ_Expected,Kequ_Expected_Real
integer Action_Type
logical Yes_Delete_Coh_Ele
logical Yes_Delete_Coh_Ele_Last_Growth
real(kind=FT) Kequ(Max_Num_Cr,2)
integer tem_loc(2)
integer c_Propagate_Crack,c_Propagate_Tip
integer c_Tip_Elem,c_mat_num 
real(kind=FT) c_Cr_Tip(2)
logical Yes_Generated
real(kind=FT),ALLOCATABLE:: Coupled_Q(:,:)
real(kind=FT),ALLOCATABLE:: CalP_Pres(:)
integer num_Count,i_CalP      
real(kind=FT) a,b,c
real(kind=FT) Max_StrX_Node,Min_StrX_Node,Max_StrY_Node,Min_StrY_Node      
!2022-10-15.
real(kind=FT),ALLOCATABLE::storK(:,:,:),diag_precon_no_invert(:)
integer,ALLOCATABLE::size_local(:)
integer,ALLOCATABLE::all_local(:,:)    
! NR iteration related. 2022-10-17.
integer i_NR,Max_NR
logical Yes_NR_Convergent
real(kind=FT) Lambda_Initial_Guess,NR_delta_Lambda
real(kind=FT) c_Lambda,Applied_Lambda,Old_Lambda
real(kind=FT) f_K_Tol,f_K(2),d_f_K,Max_KI_eq
logical Tem_Condition_1,Tem_Condition_2
real(kind=FT) Conv_Lambda,Lambda_Tol    

! The following are test variables
real(kind=FT) ttt_Disp(2),ttt_x,ttt_y
real(kind=FT),ALLOCATABLE:: rrr(:)

!-----------------------------       
! Variable Initialization
!-----------------------------  
Yes_Last_Growth = .False.
Num_Frac = Num_Substeps


!-----------------------------        
! Formatted output statement
!-----------------------------  
1001 FORMAT(' >>  Propagation step ',I5,' of ',I5,' started:')   
! 1003 FORMAT(5X,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')   
! 1005 FORMAT('       FORCE loop ',I5,' of',I5,' of growth',I5,':')
2001 FORMAT('     Iteration step ',I5,' of ',I5,' started:') 
2002 FORMAT(5X,'PhiPsi failed to converge for plastic iteration!')   
1002 FORMAT('     Force factor is ',F5.3)   
1021 FORMAT(5X,'Range of displacement x:   ',F18.10,' m to ',F18.10,' m')
1022 FORMAT(5X,'Range of displacement y:   ',F18.10,' m to ',F18.10,' m')  
1131 FORMAT(5X,'KI and KII of crack ',I5,' tip 1 are ', F14.8,' and ',F14.8,' MPa.m^(1/2)')  
1132 FORMAT(5X,'KI and KII of crack ',I5,' tip 2 are ', F14.8,' and ',F14.8,' MPa.m^(1/2)') 
!C1331 FORMAT(7X,'KI and KII of crack ',I5,' tip 1 are ', F14.8,' and ',F14.8,' MPa.m^(1/2)')  
!C1332 FORMAT(7X,'KI and KII of crack ',I5,' tip 2 are ', F14.8,' and ',F14.8,' MPa.m^(1/2)') 
!C1131 FORMAT(5X,'KI and KII of crack ',I5,' tip 1 are ', F14.8,' and ',F14.8,' MPa.m^(1/2)')  
!C1132 FORMAT(5X,'KI and KII of crack ',I5,' tip 2 are ', F14.8,' and ',F14.8,' MPa.m^(1/2)') 
1041 FORMAT(5X,'Max Von-mises stress of all nodes:   ',E17.8,' MPa')
1042 FORMAT(5X,'Max Von-mises stress of Gauss points:   ',E17.8,' MPa')
4022 FORMAT(5X,'Convergence factor is', E16.7 ,'.')  
4032 FORMAT(5X,'NR iteration done after ',I3,' tries!') 
! 1023 FORMAT('     Sum of globalK: ',E15.8)
! 2001 FORMAT(5X,' +++  Contact iteration ',I3,' of ',I3,' started:')   
! 2002 FORMAT(5X,'PhiPsi failed to converge for contact iteration')   
1051 FORMAT(5X,'Max stress x of all nodes:   ',F15.5,' MPa')   
1052 FORMAT(5X,'Min stress x of all nodes:   ',F15.5,' MPa')   
1053 FORMAT(5X,'Max stress y of all nodes:   ',F15.5,' MPa')   
1054 FORMAT(5X,'Min stress y of all nodes:   ',F15.5,' MPa')           

4009 FORMAT('     >> N-R iteration step ',I5,' of ',I5,' started <<')   
4007 FORMAT(5X,'......................................................')   
4008 FORMAT(5X,'N-R iteration failed after ',I3,' tries!')  
5005 FORMAT(5x,'>> N-R iteration done after ',I3,' tries <<')
5000 FORMAT(5X,'Converged force factor lambda :',F11.3)    
5003 FORMAT(5X,'Converged Max_KIeq (MPa*m^1/2):',F11.3) 
5004 FORMAT(5X,'Material Para KI_c (MPa*m^1/2):',F11.3)  
! Total number of holes
num_Hole =  num_Hole +num_Circ_Hole + num_Ellip_Hole  

!-----------------------------        
! Check for water pressure inside the seam
!-----------------------------       
if(Key_Crack_Inner_Pressure==1) then
    print *,'    Inner crack pressure:  true'
else
    print *,'    Inner crack pressure:  false'
endif

!-----------------------------        
! Calculate the material matrix D for each material
!----------------------------- 
call Get_Material_Matrix

!---------------------------------------        
! Allocate memory for 8 displacement components of 4 nodes in each unit
!---------------------------------------   
ALLOCATE(Last_U_of_Ele(Num_Elem,8))
Last_U_of_Ele(1:Num_Elem,1:8) = ZR

!---------------------------------------        
! NR iteration-related parameter settings. NEWFTU2022101701.
!---------------------------------------
Max_NR = 10
Lambda_Initial_Guess =  1.0D0
NR_delta_Lambda      =  0.01D0
f_K_Tol              =  0.01D6
Lambda_Tol           =  1.0D-2

!------------------  
! Crack propagation step cycle
!------------------
Counter_Iter = 0
do igrow = 1,Num_Substeps
    print *, "  " 
    ! Initially, mark it as a non-XFEM analysis
    Yes_XFEM = .False.
    WRITE(*,1001) igrow,Num_Substeps
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !                                                         %
    ! If there are cracks, holes, or inclusions (XFEM)           %
    !                                                         %
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(num_Crack.ne.0 .or. num_Hole.ne.0  .or. num_Inclusion.ne.0 )then
        ! Marked as XFEM analysis, not FEM analysis
        Yes_XFEM = .True.
        !*****************************
        ! Confirm the enhancement node
        !*****************************
        ifra = igrow
        call Determine_Enriched_Nodes(ifra,igrow)

        !*****************************
        ! Assign numbers to enhancement nodes
        !*****************************
        call Number_Enriched_Nodes(igrow)     
        print *,'    Total_FD:',Total_FD  
        print *,'    Enrich_Freedom:',Enrich_Freedom  

        !****************************************
        ! Save crack-related files (including enhanced node numbering c_POS)
        !****************************************
        call Save_Files_Crack(igrow)
        ! If there are holes, save the coordinates of the holes
        if(num_Hole.ne.0)then
            call Save_Files_Holes(igrow)
        endif
        ! If it contains Cross, save the related files.
        if(num_Cross.ne.0)then
            call Save_Files_Cross(igrow)
        endif
        ! If inclusions are present, save the coordinates of the inclusions
        if(num_Inclusion.ne.0)then
            call Save_Files_Inclusions(igrow)
        endif

        !*****************************
        ! Calculate data related to crack calculation points
        !*****************************
        !call Cal_Crack_Points_Info(igrow)
        ! The following processing (related to water pressure) is aimed at obtaining calculation point coordinates consistent with hydraulic fracturing analysis, in order to examine the distribution of fracture openings.
        if(num_Crack.ne.0)then
            Cracks_HF_State(1:Max_Num_Cr) = 1

            call Stat_Crack_Connection(igrow)

            call Cal_HF_Crack_Points_Info_Linear(igrow)
        endif

        !*****************************
        ! Consider boundary conditions
        !*****************************
        ALLOCATE(freeDOF(Total_FD))
        call Boundary_Cond(Total_FD,igrow,freeDOF,num_FreeD)

        !*****************************
        ! Variable memory allocation
        !*****************************
        ALLOCATE(globalK(Total_FD,Total_FD))
        ALLOCATE(DISP(Total_FD))
        DISP =ZR 
        ALLOCATE(tem_DISP(num_FreeD))
        tem_DISP = ZR
        !****************************************
        ! Allocate other relevant dynamic data space
        !****************************************    
        ! ALLOCATE(DISP(Total_FD))                  !Displacements obtained in the current step, including active and fixed degrees of freedom
        ALLOCATE(Coupled_Q(Total_FD,num_Tol_CalP_Water))
        ! ALLOCATE(F(Total_FD))                     !External load vector considering fluid pressure
        ALLOCATE(CalP_Pres(num_Tol_CalP_Water))        
        ! ALLOCATE(delta_x(num_FreeD))    !Newton Raphson iteration displacement and pressure increment, note the dimensions: num_freeDOF   num_free_CalP

        !*******************************************
        ! Calculate the coupling matrix Q (if considering fluid pressure inside the seam)
        !*******************************************
        if(Key_Crack_Inner_Pressure==1) then
            call Cal_HF_Matrix_Q_Linear(igrow,Coupled_Q,Total_FD)     
        endif

        if(.not. allocated(globalF)) ALLOCATE(globalF(Total_FD))

        !*******************************************
        !
        ! NR Iteration
        !
        !*******************************************
        !NEWFTU2022101701.
        print *,'    Start NR iteration loop...' 
        Yes_NR_Convergent = .False.
        do i_NR = 1,Max_NR
            ! Screen output.
            print *,' '
            write(*,4009) i_NR,Max_NR

            ! Lambda_Initial_Guess =  1.0D0      !Initial load factor
            ! NR_delta_Lambda      =  0.01D0     !Lambda increment

            !......................................
            ! If it is the first rupture step and the first NR step
            !......................................
            if(igrow==1 .and. i_NR ==1) then
                c_Lambda = Lambda_Initial_Guess
            endif

            !##############################
            !                             #
            !                             #
            ! Initial load factor          #
            !                             #
            !                             #
            !##############################

            !.......................
            ! Current load factor.
            !.......................
            Applied_Lambda = c_Lambda
            print *,'    Applied_Lambda:',Applied_Lambda

            !.......................
            ! Cluster Payload Vector
            !.......................
            call Force_Vector(Total_FD,igrow,.False.,Applied_Lambda,globalF)

            !....................................................................
            ! If there is water pressure within the crack: convert the water pressure into a pressure vector at all calculation points (if considering fluid pressure within the crack)
            ! Note: The water pressure inside the seam is also applied proportionally. It can be changed to be applied fully. 2022-10-17.
            !....................................................................
            if(Key_Crack_Inner_Pressure==1) then
                num_Count = 0
                do i_C = 1,num_Crack                   
                    do i_CalP=1,Cracks_CalP_Num(i_C)   
                        num_Count = num_Count + 1
                        CalP_Pres(num_Count) = Crack_Pressure(i_C)*Applied_Lambda
                        Cracks_CalP_Pres(i_C,i_CalP)=Crack_Pressure(i_C)*Applied_Lambda
                    end do
                end do  
                globalF(1:Total_FD)= ZR
                call Force_Vector(Total_FD,igrow,.False.,ONE,globalF)
                globalF(freeDOF(1:num_FreeD))= globalF(freeDOF(1:num_FreeD))+ &
                MATMUL(Coupled_Q(freeDOF(1:num_FreeD),1:num_Tol_CalP_Water), CalP_Pres(1:num_Tol_CalP_Water))
            endif

            !.........................................................
            !       If solver 11: Element by Element, diagonally
            !       preconditioned conjugate gradient solver
            !.........................................................
            if(Key_SLOE==11)then
                !ALLOCATE(DISP(Total_FD))
                ! Only the first NR iteration step requires assembling the stiffness matrix.
                if(i_NR==1) then
                    ALLOCATE(storK(MDOF_2D,MDOF_2D,Num_Elem))
                    ALLOCATE(diag_precon_no_invert(0:num_FreeD))
                    ALLOCATE(size_local(Num_Elem))
                    ALLOCATE(all_local(MDOF_2D,Num_Elem))                  
                    call EBE_XFEM_PCG_Get_K(igrow,Applied_Lambda,num_FreeD,freeDOF, Total_Num_G_P,storK,size_local,all_local, &
                    diag_precon_no_invert)
                endif                  
                call EBE_XFEM_PCG_with_K(igrow,Applied_Lambda,cg_tol,max_cg, &
                num_FreeD,freeDOF(1:num_FreeD),globalF(freeDOF(1:num_FreeD)), &
                storK,size_local,all_local,diag_precon_no_invert(0:num_FreeD), DISP)
                !...............................
                ! Other solvers
                !...............................
            elseif(Key_SLOE/=11)then
                ! Assemble the global stiffness matrix and solve for the initial displacement, used for cohesive crack iteration
                call Assemble_Stiffness_Matrix_XFEM(igrow,globalK,Total_FD,Total_Num_G_P)
                print *,'    Sum of globalK: ',  sum(globalK) 
                print *,'    Solving displacements......' 
                call Matrix_Solve_LSOE(0,1,Key_SLOE,globalK(freeDOF(1:num_FreeD),  &
                freeDOF(1:num_FreeD)),globalF(freeDOF(1:num_FreeD)), tem_DISP,num_FreeD)
                DISP(1:Total_FD) =ZR
                DISP(freeDOF(1:num_FreeD)) = tem_DISP     
            endif

            !.......................
            ! Calculate crack width
            !.......................
            call Cal_Crack_Aperture(igrow,DISP)

            !.......................
            ! Calculate stress intensity factor
            !.......................
            if (Key_SIFs_Method==1) then
                call Cal_SIFs_DIM(igrow,DISP)              
            elseif (Key_SIFs_Method==2) then
                call Cal_SIFs_IIM(igrow,.True.,DISP)     
            end if
            do i_C = 1,num_Crack
                WRITE(*,1131)i_C,KI(i_C,1)/1.0D6,KII(i_C,1)/1.0D6
                WRITE(*,1132)i_C,KI(i_C,2)/1.0D6,KII(i_C,2)/1.0D6
            end do

            !..............................................
            ! Determine the crack number and crack tip number with the maximum equivalent stress intensity factor
            !..............................................
            ! Calculate the equivalent stress intensity factor
            Kequ(1:num_Crack,1:2) =ZR
            do i_C = 1,num_Crack
                do i_Tip =1,2
                    call Cal_Equivalent_K(i_C,i_Tip,KI(i_C,i_Tip),KII(i_C,i_Tip),Kequ(i_C,i_Tip))
                enddo
            enddo

            Max_KI_eq    = Maxval(Kequ(1:num_Crack,1:2))
            tem_loc(1:2) = MaxLoc(Kequ(1:num_Crack,1:2))
            c_Propagate_Crack = tem_loc(1)
            c_Propagate_Tip   = tem_loc(2)
            ! Specify the crack that can propagate (only one crack is allowed to propagate)
            Cracks_Allow_Propa(1:num_Crack)       = 0
            Cracks_Allow_Propa(c_Propagate_Crack) = 1



            !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            ! STEP 4: Determine Lambda_3 (based on KI) using linear interpolation so that the stress intensity factor equals the desired value
            !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            ! Obtain crack tip coordinates
            if(c_Propagate_Tip==1)then
                c_Cr_Tip(1:2)=Cr_First_Tip(c_Propagate_Crack,1:2)
            elseif(c_Propagate_Tip==2) then
                c_Cr_Tip(1:2)=Cr_Second_Tip(c_Propagate_Crack,1:2)
            endif
            ! Obtain the unit number where the crack tip is located
            call Cal_Ele_Num_by_Coors(c_Cr_Tip(1),c_Cr_Tip(2),c_Tip_Elem)
            ! Material number of the crack tip unit
            c_mat_num = Elem_Mat(c_Tip_Elem)
            ! Expected fracture toughness
            ! Kequ_Expected = Material_Para(c_mat_num,6)*2.0D0   ! Fracture toughness, note that a magnification factor needs to be set here!!!!!!!!!!!!
            Kequ_Expected_Real  = Material_Para(c_mat_num,6)
            !Kequ_Expected  = Material_Para(c_mat_num,6)*1.05D0   ! Fracture toughness. Note that you need to set a magnification factor here, otherwise the crack may not propagate!!!!!!!!!!!!
            Kequ_Expected  = Material_Para(c_mat_num,6)*1.0001D0


            ! Calculate f_K(1).
            f_K(1) = Max_KI_eq - Kequ_Expected

            !##############################
            !                             #
            !                             #
            ! Increase load factor        #
            !                             #
            !                             #
            !##############################

            !.......................
            ! Current load factor.
            !.......................
            Applied_Lambda = c_Lambda + NR_delta_Lambda
            print *,'    Applied_Lambda:',Applied_Lambda

            !.......................
            ! Cluster Payload Vector
            !.......................
            call Force_Vector(Total_FD,igrow,.False.,Applied_Lambda,globalF)

            !....................................................................
            ! If there is water pressure within the crack: convert the water pressure into a pressure vector at all calculation points (if considering fluid pressure within the crack)
            ! Note: The water pressure inside the seam is also applied proportionally. It can be changed to full application. 2022-10-17.
            !....................................................................
            if(Key_Crack_Inner_Pressure==1) then
                num_Count = 0
                do i_C = 1,num_Crack                   
                    do i_CalP=1,Cracks_CalP_Num(i_C)   
                        num_Count = num_Count + 1
                        CalP_Pres(num_Count) = Crack_Pressure(i_C)*Applied_Lambda
                        Cracks_CalP_Pres(i_C,i_CalP)=Crack_Pressure(i_C)*Applied_Lambda
                    end do
                end do  
                globalF(1:Total_FD)= ZR
                call Force_Vector(Total_FD,igrow,.False.,ONE,globalF)
                globalF(freeDOF(1:num_FreeD))= globalF(freeDOF(1:num_FreeD))+ &
                MATMUL(Coupled_Q(freeDOF(1:num_FreeD),1:num_Tol_CalP_Water), CalP_Pres(1:num_Tol_CalP_Water))
            endif

            !.........................................................
            !       If solver 11: Element by Element, diagonally
            !       preconditioned conjugate gradient solver
            !.........................................................
            if(Key_SLOE==11)then
                !ALLOCATE(DISP(Total_FD))
                !ALLOCATE(storK(MDOF_2D,MDOF_2D,Num_Elem))
                !ALLOCATE(diag_precon_no_invert(0:num_FreeD))
                !ALLOCATE(size_local(Num_Elem))
                !ALLOCATE(all_local(MDOF_2D,Num_Elem))                  
                !call EBE_XFEM_PCG_Get_K(igrow,Lambda,num_FreeD,freeDOF, &
                !                  Total_Num_G_P,storK,size_local,all_local, &
                !                  diag_precon_no_invert)
                call EBE_XFEM_PCG_with_K(igrow,Applied_Lambda,cg_tol,max_cg, &
                num_FreeD,freeDOF(1:num_FreeD),globalF(freeDOF(1:num_FreeD)), &
                storK,size_local,all_local,diag_precon_no_invert(0:num_FreeD), DISP)
                !...............................
                ! Other solvers
                !...............................
            elseif(Key_SLOE/=11)then
                ! Assemble the global stiffness matrix and solve for the initial displacement, used for cohesive crack iteration
                call Assemble_Stiffness_Matrix_XFEM(igrow,globalK,Total_FD,Total_Num_G_P)
                print *,'    Sum of globalK: ',  sum(globalK) 
                print *,'    Solving displacements......' 
                call Matrix_Solve_LSOE(0,1,Key_SLOE,globalK(freeDOF(1:num_FreeD),  &
                freeDOF(1:num_FreeD)),globalF(freeDOF(1:num_FreeD)), tem_DISP,num_FreeD)
                DISP(1:Total_FD) =ZR
                DISP(freeDOF(1:num_FreeD)) = tem_DISP     
            endif

            !.......................
            ! Calculate crack width
            !.......................
            call Cal_Crack_Aperture(igrow,DISP)

            !.......................
            ! Calculate stress intensity factor
            !.......................
            if (Key_SIFs_Method==1) then
                call Cal_SIFs_DIM(igrow,DISP)              
            elseif (Key_SIFs_Method==2) then
                call Cal_SIFs_IIM(igrow,.True.,DISP)     
            end if
            do i_C = 1,num_Crack
                WRITE(*,1131)i_C,KI(i_C,1)/1.0D6,KII(i_C,1)/1.0D6
                WRITE(*,1132)i_C,KI(i_C,2)/1.0D6,KII(i_C,2)/1.0D6
            end do

            !..............................................
            ! Determine the crack number and crack tip number with the maximum equivalent stress intensity factor
            !..............................................
            ! Calculate the equivalent stress intensity factor
            Kequ(1:num_Crack,1:2) =ZR
            do i_C = 1,num_Crack
                do i_Tip =1,2
                    call Cal_Equivalent_K(i_C,i_Tip,KI(i_C,i_Tip),KII(i_C,i_Tip),Kequ(i_C,i_Tip))
                enddo
            enddo


            ! Therefore, the maximum equivalent stress intensity
            Max_KI_eq    = Maxval(Kequ(1:num_Crack,1:2))



            ! Calculate f_K(1).
            f_K(2) = Max_KI_eq - Kequ_Expected

            !##############################
            !                             #
            !                             #
            ! Derivative Calculation and NR Update        #
            !                             #
            !                             #
            !##############################

            ! Calculate the derivative d_f_K of f_K.
            print *,'    Calculate the derivative of f_K...'
            d_f_K = (f_K(2)-f_K(1))/NR_delta_Lambda

            ! Update the time step according to the Newton-Raphson method.
            Old_Lambda = c_Lambda
            c_Lambda   = c_Lambda - f_K(1)/d_f_K
            print *,'    Updated Lambda:',c_Lambda

            ! Convergence check.
            Tem_Condition_1 = .False.
            Tem_Condition_2 = .False.
            if(abs(f_K(1))<= f_K_Tol) then                
                Tem_Condition_1 =  .True.
            endif
            print *, '    Abs(f_K(1)),f_K_Tol:', abs(f_K(1))/1.0D6,f_K_Tol/1.0D6

            Conv_Lambda = abs((c_Lambda-Old_Lambda)/Old_Lambda)
            print *, '    Conv_Lambda,Lambda_Tol:', Conv_Lambda,Lambda_Tol 

            if(Conv_Lambda <= Lambda_Tol) then                         
                Tem_Condition_2 =  .True.
            endif  
            if(Tem_Condition_1 .and. Tem_Condition_2) then
                print *,' '
                write(*,5005) i_NR
                !write(*,4007) 
                write(*,*) '    .........................................'
                Yes_NR_Convergent = .True.
                write(*,5000) c_Lambda
                write(*,5003) Max_KI_eq/1.0D6
                !write(*,5003) Ave_KI_eq_3D/1.0D6
                write(*,5004) Kequ_Expected_Real/1.0D6
                !write(*,4007) 
                write(*,*) '    .........................................'
                print *,' '
                ! Save pressure and time data to a *.wbpt file
                !call Save_HF_wbpt_File(i_WB,i_Stage,i_Prop,c_Pres,c_Time)
                exit
            endif
            if(i_NR ==Max_NR) then
                write(*,4007) 
                write(*,4008) Max_NR
                write(*,4007) 
                call Warning_Message('S',Keywords_Blank) 
            endif

        enddo


        !*******************************************
        !
        ! Results corresponding to the load factor after convergence calculation
        !
        !*******************************************

        !.......................
        ! Current load factor.
        !.......................
        Applied_Lambda = c_Lambda
        !Applied_Lambda = Applied_Lambda*0.98D0
        print *,'    Real applied_Lambda:',Applied_Lambda

        !.......................
        ! Cluster Payload Vector
        !.......................
        call Force_Vector(Total_FD,igrow,.False.,Applied_Lambda,globalF)

        !....................................................................
        ! If there is water pressure within the crack: convert the water pressure into a pressure vector at all calculation points (if considering fluid pressure within the crack)
        ! Note: The water pressure inside the seam is also applied proportionally. It can be changed to be applied fully. 2022-10-17.
        !....................................................................
        if(Key_Crack_Inner_Pressure==1) then
            num_Count = 0
            do i_C = 1,num_Crack                   
                do i_CalP=1,Cracks_CalP_Num(i_C)   
                    num_Count = num_Count + 1
                    CalP_Pres(num_Count) = Crack_Pressure(i_C)*Applied_Lambda
                    Cracks_CalP_Pres(i_C,i_CalP)=Crack_Pressure(i_C)*Applied_Lambda
                end do
            end do  
            globalF(1:Total_FD)= ZR
            call Force_Vector(Total_FD,igrow,.False.,ONE,globalF)
            globalF(freeDOF(1:num_FreeD))= globalF(freeDOF(1:num_FreeD))+ &
            MATMUL(Coupled_Q(freeDOF(1:num_FreeD),1:num_Tol_CalP_Water), CalP_Pres(1:num_Tol_CalP_Water))
        endif

        !.........................................................
        !       If solver 11: Element by Element, diagonally
        !       preconditioned conjugate gradient solver
        !.........................................................
        if(Key_SLOE==11)then
            !ALLOCATE(DISP(Total_FD))
            !ALLOCATE(storK(MDOF_2D,MDOF_2D,Num_Elem))
            !ALLOCATE(diag_precon_no_invert(0:num_FreeD))
            !ALLOCATE(size_local(Num_Elem))
            !ALLOCATE(all_local(MDOF_2D,Num_Elem))                  
            !call EBE_XFEM_PCG_Get_K(igrow,Lambda,num_FreeD,freeDOF, &
            !                  Total_Num_G_P,storK,size_local,all_local, &
            !                  diag_precon_no_invert)
            call EBE_XFEM_PCG_with_K(igrow,Applied_Lambda,cg_tol,max_cg, &
            num_FreeD,freeDOF(1:num_FreeD),globalF(freeDOF(1:num_FreeD)), &
            storK,size_local,all_local,diag_precon_no_invert(0:num_FreeD), DISP)
            !...............................
            ! Other solvers
            !...............................
        elseif(Key_SLOE/=11)then
            ! Assemble the global stiffness matrix and solve for the initial displacement, used for cohesive crack iteration
            call Assemble_Stiffness_Matrix_XFEM(igrow,globalK,Total_FD,Total_Num_G_P)
            print *,'    Sum of globalK: ',  sum(globalK) 
            print *,'    Solving displacements......' 
            call Matrix_Solve_LSOE(0,1,Key_SLOE,globalK(freeDOF(1:num_FreeD), freeDOF(1:num_FreeD)),globalF(freeDOF(1:num_FreeD)), &
            tem_DISP,num_FreeD)
            DISP(1:Total_FD) =ZR
            DISP(freeDOF(1:num_FreeD)) = tem_DISP     
        endif


        Max_Disp_x = maxval(DISP(1:2*Num_Node:2))
        Min_Disp_x = minval(DISP(1:2*Num_Node:2))
        Max_Disp_y = maxval(DISP(2:2*Num_Node:2))
        Min_Disp_y = minval(DISP(2:2*Num_Node:2))
        WRITE(*,1021) Min_Disp_x,Max_Disp_x
        WRITE(*,1022) Min_Disp_y,Max_Disp_y  
        !.......................
        ! Calculate crack width
        !.......................
        call Cal_Crack_Aperture(igrow,DISP)

        !.......................
        ! Calculate stress intensity factor
        !.......................
        if (Key_SIFs_Method==1) then
            call Cal_SIFs_DIM(igrow,DISP)              
        elseif (Key_SIFs_Method==2) then
            call Cal_SIFs_IIM(igrow,.True.,DISP)     
        end if
        do i_C = 1,num_Crack
            WRITE(*,1131)i_C,KI(i_C,1)/1.0D6,KII(i_C,1)/1.0D6
            WRITE(*,1132)i_C,KI(i_C,2)/1.0D6,KII(i_C,2)/1.0D6
        end do


        ! Clear memory variables
        !if(allocated(globalF)) DEALLOCATE(globalF)



        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !                                                 %
        ! If no cracks, holes, or inclusions (FEM)        %
        !                                                 %
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif (num_Crack.eq.0 .and. num_Hole.eq.0 .and. num_Inclusion.eq.0 )then 
        !*************************
        ! Total degrees of freedom
        !*************************
        Total_FD = 2*Num_Node
        print *,'    Total_FD:',Total_FD  
        
        !*************************
        ! Cluster payload
        !*************************
        if(.not. allocated(globalF)) ALLOCATE(globalF(Total_FD))
        call Force_Vector(Total_FD,igrow,.False.,Lambda,globalF)
        print *,'    Sum of globalF:',sum(globalF)  

        !*************************
        ! Consider boundary conditions
        !*************************
        ALLOCATE(freeDOF(Total_FD))
        call Boundary_Cond(Total_FD,igrow,freeDOF,num_FreeD)   

        ALLOCATE(DISP(Total_FD))
        DISP(1:Total_FD) = ZR

        !*************************
        ! Assembly stiffness matrix
        !*************************

        !*********************************************************
        !       If solver 11: Element by Element, diagonally
        !       preconditioned conjugate gradient solver
        !*********************************************************
        !NEWFTU2022101501.
        if(Key_SLOE==11)then                
            call Ele_by_Ele_FEM_PCG(igrow,Lambda,cg_tol,max_cg, &
            num_FreeD,freeDOF(1:num_FreeD), globalF(freeDOF(1:num_FreeD)),DISP,Total_Num_G_P)
            !******************************
            ! Other solvers
            !******************************
        elseif(Key_SLOE/=11)then
            print *,'    Assembling K......'    
            ALLOCATE(globalK(Total_FD,Total_FD))
            call Assemble_Stiffness_Matrix_FEM(igrow,globalK,Total_FD,Total_Num_G_P)
            !*************************
            ! Solve for displacement
            !*************************
            ALLOCATE(tem_DISP(num_FreeD))
            print *,'    Solving displacements......'           
            call Matrix_Solve_LSOE(0,1,Key_SLOE,globalK(freeDOF(1:num_FreeD), freeDOF(1:num_FreeD)),globalF(freeDOF(1:num_FreeD)), &
            tem_DISP,num_FreeD)

            DISP(freeDOF(1:num_FreeD)) = tem_DISP   
        endif

        Max_Disp_x = maxval(DISP(1:Total_FD:2))
        Min_Disp_x = minval(DISP(1:Total_FD:2))
        Max_Disp_y = maxval(DISP(2:Total_FD:2))
        Min_Disp_y = minval(DISP(2:Total_FD:2))
        WRITE(*,1021) Min_Disp_x,Max_Disp_x
        WRITE(*,1022) Min_Disp_y,Max_Disp_y
    end if


    !-------------
    ! Save displacement
    !-------------
    call Save_Disp(igrow,1)

    !*********************************
    ! Save VTK file, 2021-07-16
    !*********************************
    call Save_vtk_file(igrow)

    !***********************
    ! Save VTK CRACK file.
    !***********************
    call Save_vtk_file_for_Crack(igrow)

    !-----------------------------
    ! Save load-displacement curve file *.fdcu
    !-----------------------------
    if(Key_Save_f_d_Curve==1)then
        call Save_File_F_D_Curve(igrow,c_Lambda,f_d_Curve_node_num)
    endif

    !-----------------------------
    ! Save load-COD curve file *.fccu
    !-----------------------------
    !NEWFTU2022101502.
    if(Key_Save_f_COD_Curve==1)then
        call Save_File_F_COD_Curve(igrow,c_Lambda,f_COD_Curve_Crack_Num)
    endif

    !****************************************************
    ! If the contact condition of the crack surface is considered, save the element contact state elcs file.
    !**************************************************** 
    if((num_Crack.ne.0) .and. (Key_Contact/=0)) then 
        call Save_Ele_Contac_State(igrow)
    endif

    !***********************************************
    ! Save the corresponding crack-related files (including the enhanced node number c_POS)
    !*********************************************** 
    if(num_Crack.ne.0) then 
        call Save_Files_Crack(igrow)  
    endif    

    !***********************************************
    ! Save the crack width and other relevant information, including the coordinates of calculation points,
    ! Cohesive crack calculation point cohesive strength (2017-04-24), etc.
    !***********************************************
    if(num_Crack.ne.0) then 
        call Save_Files_Cr_CalP(igrow)   
    endif

    !-----------------------------------
    ! Calculate and save stress if necessary
    !----------------------------------- 
    if (Key_Post_CS_N_Strs ==1) then
        ! Allocate memory space for node stress
        ALLOCATE(Stress_xx_Node(num_Node))
        ALLOCATE(Stress_yy_Node(num_Node))
        ALLOCATE(Stress_xy_Node(num_Node))
        ALLOCATE(Stress_vm_Node(num_Node))
        ! If thermal stress is considered, allocate memory for nodal thermal stress.
        if(Key_Thermal_Stress==1)then
            ALLOCATE(TStress_xx_Node(num_Node))
            ALLOCATE(TStress_yy_Node(num_Node))
            ALLOCATE(TStress_xy_Node(num_Node))
            ALLOCATE(TStress_vm_Node(num_Node))
        endif
        ! Calculate the nodal stress and store it in the global variables Stress_xx_Node, Stress_yy_Node, Stress_xy_Node, Stress_vm_Node
        if(Yes_XFEM.eqv..True.)then
            call Get_Node_Stress_XFEM(igrow,DISP)
        else
            call Get_Node_Stress_FEM(igrow)
        endif
        ! Maximum Mises stress at the screen output node
        Max_vm_Node=maxval(Stress_vm_Node(1:num_Node))/1.0D6
        WRITE(*,1041) Max_vm_Node       
        ! Save node stress
        call Save_Stress_Node(igrow,1)
    end if

    !-----------------------------------------------------
    ! If needed, calculate and save the Gauss point coordinates for post-processing.
    !-----------------------------------------------------      
    if (Key_Post_CS_G_Coor==1) then 
        ALLOCATE(Gauss_CoorX(Total_Num_G_P))        
        ALLOCATE(Gauss_CoorY(Total_Num_G_P)) 
        ! Obtain Gauss point coordinates (Note: Actually, during the process of assembling the stiffness matrix
        ! The Gaussian point coordinates have already been indirectly calculated, recalculated here)
        call Cal_Gauss_Coors(igrow,Total_Num_G_P,Gauss_CoorX,Gauss_CoorY)         
        ! Save Gauss point coordinates
        call Save_Gauss_Coors(igrow,Total_Num_G_P,Gauss_CoorX,Gauss_CoorY)
    end if   

    !***************************************************
    ! If necessary, calculate and save the Gauss point displacements for post-processing.
    !***************************************************     
    if (Key_Post_CS_G_Disp==1) then 
        ALLOCATE(DISP_x_Gauss(Total_Num_G_P))        
        ALLOCATE(DISP_y_Gauss(Total_Num_G_P)) 
        ! Obtain Gauss point displacement
        call Get_Gauss_Disps(igrow,DISP,Total_Num_G_P)         
        ! Save Gauss point displacement
        call Save_Gauss_Disps(igrow,Total_Num_G_P)
    end if   

    !******************************************************
    ! If needed, calculate and save the Gauss point stress for post-processing display.
    !******************************************************     
    if (Key_Post_CS_G_Strs==1) then 
        ALLOCATE(Stress_xx_Gauss(Total_Num_G_P))        
        ALLOCATE(Stress_yy_Gauss(Total_Num_G_P))  
        ALLOCATE(Stress_xy_Gauss(Total_Num_G_P))  
        ALLOCATE(Stress_vm_Gauss(Total_Num_G_P))  
        ! If thermal stress is considered, allocate memory for thermal stress at Gauss points
        if(Key_Thermal_Stress==1)then
            ALLOCATE(TStress_xx_Gauss(Total_Num_G_P))        
            ALLOCATE(TStress_yy_Gauss(Total_Num_G_P))  
            ALLOCATE(TStress_xy_Gauss(Total_Num_G_P))  
            ALLOCATE(TStress_vm_Gauss(Total_Num_G_P))  
        endif
        ! Calculate the nodal stress and store it in the global variables Stress_xx_Node, Stress_yy_Node, Stress_xy_Node, Stress_vm_Node
        if(Yes_XFEM.eqv..True.)then
            call Get_Gauss_Stress_XFEM(igrow,DISP)
        else
            call Get_Gauss_Stress_FEM(igrow)
        endif
        ! Display the maximum Mises stress at Gauss points on the screen
        Max_vm_Gauss=maxval(Stress_vm_Gauss(1:Total_Num_G_P))/1.0D6
        WRITE(*,1042) Max_vm_Gauss 
        ! Save Gauss point stress
        call Save_Gauss_Stress(igrow,Total_Num_G_P)
    end if

    !***********************************************************************************
    ! Calculate and save shaped cracks for post-processing display. 
    ! Save to sccx, sccy and scdx, scdy files.
    ! sccx: store the x coodinates of shaped crack points. 
    !       Each line store 1 crack x cooridnate
    ! sccy: store the y coodinates of shaped crack points. 
    !       Each line store 1 crack y cooridnate
    ! sccx: store the x coodinates of shaped crack points. 
    !       Each line store 1 crack x cooridnate
    ! sccy: store the y coodinates of shaped crack points. 
    !       Each line store 1 crack y cooridnate
    ! NEWFTU-2025122602.
    !***********************************************************************************
    call Cal_Shaped_Cracks_2D(igrow,DISP)

    !****************************************************
    ! If the contact condition of the crack surface is considered, save the element contact state elcs file.
    !**************************************************** 
    if((num_Crack.ne.0) .and. (Key_Contact/=0)) then 
        call Save_Ele_Contac_State(igrow)
    end if

    !**************************************
    ! If the crack is allowed to propagate and is under maximum circumferential tensile stress
    ! According to the guidelines, calculate the stress intensity factor
    !************************************** 
    if(num_Crack.ne.0)then
        call Save_SIFs_KI_and_KII(igrow)
    endif

    ! Save the number of iterations for each step
    call Save_HF_time(1,igrow,igrow,0.0D0)   

    !---------------------------------------------------- 
    ! If the crack is allowed to extend, then determine whether it extends. If an extension occurs, calculate 
    ! the new crack tip coordinates, add new crack segments, and update crack coordinates
    !---------------------------------------------------- 
    if(num_Crack.ne.0)then
        if (Key_Propagation ==1) then
            ! Determine the direction and extent of crack propagation
            call Check_Crack_Grows(1,igrow,Yes_Growth) 

            if(igrow == Num_Substeps)then
                print *,'    |<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|'
                print *,'    |    All propagation steps done!  |'
                print *,'    |<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|'
            endif    
        end if
    endif

    !199   continue   

    !-----------------------------------------------------
    ! If necessary, calculate and save the tangential relative displacement at each crack calculation point.
    !-----------------------------------------------------   
    if (Key_Post_S_TanDisp==1)then
        print *,'    Cal and save tangent relative disp of cracks.'
        allocate(Cracks_CalP_Tan_Aper(num_Crack,Max_Num_Cr_CalP))
        call Cal_Crack_Tan_Relative_Disp(igrow,DISP,Cracks_CalP_Tan_Aper)
        call Save_Crack_Tan_Relative_Disp(igrow,Cracks_CalP_Tan_Aper)
        deallocate(Cracks_CalP_Tan_Aper)
    endif
    
        
    !-------------------------------------------------------------
    ! Calculate Energy. NEWFTU-2026062801.
    !-------------------------------------------------------------
    if (Key_Cal_Energy==1) then
        call Cal_Energy_2D(igrow, DISP, globalF)
        write(*,'(A,F16.8,A)')   '     Elastic Strain Energy  = ', Elastic_Strain_Energy,' J'
        write(*,'(A,F16.8,A)')   '     Fracture Energy        = ', Fracture_Energy,' J'
        write(*,'(A,F16.8,A)')   '     External_Work          = ', External_Work,' J'
        write(*,'(A,F16.8,A)')   '     Residual               = ', Residual_Energy,' J'
        write(*,'(A,F13.5,A)')   '     Normalized Residual    = ', Normalized_Residual_Energy,'%'
        call Save_Energy(igrow,Elastic_Strain_Energy,Fracture_Energy,External_Work,Residual_Energy,Normalized_Residual_Energy)
    endif

    !----------------
    ! Clear dynamic array
    !---------------- 
    if (allocated(Coupled_Q)) DEALLOCATE(Coupled_Q)
    if (allocated(CalP_Pres)) DEALLOCATE(CalP_Pres)            
    if(allocated(globalK)) DEALLOCATE(globalK)
    if(allocated(freeDOF)) DEALLOCATE(freeDOF)
    if(allocated(globalF)) DEALLOCATE(globalF)
    if(allocated(DISP)) DEALLOCATE(DISP) 
    if(allocated(tem_disp)) DEALLOCATE(tem_disp) 
    if((num_Crack.ne.0) .and. (Key_Contact/=0)) then 
        DEALLOCATE(Coh_DISP) 
    endif
    ! Clear the memory space for node stress
    if (Key_Post_CS_N_Strs ==1) then
        DEALLOCATE(Stress_xx_Node)
        DEALLOCATE(Stress_yy_Node)
        DEALLOCATE(Stress_xy_Node)
        DEALLOCATE(Stress_vm_Node)
        if(Key_Thermal_Stress==1)then
            DEALLOCATE(TStress_xx_Node)
            DEALLOCATE(TStress_yy_Node)
            DEALLOCATE(TStress_xy_Node)
            DEALLOCATE(TStress_vm_Node)
        endif
    end if     
    ! Clear Gauss point coordinate memory
    if (Key_Post_CS_G_Coor==1) then 
        DEALLOCATE(Gauss_CoorX)
        DEALLOCATE(Gauss_CoorY)          
    endif
    ! Clear Gauss point displacement memory
    if (Key_Post_CS_G_Disp==1) then 
        DEALLOCATE(DISP_x_Gauss)
        DEALLOCATE(DISP_y_Gauss)          
    endif
    ! Clear Gauss point stress memory space
    if (Key_Post_CS_G_Strs==1) then 
        DEALLOCATE(Stress_xx_Gauss)
        DEALLOCATE(Stress_yy_Gauss)      
        DEALLOCATE(Stress_xy_Gauss)
        DEALLOCATE(Stress_vm_Gauss) 
        if(Key_Thermal_Stress==1)then
            DEALLOCATE(TStress_xx_Gauss)
            DEALLOCATE(TStress_yy_Gauss)      
            DEALLOCATE(TStress_xy_Gauss)
            DEALLOCATE(TStress_vm_Gauss) 
        endif
    endif    
    if(Key_Contact/=0 .and. num_Crack.ne.0)then
        DEALLOCATE(ori_globalK)
    endif

    !2022-10-15. NEWFTU2022101501.
    if (allocated(storK)) DEALLOCATE(storK)    
    if (allocated(diag_precon_no_invert)) then 
        DEALLOCATE(diag_precon_no_invert)    
    endif
    if (allocated(size_local)) DEALLOCATE(size_local)    
    if (allocated(all_local)) DEALLOCATE(all_local)  


end do

!200 continue

!---------------------------------------       
! Release memory after the calculation to prevent unforeseen problems.
!--------------------------------------- 
if (Key_Contact /= 0) then
    DEALLOCATE(Elem_Conta_Sta)
    DEALLOCATE(Elem_Conta_Sta_Last)
    DEALLOCATE(Ele_NumCalP)
    DEALLOCATE(Ele_CalPNum)
end if


RETURN
END SUBROUTINE PhiPsi2D_Static_FORCE_Control_CFCP1


