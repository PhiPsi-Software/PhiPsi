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
 
SUBROUTINE D3_HF_Get_Pres_by_NR_Method(i_WB,i_Stage,i_Prop,                    &
              isub,i_Time,c_Stage_Q,c_Time,Max_Pres_Steps,                     &
              num_FreeD,c_Total_FD,c_num_Tol_CalP_Water,freeDOF,c_F_U,F,Lambda,&
              diag_precon_no_invert,DISP,                                      &
              NR_delta_Pres,f_Vol_Tol,Pres_Tol,c_Pres,Output_Pres,             &
              K_CSR_NNZ,K_CSR_aa,K_CSR_ja,K_CSR_ia)     
       

! 2022-07-01.    

!#############################
! Read public variable module
!#############################
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Global_XFEM_Elements

!#####################################################
! Read subroutine interface module (activate compiler
! parameter consistency check)
!#####################################################
use Global_EBE_XFEM_PCG_3D_with_K
use Global_EBE_Determine_Contact_State_by_Iteration_3D
use Global_D3_HF_Const_Pres_to_F_Vector

!###########################
! Variable Type Declaration
!###########################
implicit none
integer,intent(in)::i_WB,i_Stage,i_Prop,isub,i_Time,Max_Pres_Steps
integer,intent(in)::num_FreeD,c_Total_FD,c_num_Tol_CalP_Water
real(kind=FT),intent(in)::c_Stage_Q,c_Time
integer,intent(in)::freeDOF(num_FreeD)
real(kind=FT),intent(in)::diag_precon_no_invert(0:num_FreeD)
real(kind=FT),intent(in)::c_F_U(c_Total_FD)
!2025-12-02.
integer,intent(in)::K_CSR_NNZ
real(kind=FT),intent(in)::K_CSR_aa(K_CSR_NNZ)
integer,intent(in)::K_CSR_ja(K_CSR_NNZ)
integer,intent(in)::K_CSR_ia(num_FreeD+1)

real(kind=FT),intent(out)::F(c_Total_FD)
real(kind=FT),intent(in)::c_Pres,NR_delta_Pres
real(kind=FT),intent(in)::Lambda
real(kind=FT),intent(in)::f_Vol_Tol,Pres_Tol
real(kind=FT),intent(out)::Output_Pres
real(kind=FT),intent(out)::DISP(Total_FD)
!2023-09-26. NEWFTU-2025120201.
real(kind=FT),ALLOCATABLE:: F_each_NR(:,:)
real(kind=FT),ALLOCATABLE:: Output_Pres_each_NR(:)
real(kind=FT),ALLOCATABLE:: DISP_each_NR(:,:)
real(kind=FT),ALLOCATABLE:: Criterion_1_NR(:),Criterion_2_NR(:)
logical,ALLOCATABLE:: Flag_1_NR(:),Flag_2_NR(:)
integer,ALLOCATABLE::sort_INDEX(:)
integer i_Try,Current_Try_Index
real(kind=FT):: DISP_temp(num_FreeD)

integer i_Pres
real(kind=FT) f_Vol(2)
real(kind=FT) Old_Pres
! real(kind=FT) NR_delta_Pres, NR_delta_Time   ! Small increments in N-R iteration
real(kind=FT) Applied_Pres,d_f_Vol
logical Yes_P_Convergent
real(kind=FT) Max_Cr_Vol,Min_Cr_Vol
real(kind=FT) DISP2(Total_FD),tem_Pres
integer i_C
       
2001 FORMAT(8X,'Max volume of crack:   ',E12.5,' m^3')       
2002 FORMAT(8X,'Min volume of crack:   ',E12.5,' m^3')  
3003 FORMAT(8X,'Sum of abs(Force vector):   ',E16.4)      
4002 FORMAT(8X,'Pressure-step ',I3,' (',I3,') | Time-Step ',I3,' | Prop-Step ',I3,&
               ' | Stage ',I3,' | WB ',I3,' started...')   
3100 FORMAT(8X,'Mass consrevation value:      ',E12.5,' m^3 / ', E12.5,' m^3')    
3111 FORMAT(8X,'Mass consrevation value:      ',E12.5,' m^3 / ', E12.5,' m^3 (relaxed)')    
3101 FORMAT(8X,'Pressure convergence value:   ',E12.5,' / ', E12.5)   
3211 FORMAT(8X,'Pressure convergence value:   ',E12.5,' / ', E12.5,' (relaxed)')   
4006 FORMAT(8X,'N-R Pressure iteration converged after ',I3,' tries!')   
4007 FORMAT(8X,'..................................................')   
4008 FORMAT(8X,'N-R Pressure iteration failed after ',I3,' tries!')  
4009 FORMAT(8X,'N-R Pressure iteration failed after ',I3,' tries and check!')  
2004 FORMAT(5X,'Max aperture of crack ',I7,' is ',E12.5,' mm')  
2005 FORMAT(5X,'Min aperture of crack ',I7,' is ',E12.5,' mm')   
      
tem_Pres = c_Pres


!2023-09-26. IMPROV2023092601.
ALLOCATE(F_each_NR(Max_Pres_Steps,c_Total_FD))
ALLOCATE(Output_Pres_each_NR(Max_Pres_Steps))
ALLOCATE(DISP_each_NR(Max_Pres_Steps,Total_FD))
ALLOCATE(Criterion_1_NR(Max_Pres_Steps))
ALLOCATE(Criterion_2_NR(Max_Pres_Steps))
ALLOCATE(Flag_1_NR(Max_Pres_Steps))
ALLOCATE(Flag_2_NR(Max_Pres_Steps))
ALLOCATE(sort_INDEX(Max_Pres_Steps))
Flag_1_NR(1:Max_Pres_Steps) = .False.
Flag_2_NR(1:Max_Pres_Steps) = .False.

!######################
! Pressure step cycle.
!######################
do i_Pres = 1,Max_Pres_Steps
    ! Information output
    print *,' '
    write(*,4002) i_Pres,Max_Pres_Steps,i_Time,i_Prop, i_Stage,i_WB  
    
    ! If it's the first stress step.
    if(i_Pres==1)then
          Yes_P_Convergent = .False.
    endif
              
    !***********************************************************************
    !*                                                                    *
    !*                                                                    *
    !*     STEP 1: Apply normal pressure and calculate f_Vol(1).          *
    !*                                                                    *
    !*                                                                    *
    !***********************************************************************
    Applied_Pres =  tem_Pres
    ! Applied_Pres = tem_Pres - NR_delta_Pres ! Central difference calculation of derivative.
    ! 2022-08-12. IMPROV2022081201.
    
    
    !///////////////////////////////////////////////////////////////////////////////////////////////
    ! Convert the water pressure into a pressure vector for all calculation points and output the F
    ! vector.
    !///////////////////////////////////////////////////////////////////////////////////////////////
    call D3_HF_Const_Pres_to_F_Vector(isub,Applied_Pres,num_FreeD,Total_FD,num_Tol_CalP_Water,  &
                                      freeDOF(1:num_FreeD),c_F_U,F)     
    WRITE(*,3003) sum(abs(F))
    
    !2023-09-26.
    F_each_NR(i_Pres,:) = F
    
    !/////////////////////////////////////
    ! Ignoring contact: EBE-PCG solution.
    !/////////////////////////////////////
    if(Key_Contact == 0 .or. Key_Contact==5) then
        if(Key_SLOE==11)then
            print *,'       PCG-EBE: solving KU=F...'
            call EBE_XFEM_PCG_3D_with_K(isub,Lambda,cg_tol,                                &
                    max_cg,num_FreeD,freeDOF(1:num_FreeD),F(freeDOF(1:num_FreeD)),DISP,    &
                    diag_precon_no_invert(0:num_FreeD))
        !NEWFTU-2025120201.
        else
            print *,'       Solving KU=F...'
#ifndef Silverfrost
            call Matrix_Solve_LSOE_Sparse(1,1,Key_SLOE,&
                    K_CSR_NNZ,&
                    K_CSR_aa(1:K_CSR_NNZ),&
                    K_CSR_ja(1:K_CSR_NNZ),&
                    K_CSR_ia(1:num_FreeD+1),&
                    F(freeDOF(1:num_FreeD)),&
                    DISP_temp(1:num_FreeD),num_FreeD)
            DISP(freeDOF(1:num_FreeD)) =DISP_temp(1:num_FreeD)      
#endif
        endif
        
    endif
    
    !/////////////////////////////////////////////////////////////////////////
    ! Consider contact: EBE-PCG contact iterative solution. NEWFTU2022081001.
    !/////////////////////////////////////////////////////////////////////////
    if(Key_Contact/=0 .and. Key_Contact/=5)then
        if(Key_SLOE==11)then
            print *,'       PCG-EBE: solving KU=F...'
            call EBE_XFEM_PCG_3D_with_K(isub,Lambda,cg_tol,                                &
                max_cg,num_FreeD,freeDOF(1:num_FreeD),F(freeDOF(1:num_FreeD)),DISP,    &
                diag_precon_no_invert(0:num_FreeD))
        !NEWFTU-2025120201.
        else
            print *,'       Solving KU=F...'
#ifndef Silverfrost
            call Matrix_Solve_LSOE_Sparse(1,1,Key_SLOE,&
                    K_CSR_NNZ,&
                    K_CSR_aa(1:K_CSR_NNZ),&
                    K_CSR_ja(1:K_CSR_NNZ),&
                    K_CSR_ia(1:num_FreeD+1),&
                    F(freeDOF(1:num_FreeD)),&
                    DISP_temp(1:num_FreeD),num_FreeD)
            DISP(freeDOF(1:num_FreeD)) =DISP_temp(1:num_FreeD)   
#endif
        endif
        
        ! Contact iteration
        print *,'       Determine contact states by iteration...' 
        if(Key_SLOE==11)then 
            call EBE_Determine_Contact_State_by_Iteration_3D(isub,cg_tol,max_cg,&
                   num_FreeD,freeDOF(1:num_FreeD),                              &      
                   F,DISP,                                                      &
                   diag_precon_no_invert(0:num_FreeD),8)
        else
            !TBD
        endif
    endif 
    
    !2023-09-26.
    DISP_each_NR(i_Pres,:) = DISP

    !///////////////////////////////////
    ! Calculate crack width and volume.
    !///////////////////////////////////
    print *,'       Calculating crack apertures...'  
    Cracks_Volume_Old  = Cracks_Volume
    call Cal_Crack_Aperture_3D(isub,DISP)
    Max_Cr_Vol = maxval(Cracks_Volume(1:num_Crack))
    Min_Cr_Vol = minval(Cracks_Volume(1:num_Crack))
    WRITE(*,2001) Max_Cr_Vol
    WRITE(*,2002) Min_Cr_Vol
    
    !/////////////////////
    ! Calculate f_Vol(1).
    !/////////////////////
    print *,'       Calculate f_Vol(p)...'
    f_Vol(1) = ZR
    do i_C = 1,num_Crack
          if(Crack_Type_Status_3D(i_C,1)==1) then
              f_Vol(1) = f_Vol(1) + Cracks_Volume(i_C)
          endif
    enddo
    f_Vol(1) = f_Vol(1) - c_Stage_Q*c_Time           
    

    !******************************************************
    !*                                                   *
    !*                                                   *
    !*   STEP 2: Increase the pressure by a small amount *
    !*           and calculate f_Vol(2).                 *
    !*                                                   *
    !******************************************************
    Applied_Pres =  tem_Pres + NR_delta_Pres
    
    !/////////////////////////////////////////////////////
    ! Convert the water pressure into a pressure vector
    ! for all calculation points and output the F vector.
    !/////////////////////////////////////////////////////
    call D3_HF_Const_Pres_to_F_Vector(isub,Applied_Pres,num_FreeD,Total_FD,num_Tol_CalP_Water,  &
                                freeDOF(1:num_FreeD),c_F_U,F)     
    WRITE(*,3003) sum(abs(F))
    
    !/////////////////////////////////////
    ! Ignoring contact: EBE-PCG solution.
    !/////////////////////////////////////
    if(Key_Contact == 0 .or. Key_Contact==5) then
        if(Key_SLOE==11)then
            print *,'       PCG-EBE: solving KU=F...'
            call EBE_XFEM_PCG_3D_with_K(isub,Lambda,cg_tol,                              &
                max_cg,num_FreeD,freeDOF(1:num_FreeD),F(freeDOF(1:num_FreeD)),DISP2, & 
                diag_precon_no_invert(0:num_FreeD))
        !NEWFTU-2025120201.
        else
            print *,'       Solving KU=F...'
#ifndef Silverfrost
            call Matrix_Solve_LSOE_Sparse(1,1,Key_SLOE,&
                    K_CSR_NNZ,&
                    K_CSR_aa(1:K_CSR_NNZ),&
                    K_CSR_ja(1:K_CSR_NNZ),&
                    K_CSR_ia(1:num_FreeD+1),&
                    F(freeDOF(1:num_FreeD)),&
                    DISP_temp(1:num_FreeD),num_FreeD)
            DISP2(freeDOF(1:num_FreeD)) =DISP_temp(1:num_FreeD)   
#endif
        endif
    endif
    
    !/////////////////////////////////////////////////////////////////////////
    ! Consider contact: EBE-PCG contact iterative solution. NEWFTU2022081001.
    !/////////////////////////////////////////////////////////////////////////
    !if(Key_Contact /= 0) then
    if(Key_Contact/=0 .and. Key_Contact/=5)then 
        if(Key_SLOE==11)then 
            print *,'       PCG-EBE: solving KU=F...'
            call EBE_XFEM_PCG_3D_with_K(isub,Lambda,cg_tol,                              &
                max_cg,num_FreeD,freeDOF(1:num_FreeD),F(freeDOF(1:num_FreeD)),DISP2, & 
                diag_precon_no_invert(0:num_FreeD))
        !NEWFTU-2025120201.
        else
            print *,'       Solving KU=F...'
#ifndef Silverfrost
            call Matrix_Solve_LSOE_Sparse(1,1,Key_SLOE,&
                    K_CSR_NNZ,&
                    K_CSR_aa(1:K_CSR_NNZ),&
                    K_CSR_ja(1:K_CSR_NNZ),&
                    K_CSR_ia(1:num_FreeD+1),&
                    F(freeDOF(1:num_FreeD)),&
                    DISP_temp(1:num_FreeD),num_FreeD)
            DISP2(freeDOF(1:num_FreeD)) =DISP_temp(1:num_FreeD)   
#endif
        endif
        ! Contact iteration.
        if(Key_SLOE==11)then 
            call EBE_Determine_Contact_State_by_Iteration_3D(isub,cg_tol,max_cg,&
               num_FreeD,freeDOF(1:num_FreeD),                              &      
               F,DISP2,                                                     &
               diag_precon_no_invert(0:num_FreeD),8)
        
        else

        endif
    endif
    
    !///////////////////////////////////
    ! Calculate crack width and volume.
    !///////////////////////////////////
    print *,'       Calculating crack apertures...'  
    Cracks_Volume_Old  = Cracks_Volume
    call Cal_Crack_Aperture_3D(isub,DISP2)
    Max_Cr_Vol = maxval(Cracks_Volume(1:num_Crack))
    Min_Cr_Vol = minval(Cracks_Volume(1:num_Crack))
    WRITE(*,2001) Max_Cr_Vol
    WRITE(*,2002) Min_Cr_Vol
    
    !//////////////////////////////////////////////////////////////
    ! If the crack volume is negative, adjust it to 0. 2022-08-05.
    !//////////////////////////////////////////////////////////////
    ! do i_C = 1,num_Crack
    !    if (Cracks_Volume(i_C)<=ZR) then
    !        Cracks_Volume(i_C)=ZR
    !    endif
    ! enddo
    
    !/////////////////////
    ! Calculate f_Vol(2).
    !/////////////////////
    print *,'       Calculate f_Vol(p+delta_p)...'
    f_Vol(2) = ZR
    do i_C = 1,num_Crack
          if(Crack_Type_Status_3D(i_C,1)==1) then
              f_Vol(2) = f_Vol(2) + Cracks_Volume(i_C)
          endif
    enddo
    f_Vol(2) = f_Vol(2) - c_Stage_Q*c_Time               
    
    !*****************************************************
    !*                                                  *
    !*                                                  *
    !* STEP3: Calculate the derivative of f_Vol,d_f_Vol.*
    !*                                                  *
    !*                                                  *
    !*****************************************************
    print *,'       Calculate the derivative of f_Vol...'
    d_f_Vol = (f_Vol(2)-f_Vol(1))/NR_delta_Pres
    
    ! d_f_Vol = (f_Vol(2) - f_Vol(1)) / (TWO * NR_delta_Pres) ! Central difference calculation of
    ! derivative. 2022-08-12. IMPROV2022081201.
    
    !*****************************************************
    !*                                                  *
    !*                                                  *
    !* STEP4: Update the pressure according to the      *
    !*        Newton-Raphson method.                    *
    !*                                                  *
    !*****************************************************
    Old_Pres = tem_Pres
    tem_Pres = tem_Pres - f_Vol(1)/d_f_Vol
    
    Output_Pres_each_NR(i_Pres) = tem_Pres
    
    ! Explanation: For calculating derivatives using central differences, it is necessary to compute the
    ! crack opening at position f_0, which increases the computational workload.
    !      To be done
    
    print *,'       Updated pressure (MPa):',tem_Pres/1.0D6
    !*****************************************************
    !*                                                  *
    !*                                                  *
    !*              STEP5: Convergence Check.           *
    !*                                                  *
    !*                                                  *
    !*****************************************************
    !Tem_Condition_1 = .False.
    !Tem_Condition_2 = .False.
    
    Criterion_1_NR(i_Pres) = abs(f_Vol(1))
    
    if(Criterion_1_NR(i_Pres) <= f_Vol_Tol) then
        !Tem_Condition_1 =  .True.
        Flag_1_NR(i_Pres) = .True.
    endif
    write(*,3100)  Criterion_1_NR(i_Pres) ,f_Vol_Tol
    
    
    Criterion_2_NR(i_Pres) = abs((tem_Pres-Old_Pres)/Old_Pres)
    
    if(Criterion_2_NR(i_Pres) <= Pres_Tol) then                            
        !Tem_Condition_2 =  .True.
        Flag_2_NR(i_Pres) = .True.
    endif  
    write(*,3101)  Criterion_2_NR(i_Pres),Pres_Tol
    
    ! If both conditions are met.
    if((Flag_1_NR(i_Pres).eqv. .True.) .and. (Flag_2_NR(i_Pres).eqv. .True.)) then   
          write(*,*) '       ..................................................'
          write(*,4006) i_Pres
          write(*,*) '       ..................................................'
          exit
    endif

    ! ==========
    ! OPTION 1.   
    ! ==========
    !    If the calculation still does not converge after all steps are completed, 
    !    an error message will be displayed and the computation will be terminated.
    !    if(i_Pres == Max_Pres_Steps) then
    !          write(*,4007) 
    !          write(*,4008) Max_Pres_Steps
    !          write(*,4007) 
    !          call Warning_Message('S',Keywords_Blank) 
    !    endif
    
    ! ==========
    ! OPTION 2.
    ! ==========
    !2023-09-26. IMPROV2023092601.
    ! If after completing all steps it still does not converge, then find the best result from the
    ! previous calculations. 2023-09-26.
    if(i_Pres == Max_Pres_Steps .and. SlipWater_Pres_Step_Conv_Check==1) then
          ! Sort Criterion_1_NR(:) from small to large.
          call Vector_Sort_Dou_with_Index_v_HSL_KB07(Criterion_1_NR(1:Max_Pres_Steps),Max_Pres_Steps,sort_INDEX(1:Max_Pres_Steps))
          ! Find the best result that meets the conditions from small to large.
          do i_Try=1,Max_Pres_Steps
            Current_Try_Index = sort_INDEX(i_Try)
            ! If the criterion 2 corresponding to the minimum value of Criterion_1_NR(:) is satisfied.
            if(Flag_2_NR(Current_Try_Index) .eqv. .True.)then
                F           = F_each_NR(Current_Try_Index,:)
                Output_Pres = Output_Pres_each_NR(Current_Try_Index)
                DISP        = DISP_each_NR(Current_Try_Index,:)
                exit
            endif
            
            if (i_Try == Max_Pres_Steps) then
                ! Pop up an error message.
                write(*,4007) 
                write(*,4009) Max_Pres_Steps
                write(*,4007) 
                call Warning_Message('S',Keywords_Blank) 
            endif
            
          enddo
    endif
  
enddo

Output_Pres = tem_Pres

! 2024-02-27. BUGFIX2024022701. Added the following code to update the hydraulic fracturing crack
! internal water pressure.
do i_C=1,num_Crack
    ! If it is an HF crack.
    if(Crack_Type_Status_3D(i_C,1) == 1) then 
        Crack_Pressure(i_C) = Output_Pres
    endif
enddo


!2023-09-26. IMPROV2023092601.
DEALLOCATE(F_each_NR)
DEALLOCATE(Output_Pres_each_NR)
DEALLOCATE(DISP_each_NR)
DEALLOCATE(Criterion_1_NR)
DEALLOCATE(Criterion_2_NR)
DEALLOCATE(Flag_1_NR)
DEALLOCATE(Flag_2_NR)
DEALLOCATE(sort_INDEX)

RETURN
END SUBROUTINE D3_HF_Get_Pres_by_NR_Method
