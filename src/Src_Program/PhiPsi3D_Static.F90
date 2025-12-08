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
 
SUBROUTINE PhiPsi3D_Static
! Three-dimensional static analysis.

!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common   
use Global_Filename
use Global_Model
use Global_XFEM_Elements
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
use Global_DISP
use Global_HF
use Global_Stress
use Global_Strain
use Global_POST
use Global_Ragged_Array_Real_Classs
      
!----------------------------------------------------------------------------------
! Read subroutine interface module (activate compiler parameter consistency check)
!----------------------------------------------------------------------------------
use Global_Inter_Matrix_Solve_LSOE
use Global_Ele_by_Ele_XFEM_PCG_3D
use Global_EBE_Determine_Contact_State_by_Iteration_3D

!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
integer isub
logical Yes_Last_Growth
real(kind=FT) Lambda
integer,ALLOCATABLE::freeDOF(:),fixedDOF(:)
real(kind=FT),ALLOCATABLE:: globalK(:,:)
real(kind=FT),ALLOCATABLE:: ori_globalK(:,:)
integer Total_Num_G_P
real(kind=FT),ALLOCATABLE::tem_DISP(:)
integer num_FreeD
integer i_C
logical Yes_Growth(Max_Num_Cr_3D)
integer ifra
integer num_FixedD
real(kind=FT),ALLOCATABLE:: F_U(:),F(:)
real(kind=FT),ALLOCATABLE:: CalP_Pres(:)
real(kind=FT),ALLOCATABLE:: delta_x(:)
integer date_time(8)
integer(LIT) F_time
character*10  current_data
character(200) c_File_name_1
real(kind=FT) Max_Cr_Vol,Min_Cr_Vol
integer i_CalP,num_Count
real(kind=FT),ALLOCATABLE:: Contact_DISP(:)
real(kind=FT),ALLOCATABLE:: Vector_Dc(:)
integer i_FD
logical New_Crack_Flag
real(kind=FT),ALLOCATABLE::diag_precon(:) 
logical Inactive_Key_Multi_TipEnrNode
integer i_WB,i_Stage,i_Prop
integer i_Dof      
integer,ALLOCATABLE::tem_all_local(:,:)
real(kind=FT),ALLOCATABLE:: K_CSR_aa(:)
integer,ALLOCATABLE:: K_CSR_ja(:)
integer,ALLOCATABLE:: K_CSR_ia(:)
integer(kind=LIT) K_CSR_NNZ_Max
integer K_CSR_NNZ
      
!-------------------------------------
! Initial value of temporary variable
!-------------------------------------
Yes_Last_Growth = .False.
New_Crack_Flag = .False.

!----------------------------
! Formatted output statement
!----------------------------
1001 FORMAT('  >> Substep ',I5,' of ',I5,' started:')   
1002 FORMAT('     Force factor is ',F5.3)   
2001 FORMAT(5X,'Max volume of crack:   ',E13.6,' m^3')       
2002 FORMAT(5X,'Min volume of crack:   ',E13.6,' m^3')  
2003 FORMAT(5X,'Volume of crack:   ',E13.6,' m^3') 
2004 FORMAT(5X,'Max aperture of crack ',I7,' is ',E13.6,' mm')  
2005 FORMAT(5X,'Min aperture of crack ',I7,' is ',E13.6,' mm')   
3001 FORMAT(5X,'Elapsed CPU time: ',I7,' s, about ',F7.2,' mins')  
3002 FORMAT(5X,'Sum of global K:   ',E16.5)      
3003 FORMAT(5X,'Sum of abs(Force vector):   ',E16.4)      
4001 FORMAT(5X,'Number of crack is ',I3,' / ',I5)   
1171 FORMAT('  >> Substep ',I3,' of stage ',I3,' /',I3,' of WB ',I3,' /',I3,' started:')   
1172 FORMAT('     Step count:',I3) 

!------------------------------------------
! Check for water pressure inside the seam
!------------------------------------------
if(Key_Crack_Inner_Pressure==1) then
  print *,'    Inner crack pressure:  true'
else
  print *,'    Inner crack pressure:  false'
endif

!---------------------------------------------------
! Calculate the material matrix D for each material
!---------------------------------------------------
call Get_Material_Matrix_3D
!If has composite material.
if (any(Material_Type == 5))then   
  !Get the rotation matrix for the composite material for each element.
  call Get_Ele_Composite_Mat_Rot_Matrix 
endif
      
!-------------------------------------------------------
! Generate initial natural fractures. NEWFTU2022061001.
!-------------------------------------------------------
if (Key_Random_NaCr==1) then 
    call Tool_Generate_Natural_Fractures_3D
endif

!----------------------------------------------------------------------------------------
! Set natural fractures based on the data in the kpp file. 2023-08-24. NEWFTU2023082401.
!----------------------------------------------------------------------------------------
if (Key_Random_NaCr==3) then 
  call Tool_Set_Natural_Fractures_by_kpp_3D
endif 
      
!-----------------------------------------------------------------
! Save wellbore coordinate data (if available). NEWFTU2022110901.
!-----------------------------------------------------------------
if(num_Wellbore>=1)then
    print *, "    Saving *.wbfp file of wellbore..." 
    call Save_HF_wbfp_File
endif

!------------------------------------------------------------------------------------------------
! Calculate the traditional degree-of-freedom displacement field under in-situ stress (far-field
! stress), and calculate the in-situ stress level, then save it.
! Node Stress and Gauss Point Stress
! Note: There are no enhancement nodes in the model at this time
!------------------------------------------------------------------------------------------------
if(Key_InSitu_Strategy /=0)then
  ! Key_InSitu_Strategy=4 does not require calculating the initial stress field, but is directly
  ! defined or read from a file, 2022-06-03.
  if(Key_InSitu_Strategy /=4)then
      call Cal_InSitu_Stress_3D
  endif
endif

!---------------------------------------------------------------------
! Initial ground stress treatment. Key_InSitu_Strategy=4. 2022-06-03.
! Ref: Cruz_2019_An XFEM implementation in Abaqus to model 
! intersections between fractures in porous rocks.pdf, Section 3.5
! NEWFTU2022060301.
!---------------------------------------------------------------------
if(Key_InSitu_Strategy ==4)then
  ! Obtain internal force (stress at each conventional Gauss point): obtained according to user custom
  ! settings or input file.
  ! Store variables such as InSitu_Strs_Gaus_xx(num_Elem, Num_Gauss_P_FEM_3D).
  print *,'    Getting internal stress/strain of Gauss points.'  
  call D3_Get_Internal_Stress_and_Strain
endif
      
!----------------------------------------------------------------------------------------------
! Make random minor adjustments to each initial crack. Note that these adjustments are random.
!----------------------------------------------------------------------------------------------
if(num_Crack.ne.0 .and. Key_Adj_Ini_Crack_3D==1)then
  print *,'    Adjust initial 3D cracks...'  
  call D3_Adjust_Initial_Crack
endif

!---------------------------------------------------------------------------
! Segment each initial crack (generate triangular geometric crack elements)
!---------------------------------------------------------------------------
if(num_Crack.ne.0)then
  print *,'    Meshing initial 3D cracks...'  
  call D3_Mesh_Initial_Crack
endif

!----------------------------------------------------
! Inspect and adjust each initial crack. 2022-08-01.
!----------------------------------------------------
if(num_Crack.ne.0 .and. Key_Check_and_Adjust_Cracks_3D>=1)then
  print *,'    Check and adjust initial 3D cracks...'  
  call D3_Check_and_Adjust_Cracks
endif    

!----------------------------------------------------------------------------------------------------
! Logic setting: If it is not a staged hydraulic fracturing analysis, set the corresponding cycle to
! 1, 2022-06-02.
!----------------------------------------------------------------------------------------------------
if(Key_HF_Multistage_3D==0) then
    num_Wellbore=1
    num_Stages_Wellbores(1) = 1      
    if(Key_Crack_Inner_Pressure==1)then
        do i_C =1,num_Crack
              if (abs(Crack_Pressure(i_C)) > Tol_10) then
                  Flag_HF_3D = 1
                  Crack_Type_Status_3D(i_C,1) = 1
              else
                  Crack_Type_Status_3D(i_C,1) = 0 
              endif
        enddo
    endif
elseif(Key_HF_Multistage_3D==1) then
    Flag_HF_3D = 1
endif
      
!------------------------------        
! Cycle between each load step
!------------------------------  
isub = 0
do i_WB = 1,num_Wellbore
    do i_Stage = 1,num_Stages_Wellbores(i_WB) 
        do i_Prop = 1,Num_Substeps
            print *, "  " 
            isub = isub + 1
            if(Key_HF_Multistage_3D==0) then
              WRITE(*,1001) i_Prop,Num_Substeps
            elseif (Key_HF_Multistage_3D==1) then
              WRITE(*,1171) i_Prop,i_Stage,num_Stages_Wellbores(i_WB),i_WB,num_Wellbore
              WRITE(*,1172) isub
            endif          

        !*****************************************
        ! If it is staged fracturing, 2022-05-31.
        !*****************************************
        if(Key_HF_Multistage_3D==1) then
            ! Generate initial fractures along the wellbore (WB).
            call D3_HF_Generate_Initial_Cracks_of_WB(i_Prop,i_WB,i_Stage)
                
            ! Make slight adjustments to each initial crack.
            if(num_Crack.ne.0 .and. Key_Adj_Ini_Crack_3D==1)then
                print *,'    Adjust initial 3D cracks...'  
                call D3_Adjust_Initial_Crack
            endif
                     
            ! Segment each initial crack (to generate triangular geometric crack elements), only for cracks that
            ! have not been discretized yet.
            if(num_Crack.ne.0)then
                print *,'    Meshing initial 3D cracks...'  
                call D3_Mesh_Initial_Crack
            endif
        endif

          
        !*********************************************************************
        ! Check the connectivity of the cracks. 2022-06-10. NEWFTU2022061001.
        ! And update the crack type based on the connectivity relationship.
        !*********************************************************************
        !call D3_Check_Crack_Connections(i_WB,i_Stage,i_Prop,isub)
      
        !******************************************
        ! Mesh the newly added cracks. 2023-01-09.
        !******************************************
        print *,'    Meshing for 3D new cracks...'  
        call D3_Mesh_Initial_Crack
      
        !********************
        ! Obtain load factor
        !********************
        call Force_Factor(Lambda,isub,Yes_Last_Growth)
        WRITE(*,1002) Lambda
      
      
        !***********************************************
        !if new crack emerged in the previous FEM step, 
        !then add crack number.
        !***********************************************
        if(New_Crack_Flag .eqv..True.)then
            num_Crack =  num_Crack +1
            New_Crack_Flag = .False.
            !Mesh the initial crack.
            print *,'    Meshing updated 3D cracks...'
            call D3_Mesh_Initial_Crack
        endif
      
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !                                             %
      !                                             %
      !                                             %
      !                                             %
      !                                             %
      !                                             %
      !      If there are cracks (XFEM)             %
      !                                             %
      !                                             %
      !                                             %
      !                                             %
      !                                             %
      !                                             %
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(num_Crack.ne.0)then
          !*************************************************************************************
          ! Check whether the maximum allowable number of cracks has been exceeded (2021-08-21)
          !*************************************************************************************
          write (*,4001) num_Crack,Max_Num_Cr_3D
          if(num_Crack > Max_Num_Cr_3D)then
              print *, '    Error :: num_Crack > Max_Num_Cr_3D!'
              print *, '             Error in PhiPsi3D_Static.f.'
              call Warning_Message('S',Keywords_Blank)  
          endif
          
          !*********************************************************************************
          ! Adjust the crack edge nodes (to prevent them from being on the surface of solid
          ! elements)
          !*********************************************************************************
          print *,'    Adjusting crack vertexs...'   
          call D3_Crack_Vertex_Adjust(isub)
          
          !***************************************************************
          ! Calculate the local coordinate system of the crack edge nodes
          !***************************************************************
          print *,'    Calculating local coordinate system of crack vertexs...'   
          call D3_Crack_Vertex_Local_Coor_Sys(isub)

          !***********************************************
          ! Prepare crack front segmentation, 2021-08-20.
          !***********************************************
          if(key_front_segmentation==1)then
              print *,'    Prepare for crack frontsegmentation...'  
              call D3_Prepare_Front_Segmentation(isub)
          endif
          
          !******************************
          ! Confirm the enhancement node
          !******************************
          ifra = 1
          if (Key_Local_Mesh_Refine == 0)then
              print *,'    Determine enriched node...'
              call Determine_Enriched_Nodes_3D(ifra,isub)
          endif
          
          !*************************************************
          ! Enhanced Unit Local Optimization (2021-08-06)
          ! Attention: This operation will change the total 
          ! number of elements, nodes, and coordinates of 
          ! added nodes. The enriched scheme is referred to
          ! my paper.
          !*************************************************
#ifndef github
          if (Key_Local_Mesh_Refine >=1)then
              Flag_Local_Refined =.False.
              if (isub>=2)then
                  ! Reload the original geometric model file, otherwise the encryption will be applied on top of the
                  ! previously encrypted mesh, resulting in the mesh becoming increasingly dense.
                  print *,'    Reload initial mesh...'
                  call Read_Geo_3D
                  print *,'    Update Mesh info of 3D cracks...'  
                  call D3_Update_Mesh_Info_Crack(isub)
              endif
              !Determine enriched nodes and elements for the first time.
              if (Key_Multi_TipEnrNode   ==  1) then
                  Key_Multi_TipEnrNode   =  0
                  Inactive_Key_Multi_TipEnrNode = .True.
              endif
              print *,'    Get enriched node of the initial mesh...'
              call Determine_Enriched_Nodes_3D(ifra,isub)
              print *,'    Local mesh refinement procedure...'
              print *,'    > Refine enriched elements...'
              call Refine_Enriched_Elements_3D(ifra,isub)   
              Flag_Local_Refined = .True.
              
              !Determine enriched nodes and elements again.
              if(Inactive_Key_Multi_TipEnrNode)then
                  Key_Multi_TipEnrNode          =  1
                  Inactive_Key_Multi_TipEnrNode = .False.
              endif
              print *,'    > Update Mesh info of 3D cracks...'  
              call D3_Update_Mesh_Info_Crack(isub)
              !endif
              print *,'    > Adjusting crack vertexs...'
              call D3_Crack_Vertex_Adjust(isub)
              print *,'    > Compute local cs of crack vertexs...'   
              call D3_Crack_Vertex_Local_Coor_Sys(isub)   
              print *,'    > Determine enriched node again...'
              call Determine_Enriched_Nodes_3D(ifra,isub)  
              print *,'    End of Local mesh refinement procedure'
          endif
#endif 
          !*************************************
          ! Assign numbers to enhancement nodes
          !*************************************
          call Number_Enriched_Nodes_3D(isub)     
          print *,'    Total_FD:',Total_FD  
          print *,'    Enrich_Freedom:',Enrich_Freedom  

          !******************************
          ! Consider boundary conditions
          !******************************
          ALLOCATE(freeDOF(Total_FD))
          ALLOCATE(fixedDOF(Total_FD))
          
          call Boundary_Cond_3D(Total_FD,isub,freeDOF,num_FreeD,fixedDOF,num_FixedD)    
 
          ! 2022-09-25. Used for determinant calculation and analysis.
          if(Key_Determinant==1)then
              ALLOCATE(freeDOF_for_Det(Total_FD))
              freeDOF_for_Det = freeDOF
          endif
          
          
          !*************************************
          ! Show degrees of freedom information
          !*************************************
1101      FORMAT(5X,'Total DOFs of solid:       ',I8)
1102      FORMAT(5X,'Total Enriched DOFs:       ',I8) 
1103      FORMAT(5X,'Total active DOFs of solid:',I8)   
          WRITE(*,1101) Total_FD 
          WRITE(*,1102) Enrich_Freedom  
          WRITE(*,1103) num_FreeD


          !*****************************************************************************************
          ! Obtain the list of newly added XFEM elements and the list of XFEM elements that require
          ! an updated stiffness matrix.
          !2022-06-24. NEWFTU2022062403.
          !*****************************************************************************************
          call D3_Get_New_and_Updated_XFEM_Elements(isub)
          
          !***************************************************************************
          ! Calculate related data for crack calculation points (fluid element nodes)
          !***************************************************************************
1104      FORMAT(5X,'Total DOFs of fluid:',I8)  
          call D3_Cal_HF_Crack_Points_Info_Linear(isub)
          num_Tol_CalP_Water = 0
          do i_C=1,Num_Crack
            ! If it is an HF fracture
            if(Crack_Type_Status_3D(i_C,1) == 1) then
                num_Tol_CalP_Water= num_Tol_CalP_Water +Cracks_Real_CalP_Num_3D(i_C)        
            endif
          end do 
          
          WRITE(*,1104) num_Tol_CalP_Water        

          !*********************************************************************************
          ! If it is a block integration algorithm, then perform the relevant calculations.
          ! NEWFTU2022072701.
          !*********************************************************************************
          if (Key_Integral_Sol==4) then
              call D3_Prepare_Subdivision_Integration(isub)
          endif
          
          !********************************************
          ! Allocate other relevant dynamic data space
          !********************************************
          ALLOCATE(F_U(Total_FD))
          
          ! ALLOCATE(Coupled_Q_3D(Total_FD, num_Tol_CalP_Water))! Fluid-structure coupling matrix
          ! Using an uneven array. 2022-09-16
          if(allocated(Coupled_Q_3D)) deallocate(Coupled_Q_3D)  
          if(allocated(Coupled_Q_3D_Index)) deallocate(Coupled_Q_3D_Index)  
          ALLOCATE(Coupled_Q_3D(Total_FD))   
          ALLOCATE(Coupled_Q_3D_Index(Total_FD))
          
          ALLOCATE(F(Total_FD))
          ALLOCATE(CalP_Pres(num_Tol_CalP_Water))        
          ALLOCATE(delta_x(num_FreeD))
          
          !***************************
          ! Initialize some variables
          !***************************
          F(1:Total_FD)             = ZR
          F_U(1:Total_FD)           = ZR

          !*********************************
          ! Calculate the coupling matrix Q
          !*********************************
          if(Flag_HF_3D==1) then
#ifndef Silverfrost
            call Cal_HF_Matrix_Q_Linear_3D(isub,Total_FD,num_Tol_CalP_Water) 
#endif
#ifdef Silverfrost
            print *,'    ERROR :: Silverfrost compiler failed to compile Cal_HF_Matrix_Q_Linear_3D.f90!'
            print *,'             In PhiPsi3D_Static.F90.'
            call Warning_Message('S',Keywords_Blank)
#endif
          endif               
              
          !**************************************************************************
          ! Convert water pressure into a pressure vector for all calculation points
          !**************************************************************************
          ! If it is staged fracturing, June 2, 2022.
          if(Flag_HF_3D==1) then
              num_Count = 0
              do i_C = 1,num_Crack    
                ! If it is an HF fracture
                if(Crack_Type_Status_3D(i_C,1) == 1) then
                  ! If it weren't for staged fracturing
                  if (Key_HF_Multistage_3D==0) then
                      do i_CalP=1,Cracks_Real_CalP_Num_3D(i_C)   
                        num_Count = num_Count + 1
                        CalP_Pres(num_Count)=Crack_Pressure(i_C) - Cracks_FluidEle_CalP_Glo_Insitu(num_Count)
                         
                      end do
                  endif       
                  ! If it is staged fracturing
                  if (Key_HF_Multistage_3D==1) then
                      do i_CalP=1,Cracks_Real_CalP_Num_3D(i_C)   
                        num_Count = num_Count + 1
                        CalP_Pres(num_Count) = 10.0D6 -Cracks_FluidEle_CalP_Glo_Insitu(num_Count)
                        if(mod(i_CalP,100)==0)then
                          print *,'     ===    10.0 MPa pressure applied here!   ==='
                        endif
                      end do
                  endif                      
                endif
              end do    
          endif              
     
          !***********************************************
          ! Assembly load vector (without fluid pressure)
          !***********************************************
          call Force_Vector_3D(Total_FD,isub,Lambda,F_U)

          
          !*****************************************************
          ! Assemble the load vector considering fluid pressure
          !*****************************************************
          if(Flag_HF_3D==1)then
              do i_Dof =1,num_FreeD
                if(allocated(Coupled_Q_3D(freeDOF(i_Dof))%row))then 
                    ! Each fluid element has 3 fluid nodes. BUGFIX2022111101.
                    F(freeDOF(i_Dof)) =  F_U(freeDOF(i_Dof))+   &
                            dot_product(Coupled_Q_3D(freeDOF(i_Dof))%row(1:3), &
                                        CalP_Pres(Coupled_Q_3D_Index(freeDOF(i_Dof))%row(1:3)))
                endif
              enddo
          else
              F(freeDOF(1:num_FreeD)) =  F_U(freeDOF(1:num_FreeD))
          endif
          WRITE(*,3003) sum(abs(F))
          
          ! Save node loads (including enhanced nodes) for Matlab post-processing, 2022-04-23,
          ! NEWFTU2022042301.
          call Save_Dof_Force_3D(isub,F,Total_FD)
          
          !****************************************************
          !       If solver 11: Element by Element, diagonally
          !       preconditioned conjugate gradient solver
          !****************************************************
          if(Key_SLOE==11)then           
              ALLOCATE(DISP(Total_FD))
              !ALLOCATE(storK_XFEM(MDOF_3D,MDOF_3D,num_XFEM_Elem))   !2022-04-16
              ALLOCATE(storK_XFEM(num_XFEM_Elem))     
              ! storK_FEM is calculated only once because it will not change afterward and will only decrease.
              if(Key_EBE_Sym_Storage_K==0)then
                  if (allocated(storK_FEM).eqv. .False.) then
                      ALLOCATE(storK_FEM(24,24,num_FEM_Elem))  
                      ALLOCATE(Size_Local_3D(Num_Elem))
                      ALLOCATE(All_Local_3D(MDOF_3D,Num_Elem))   
                  endif
              endif
              if(Key_EBE_Sym_Storage_K==1)then
                  if (allocated(storK_FEM_Sym).eqv. .False.) then
                      ALLOCATE(storK_FEM_Sym(300,num_FEM_Elem))  
                      ALLOCATE(Size_Local_3D(Num_Elem))
                      ALLOCATE(All_Local_3D(MDOF_3D,Num_Elem))   
                  endif
              endif
       
              ! Expand the all_local array. !IMPROV2023011101. 2023-01-11. Avoid the problem of insufficient
              ! MDOF_3D.
              if(MDOF_3D > Old_MDOF_3D)then
#ifndef Silverfrost 
                  print *,'    Extending matrix All_Local_3D...'  
                  allocate(tem_all_local(MDOF_3D,Num_Elem))  
                  tem_all_local(1:MDOF_3D,1:Num_Elem)  = 0
                  tem_all_local(1:Old_MDOF_3D,1:Num_Elem) = All_Local_3D(1:Old_MDOF_3D,1:Num_Elem)
                  deallocate(All_Local_3D)
                  call move_alloc(tem_all_local,All_Local_3D)
#endif
#ifdef Silverfrost  
                  print *, '    Error :: Extending matrix All_Local_3D is not valid for Silverfrost compiler!'
                  print *, '             move_alloc function in Silverfrost is faulty, thus code is commented out!'
                  print *, '             In PhiPsi3D_Static.F90.'
                  call Warning_Message('S',Keywords_Blank)
#endif
              endif
              
              !ALLOCATE(storK_FEM(MDOF_3D,MDOF_3D,num_FEM_Elem))               !2022-04-16
              ALLOCATE(diag_precon(0:num_FreeD))
              
              call Ele_by_Ele_XFEM_PCG_3D(isub,Lambda,cg_tol,max_cg,             &
                    num_FreeD,freeDOF(1:num_FreeD),F(freeDOF(1:num_FreeD)),DISP, & 
                    diag_precon(0:num_FreeD))

              ! Save the number of Gaussian integration points for each element to a *.elgn file for
              ! post-processing. NEWFTU2022071601.
              if (Key_Post_Elements_Gauss_Num== 1) then
                  call Save_Gauss_Num_of_Elements(isub)
              endif
 
              ! Read test from binary file
              !CALL Read_storK_XFEM_from_Binary_file(MDOF_3D,31,Ele_storK_XFEM(1:MDOF_3D,1:MDOF_3D))   
              goto 555
          endif
              
          !***************************
          ! Assembly stiffness matrix
          !***************************
          if(Key_K_Sparse==0)then
              ALLOCATE(globalK(Total_FD,Total_FD))
              print *,'    Assembling the stiffness matrix K...'
              call Assemble_Stiffness_Matrix_XFEM_3D(isub,globalK,Total_FD,Total_Num_G_P)   
              write(*,3002) sum(globalK)    
          elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
!#ifndef github
              print *,'    Get max NNZ before assembling K......'
              call Assemble_Stiffness_Matrix_SPARS_XFEM_3D_Get_MaxNNZ(isub,&
                          freeDOF(1:num_FreeD),num_FreeD,Total_FD,K_CSR_NNZ_Max)
              print *,'    K_CSR_NNZ_Max:',K_CSR_NNZ_Max    
              
              if (allocated(K_CSR_aa)) DEALLOCATE(K_CSR_aa)
              if (allocated(K_CSR_ja)) DEALLOCATE(K_CSR_ja)
              if (allocated(K_CSR_ia)) DEALLOCATE(K_CSR_ia)
              ALLOCATE(K_CSR_aa(K_CSR_NNZ_Max))
              ALLOCATE(K_CSR_ja(K_CSR_NNZ_Max))
              ALLOCATE(K_CSR_ia(num_FreeD+1))
              print *,'    Assemble the sparse K...'
              ! Assemble the sparse stiffness matrix (automatically considers boundary conditions and removes
              ! sparse matrix elements corresponding to constrained degrees of freedom)
              call Assemble_Stiffness_Matrix_SPARS_XFEM_3D(isub, &
                   freeDOF(1:num_FreeD),num_FreeD,&
                   fixedDOF(1:num_FixedD),num_FixedD,     &
                   K_CSR_aa,K_CSR_ja,K_CSR_ia,&
                   K_CSR_NNZ_Max,K_CSR_NNZ,&
                   Total_FD,Total_Num_G_P)
!#endif
#endif
#endif
          endif

          ! Save the number of Gaussian integration points for each element to a *.elgn file for
          ! post-processing. NEWFTU2022071601.
          if (Key_Post_Elements_Gauss_Num== 1) then
              call Save_Gauss_Num_of_Elements(isub)
          endif
              
          !******************************************************************************************
          !       Pre_Conditioner     DATE: 2020-02-11
          ! Condition number reduction, used to reduce the condition number of the stiffness matrix,
          ! thereby lowering the system of equations
          ! Pathological characteristics, reduced computational load
          !******************************************************************************************
          if(Key_K_Sparse==0 .and. Key_Pre_Conditioner==1)then
            !====================
            ! Global Formation K
            !====================
            ALLOCATE(Vector_Dc(Total_FD))
            print *,'    Applying the pre-conditioner to K...'  
            call Apply_Pre_Conditioner_to_K(Total_FD,globalK,Vector_Dc)
            print *,'    Applying the pre-conditioner to F...'  
            do i_FD =1,Total_FD
                F(i_FD) = F(i_FD)*Vector_Dc(i_FD)                     
            enddo 
          endif
              
          !***************************************************************
          ! Solve for displacement (solve the linear system of equations)
          !***************************************************************
          ALLOCATE(DISP(Total_FD))
          ALLOCATE(tem_DISP(num_FreeD))
          print *,'    Solving displacements......'     
          if(Key_K_Sparse==0)then
              call Matrix_Solve_LSOE(0,1,Key_SLOE,globalK(freeDOF(1:num_FreeD),    & 
                    freeDOF(1:num_FreeD)),F(freeDOF(1:num_FreeD)),delta_x(1:num_FreeD),num_FreeD) 
          elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
#ifndef Silverfrost
!#ifndef github
              call Matrix_Solve_LSOE_Sparse(0,1,Key_SLOE,&
                        K_CSR_NNZ,&
                        K_CSR_aa(1:K_CSR_NNZ),&
                        K_CSR_ja(1:K_CSR_NNZ),&
                        K_CSR_ia(1:num_FreeD+1),&
                        F(freeDOF(1:num_FreeD)),&
                        delta_x(1:num_FreeD),num_FreeD)
!#endif
#endif
#endif
#endif
          endif
          
          DISP(1:Total_FD) = ZR
          DISP(freeDOF(1:num_FreeD)) = delta_x(1:num_FreeD)
          
          ! Full Matrix Pre_Conditioner DISP Processing, 2020-02-12
          if(Key_K_Sparse==0 .and. Key_Pre_Conditioner==1) then       
              if(Key_Pre_Conditioner==1)then
                  !Apply Vector_Dc to D.
                  print *,'    Applying the pre-conditioner to D...'  
                  do i_FD =1,Total_FD
                      DISP(i_FD) = DISP(i_FD)*Vector_Dc(i_FD)                     
                  enddo     
                  DEALLOCATE(Vector_Dc)
              endif
          endif
     
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! If considering contact, January 16, 2020
          ! Note: 3D contact analysis does not support sparse matrices
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          if(Key_Contact==1) then
            ALLOCATE(Contact_DISP(Total_FD))
            Contact_DISP(1:Total_FD) = ZR
            ! Given the initial value of Contact_DISP (considering the displacement at contact as the initial
            ! value)
            Contact_DISP(freeDOF(1:num_FreeD))=delta_x(1:num_FreeD)
            ! ------Used for contacting iterative ori_globalK
            ALLOCATE(ori_globalK(Total_FD,Total_FD))
            ori_globalK = globalK          
            ! call Determine_Contact_State_by_Iteration_3D(isub, isub, isub, Contact_DISP, Total_FD, & !
            ! Contact_DISP is both input and output
            ! Usual_Freedom, Enrich_Freedom, freeDOF, num_FreeD, F_U, ori_globalK, globalK) ! Number of
            ! traditional degrees of freedom and number of enriched degrees of freedom (note that both of these
            ! are global variables), globalK is the output
            ! BUGFIX2022080601. Change F_U to F.
            call Determine_Contact_State_by_Iteration_3D(isub,isub,isub,Contact_DISP,Total_FD,   &
                   Usual_Freedom,Enrich_Freedom,&
                   freeDOF,num_FreeD,F,ori_globalK,globalK)
            ! Return the displacement calculated after the new contact consideration
            DISP = Contact_DISP
          endif
          
          555 CONTINUE
  
          !====================================================
          !Element by Element Contact Calculation, 2020-03-31.
          !====================================================
          if(Key_SLOE==11 .and. Key_Contact /= 0) then
            !call EBE_Determine_Contact_State_by_Iteration_3D(isub,cg_tol,max_cg,num_FreeD,freeDOF(1:num_FreeD),
            !& !cg_tol,max_cg: see module Global_Common
            ! F_U, DISP, storK_XFEM, storK_FEM, Size_Local_3D, All_Local_3D, diag_precon(0:num_FreeD)) ! DISP is
            ! both input and output
            ! BUGFIX2022080601. Change F_U to F.
            print *,'    Determine contact states by iteration...' 
            call EBE_Determine_Contact_State_by_Iteration_3D(isub,cg_tol,max_cg,num_FreeD,&
                  freeDOF(1:num_FreeD), F,DISP, diag_precon(0:num_FreeD),5)
          endif   
     
          !*******************
          ! Save displacement
          !*******************
          call Tool_Print_to_Screen(2,5)
          call Save_Disp(isub,1)
          
          !***********************
          ! Calculate crack width
          !***********************
          print *,'    Calculating crack apertures......'   
          call Cal_Crack_Aperture_3D(isub,DISP)
          Max_Cr_Vol = maxval(Cracks_Volume(1:num_Crack))
          Min_Cr_Vol = minval(Cracks_Volume(1:num_Crack))
          
          do i_C = 1,num_Crack
              WRITE(*,2004) i_C,maxval(Cracks_CalP_Aper_3D(i_C)%row(1:Cracks_CalP_Num_3D(i_C)))*1.0D3  
              WRITE(*,2005) i_C,minval(Cracks_CalP_Aper_3D(i_C)%row(1:Cracks_CalP_Num_3D(i_C)))*1.0D3     
          enddo
          
          if (num_Crack >1)then
              WRITE(*,2001) Max_Cr_Vol
              WRITE(*,2002) Min_Cr_Vol
          elseif (num_Crack ==1) then
              WRITE(*,2003) Cracks_Volume(1)  
          endif

          !*********************************************
          ! Calculate fracture permeability. 2022-11-26
          !*********************************************
          print *,'    Calculating crack and element permeability......'   
          call Cal_Crack_Permeability_3D(isub)
             
          !********************************************************************
          ! Save crack-related files (including enhanced node numbering c_POS)
          !********************************************************************
          call Save_Files_Crack_3D(isub)    
          
          !**************************************
          ! Save crack width related information
          !**************************************          
          call Save_Files_Cr_CalP_3D(isub)          
          
          !***********************************
          ! Calculate stress intensity factor
          !***********************************
          print *,'    Calculating stress intensity factors......'   
          if (CFCP==3 .and. Key_SIFs_Method==1) then
              call Cal_SIFs_DIM_3D(isub,DISP,5)
              call Tool_Print_to_Screen(1,5)
          end if

          !****************************************************
          ! Save VTK CRACK file, 2023-07-28. NEWFTU2023072801.
          !****************************************************
          call Save_vtk_file_for_Crack(isub)

         !******************************************************************************************
         ! If the crack is allowed to extend, then determine whether it extends.
         ! If an extension occurs, calculate the new crack tip coordinates, add new crack segments,
         ! and update.
         ! Crack coordinates.
         !******************************************************************************************
         if (Key_Propagation ==1) then
           ! Determine the direction and size of crack propagation
           call Check_Crack_Grows_3D(1,isub,Yes_Growth)
           ! Smooth treatment of crack front edge, 2021-09-08
           if(Key_Smooth_Front>=1 .and. Key_Smooth_Front/=6)then
               call Crack_Front_Smooth_3D(1,isub)
           endif
           ! Save the maximum principal stress direction vectors of the discrete crack edge nodes
           call Save_Files_Crack_3D_Extra(isub)
           
           ! Check for any new cracks, 2020-03-12.
           New_Crack_Flag =.false.
           if (Key_Initiation ==1 .and. Key_Max_Num_Initiation_Cracks>=2 &
              .and. (Num_Initiation_Cracks < Key_Max_Num_Initiation_Cracks)) then
                call Check_Crack_initialize_3D(1,isub,New_Crack_Flag)  
           endif
         
           ! If a crack has propagated, update the crack propagation flag: Yes_Last_Growth
           if (any(Yes_Growth).eqv..True.)then
               Yes_Last_Growth = .True.
           ! If no cracks have propagated, exit the program.
           else
               Yes_Last_Growth = .False.
               if(isub < Num_Substeps .and. New_Crack_Flag .eqv. .false.)then
                   print *,'    ---$---$---$---$---$---$---$---$---$---S---S---S---S---'
                   print *,'    Warning :: No crack propapated, program will be ended! '
                   print *,'    ---$---$---$---$---$---$---$---$---$---S---S---S---S---'
               elseif(isub == Num_Substeps)then
                   print *,'    |<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|'
                   print *,'    |    All propagation steps done!  |'
                   print *,'    |<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|'
               endif
           endif
         end if
                         
         !*************************************************************
         ! Determine the new number of cracks and update the next step
         !num_Crack_Log(isub+1)
         !*************************************************************
         num_Crack_Log(isub+1)  = num_Crack
              
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !                                             %
      !                                             %
      !                                             %
      !                                             %
      !                                             %
      !                                             %
      !   If there are no fractures (FEM)           %
      !                                             %
      !                                             %
      !                                             %
      !                                             %
      !                                             %
      !                                             %
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      elseif (num_Crack.eq.0)then 
          !**************************
          ! Total degrees of freedom
          !**************************
          Total_FD = 3*Num_Node
          print *,'    Total_FD:',Total_FD
          
          !************************
          !  Get the force vector.
          !************************
          ALLOCATE(F(Total_FD))
          call Force_Vector_3D(Total_FD,isub,Lambda,F)
          
          !******************************
          ! Consider boundary conditions
          !******************************
          ALLOCATE(freeDOF(Total_FD))
          ALLOCATE(fixedDOF(Total_FD))
          call Boundary_Cond_3D(Total_FD,isub,freeDOF,num_FreeD,fixedDOF,num_FixedD) 
 
          !****************************************************
          !       If solver 11: Element by Element, diagonally
          !       preconditioned conjugate gradient solver
          !****************************************************
          if(Key_SLOE==11)then
              ALLOCATE(DISP(Total_FD))
              call Ele_by_Ele_FEM_PCG_3D(isub,Lambda,cg_tol,max_cg,    &
                     num_FreeD,freeDOF(1:num_FreeD),F(freeDOF(1:num_FreeD)),DISP)       
              goto 666
          endif
          
          !***************************
          ! Assembly stiffness matrix
          !***************************
          if(Key_K_Sparse==0)then
              print *,'    Assembling the stiffness matrix K......'
              ALLOCATE(globalK(Total_FD,Total_FD))
              call Assemble_Stiffness_Matrix_FEM_3D(isub,globalK,Total_FD,Total_Num_G_P)
              write(*,3002) sum(globalK)
          elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
#ifndef Silverfrost
!#ifndef github
              ! Get the number of maximum non-zero elements.
              print *,'    Get max NNZ before assembling K......'
              call Assemble_Stiffness_Matrix_SPARS_FEM_3D_Get_MaxNNZ(isub,&
                          freeDOF(1:num_FreeD),num_FreeD,Total_FD,K_CSR_NNZ_Max)
              print *,'    K_CSR_NNZ_Max:',K_CSR_NNZ_Max           
              if (allocated(K_CSR_aa)) DEALLOCATE(K_CSR_aa)
              if (allocated(K_CSR_ja)) DEALLOCATE(K_CSR_ja)
              if (allocated(K_CSR_ia)) DEALLOCATE(K_CSR_ia)
              ALLOCATE(K_CSR_aa(K_CSR_NNZ_Max))
              ALLOCATE(K_CSR_ja(K_CSR_NNZ_Max))
              ALLOCATE(K_CSR_ia(num_FreeD+1))
              ! Assemble the sparse stiffness matrix (automatically considers boundary conditions and removes
              ! sparse matrix elements corresponding to constrained degrees of freedom)
              print *,'    Assembling the sparse K......'
              call Assemble_Stiffness_Matrix_SPARS_FEM_3D(isub,&
                      freeDOF(1:num_FreeD),num_FreeD,&
                      K_CSR_aa,K_CSR_ja,K_CSR_ia,&
                      K_CSR_NNZ_Max,K_CSR_NNZ,Total_FD,Total_Num_G_P)
!#endif
#endif
#endif
#endif
          endif
          
          !************************
          ! Solve for displacement
          !************************
          ALLOCATE(DISP(Total_FD))
          DISP(1:Total_FD) = ZR
          ALLOCATE(tem_DISP(num_FreeD))
          
          print *,'    Solving displacements......'   
          if(Key_K_Sparse==0)then
              call Matrix_Solve_LSOE(0,1,Key_SLOE,globalK(freeDOF(1:num_FreeD),    & 
                    freeDOF(1:num_FreeD)),F(freeDOF(1:num_FreeD)),tem_DISP,num_FreeD) 
          elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
!#ifndef github
              call Matrix_Solve_LSOE_Sparse(0,1,Key_SLOE,&
                        K_CSR_NNZ,&
                        K_CSR_aa(1:K_CSR_NNZ),&
                        K_CSR_ja(1:K_CSR_NNZ),&
                        K_CSR_ia(1:num_FreeD+1),&
                        F(freeDOF(1:num_FreeD)),&
                        tem_DISP,num_FreeD)
!#endif
#endif
#endif
          endif
          
          DISP(freeDOF(1:num_FreeD)) = tem_DISP
          
666       continue    
          
          call Tool_Print_to_Screen(2,5)
    
          !*******************
          ! Save displacement
          !*******************
          call Save_Disp(isub,1)
          
          !***************************************************
          ! Check if any new cracks have appeared, 2020-03-12
          !***************************************************
          if (Key_Initiation ==1) then
              call Check_Crack_initialize_3D(1,isub,New_Crack_Flag)  
              ! First activation of XFEM. 2023-01-23. BUGFIX2023012301.
              if(First_XFEM_Step == 0) then
                   First_XFEM_Step = 1  
              endif
          endif              
      end if
          
      !#################################################################################
      ! Displacement in the Column Coordinate System of the Computing Node (2021-09-11)
      !#################################################################################
      if(Key_Post_CS_N_Stra_Cylindr==1)then
          ! Calculate the transformation matrix from Cartesian coordinates to cylindrical coordinates
          ALLOCATE(Theta_Cartesian_to_Cylinder(Num_Elem,8))
          ALLOCATE(Theta_Cartesian_to_Cylinder_Node(Num_Node))
          call Get_Theta_Cartesian_to_Cylinder
          ! Allocate memory
          ALLOCATE(DISP_Cylinder(Total_FD))             
          !
          call Convert_DISP_to_Cylinder_3D(Total_FD,DISP(1:Total_FD),DISP_Cylinder(1:Total_FD))  
          call Tool_Print_to_Screen(3,5)
          !Save nodal DISP.
          call Save_Disp(isub,2)
          
          DEALLOCATE(Theta_Cartesian_to_Cylinder)
          DEALLOCATE(Theta_Cartesian_to_Cylinder_Node)
          DEALLOCATE(DISP_Cylinder)
      endif 
          
      !###########################
      ! Computational Node Stress
      !###########################
      if(Key_Post_CS_N_Strs==1)then
          ALLOCATE(Stress_xx_Node(num_Node),Stress_yy_Node(num_Node),Stress_zz_Node(num_Node))
          ALLOCATE(Stress_xy_Node(num_Node),Stress_yz_Node(num_Node),Stress_xz_Node(num_Node),Stress_vm_Node(num_Node))      
          ! No XFEM.
          if(Enrich_Freedom==0) then
              call Get_Node_Stress_FEM_IN_OUT_3D(.True.,isub,DISP,   &
                    Stress_xx_Node,Stress_yy_Node,Stress_zz_Node,    & 
                    Stress_xy_Node,Stress_yz_Node,Stress_xz_Node,Stress_vm_Node)
          !XFEM.
          else
              ! Took a very long time (issue fixed on 2020-12-27)
              call Get_Node_Str_XFEM_IN_OUT_3D(.True., isub,             &
                   DISP,1,1,Stress_xx_Node,Stress_yy_Node,Stress_zz_Node,&
                   Stress_xy_Node,Stress_yz_Node,Stress_xz_Node,Stress_vm_Node)
          endif
          call Tool_Print_to_Screen(4,5)
          !Save nodal stress.
          call Save_Stress_Node(isub,1)
      endif    

      !######################################
      ! Calculation Node Strain (2021-09-10)
      !######################################
      if(Key_Post_CS_N_Stra==1)then 
          ALLOCATE(Strain_xx_Node(num_Node),Strain_yy_Node(num_Node),Strain_zz_Node(num_Node))
          ALLOCATE(Strain_xy_Node(num_Node),Strain_yz_Node(num_Node),Strain_xz_Node(num_Node),Strain_vm_Node(num_Node))  
          ! Accurate calculation, very time-consuming (issue fixed on 2020-12-27)
          call Get_Node_Str_XFEM_IN_OUT_3D(.True., isub,            &
                 DISP,2,1,Strain_xx_Node,Strain_yy_Node,Strain_zz_Node,   &
                 Strain_xy_Node,Strain_yz_Node,Strain_xz_Node,Strain_vm_Node)
          !Save nodal strain.
          call Tool_Print_to_Screen(5,5)
          call Save_Strain_Node(isub,1)
      endif    

      !########################################################################################
      ! Compute node stress, results in cylindrical coordinate system, such as circumferential
      ! stress (2021-09-11)
      !########################################################################################
      if(Key_Post_CS_N_Strs_Cylindr==1)then
          ! Calculate the transformation matrix from Cartesian coordinates to cylindrical coordinates
          ALLOCATE(Theta_Cartesian_to_Cylinder(Num_Elem,8))
          ALLOCATE(Theta_Cartesian_to_Cylinder_Node(Num_Node))
          call Get_Theta_Cartesian_to_Cylinder
          
          ! Allocate memory (tt represents theta_theta).
          ALLOCATE(Stress_Crr_Node(num_Node),Stress_Ctt_Node(num_Node),Stress_Czz_Node(num_Node))
          ALLOCATE(Stress_Crt_Node(num_Node),Stress_Ctz_Node(num_Node),Stress_Crz_Node(num_Node))
          ALLOCATE(Stress_Cvm_Node(num_Node))  
          ! No XFEM.
          if(Enrich_Freedom==0) then
              call Get_Node_Stress_FEM_IN_OUT_3D(.True.,isub,DISP,        &
                         Stress_Crr_Node,Stress_Ctt_Node,                 &  
                         Stress_Czz_Node, Stress_Crt_Node,Stress_Ctz_Node,&
                         Stress_Crz_Node,Stress_Cvm_Node)
          !XFEM.
          else
              ! Accurate calculation, extremely time-consuming (issue fixed on 2020-12-27)
              call Get_Node_Str_XFEM_IN_OUT_3D(.True., isub,              &
                         DISP,1,2,Stress_Crr_Node,Stress_Ctt_Node,        &
                         Stress_Czz_Node, Stress_Crt_Node,Stress_Ctz_Node,&
                         Stress_Crz_Node,Stress_Cvm_Node)
          endif 
          call Tool_Print_to_Screen(6,5)
          !Save nodal stress.
          call Save_Stress_Node(isub,2)
          DEALLOCATE(Theta_Cartesian_to_Cylinder)
          DEALLOCATE(Theta_Cartesian_to_Cylinder_Node)
      endif        
          
      !###########################################################################################
      ! Calculated node strain, results in cylindrical coordinate system, such as circumferential
      ! strain (2021-09-10)
      !###########################################################################################
      if(Key_Post_CS_N_Stra_Cylindr==1)then
          ! Calculate the transformation matrix from Cartesian coordinates to cylindrical coordinates
          ALLOCATE(Theta_Cartesian_to_Cylinder(Num_Elem,8))
          ALLOCATE(Theta_Cartesian_to_Cylinder_Node(Num_Node))
          call Get_Theta_Cartesian_to_Cylinder
          
          ! Allocate memory (tt represents theta_theta).
          ALLOCATE(Strain_Crr_Node(num_Node),Strain_Ctt_Node(num_Node),Strain_Czz_Node(num_Node))
          ALLOCATE(Strain_Crt_Node(num_Node),Strain_Ctz_Node(num_Node),Strain_Crz_Node(num_Node))
          ALLOCATE(Strain_Cvm_Node(num_Node))              
          ! Accurate calculation, extremely time-consuming (issue fixed on 2020-12-27)
          call Get_Node_Str_XFEM_IN_OUT_3D(.True., isub,             &
                 DISP,2,2,Strain_Crr_Node,Strain_Ctt_Node,Strain_Czz_Node,  &
                 Strain_Crt_Node,Strain_Ctz_Node,Strain_Crz_Node,Strain_Cvm_Node)
          call Tool_Print_to_Screen(6,5)
          !Save nodal strain.
          call Save_Strain_Node(isub,2)
          DEALLOCATE(Theta_Cartesian_to_Cylinder)
          DEALLOCATE(Theta_Cartesian_to_Cylinder_Node)
          !----------------------------------------
          !----------------------------------------
      endif            
      
      !***************************
      ! Save VTK file, 2021-07-16
      !***************************
      call Save_vtk_file(isub)
      
      !*********************************************************************
      ! Check the connectivity of the cracks. 2022-06-10. NEWFTU2022061001.
      ! And update the crack type based on the connectivity relationship.
      ! IMPROV2024022202.
      !*********************************************************************
      call D3_Check_Crack_Connections(i_WB,i_Stage,i_Prop,isub)
          
      !###############################
      ! Clearing dynamic arrays, etc.
      !###############################
      if(Key_K_Sparse==0)then
          if (allocated(globalK)) DEALLOCATE(globalK)
      elseif(Key_K_Sparse==1)then
#ifdef gfortran
#ifdef notcbfortran
          if (allocated(K_CSR_aa)) DEALLOCATE(K_CSR_aa)
          if (allocated(K_CSR_ja)) DEALLOCATE(K_CSR_ja)
          if (allocated(K_CSR_ia)) DEALLOCATE(K_CSR_ia)
#endif
#endif
      endif
      if(Key_Contact/=0 .and. num_Crack.ne.0)then
          if (allocated(Contact_DISP)) DEALLOCATE(Contact_DISP)          
          if (allocated(ori_globalK)) DEALLOCATE(ori_globalK)
      endif          
      !DEALLOCATE(globalF)
      if (allocated(freeDOF)) DEALLOCATE(freeDOF)
      if (allocated(freeDOF_for_Det)) DEALLOCATE(freeDOF_for_Det)
      if (allocated(fixedDOF)) DEALLOCATE(fixedDOF)
      if (allocated(DISP)) DEALLOCATE(DISP)
      if (allocated(tem_DISP)) DEALLOCATE(tem_DISP)
      if (allocated(F)) DEALLOCATE(F)
      
      if(num_Crack.ne.0)then
          if (allocated(F_U)) DEALLOCATE(F_U)
          !if (allocated(Coupled_Q_3D)) DEALLOCATE(Coupled_Q_3D)
          if (allocated(CalP_Pres)) DEALLOCATE(CalP_Pres)
          if (allocated(delta_x)) DEALLOCATE(delta_x)   
      endif
      
      if(Key_Post_CS_N_Strs==1)then
          if (allocated(Stress_xx_Node))DEALLOCATE(Stress_xx_Node)
          if (allocated(Stress_yy_Node))DEALLOCATE(Stress_yy_Node)
          if (allocated(Stress_zz_Node))DEALLOCATE(Stress_zz_Node)
          if (allocated(Stress_xy_Node))DEALLOCATE(Stress_xy_Node)
          if (allocated(Stress_yz_Node))DEALLOCATE(Stress_yz_Node)
          if (allocated(Stress_xz_Node))DEALLOCATE(Stress_xz_Node)
          if (allocated(Stress_vm_Node))DEALLOCATE(Stress_vm_Node)              
      endif
      
      if(Key_Post_CS_N_Stra==1)then
          if (allocated(Strain_xx_Node))DEALLOCATE(Strain_xx_Node)
          if (allocated(Strain_yy_Node))DEALLOCATE(Strain_yy_Node)
          if (allocated(Strain_zz_Node))DEALLOCATE(Strain_zz_Node)
          if (allocated(Strain_xy_Node))DEALLOCATE(Strain_xy_Node)
          if (allocated(Strain_yz_Node))DEALLOCATE(Strain_yz_Node)
          if (allocated(Strain_xz_Node))DEALLOCATE(Strain_xz_Node)
          if (allocated(Strain_vm_Node))DEALLOCATE(Strain_vm_Node)              
      endif

      if(Key_Post_CS_N_Stra_Cylindr==1)then
          if (allocated(Strain_Crr_Node))DEALLOCATE(Strain_Crr_Node)
          if (allocated(Strain_Ctt_Node))DEALLOCATE(Strain_Ctt_Node)
          if (allocated(Strain_Czz_Node))DEALLOCATE(Strain_Czz_Node)
          if (allocated(Strain_Crt_Node))DEALLOCATE(Strain_Crt_Node)
          if (allocated(Strain_Ctz_Node))DEALLOCATE(Strain_Ctz_Node)
          if (allocated(Strain_Crz_Node))DEALLOCATE(Strain_Crz_Node)
          if (allocated(Strain_Cvm_Node))DEALLOCATE(Strain_Cvm_Node)              
      endif
      
      !Clear solver 11 (EBE) related temporary variables.
      if (allocated(storK_XFEM)) DEALLOCATE(storK_XFEM)       
      if (allocated(diag_precon)) DEALLOCATE(diag_precon)    
      
      
      ! If all cracks have stopped propagating and no new cracks have formed, terminate the program.
      if(num_Crack.ne.0)then
          if((Yes_Last_Growth .eqv. .False.) .and. (New_Crack_Flag .eqv. .false.))then
              goto 200
          endif
      endif
      
      ! CPU time consumption.
      call Tool_Get_Current_Time(current_data,date_time,F_time)
      WRITE(*,3001) F_time-S_time,(dble(F_time)-dble(S_time))/Con_60

         
        enddo
    enddo
enddo
      
200 continue
      
!-----------------------------------------------------------
! Delete temporary binary files for the group sparse matrix
!-----------------------------------------------------------
if (file_Sparse_K_Location_COO_bin .eqv. .True.)then
    c_File_name_1=trim(Full_Pathname)//'_Sparse_K_Location_COO.bin'
    ! call system('del c_File_name_1') !Sparse_K_Location_COO.bin file, stores the indices of the
    ! non-zero elements corresponding to the stiffness matrix elements
    print *,'    Deleting temporary binary files......'  
    Open(1234,File =c_File_name_1,Access='Direct',Form = 'Unformatted', RecL = 4)
    close(1234, status='delete')
endif
if (file_Sparse_K_Mask_COO_bin .eqv. .True.)then
    c_File_name_1=trim(Full_Pathname)//'_Sparse_K_Mask_COO.bin'
    ! call system('del c_File_name_1') !Sparse_K_Location_COO.bin file, stores the indices of the
    ! non-zero elements corresponding to the stiffness matrix elements
    print *,'    Deleting temporary binary files......'  
    Open(1234,File =c_File_name_1,Access='Direct',Form = 'Unformatted', RecL = 4)
    close(1234, status='delete')    
endif      
      
      
RETURN
END SUBROUTINE PhiPsi3D_Static
