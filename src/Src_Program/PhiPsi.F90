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
 
PROGRAM PhiPsi
!This is the main program of PhiPsi.
!
!---------------------
!DESCRIPTION of files
!---------------------
! Jzcr: Coordinate file of circular inclusions
! Jzpx: Polygon Mixed X Coordinate File
! Jzpy: Y-coordinate file of polygon inclusions
! strn: Stress at the node
! sttn: Thermal stress of the node
! strg: Stress at Gauss points
! sttg: Thermal stress at Gauss points
! eqac: earthquake acceleration value
! fbvl: Nodal numbers of boundary value problems and their corresponding 
!       scalar values (calculation results)
! fbfx: Flow in the x-direction at the field point (calculated result)
! fbfy: Flow in the y-direction at the field point (calculated result)
! fbqn: Node numbers and flow values at the boundary of normal flux in 
!       field problems
! fbiv: Field problem boundary node numbers and their corresponding 
!       scalar
!       values (initial values, used for transient analysis)
! gasp: Evaluating shale gas flow rates and production at different 
!       times
! bhpc: Time and pressure of bottom-hole pressure extraction curve
! post: Summary of other information, such as analysis type, etc.
! cape: crack width w
! cond: Conductivity of proppant supporting fractures Cf
! pokf: The permeability of fractures supported by proppant kf = Cf / w
! wppg: Used for shale gas production evaluation and analysis, it stores
!       information on the aperture of supporting fractures
! plep: Equivalent plastic strain at Gauss points in plastic analysis
! fdcu: Load-displacement curve of a certain node, F-F Curve
! cohx: Cohesive crack calculation point cohesive strength X-direction 
!       component
! cohy: Cohesive crack calculation point cohesive force component in the
!       y-direction
! elcs: Contact status of the element relative to each crack
! elco: Cohesive fracture state of the element relative to each crack
! njel: The Junction node corresponds to the Junction element number
! ncel: The Cross element number corresponding to the Cross-enhanced node

use Global_Float_Type
use Global_Common
use Global_Filename
use Global_HF
use Global_Crack
use Global_Model
use Module_Tool_Logger
use omp_lib
use iso_fortran_env
use ISO_C_BINDING
use Function_MATMUL_LP
use omp_lib
#ifndef Silverfrost
use fpm_environment
#endif 

implicit none
!include 'mpif.h'
integer  date_time(8)
integer(LIT) F_time
character*10  current_data
integer nthreads
CHARACTER(32) :: login_user
character(200) temp_log
logical Flag_Read_Geo
#ifdef Silverfrost
character(len=256) :: command
character(len=256) :: cwd_output_file
#endif  
1000 FORMAT('  ')
1001 FORMAT('     Strat time: ',A4,'-',A2,'-',A2,', ',I2,':',I2,':',I2)
1002 FORMAT('     End time: ',A4,'-',A2,'-',A2,', ',I2,':',I2,':',I2)
1003 FORMAT('     Total elapsed CPU time: ',I8,' s, about ',F10.4,' mins')
2001 FORMAT('Start time: ',A8,', ',I2,':',I2,':',I2)
2002 FORMAT('End time: ',A8,', ',I2,':',I2,':',I2)
2003 FORMAT('Total elapsed CPU time: ',I8,' s, about ',F10.4,' mins')

!******************************************
!                                        *
! Generating log file PhiPsi_Win64.log   *
!                                        *
!******************************************
#ifdef sffortran
call log_startup ( 'PhiPsi_Win64.log' )
#endif

#ifdef sffortranmpi
call log_startup ( 'PhiPsi_Win64_MPI.log' )
#endif

#ifdef cbfortran
call log_startup ( 'PhiPsi_CBFortran.log' )
#endif

#ifdef linux
call log_startup ( 'PhiPsi_Win64.log' )
#endif

#ifdef ifort
call log_startup ( 'PhiPsi_Intel_Compiler.log' )
#endif

#ifdef ifx
call log_startup ( 'PhiPsi_Intel_Compiler_ifx.log' )
#endif

#ifdef Silverfrost
call log_startup ( 'PhiPsi_Silverfrost_Compiler.log' )
#endif


call log_configure ( "writeonstdout" , .false. )
call log_msg('------------------------------------------')
call log_msg('/                                        /')
call log_msg('/           PhiPsi Runing Log            /') 
call log_msg('/                                        /')
call log_msg('------------------------------------------')


!************************************
!                                  *
!             Welcome              *
!                                  *
!************************************
CALL Welcome

!************************************
!                                  *
!  Output PhiPsi_Win64.INFO file.  *
!  NEWFTU2023072901.               *
!                                  *
!************************************
open(101,file='PhiPsi_Win64.INFO',status='unknown')     
    write(101,*)  trim(PhiPsi_Version)
    write(101,*)  trim(PhiPsi_Release_Date)
close(101)    

!************************************
!                                  *
!    System Type. 2023-08-10.      *  
!    NEWFTU2023081001.             *
!                                  *
!************************************
#ifndef Silverfrost
    Operation_System_Type = get_os_type()
    if(Operation_System_Type==1)then
        write(*,'(*(g0))')  '     Operation system: Linux ',bit_size(0_C_INTPTR_T),'-bit.'
        String_Connector ='/'
    elseif(Operation_System_Type==2)then
        write(*,'(*(g0))')  '     Operation system: Mac OS ',bit_size(0_C_INTPTR_T),'-bit.' 
        String_Connector ='/'
    elseif(Operation_System_Type==3)then
        write(*,'(*(g0))')  '     Operation system: Windows ',bit_size(0_C_INTPTR_T),'-bit.'
        String_Connector ='\'
    elseif(Operation_System_Type==0)then
        write(*,'(*(g0))')  '     Operation system: Unknown.'   
    endif
#endif

#ifdef Silverfrost
    print *,'    Operation system: Windows.'
    Operation_System_Type=3
    String_Connector ='\'
#endif


!************************************
!                                  *
!           Login user.            *
!                                  *
!************************************
#ifndef Silverfrost
CALL GETLOG(login_user)
print *,'    Username: ',trim(login_user),'.'
#endif

!*****************************************
!                                       *
!Get current work directory of PhiPsi.  *
!                                       *
!*****************************************
#ifndef Silverfrost
CALL getcwd(PhiPsi_Current_Directory) 
#endif
#ifdef Silverfrost
!For the Silverfrost compiler, since getcwd is not supported,
!the current working directory is obtained using system.
cwd_output_file = 'cwd_output.txt'
command = 'cd > ' // trim(cwd_output_file)
!Call the Windows cd command using the system command and save the 
!result to a file
call system(trim(command))
!Open file to read the working directory
open(unit=10, file=trim(cwd_output_file), status='old', action='read') 
read(10, '(A)') PhiPsi_Current_Directory
close(10)
command = 'del ' // trim(cwd_output_file)
call system(trim(command))
#endif
print *,'    PhiPsi directory: ',trim(PhiPsi_Current_Directory),'.'

!************************************
!                                  *
!      Get initial timestamp.      *
!                                  *
!************************************
call Tool_Get_Current_Time(current_data,date_time,S_time)
WRITE(*,1000)
WRITE(*,1001) current_data(1:4),current_data(5:6),current_data(7:8),date_time(5),date_time(6),date_time(7)
WRITE(*,1000)
!log file
WRITE(temp_log,2001)current_data,date_time(5),date_time(6),date_time(7)
call Tool_chrpak_s_blank_delete(temp_log(13:))
call log_msg(temp_log)



!************************************
!                                  *
!           Inilization.           *
!                                  *
!************************************
call Initialize

!************************************
!                                  *
!          Read information.       *
!                                  *
!************************************
! OPTION 1: Directly set the relevant control parameters 
!           in the Fortran source code.
!           IMPROV2022101501.
!call PhiPsi_Input_Keywords  !Read keywords.
!call PhiPsi_Input_Filename  !Read filename.

!OPTION 2: Read from the keywords file. 
! Read relevant control parameters from the keyword file *.kpp
call PhiPsi_Read_Input

!*****************************************
!                                       *
! Save the current working path to the  *
! directory where PhiPsi is located.    *
! Saved in the current_folder.dat file. *
!                                       *
!*****************************************
call Save_current_folder_dat_file

!*****************************************
!                                       *
!Clear all results file if necessary.   *
!                                       *
!*****************************************
if (Key_Clear_All==1)  then
  ! OPTION 1: Use the filesys open-source library in the Simply Fortran compiler
  !           Only applicable to lower versions of the gfortran compiler.
  !call Clear_All_Results_Files

  ! OPTION 2: Using Python 3, suitable for all kinds of compilers. 2022-07-25.
  !          NEWFTU2022072501.
#ifndef moacos
  call Clear_All_Results_Files_v2
#endif
endif

!*****************************************
!                                       *
!    Visit PhiPsi.top if necessary.     *
!                                       *
!*****************************************
if(Key_Visit_PhiPsi_top==1) then
    call Tool_Internet_by_using_CURL
endif


!************************************************
!                                              *
!  Preparation for openMP parallel computing.  *
!                                              *
!************************************************
print *,' '
!$omp parallel shared(max_nthreads)
nthreads = omp_get_thread_num( )
if (nthreads .eq. 0 ) then
  write ( *, '(a,i4)' ) '     OpenMP, the number of processors available is',  omp_get_num_procs( )
  write ( *, '(a,i4)' ) '     OpenMP, the number of threads available is',  omp_get_num_threads( )
  max_nthreads = omp_get_num_threads( )
end if
!$omp end parallel

#ifndef Silverfrost
if (Key_Num_Process==99) then
  call omp_set_num_threads(max_nthreads)
  write ( *, '(a,i4)' ) '     OpenMP, number of threads has been set as',  max_nthreads
else
  if(Key_Num_Process<=max_nthreads) then
      call omp_set_num_threads(Key_Num_Process)
      write ( *, '(a,i4)' ) '     OpenMP, number of threads has been set as',  Key_Num_Process
  else
      call omp_set_num_threads(max_nthreads)
      write ( *, '(a,i4)' ) '     OpenMP, number of threads has been set as',  max_nthreads
  endif
endif
#endif

print *,' '



!************************************************
!                                              *
!          Read model information.             *
!                                              *
!************************************************
Flag_Read_Geo = .True.
! If it is rigid body dynamics analysis of small balls or molecular
! dynamics analysis, there is no need to input any information.
if(Key_Analysis_Type==11 .or.Key_Analysis_Type==21)then
  Flag_Read_Geo = .False.
endif

#ifdef gfortran
#ifndef github
!If PD analysis.
if(Key_Analysis_Type==31)then
  CALL PD_Read_Geo
  Flag_Read_Geo = .False.
endif
#endif
#endif

!If Impact analysis, 2021-07-23.
if(Key_Analysis_Type==12)then
  Flag_Read_Geo = .False.
endif

!If PFEM solver, 2021-12-19.
if(Key_Analysis_Type==51 .or. Key_Analysis_Type==52 )then
  Flag_Read_Geo = .False.
endif

!************************************************
!                                              *
!    If FEM or XFEM analysis. call Read_Geo    *
!    or Read_Geo_3D.                           *
!                                              *
!************************************************
if (Flag_Read_Geo .eqv. .True.)then
    if (Key_Dimension == 2) then
        CALL Read_Geo
    elseif(Key_Dimension == 3) then
        CALL Read_Geo_3D
    end if
endif

!************************************************
!                                              *
!       Unit conversion if necessary.          *
!                                              *
!************************************************
call Unit_Conversion

!************************************************
!                                              *
!    Check and display input information.      *
!                                              *
!************************************************
CALL Input_Check_Display

!************************************************
!                                              *
!    Save *.post file for post-processing.     *
!                                              *
!************************************************
call Save_Info_for_Matlab_Post


!************************************************
!                                              *
!                                              *
!                                              *
!                 2D problems                  *
!                                              *
!                                              *
!                                              *
!************************************************
print *, " "
print *, " >> SOLUTION STAGE..."
if (Key_Dimension == 2) then
  !~~~~~~~~~~~~~~~~~~~~~~~~~
  !(1)Quasi-static Analysis
  !~~~~~~~~~~~~~~~~~~~~~~~~~
  if (Key_Analysis_Type==1) then
      ! If it is a cohesive crack
      if(Key_TipEnrich ==4)then
#ifndef github
          !If it is a static load analysis (all loads are applied
          !at once and then remain unchanged)
          if(Key_Force_Control/=4)then
            print *,'    Main program: PhiPsi2D_Static_Cohesive'
            call log_msg('Main program: PhiPsi2D_Static_Cohesive')
            call log_msg('Performing simulation...')
            CALL PhiPsi2D_Static_Cohesive
            
          !If it is the cohesive crack load control algorithm 1
          elseif(Key_Force_Control==4)then
            print *,'    Main program: PhiPsi2D_Static_Cohesive_FORCE_Control'
            call log_msg('Main program: PhiPsi2D_Static_Cohesive_FORCE_Control')
            call log_msg('Performing simulation...')
            CALL PhiPsi2D_Static_Cohesive_FORCE_Control
          endif
#endif            
      !Normal Quasi-static Analysis
      else
          !if *Key_Force_Control==5(make sure that only one crack 
          !propagates druing one step)
          if(Key_Force_Control==5) then
#ifndef github
              if (CFCP==1)then
                  print *,'    Main program: PhiPsi2D_Static_FORCE_Control_CFCP1'
                  call log_msg('Main program: PhiPsi2D_Static_FORCE_Control_CFCP1')
                  call log_msg('Performing simulation...')
                  CALL PhiPsi2D_Static_FORCE_Control_CFCP1
              elseif(CFCP==2)then
                  print *,'    Main program: PhiPsi2D_Static_FORCE_Control_CFCP2'
                  call log_msg('Main program: PhiPsi2D_Static_FORCE_Control_CFCP2')
                  call log_msg('Performing simulation...')
                  CALL PhiPsi2D_Static_FORCE_Control_CFCP2
              endif
#endif              
          else
              print *,'    Main program: PhiPsi2D_Static'
              call log_msg('Main program: PhiPsi2D_Static')
              call log_msg('Performing simulation...')
              CALL PhiPsi2D_Static
          endif
      endif

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (2) Implicit Dynamic Analysis
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==2) then
#ifndef github
      print *,'    Main program: PhiPsi2D_I_Dynamic'
      call log_msg('Main program: PhiPsi2D_I_Dynamic')
      call log_msg('Performing simulation...')
      CALL PhiPsi2D_I_Dynamic
#endif  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (3) Hydraulic Fracturing Static Analysis
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==3) then
      !///////////////////////////////////////////////////////////////////////////////
      ! Given the extended step size, mass conservation is implemented in the Newton 
      ! iteration step, with delta_time constantly changing. Overall, this scheme
      ! The computational load is small, but the stability is relatively low, and the
      ! solution is very sensitive to the initial guess.
      ! Besides, line searching scheme is introduced to improve staility.
      !///////////////////////////////////////////////////////////////////////////////
      if (Key_HF_Secant_TS==0)then
          print *,'    Main program: PhiPsi2D_Static_HF'
          call log_msg('Main program: PhiPsi2D_Static_HF')
          call log_msg('Performing simulation...')
          call PhiPsi2D_Static_HF

      !///////////////////////////////////////////////////////////////////////////////
      ! Given an extended step size, add a layer of secant method iteration to find 
      ! the time step that best satisfies mass conservation, namely the N-R iteration
      ! process. The time step remains unchanged, and the overall calculation in this
      ! scheme may be slightly larger, but the stability is relatively high.
      ! and the solution is not sensitive to the initial guess.
      !///////////////////////////////////////////////////////////////////////////////
      elseif (Key_HF_Secant_TS==1)then
#ifndef github
         print *,'    Main program: PhiPsi2D_Static_HF_with_SIM'
         call log_msg('Main program: PhiPsi2D_Static_HF_with_SIM')
         call log_msg('Performing simulation...')
         call PhiPsi2D_Static_HF_with_SIM
#endif  
      endif
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (4) Test procedure applying only uniformly 
  !     distributed water pressure inside the crack
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==4) then
#ifndef github
      print *,'    Main program: PhiPsi2D_Static_HF_SlipWater'
      call log_msg('Main program: PhiPsi2D_Static_HF_SlipWater')
      call log_msg('Performing simulation...')
      call PhiPsi2D_Static_HF_SlipWater
#endif  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (5) Theoretical water pressure of 
  !     viscosity-dominated problems applied within cracks
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==5) then
#ifndef github
     print *,'    Main program: PhiPsi2D_Static_WP_TheoVisDomi'
     call log_msg('Main program: PhiPsi2D_Static_WP_TheoVisDomi')
     call log_msg('Performing simulation...')
     call PhiPsi2D_Static_WP_TheoVisDomi
#endif  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (6) Explicit Dynamic Analysis
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==6) then
#ifndef github
      print *,'    Main program: PhiPsi2D_E_Dynamic'
      call log_msg('Main program: PhiPsi2D_E_Dynamic')
      CALL PhiPsi2D_E_Dynamic
#endif 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (7) Nonlinear analysis (plastic materials, etc.)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==7) then
#ifndef github
      print *,'    Main program: PhiPsi2D_Static_NonLinear'
      call log_msg('Main program: PhiPsi2D_Static_NonLinear')
      call log_msg('Performing simulation...')
      call PhiPsi2D_Static_NonLinear
#endif
  !~~~~~~~~~~
  !(8) blank
  !~~~~~~~~~~
  elseif (Key_Analysis_Type==8) then


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (11) Multi-body dynamics simulation of randomly generated balls
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==11) then
#ifndef github
      print *,'    Main program: PhiPsi2D_Rigid_Balls_Contact'
      call log_msg('Main program: PhiPsi2D_Rigid_Balls_Contact')
      call log_msg('Performing simulation...')
      call PhiPsi2D_Rigid_Balls_Contact
#endif
  !~~~~~~~~~~~~~~~~~~~~~~~~
  !(12) Impact, 2021-07-23
  !~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==12) then
#ifdef gfortran  
#ifndef github
      print *,'    Main program: PhiPsi_Impact'
      call log_msg('Main program: PhiPsi_Impact')
      call log_msg('Performing simulation...')
      call PhiPsi_Impact
#endif
#endif
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (15) Quasi-static field problems.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==15) then
#ifndef github
      print *,'    Main program: PhiPsi2D_Static_Field_Problem'
      call log_msg('Main program: PhiPsi2D_Static_Field_Problem')
      call log_msg('Performing simulation...')
      call PhiPsi2D_Static_Field_Problem
#endif
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (16) Basic Transient Field Problems, A First Course in the Finite 
  !      Element Method, 4th Edition, p. 686
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==16) then
#ifndef github
      print *,'    Main program: PhiPsi2D_I_Dynamic_Field_Problem'
      call log_msg('Main program: PhiPsi2D_I_Dynamic_Field_Problem')
      call log_msg('Performing simulation...')
      call PhiPsi2D_I_Dynamic_Field_Problem
#endif
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (17) Transient field problems, shale gas and other output 
  !      assessments, coupled with deformation field analysis.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==17) then
#ifndef github
      print *,'    Main program: PhiPsi2D_I_Dynamic_Field_Problem_Coupled_Deform'
      call log_msg ('Main program: PhiPsi2D_I_Dynamic_Field_Problem_Coupled_Deform')
      call log_msg('Performing simulation...')
      call PhiPsi2D_I_Dynamic_Field_Problem_Coupled_Deform
#endif
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (21) 2D Molecular Dynamics Simulation
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==21) then
#ifdef gfortran
#ifndef github
      print *,'    Main program: PhiPsi_Molecular_Dynamics'
      call log_msg('Main program: PhiPsi_Molecular_Dynamics')
      call log_msg('Performing simulation...')
      call PhiPsi_Molecular_Dynamics
#endif
#endif
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !(31) Peridynamic simulation
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==31) then
#ifdef gfortran
#ifndef github  
      print *,'    Main program: PhiPsi_Peridynamic'
      call log_msg('Main program: PhiPsi_Peridynamic')
      call log_msg('Performing simulation...')
      call PhiPsi_Peridynamic
#endif      
#endif 
  !~~~~~~~~~~~~~
  ! (41) blank.
  !~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (51) PFEM Solver -- Solid Element Structural Analysis, 2021-12-19.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==51) then
#ifndef Silverfrost
#ifndef github
      print *,'    Main program: PhiPsi_PFEM_p56'
      call log_msg('Main program: PhiPsi_PFEM_p56')
      call log_msg('Performing simulation...')
      call PhiPsi_PFEM_p56
#endif 
#endif 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (52) PFEM Solver -- Solid Element Structural Analysis with Abaqus
  !      Umat, 2021-12-19.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==52) then
#ifndef Silverfrost
#ifndef github
      print *,'    Main program: PhiPsi_PFEM_p57'
      call log_msg('Main program: PhiPsi_PFEM_p57')
      call log_msg('Performing simulation...')
      call PhiPsi_PFEM_p57
#endif 
#endif 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (53) PFEM Solver -- Structure-Flow Direct Coupling, 2021-12-19.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==53) then
#ifndef Silverfrost
#ifndef github
      print *,'    Main program: PhiPsi_PFEM_p94'
      call log_msg('Main program: PhiPsi_PFEM_p94')
      call log_msg('Performing simulation...')
      call PhiPsi_PFEM_p94
#endif 
#endif 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (54) PFEM Solver -- Structure-Flow Direct Coupling, EBE PCG, 2021-12-19.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==54) then
#ifndef Silverfrost
#ifndef github
      print *,'    Main program: PhiPsi_PFEM_p95'
      call log_msg('Main program: PhiPsi_PFEM_p95')
      call log_msg('Performing simulation...')
      call PhiPsi_PFEM_p95
#endif 
#endif 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (55) PFEM Solver -- Structure-Flow Direct Coupling - Mohr-Coulomb Plasticity, 
  !                     Global Stiffness Matrix, 2021-12-20
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==55) then
#ifndef Silverfrost
#ifndef github
      print *,'    Main program: PhiPsi_PFEM_p96'
      call log_msg('Main program: PhiPsi_PFEM_p96')
      call log_msg('Performing simulation...')
      call PhiPsi_PFEM_p96
#endif 
#endif 
  end if

!************************************************
!                                              *
!                                              *
!                                              *
!                 3D problems                  *
!                                              *
!                                              *
!                                              *
!************************************************
elseif (Key_Dimension == 3) then
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Conventional 3D Static Analysis.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (Key_Analysis_Type==1) then
      ! For New Australia Use. 2022-09-10.
      if(Key_XA ==1)then
#ifndef github
          print *,'    Main program: PhiPsi3D_Static_XA'
          call log_msg('Main program: PhiPsi3D_Static_XA')
          call log_msg('Performing simulation...')
          call PhiPsi3D_Static_XA
#endif 
      ! Exclusive for Xinao. 2023-03-13. Used for preparing Fortran library.
      elseif(Key_XA ==2)then
#ifndef github
          print *,'    Main program: PhiPsi3D_Static_XA_v2'
          call log_msg('Main program: PhiPsi3D_Static_XA_v2')
          call log_msg('Performing simulation...')
          call PhiPsi3D_Static_XA_v2
#endif 
      else
          print *,'    Main program: PhiPsi3D_Static'
          call log_msg('Main program: PhiPsi3D_Static')
          call log_msg('Performing simulation...')
          CALL PhiPsi3D_Static
      endif
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (3) Hydraulic Fracturing Static Analysis, 2021-12-03.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==3) then
#ifndef github
      ! Picard Iterative 3D Hydraulic Fracturing Analysis
      print *,'    Main program: PhiPsi3D_Static_HF_Picard'
      call log_msg('Main program: PhiPsi3D_Static_HF_Picard')
      call log_msg('Performing simulation...')
      call PhiPsi3D_Static_HF_Picard
#endif 

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (4) 3D water slide pressure. Normal water pressure. 
  !     Dual loop and bisection method. 2022-06-28.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==4) then
      print *,'    Main program: PhiPsi3D_Static_HF_SlipWater'
      call log_msg('Main program: PhiPsi3D_Static_HF_SlipWater')
      call log_msg('Performing simulation...')
      call PhiPsi3D_Static_HF_SlipWater

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (6) Explicit dynamic analysis, 2021-08-23
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==6) then
#ifndef github
      print *,'    Main program: PhiPsi3D_E_Dynamic'
      call log_msg('Main program: PhiPsi3D_E_Dynamic')
      call log_msg('Performing simulation...')
      CALL PhiPsi3D_E_Dynamic
#endif 

  !~~~~~~~~~~~
  !(8) Blank.
  !~~~~~~~~~~~
  elseif (Key_Analysis_Type==8) then

  !~~~~~~~~~~~~~~~~~~~~~~~~
  !(12) Impact, 2021-07-23
  !~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==12) then
#ifdef gfortran  
#ifndef github
      print *,'    Main program: PhiPsi_Impact'
      call log_msg('Main program: PhiPsi_Impact')
      call log_msg('Performing simulation...')
      call PhiPsi_Impact
#endif
#endif

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (21) 3D Molecular Dynamics Simulation
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==21) then
#ifdef gfortran 
#ifndef github 
      print *,'    Main program: PhiPsi_Molecular_Dynamics'
      call log_msg('Main program: PhiPsi_Molecular_Dynamics')
      call log_msg('Performing simulation...')
      call PhiPsi_Molecular_Dynamics
#endif
#endif

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (51) PFEM Solver -- Solid Element Structural Analysis, 2021-12-19
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==51) then
#ifndef Silverfrost
#ifndef github 
      print *,'    Main program: PhiPsi_PFEM_p56'
      call log_msg('Main program: PhiPsi_PFEM_p56')
      call log_msg('Performing simulation...')
      call PhiPsi_PFEM_p56
#endif
#endif
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (52) PFEM Solver -- Solid Element Structural Analysis with Abaqus 
  !      Umat, 2021-12-19
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==52) then
#ifndef Silverfrost
#ifndef github 
      print *,'    Main program: PhiPsi_PFEM_p57'
      call log_msg('Main program: PhiPsi_PFEM_p57')
      call log_msg('Performing simulation...')
      call PhiPsi_PFEM_p57
#endif
#endif
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (53) PFEM Solver -- Structure-Flow Direct Coupling, 2021-12-19
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==53) then
#ifndef Silverfrost
#ifndef github 
      print *,'    Main program: PhiPsi_PFEM_p94'
      call log_msg('Main program: PhiPsi_PFEM_p94')
      call log_msg('Performing simulation...')
      call PhiPsi_PFEM_p94
#endif
#endif
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (54) PFEM Solver -- Structure-Flow Direct Coupling, 2021-12-19.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==54) then
#ifndef Silverfrost
#ifndef github 
      print *,'    Main program: PhiPsi_PFEM_p95'
      call log_msg('Main program: PhiPsi_PFEM_p95')
      call log_msg('Performing simulation...')
      call PhiPsi_PFEM_p95
#endif
#endif
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! (55) PFEM Solver -- Structure-Flow Direct Coupling - Mohr-Coulomb Plasticity, 
  !      Global Stiffness Matrix, 2021-12-20
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==55) then
#ifndef Silverfrost
#ifndef github 
      print *,'    Main program: PhiPsi_PFEM_p96'
      call log_msg('Main program: PhiPsi_PFEM_p96')
      call log_msg('Performing simulation...')
      call PhiPsi_PFEM_p96
#endif
#endif
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !(61) 3D Fracturing experiment. 2023-01-23.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  elseif (Key_Analysis_Type==61) then
#ifndef github 
      print *,'    Main program: PhiPsi3D_Static_Fracturing_Experiment'
      call log_msg('Main program: PhiPsi3D_Static_Fracturing_Experiment')
      call log_msg('Performing simulation...')
      call PhiPsi3D_Static_Fracturing_Experiment
#endif 
  end if
end if

!************************************
!                                  *
!      Get endding timestamp.      *
!                                  *
!************************************
call Tool_Get_Current_Time(current_data,date_time,F_time)
print *,' '
WRITE(*,1002) current_data(1:4),current_data(5:6),current_data(7:8),date_time(5),date_time(6),date_time(7)
WRITE(*,1003) F_time-S_time,(dble(F_time)-dble(S_time))/Con_60
!log file
WRITE(temp_log,2002)current_data,date_time(5),date_time(6),date_time(7)
call Tool_chrpak_s_blank_delete(temp_log(13:))
call log_msg(temp_log)
WRITE(temp_log,2003) F_time-S_time,(dble(F_time)-dble(S_time))/Con_60
!call Tool_chrpak_s_blank_delete(temp_log)  !delete blanks
call log_msg(temp_log)
call log_msg(' ')

!***************************************
!                                     *
! Clear memory variables. 2022-08-04. *
!                                     *
!***************************************
call Clear_Memory

!************************************
!                                  *
!        Play ending sound.        *
!                                  *
!************************************
if(Key_Play_Sounds==1)then
  !----------
  ! option 1
  !----------
#ifndef macos  
#ifndef github 
  call system('cd .\Python_Tools\ && PhiPsi_Play_Done.py')
#endif
#endif
  !----------
  ! option 2
  !----------
  !beep = char(7)
  !write (*, 10100) beep
endif

!***************************************
!                                     *
! Closing log file PhiPsi_Win64.log.  *
!                                     *
!***************************************
call log_shutdown()


!******************************************
!                                        *
! Interaction at the end of the program. *
!                                        *
!******************************************
! Do not close the window after the program ends; wait for the user's response.
if (Key_Close_Window==0) then
  print *,' '
  write( *, * ) '    Press any key to exit PhiPsi Kernal.'
  read( *, * )
! Close the window immediately after the program ends
elseif(Key_Close_Window==1)then
  !stop
  !IMPROV2022101401. Comment out stop, otherwise it will produce abnormal output:
  !                  Note: The following floating-point exceptions are 
  !                  signalling: IEEE_INVALID_FLAG IEEE_DENORMAL
  !Ref: https://stackoverflow.com/questions/44308577/ieee-underflow-flag-ieee-denormal-in-fortran-77
endif


END PROGRAM PhiPsi



