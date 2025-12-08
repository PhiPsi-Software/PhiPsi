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
 
SUBROUTINE Welcome
!This subroutine displays information about the program.

use Global_Common
use iso_fortran_env

implicit none
character(120) :: print_compiler,print_release_date,print_phipsi_version

!=============Modify the following part=============
PhiPsi_Version      = '1.56.16'
!PhiPsi_Release_Date = 'December 24, 2024'
!PhiPsi_Release_Date = 'January 27, 2024 '
!PhiPsi_Release_Date = 'Feburary 26, 2024'     
!PhiPsi_Release_Date = 'March 21, 2024   '
!PhiPsi_Release_Date = 'April 11, 2024   '
!PhiPsi_Release_Date = 'July 25, 2024    '
!PhiPsi_Release_Date = 'September 7, 2024'
!PhiPsi_Release_Date = 'October 4, 2024  '
!PhiPsi_Release_Date = 'November 7, 2024' 
!PhiPsi_Release_Date = 'January 26, 2025' 
PhiPsi_Release_Date = 'November 21, 2024' 
!=======================END=========================

print_phipsi_version=repeat(' ',4)//'****'//repeat(' ',16)//'Version '//&
                     trim(PhiPsi_Version)//&
                     repeat(' ',22-len(trim(PhiPsi_Version)))//'****'
                     
print_release_date=repeat(' ',4)//'>'//repeat(' ',1)//'Release:  '//&
                   trim(PhiPsi_Release_Date)//&
                   repeat(' ',41-len(trim(PhiPsi_Release_Date)))//'<'

#ifndef Silverfrost
compiler            = compiler_version()
if(compiler(1:3) == 'GCC')then
    print_compiler  = '    > Compiler: '//trim(compiler(1:19))//'                       <'
    Compiler_Type   = 1
elseif(compiler(1:25) == 'Intel(R) Fortran Intel(R)') then
    print_compiler  = '    > Compiler: '//trim(compiler(1:16))//' classical '//compiler(96:105)//'    <'
    Compiler_Type   = 2
elseif(compiler(1:25) == 'Intel(R) Fortran Compiler') then
    print_compiler  = '    > Compiler: '//trim(compiler(1:16))//' IFX '//compiler(76:84)//'           <'
    Compiler_Type   = 3
endif
#endif

#ifdef Silverfrost  
    compiler            = COMPILER_VERSION() 
    print_compiler  = '    > Compiler: Silverfrost Fortran '//trim(compiler(1:11))//'          <'
    Compiler_Type   = 4
#endif


print *,"                                                          "
print *,"    ======================================================"
print *,"    ****                                              ****"
print *,"    ****          ___  _      _    ___       _        ****"
print *,"    ****         / _ \| |__  (_)  / _ \ ___ (_)       ****"
print *,"    ****        / /_)/| '_ \ | | / /_)// __|| |       ****"
print *,"    ****       / ___/ | | | || |/ ___/ \__ \| |       ****"
print *,"    ****       \/     |_| |_||_|\/     |___/|_|       ****"
print *,"    ****                                              ****"
print *, trim(print_phipsi_version)
print *,"    ****                                              ****"
print *,"    ======================================================"
print *,"    > About:    A Fortran Program to Simulate 2D and 3D  <"
print *,"    >           Crack Propagation Problems, Nonlinear    <"
print *,"    >           Material Problems, and Explicit Dynamic  <"
print *,"    >           Problems based Mainly on FEM and XFEM.   <"
print *, trim(print_release_date)
print *, trim(print_compiler)
print *,"    > Author:   Fang Shi, shifang@hyit.edu.cn            <"
print *,"    > Website:  phipsi.top, copyright(c) 2016-2025       <"
print *,"    ======================================================"
print *,"                                                          "

RETURN
END SUBROUTINE Welcome
