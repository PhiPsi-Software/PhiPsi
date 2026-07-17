!-----------------------------------------------------------
! Brief: Print the program banner, version, and compiler info.
!
! Notes:   Sets the PhiPsi version and release-date globals and writes
!          an ASCII banner plus compiler identification to the console.
!-----------------------------------------------------------

SUBROUTINE Welcome
!This subroutine displays information about the program.

use Global_Common
use iso_fortran_env

implicit none
character(120) :: print_compiler,print_release_date,print_phipsi_version

!=============Modify the following part=============
PhiPsi_Version      = '1.59'
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
!PhiPsi_Release_Date = 'November 21, 2025' 
!PhiPsi_Release_Date = 'December 28, 2025'
!PhiPsi_Release_Date = 'January 23, 2026'
!PhiPsi_Release_Date = 'Feburary 15, 2026' 
!PhiPsi_Release_Date = 'May 21, 2026     ' 
PhiPsi_Release_Date = 'July 15, 2026    '  
!=======================END=========================

print_phipsi_version=repeat(' ',4)//'****'//repeat(' ',16)//'Version '// trim(PhiPsi_Version)// &
repeat(' ',22-len(trim(PhiPsi_Version)))//'****'
                     
print_release_date=repeat(' ',4)//'>'//repeat(' ',1)//'Release:  '// trim(PhiPsi_Release_Date)// &
repeat(' ',41-len(trim(PhiPsi_Release_Date)))//'<'

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
print *,"    > Website:  phipsi.top, copyright(c) 2016-2026       <"
print *,"    ======================================================"
print *,"                                                          "

RETURN
END SUBROUTINE Welcome
