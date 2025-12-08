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
 
SUBROUTINE Set_OpenMP
! Setting up OpenMP. 2023-03-25.

use Global_Float_Type
use Global_Common
use omp_lib
use iso_fortran_env
use ISO_C_BINDING

implicit none
integer nthreads

!-------------------------------------------
!Preparation for openMP parallel computing.
!-------------------------------------------
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
  write ( *, '(a,i4)' ) '     OpenMP, number of processors has been set as',  max_nthreads
else
  if(Key_Num_Process<=max_nthreads) then
      call omp_set_num_threads(Key_Num_Process)
      write ( *, '(a,i4)' ) '     OpenMP, number of processors has been set as',  Key_Num_Process
  else
      call omp_set_num_threads(max_nthreads)
      write ( *, '(a,i4)' ) '     OpenMP, number of processors has been set as',  max_nthreads              
  endif
endif
#endif
print *,' '      


RETURN
END SUBROUTINE Set_OpenMP
