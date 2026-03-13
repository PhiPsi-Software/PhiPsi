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
 
      SUBROUTINE Matrix_Count_Row_Int_m_x_2(m,Matrix,
     &                                      Count_Row,Num_Once,
     &                                      key_process)   
      ! Count the occurrences of each row in the matrix (m,2) and store them in Count_Row.
      ! 2020-03-27: Add openMP, disable key_process.    

      use Global_Float_Type    
      use CMD_Progress
      use Global_Common
#ifndef Silverfrost
      use omp_lib
#endif 
      implicit none
      
      integer,intent(in)::m
      integer,intent(in)::key_process
      integer,intent(in)::Matrix(m,2)
      integer,intent(out)::Count_Row(m)
      integer,intent(out)::Num_Once
      
      integer i
      
      type(CLS_CMD_Progress)::Progress

      Count_Row(1:m) = 0
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(static)  
      do i=1,m    
          Count_Row(i) =count((Matrix(1:m,1) .EQ. Matrix(i,1)) .and.
     &                        (Matrix(1:m,2) .EQ. Matrix(i,2)))
      end do      
      !$omp end parallel do 
      
      Num_Once= 0
      Num_Once=count(Count_Row .eq. 1)
      
      return
      END SUBROUTINE Matrix_Count_Row_Int_m_x_2
    


