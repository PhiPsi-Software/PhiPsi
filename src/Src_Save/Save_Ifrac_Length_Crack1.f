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
 
      SUBROUTINE Save_Ifrac_Length_Crack1(ifra_Length,Num_Frac)
c     Save the length of Fracture 1 for each fracture step
c     upon completion of all fracture calculations.

      use Global_Float_Type      
      use Global_Common
      use Global_Filename
      use Global_POST
      
      implicit none
      
      real(kind=FT),intent(in):: ifra_Length(Num_Frac)
      integer,intent(in):: Num_Frac
      integer i
      character(200) c_File_name_1
      
      if (Key_Save_Nothing==1) return
      
      print *,'    Saving length (m) of crack 1 of each time step...'
      c_File_name_1   =  trim(Full_Pathname)//'.ilth'   
      open(201,file=c_File_name_1,status='unknown') 
      write(201,*) '    Length of crack 1 of each time step' 
      do i=1,Num_Frac
          if(ifra_Length(i)>ZR)then
              write(201, '(F12.5)') ifra_Length(i)
          endif
      end do
      close(201)   
 
      RETURN
      END SUBROUTINE Save_Ifrac_Length_Crack1