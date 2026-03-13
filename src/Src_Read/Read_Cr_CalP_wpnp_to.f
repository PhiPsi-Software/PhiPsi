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
 
      SUBROUTINE Read_Cr_CalP_Wpnp_to(Result_File_Num,Variable,m,n)
      ! Read the wpnp (Width of Proppant No Pressure) file, which contains
      ! Crack opening of the support crack under an uncompressed state (obtained from fluid-solid coupling
      ! analysis)
      ! Used to calculate the opening of support cracks
      
      use Global_Float_Type      
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Filename
      
      implicit none
      
      integer,intent(in)::m,n,Result_File_Num
      real(kind=FT),intent(out)::Variable(m,n)
      
      integer i,j
      character(200) c_File_name_1
      character(5) temp
                     
      print *,'    Reading wpnp of crack points......'
      write(temp,'(I5)') Result_File_Num
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'wpnp'//'_'//ADJUSTL(temp)            
      
     
      select case(Key_Data_Format)
      case(1)
          open(101,file=c_File_name_1,status='unknown')         
          do i=1,num_Crack
              read(101,'(2000E20.12)') (Variable(i,
     &                                   j),j=1,Cracks_CalP_Num(i))
          end do
          close(101) 
      case(2)
          open(101,file=c_File_name_1,status='unknown',
     &                                form='unformatted')         
          do i=1,num_Crack
              read(101) (Variable(i,j),j=1,Cracks_CalP_Num(i))
          end do
          close(101) 
      end select
      
      RETURN
      END SUBROUTINE Read_Cr_CalP_Wpnp_to
