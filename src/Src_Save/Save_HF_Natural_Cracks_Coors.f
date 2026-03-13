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
 
      SUBROUTINE Save_HF_Natural_Cracks_Coors(c_num_Na_Cr,
     &                                        c_Na_Crack_Coor)
c     Save natural fracture coordinates (used for post-processing).

      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Model
      use Global_Filename
      use Global_POST
      
      implicit none
      integer,intent(in)::c_num_Na_Cr              
      real(kind=FT),intent(in)::c_Na_Crack_Coor(c_num_Na_Cr,2,2)      
      integer i
      character(200) c_File_name_1,c_File_name_2
      real(kind=FT) tem_1 
      
      tem_1 = 1000.0D0
      
      if (Key_Save_Nothing==1) return 
      
      print *,'    Saving coordinates of natural cracks...'
      c_File_name_1   =  trim(Full_Pathname)//'.ncrx'
      c_File_name_2   =  trim(Full_Pathname)//'.ncry'      
      open(101,file=c_File_name_1,status='unknown')    
      open(102,file=c_File_name_2,status='unknown')  
      if(Key_Unit_System==1)then     
          do i=1,c_num_Na_Cr
            write(101,'(E20.12,A,E20.12)')
     &            c_Na_Crack_Coor(i,1,1),',',c_Na_Crack_Coor(i,2,1)
            write(102,'(E20.12,A,E20.12)')
     &            c_Na_Crack_Coor(i,1,2),',',c_Na_Crack_Coor(i,2,2)
          end do
      elseif(Key_Unit_System==2)then 
          do i=1,c_num_Na_Cr
            write(101,'(E20.12,A,E20.12)')
     &      c_Na_Crack_Coor(i,1,1)/tem_1,',',
     &      c_Na_Crack_Coor(i,2,1)/tem_1
            write(102,'(E20.12,A,E20.12)')
     &      c_Na_Crack_Coor(i,1,2)/tem_1,',',
     &      c_Na_Crack_Coor(i,2,2)/tem_1
          end do
      endif
      
      close(101)
      close(102)   
      RETURN
      END SUBROUTINE Save_HF_Natural_Cracks_Coors
