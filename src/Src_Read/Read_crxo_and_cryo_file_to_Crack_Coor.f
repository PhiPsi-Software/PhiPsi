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
 
      SUBROUTINE Read_crxo_and_cryo_file_to_Crack_Coor(Result_File_Num)
      ! Read fracture coordinates and the number of fractures from crxo and cryo files
      ! Note the difference between crxo and cryo files and crax and cray files; the latter are the edge
      ! crack coordinates shifted into the model.
      
      use Global_Float_Type      
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Filename
      
      implicit none
      
      integer,intent(in)::Result_File_Num
      
      integer i,j,stat_read,i_C,i_Cr_P
      character(200) c_File_name_1,c_File_name_2
      character(5) temp
      real(kind=FT) Variable_x(100,1000),Variable_y(100,1000)
                     
      print *,'    Reading crox and croy file......'
      write(temp,'(I5)') Result_File_Num
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'crxo'//'_'//ADJUSTL(temp)       
      c_File_name_2   =  trim(Full_Pathname)//'.'
     &             //'cryo'//'_'//ADJUSTL(temp)   
     
      Variable_x(1:100,1:1000) = ZR
      Variable_y(1:100,1:1000) = ZR
      
      open(101,file=c_File_name_1,status='unknown')    
      open(102,file=c_File_name_2,status='unknown') 
      
      num_Crack = 1
      do    
          read(101,'(2000E20.12)',iostat=stat_read) (
     &              Variable_x(num_Crack,j),j=1,1000)
          if(stat_read/=0) exit
          num_Crack = num_Crack + 1
      enddo
      close(101) 
      
      num_Crack = 1
      do    
          read(102,'(2000E20.12)',iostat=stat_read) (
     &              Variable_y(num_Crack,j),j=1,1000)
          if(stat_read/=0) exit
          num_Crack = num_Crack + 1
      enddo
      close(102) 
      
      num_Crack = num_Crack -1
      
      do i_C = 1,num_Crack
          do i_Cr_P = 1,1000
              if(Variable_x(i_C,i_Cr_P)/= 0 .or.
     &           Variable_y(i_C,i_Cr_P)/= 0) then
                  Crack_Coor(i_C,i_Cr_P,1)=Variable_x(i_C,i_Cr_P)
                  Crack_Coor(i_C,i_Cr_P,2)=Variable_y(i_C,i_Cr_P)
                  Each_Cr_Poi_Num(i_C) = i_Cr_P
              endif
          enddo
      enddo
      
      
      RETURN
      END SUBROUTINE Read_crxo_and_cryo_file_to_Crack_Coor
