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
 
      SUBROUTINE Save_Gauss_Disps(isub,Total_Num_G_P)
      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_DISP
      use Global_POST
      
      implicit none
      
      integer,intent(in)::isub,Total_Num_G_P
      integer i
      character(200) c_File_name_1
      character(5) temp   
      real(kind=FT) tem_1
      
      if (Key_Save_Nothing==1) return 
      
      tem_1 = 1000.0D0
      
      print *,'    Saving displacement of Gauss points...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'disg'//'_'//ADJUSTL(temp)   

      select case(Key_Data_Format)
      case(1)  
          open(201,file=c_File_name_1,status='unknown')   
          if(Key_Unit_System==1)then     
            do i=1,Total_Num_G_P
              write(201, '(I8,2E20.12)') i,DISP_x_Gauss(i),
     &                                     DISP_y_Gauss(i)
            end do  
          elseif(Key_Unit_System==2)then
            do i=1,Total_Num_G_P
              write(201, '(I8,2E20.12)') i,DISP_x_Gauss(i)/tem_1,
     &                                     DISP_y_Gauss(i)/tem_1
            end do  
          endif    
          close(201)  
      case(2) 
          open(203,file=c_File_name_1,status='unknown',
     &             form='unformatted',access='stream')     
          write(203) (DISP_x_Gauss(i),DISP_y_Gauss(i),
     &                  i=1,Total_Num_G_P)  
          close(203)         
      end select
 
      RETURN
      END SUBROUTINE Save_Gauss_Disps
