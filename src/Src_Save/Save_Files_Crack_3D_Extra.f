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
 
      SUBROUTINE Save_Files_Crack_3D_Extra(isub)
      use Global_Float_Type
      use Global_Common
      use Global_Crack_Common
      use Global_Crack_3D
      use Global_Model
      use Global_Filename
      use Global_POST
      
      implicit none
      
      integer isub
      integer i,j
      character(200) c_File_name_1
      character(200) c_File_name_2
      character(200) c_File_name_3

      
      character(5) temp  
      real(kind=FT) tem_1  
      
      if (Key_Save_Nothing==1) return
      
      tem_1 = 1000.0D0
      
      select case(Key_Data_Format)
      case(1:2)
          write(temp,'(I5)') isub
          if(CFCP==2) then
              print *,'    Saving S1 vector of crack node...'
              c_File_name_1   =  trim(Full_Pathname)//'.cndx'
     &                                           //'_'//ADJUSTL(temp)     
              c_File_name_2   =  trim(Full_Pathname)//'.cndy'
     &                                           //'_'//ADJUSTL(temp)    
              c_File_name_3   =  trim(Full_Pathname)//'.cndz'
     &                                           //'_'//ADJUSTL(temp)  
              open(701,file=c_File_name_1,status='unknown') 
              open(702,file=c_File_name_2,status='unknown') 
              open(703,file=c_File_name_3,status='unknown')           
     &           
              do i=1,num_Crack
                  write(701, '(50000E20.12)') (Crack3D_Vector_S1(i)%row(
     &                             j,1),j=1,Crack3D_Meshed_Node_num(i))       
                  write(702, '(50000E20.12)') (Crack3D_Vector_S1(i)%row(
     &                             j,2),j=1,Crack3D_Meshed_Node_num(i))   
                  write(703, '(50000E20.12)') (Crack3D_Vector_S1(i)%row(
     &                             j,3),j=1,Crack3D_Meshed_Node_num(i))        
              end do
              close(701)          
              close(702)    
              close(703)    
          endif
      end select
      RETURN
      END SUBROUTINE Save_Files_Crack_3D_Extra
