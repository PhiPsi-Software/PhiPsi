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
 
      SUBROUTINE Save_Files_Cross(isub)

      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_Cross
      use Global_POST
      
      implicit none
      
      integer isub
      integer i,j
      character(200) c_File_name_1
      character(5) temp   
      
      if (Key_Save_Nothing==1) return
      
      write(temp,'(I5)') isub
      
      select case(Key_Data_Format)

      case(1:2)
          print *,'    Saving coordinates of crosses...'
          write(temp,'(I5)') isub
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'cscr'   
          open(101,file=c_File_name_1,status='unknown')    
          
          do i=1,num_Cross
              write(101, '(2E20.12)') Cross_Point_RABCD(i,1,1:2)
          end do
          close(101)   
          
          print *,'    Saving enns file of cross...'
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'enns'//'_'//ADJUSTL(temp)       
          open(103,file=c_File_name_1,status='unknown')         
          do i=1,Num_Node
              write(103, '(200I10)') (Enriched_Node_Type_Cross(i,
     &                                    j),j=1,num_Cross)
          end do
          close(103)     
          
          print *,'    Saving elts file of hole...'
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'elts'//'_'//ADJUSTL(temp)       
          open(104,file=c_File_name_1,status='unknown')         
          do i=1,Num_Elem
              write(104, '(200I10)') 
     &                    (Elem_Type_Cross(i,j),j=1,num_Cross)
          end do
          close(104)  
          
          print *,'    Saving poss file of hole...'
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'poss'//'_'//ADJUSTL(temp)         
          open(110,file=c_File_name_1,status='unknown')      
          do i=1,Num_Node
              write(110, '(200I10)') (c_POS_Cross(i,j),j=1,num_Cross)
          end do
          close(110) 
          
          print *,'    Saving nods file of cross...'
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'nods'//'_'//ADJUSTL(temp)       
          open(103,file=c_File_name_1,status='unknown')         
          do i=1,Num_Node
              write(103, '(200I10)') (Node_Cross_elem(i,
     &                                    j),j=1,num_Cross)
          end do
          close(103)   
      end select
      RETURN
      END SUBROUTINE Save_Files_Cross
