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
 
      SUBROUTINE Save_Ele_Contac_State(isub)
      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_Contact
      use Global_POST
      
      implicit none
      
      integer,intent(in)::isub
      integer i,c_State
      character(200) c_File_name_1
      character(5) temp   
      
      if (Key_Save_Nothing==1) return
      
      print *,'    Saving contact state of elements...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'elcs'//'_'//ADJUSTL(temp)   

      select case(Key_Data_Format)
      case(1:2)  
          open(201,file=c_File_name_1,status='unknown')     
          do i=1,num_Elem
              if (maxval(Elem_Conta_Sta(i,:)) == 1) then       
                  c_State = 1
              elseif (maxval(Elem_Conta_Sta(i,:)) == 2) then   
                  c_State = 2
              elseif (maxval(Elem_Conta_Sta(i,:)) == 0) then  
                  c_State = 0
              end if
              write(201, '(I2)') c_State
          end do      
          close(201)       
      end select
 
      RETURN
      END SUBROUTINE Save_Ele_Contac_State
