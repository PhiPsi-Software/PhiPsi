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
 
      SUBROUTINE Save_Killed_Elements(isub)
      use Global_Float_Type
      use Global_Common
      use Global_Filename
      use Global_Model
      use Global_POST
      
      implicit none
      
      integer isub
      integer i,j,i_E
      character(200) c_File_name_1
      character(5) temp   
      integer num_Ele_Killed
   
      if (Key_Save_Nothing==1) return
      
      write(temp,'(I5)') isub
      print *,'    Saving killed elements...'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'kiel'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown')         
      do i=1,isub
          num_Ele_Killed = count(Ele_Killed_Each_Load_Step(i,:)>0)
          if(num_Ele_Killed>=1)then
              do j=1,num_Ele_Killed
                write(101, '(I10)') Ele_Killed_Each_Load_Step(i,j)
              enddo
          endif        
      end do
      do i_E=1,Num_Elem
        if(Elem_Break(i_E))then
            write(101, '(I10)') i_E
        endif
      enddo        
      close(101)       
      

     

      RETURN
      END SUBROUTINE Save_Killed_Elements
