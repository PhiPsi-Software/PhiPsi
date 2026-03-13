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
 
      subroutine Tool_Get_Last_Result_Number(Result_File_Num)
C     Obtain the file number of the last save step of the hydraulic fracturing analysis
      
      use Global_Float_Type
      use Global_Filename
      use Global_Common   
      
      implicit none
      integer,intent(out)::Result_File_Num
      
      integer i_Check
      character(200) c_File_name_1,c_File_name_2,c_File_name_3
      logical exist_1,exist_2,exist_3
      character(5) temp
      
      Result_File_Num = 0
      do i_Check=1,2000
          write(temp,'(I5)') i_Check
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'crax'//'_'//ADJUSTL(temp)    
          c_File_name_2   =  trim(Full_Pathname)//'.'
     &                 //'cray'//'_'//ADJUSTL(temp)    
          c_File_name_3   =  trim(Full_Pathname)//'.'
     &                 //'wpnp'//'_'//ADJUSTL(temp)    
         inquire(file=c_File_name_1, exist= exist_1) 
         inquire(file=c_File_name_2, exist= exist_2) 
         inquire(file=c_File_name_3, exist= exist_3) 
         if(exist_1 .and. exist_2 .and. exist_3)then
             Result_File_Num = i_Check
         endif
      enddo
      if(Result_File_Num ==0)then
          print *, '    Error :: cannot find *wpnp file, which'
     &            // ' is required when Key_Propped_Width=1 or' 
     &            // ' Key_Proppant_Active=1!'
          print *, '             Make sure that HF analysis has already'
     &            // ' been performed!'
          if(Key_Analysis_Type/=17)then
              call Warning_Message('S',Keywords_Blank)
          endif
      endif
      
      RETURN
      END SUBROUTINE Tool_Get_Last_Result_Number