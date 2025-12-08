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
 
SUBROUTINE Warning_Message(Mess_Type,Keywords)
! Display error message.
use Global_Float_Type
use Global_Common
use Module_Tool_Logger

implicit none

character*1  Mess_Type
character(len=*) ::  Keywords
character*100 Mess_Content,temp_log
  
1001 FORMAT('     ?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?')    
!10100 format (1x, a1)     
  select case(Mess_Type)
      
  case('E')
      Mess_Content =  '       '//'ERROR :: Keyword *'//trim(Keywords)//' Wrong!'
      WRITE(*,1001)  
      print *,trim(Mess_Content)
      WRITE(*,1001)   
      if(Key_Play_Sounds==1)then
          !**********
          ! option 1
          !**********
          call system('cd .\Python_Tools\ && PhiPsi_Play_Error.py')

          !**********
          ! option 2
          !**********
          !beep = char(7) 
          !write (*, 10100) beep                     
      endif
      !write( *, * ) '    Close the Program and Check Input Files or Press Enter to Continue.' 
      write( *, * ) '    Close the Program and Check Input Files or Press Enter to Continue.'
      read( *, * )    
      

  case('W')
      Mess_Content =  '       '//'Warning :: Keyword *'//trim(Keywords)//' Wrong!'      
      WRITE(*,1001)  
      print *,trim(Mess_Content)
      WRITE(*,1001)  
      if(Key_Play_Sounds==1)then
       !CALL SYSTEM('python d:/PhiPsi_Python/PhiPsi_Play_Sounds.py')
       !CALL SYSTEM('python d:/PhiPsi_Python/PhiPsi_Play_Sounds.py')
      endif
  case('S')
      Mess_Content =  '                    '//'An error occurred.'  
      WRITE(*,1001)  
      print *,trim(Mess_Content)
      WRITE(*,1001)   
      if(Key_Play_Sounds==1)then
          !**********
          ! option 1
          !**********
          call system('cd .\Python_Tools\ && PhiPsi_Play_Error.py')
          !**********
          ! option 2
          !**********
          !beep = char(7) 
          !write (*, 10100) beep               
      endif
      !write( *, * ) '    Close the Program and Check Input Files or Press Enter to Continue.' 
      write( *, * ) '    Close the Program and Check Input Files or Press Enter to Continue.'
      read( *, * )     
  end select
  
  temp_log = trim('ERROR OCCUR!') 
  call log_msg(temp_log)
  
  RETURN
  END SUBROUTINE Warning_Message
