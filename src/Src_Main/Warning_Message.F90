!-----------------------------------------------------------
! Brief: Display a formatted warning or error message and pause.
!
! Parameters:
!   Input:  Mess_Type - message category ('E' error, 'W' warning, 'S' stop)
!           Keywords   - keyword string to include in the message
!
! Notes:   Optionally plays an audio cue and waits for user acknowledgement
!          before continuing, depending on the message type.
!-----------------------------------------------------------

SUBROUTINE Warning_Message(Mess_Type,Keywords)
! Display error message.
use Global_Float_Type
use Global_Common
use Module_Tool_Logger

implicit none

character(len=1)  Mess_Type
character(len=*) ::  Keywords
character(len=100) Mess_Content,temp_log

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
