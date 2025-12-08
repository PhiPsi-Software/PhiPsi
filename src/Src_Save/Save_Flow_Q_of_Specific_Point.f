 
      SUBROUTINE Save_Flow_Q_of_Specific_Point(ifra,total_time,Flow_Q,
     &                                         Saved_Filename)

      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_DISP
      use Global_POST
      
      implicit none
      integer,intent(in)::ifra
      real(kind=FT),intent(in)::total_time,Flow_Q
      character(len=*),intent(in) :: Saved_Filename
      character(200) c_File_name_1
      logical alive
      
      if (Key_Save_Nothing==1) return 
      
      c_File_name_1   =  trim(Full_Pathname)//'.'//trim(Saved_Filename)   

      if(ifra==1)then
          inquire(file=c_File_name_1, exist=alive)  
          if(alive.EQV..True.)then
              OPEN  (UNIT=105, FILE=c_File_name_1, STATUS='OLD') 
              CLOSE (UNIT=105, STATUS='DELETE')
          endif
      endif

      open(301,file=c_File_name_1,status='unknown',
     &         position='append',action='write') 
      if(ifra==1)then
          write(301,*) '  ifra   |  time   | flow quantity'   
      endif   
      write(301, '(I10,2F15.8)') ifra,total_time,Flow_Q
      close(301)   
 
      RETURN
      END SUBROUTINE Save_Flow_Q_of_Specific_Point
