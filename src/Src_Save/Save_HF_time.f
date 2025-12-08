 
      SUBROUTINE Save_HF_time(imf,ifra,Counter_Iter,c_Time)

      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_DISP
      use Global_POST
      
      implicit none
      
      integer,intent(in)::imf,ifra,Counter_Iter
      real(kind=FT),intent(in)::c_Time
      character(200) c_File_name_1
      logical alive
      
      if (Key_Save_Nothing==1) return 
      
      c_File_name_1   =  trim(Full_Pathname)//'.hftm'   
      
      if(imf==1. and. ifra==1 .and. Counter_Iter==1)then
          inquire(file=c_File_name_1, exist=alive) 
          if(alive.EQV..True.)then
              OPEN  (UNIT=105, FILE=c_File_name_1, STATUS='OLD') 
              CLOSE (UNIT=105, STATUS='DELETE')
          endif
      endif

      open(301,file=c_File_name_1,status='unknown',
     &         position='append',action='write') 
      if(imf==1.and.ifra==1.and.Counter_Iter==1)then
          write(301,*) '    imf   |   ifra   | total_ter|   time   '   
      endif   
      write(301, '(3I10,1F18.5)') imf,ifra,Counter_Iter,c_Time
      close(301)    
 
      RETURN
      END SUBROUTINE Save_HF_time
