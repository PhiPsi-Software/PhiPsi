 
      SUBROUTINE Save_HF_wbpt_File(i_WB,i_Stage,i_Prop,c_Pres,c_Time)

      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_DISP
      use Global_POST
      
      implicit none
      
      integer,intent(in)::i_WB,i_Stage,i_Prop
      real(kind=FT),intent(in)::c_Time,c_Pres
      character(200) c_File_name_1
      logical alive
      
      if (Key_Save_Nothing==1) return 
      
      c_File_name_1   =  trim(Full_Pathname)//'.wbpt'   
      

      if(i_WB ==1 .and. i_Stage==1 .and. i_Prop==1)then
          inquire(file=c_File_name_1, exist=alive)  
          if(alive.EQV..True.)then
              OPEN  (UNIT=105, FILE=c_File_name_1, STATUS='OLD') 
              CLOSE (UNIT=105, STATUS='DELETE')
          endif
      endif

      open(301,file=c_File_name_1,status='unknown',
     &         position='append',action='write') 
      if(i_WB ==1 .and. i_Stage==1 .and. i_Prop==1)then
          write(301,*) '    i_WB   |   i_Stage   | i_Prop   ' 
     &       //'|    Time (s)   |     Pressure (Pa)'   
      endif   
      write(301, '(3I10,F18.5,5X,F18.5)') i_WB,i_Stage,i_Prop,
     &                                    c_Time,c_Pres
      close(301)    
 
      RETURN
      END SUBROUTINE Save_HF_wbpt_File
