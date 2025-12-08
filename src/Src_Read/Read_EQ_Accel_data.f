 
      SUBROUTINE Read_EQ_Accel_data
      
      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_Dynamic 
      
      implicit none
      
      LOGICAL alive
      character*200 temp_name
      integer Tool_Count_Lines
      real(kind=FT),ALLOCATABLE::Temp_DATA(:,:)
      logical Flag_Blank
      
      print *, "    Trying to read eqac files...." 
      
      temp_name = trim(trim(Full_Pathname)//'.eqac')
      
      inquire(file=temp_name, exist=alive)  
      
      if(alive.EQV..FALSE.)then
          print *, "    ERROR :: Can not find eqac file," 
          print *, "             which is necessary when *Key_EQ=1!" 
          call Warning_Message('S',Keywords_Blank)
      else
          num_EQ_Accel = Tool_Count_Lines(temp_name) 
          ALLOCATE( EQ_Accel_data(num_EQ_Accel))
          ALLOCATE( Temp_DATA(num_EQ_Accel,1))
          Call Tool_Read_File(temp_name,"eqac",num_EQ_Accel,1,Temp_DATA,
     &                        Flag_Blank)
          EQ_Accel_data  = Temp_DATA(:,1)
          DEALLOCATE(Temp_DATA)
      endif  
      
      RETURN
      END SUBROUTINE Read_EQ_Accel_data
