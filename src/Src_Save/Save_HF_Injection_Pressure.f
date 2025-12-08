 
      SUBROUTINE Save_HF_Injection_Pressure(ifra,c_Time)

      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      use Global_POST
      
      implicit none
      
      integer,intent(in)::ifra
      real(kind=FT),intent(in)::c_Time
      character(200) c_File_name_1
      logical alive
      real(kind=FT) c_Inj_Press
      integer i_C
      
      if (Key_Save_Nothing==1) return 
      
      print *,'    Saving injection pressure...'
      c_File_name_1   =  trim(Full_Pathname)//'.injp'   
      
      if(ifra==1)then
          inquire(file=c_File_name_1, exist=alive)  
          if(alive.EQV..True.)then
              OPEN  (UNIT=105, FILE=c_File_name_1, STATUS='OLD') 
              CLOSE (UNIT=105, STATUS='DELETE')
          endif
      endif
      

      c_Inj_Press = ZR
      do i_C=1,num_Crack
          if (i_C == Inject_Crack_Num) then                             
              c_Inj_Press = Cracks_CalP_Pres(i_C,CalP_num_InjP_Local)   
              Last_Inj_Pres = c_Inj_Press 
              c_Inj_Press = c_Inj_Press +
     &                    Cracks_CalP_Remo_Strs(i_C,CalP_num_InjP_Local) 
              exit
          endif
      enddo

      
      open(301,file=c_File_name_1,status='unknown',
     &             position='append',action='write') 
      if(ifra==1)then
          write(301,*) '    ifra   |   time   | injection pressure'   
      endif   
      write(301, '(I10,2E20.12)') ifra,c_Time,c_Inj_Press
      close(301)   
 
      RETURN
      END SUBROUTINE Save_HF_Injection_Pressure
