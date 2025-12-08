 
      SUBROUTINE Save_Coors_of_Fracture_Zone
      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_POST
      
      implicit none
      
      character(200) c_File_name_1
      c_File_name_1   =  trim(Full_Pathname)//'.fraz'  
      
      
      if (Key_Save_Nothing==1) return
      
      open(201,file=c_File_name_1,status='unknown')   
      
      if(Key_Dimension==2) then
        write(201,*) 'min_X    |    max_X    |    min_Y    |    max_Y'
        write(201, '(4F12.5)')   Frac_Zone_MinX,Frac_Zone_MaxX,
     &                           Frac_Zone_MinY,Frac_Zone_MaxY  
      elseif(Key_Dimension==3)then
        write(201,*) 'min_X    |    max_X    '
     &        //'|    min_Y    |    max_Y    |    min_Z    |    max_Z'
        write(201, '(6F12.5)')   Frac_Zone_MinX,Frac_Zone_MaxX,
     &                           Frac_Zone_MinY,Frac_Zone_MaxY,
     &                           Frac_Zone_MinZ,Frac_Zone_MaxZ        
      endif
      close(201)   
 
      RETURN
      END SUBROUTINE Save_Coors_of_Fracture_Zone
