 
      subroutine Tool_Get_Last_Result_Number(Result_File_Num)
      
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