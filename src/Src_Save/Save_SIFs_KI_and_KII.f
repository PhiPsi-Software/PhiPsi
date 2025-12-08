 
      SUBROUTINE Save_SIFs_KI_and_KII(isub)
      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Filename
      use Global_POST
      
      implicit none
      
      integer isub
      integer i,j
      character(200) c_File_name_1
      character(5) temp  
      real(kind=FT) tem_1
      
      if (Key_Save_Nothing==1) return 
      
      tem_1 = 1.0D-6*sqrt(1000.0D0)
      print *,'    Saving SIFs...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'sifs'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown')    
          
      if(Key_Unit_System==1)then    
          do i=1,num_Crack
              write(101, '(4E20.12)')
     &             KI(i,1),KII(i,1),KI(i,2),KII(i,2)
          end do
      elseif(Key_Unit_System==2)then
          
          do i=1,num_Crack
              write(101, '(4E20.12)')
     &             KI(i,1)/tem_1,KII(i,1)/tem_1,
     &             KI(i,2)/tem_1,KII(i,2)/tem_1
          end do
      endif
          
      close(101)   
      RETURN
      END SUBROUTINE Save_SIFs_KI_and_KII
