 
      SUBROUTINE Save_Disp_Real_to_POST(isub,POST_Disp)
      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_DISP
      use Global_POST
      
      implicit none
      
      integer,intent(in)::isub
      real(kind=FT),intent(in)::POST_Disp(Total_FD)
      integer i
      character(200) c_File_name_1
      character(5) temp   
      real(kind=FT) tem_1
      
      if (Key_Save_Nothing==1) return 
      
     
      print *,'    Saving real disp to POST...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'disn'//'_'//ADJUSTL(temp)   
      tem_1 = 1000.0D0
      
      select case(Key_Data_Format)
      case(1) 
          open(202,file=c_File_name_1,status='unknown')
          if (Key_Dimension == 2) then   
              if(Key_Unit_System==1)then    
                  do i=1,int(Total_FD/2)
                      write(202, '(I8,2(A,E20.12))') 
     &                       i,',',POST_Disp(2*i-1),',',POST_Disp(2*i)
                  end do
              elseif(Key_Unit_System==2)then
                  do i=1,int(Total_FD/2)
                      write(202, '(I8,2(A,E20.12))') 
     &                     i,',',POST_Disp(2*i-1)/tem_1,','
     &                          ,POST_Disp(2*i)/tem_1
                  end do
              endif 
          endif
          close(202)  
      case(2) 
          open(203,file=c_File_name_1,status='unknown',
     &                            form='unformatted',access='stream')     
          write(203) (POST_Disp(2*i-1),POST_Disp(2*i),
     &                i=1,int(Total_FD/2))
          close(203)    
      end select
 
      RETURN
      END SUBROUTINE Save_Disp_Real_to_POST
