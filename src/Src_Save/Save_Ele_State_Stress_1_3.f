 
      SUBROUTINE Save_Ele_State_Stress_1_3(isub)

      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_Stress
      use Global_POST
      
      implicit none
      
      integer,intent(in)::isub
      integer i,c_State
      character(200) c_File_name_1
      character(5) temp   
      
      if (Key_Save_Nothing==1) return
      
      print *,'    Saving stress state of elements...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'elss'//'_'//ADJUSTL(temp)   

      select case(Key_Data_Format)
      case(1:2)  
          open(201,file=c_File_name_1,status='unknown')     
          do i=1,num_Elem
              write(201, '(I2)') Ele_State_Stress_1_3(i)
          end do      
          close(201)      
      end select
 
      RETURN
      END SUBROUTINE Save_Ele_State_Stress_1_3
