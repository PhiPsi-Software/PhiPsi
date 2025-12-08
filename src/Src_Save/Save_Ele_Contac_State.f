 
      SUBROUTINE Save_Ele_Contac_State(isub)
      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_Contact
      use Global_POST
      
      implicit none
      
      integer,intent(in)::isub
      integer i,c_State
      character(200) c_File_name_1
      character(5) temp   
      
      if (Key_Save_Nothing==1) return
      
      print *,'    Saving contact state of elements...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'elcs'//'_'//ADJUSTL(temp)   

      select case(Key_Data_Format)
      case(1:2)  
          open(201,file=c_File_name_1,status='unknown')     
          do i=1,num_Elem
              if (maxval(Elem_Conta_Sta(i,:)) == 1) then       
                  c_State = 1
              elseif (maxval(Elem_Conta_Sta(i,:)) == 2) then   
                  c_State = 2
              elseif (maxval(Elem_Conta_Sta(i,:)) == 0) then  
                  c_State = 0
              end if
              write(201, '(I2)') c_State
          end do      
          close(201)       
      end select
 
      RETURN
      END SUBROUTINE Save_Ele_Contac_State
