 
      SUBROUTINE Save_Gauss_Num_of_Elements(isub)

      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_DISP
      use Global_POST
      
      implicit none
      
      integer,intent(in)::isub

      integer i
      character(200) c_File_name_1
      character(5) temp   

      if (Key_Save_Nothing==1) return 
      
      if(Key_Simple_Post==1) return
      
      print *,'    Saving number of Gauss points of elements...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'elgn'//'_'//ADJUSTL(temp)   
      select case(Key_Data_Format)
      case(1)  
          open(201,file=c_File_name_1,status='unknown')  
            do i=1,Num_Elem
              write(201, '(I8)') Elements_Gauss_Num(i)
            end do  
          close(201)  
      case(2) 
          open(201,file=c_File_name_1,status='unknown',
     &             form='unformatted',access='stream')      
          do i=1,Num_Elem
              write(201) Elements_Gauss_Num(i)
          end do
          close(201)   
      end select
 
      RETURN
      END SUBROUTINE Save_Gauss_Num_of_Elements
