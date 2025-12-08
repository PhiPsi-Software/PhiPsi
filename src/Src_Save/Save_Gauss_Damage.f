 
      SUBROUTINE Save_Gauss_Damage(isub,Total_Num_G_P)
      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_DISP
      use Global_Stress
      use Global_POST
      
      implicit none
      integer,intent(in)::isub,Total_Num_G_P
      integer i
      character(200) c_File_name_1
      character(5) temp   
      real(kind=FT) tem_2  
      
      if (Key_Save_Nothing==1) return 
      
      tem_2 = 1.0D6
      
      print *,'    Saving damage factor of Gauss points...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'damg'//'_'//ADJUSTL(temp)   

      select case(Key_Data_Format)
      case(1)  
          open(201,file=c_File_name_1,status='unknown') 
          do i=1,Total_Num_G_P
              write(201, '(I8,E20.12)') i,Damage_Gauss(i)
          end do   
          close(201)  
      case(2) 
          open(203,file=c_File_name_1,status='unknown',
     &                            form='unformatted',access='stream')     
          write(203) (Damage_Gauss(i),i=1,Total_Num_G_P)  
          close(203)      
      end select
      
      
      RETURN
      END SUBROUTINE Save_Gauss_Damage
