 
      SUBROUTINE Save_Gauss_Stress(isub,Total_Num_G_P)
      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_DISP
      use Global_Stress
      
      implicit none
      
      integer,intent(in)::isub,Total_Num_G_P
      integer i
      character(200) c_File_name_1
      character(5) temp
      real(kind=FT) tem_2  
      
      tem_2 = 1.0D6
      
      print *,'    Saving stresses of Gauss points...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'strg'//'_'//ADJUSTL(temp)   

      select case(Key_Data_Format)
      case(1)
          open(201,file=c_File_name_1,status='unknown') 
          if(Key_Unit_System==1)then
              do i=1,Total_Num_G_P
                  write(201, '(I8,4E20.12)') i,Stress_xx_Gauss(i),
     &                                         Stress_yy_Gauss(i),
     &                                         Stress_xy_Gauss(i),
     &                                         Stress_vm_Gauss(i)
              end do   
          elseif(Key_Unit_System==2)then
              do i=1,Total_Num_G_P
                  write(201, '(I8,4E20.12)') i,Stress_xx_Gauss(i)*tem_2,
     &                                         Stress_yy_Gauss(i)*tem_2,
     &                                         Stress_xy_Gauss(i)*tem_2,
     &                                         Stress_vm_Gauss(i)*tem_2
              end do   
          endif
   
          close(201)  
      case(2)
          open(203,file=c_File_name_1,status='unknown',
     &                            form='unformatted',access='stream')     
          write(203) ( Stress_xx_Gauss(i),
     &                 Stress_yy_Gauss(i),
     &                 Stress_xy_Gauss(i),
     &                 Stress_vm_Gauss(i),i=1,Total_Num_G_P)
          close(203)      
      end select
     
      
      RETURN
      END SUBROUTINE Save_Gauss_Stress
