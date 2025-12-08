 
      SUBROUTINE Save_Gauss_Coors(isub,Total_Num_G_P,
     &                            c_G_CoorX,c_G_CoorY)
      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_DISP
      use Global_POST
      
      implicit none
      
      integer,intent(in)::isub,Total_Num_G_P
      real(kind=FT),intent(in)::c_G_CoorX(Total_Num_G_P),
     &                             c_G_CoorY(Total_Num_G_P)
      integer i
      character(200) c_File_name_1
      character(5) temp  
      real(kind=FT) tem_1
      
      if (Key_Save_Nothing==1) return 
      
      tem_1 = 1000.0D0
      print *,'    Saving coordinates of Gauss points...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'gcor'//'_'//ADJUSTL(temp)   
      select case(Key_Data_Format)
      case(1) 
          open(201,file=c_File_name_1,status='unknown')  
          if(Key_Unit_System==1)then    
            do i=1,Total_Num_G_P
              write(201, '(I8,2E20.12)') i,c_G_CoorX(i),c_G_CoorY(i)
            end do  
          elseif(Key_Unit_System==2)then 
            do i=1,Total_Num_G_P
              write(201, '(I8,2E20.12)') i,c_G_CoorX(i)/tem_1,
     &                                     c_G_CoorY(i)/tem_1
            end do  
          endif
    
          close(201)  
      case(2)  
          open(203,file=c_File_name_1,status='unknown',
     &             form='unformatted',access='stream')     
          write(203) (c_G_CoorX(i),c_G_CoorY(i),i=1,Total_Num_G_P)      
          close(203)      
      end select
 
      RETURN
      END SUBROUTINE Save_Gauss_Coors
