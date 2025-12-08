 
      SUBROUTINE Save_Dof_Force_3D(isub,F,all_Solid_Freedom)
      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_Stress
      use Global_POST
      
      implicit none
      
      integer,intent(in)::isub,all_Solid_Freedom
      real(kind=FT),intent(in)::F(all_Solid_Freedom)
      
      real(kind=FT) Force_x(all_Solid_Freedom/3)
      real(kind=FT) Force_y(all_Solid_Freedom/3)
      real(kind=FT) Force_z(all_Solid_Freedom/3)
      integer i
      character(200) c_File_name_1,c_File_name_2,c_File_name_3
      character(5) temp   
      
      if (Key_Save_Nothing==1) return
      
      if(Key_Simple_Post==1) return
      
      print *,'    Saving force of all solid dofs...'
      
      do i=1,all_Solid_Freedom/3
          Force_x(i)=F(i*3-2)
          Force_y(i)=F(i*3-1)
          Force_z(i)=F(i*3)
      end do      

      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'fxdf'//'_'//ADJUSTL(temp)        
      c_File_name_2   =  trim(Full_Pathname)//'.'
     &                 //'fydf'//'_'//ADJUSTL(temp)   
      c_File_name_3   =  trim(Full_Pathname)//'.'
     &                 //'fzdf'//'_'//ADJUSTL(temp)  
     
      
      select case(Key_Data_Format)
      case(1)  
          open(202,file=c_File_name_1,status='unknown')
          open(203,file=c_File_name_2,status='unknown')
          open(204,file=c_File_name_3,status='unknown')
          do i=1,all_Solid_Freedom/3
              write(202, '(E20.12)') Force_x(i)
              write(203, '(E20.12)') Force_y(i)
              write(204, '(E20.12)') Force_z(i)
          end do       
          close(202)  
          close(203)  
          close(204)  
      case(2)  
          open(202,file=c_File_name_1,status='unknown',    
     &             form='unformatted',access='stream')
          open(203,file=c_File_name_2,status='unknown',    
     &             form='unformatted',access='stream')
          open(204,file=c_File_name_3,status='unknown',    
     &             form='unformatted',access='stream')
          do i=1,all_Solid_Freedom/3
              write(202) Force_x(i)
              write(203) Force_y(i)
              write(204) Force_z(i)
          end do       
          close(202)  
          close(203)  
          close(204)    
      end select

      RETURN
      END SUBROUTINE Save_Dof_Force_3D
