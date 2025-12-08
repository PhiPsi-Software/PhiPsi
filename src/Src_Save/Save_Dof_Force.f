 
      SUBROUTINE Save_Dof_Force(isub,F,all_Solid_Freedom)
      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_Stress
      use Global_POST
      
      implicit none
      
      integer,intent(in)::isub,all_Solid_Freedom
      real(kind=FT),intent(in)::F(all_Solid_Freedom)
      
      real(kind=FT) Force_x(all_Solid_Freedom/2)
      real(kind=FT) Force_y(all_Solid_Freedom/2)
      integer i
      character(200) c_File_name_1,c_File_name_2
      character(5) temp  
      
      if (Key_Save_Nothing==1) return
      
      print *,'    Saving force of all solid dofs...'
      

      do i=1,all_Solid_Freedom/2
          Force_x(i)=F(i*2-1)
          Force_y(i)=F(i*2)
      end do      


      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'fxdf'//'_'//ADJUSTL(temp)        
      select case(Key_Data_Format)
      case(1:2)  
          open(202,file=c_File_name_1,status='unknown')
          if (Key_Dimension == 2) then      
              do i=1,all_Solid_Freedom/2
                  write(202, '(E20.12)') Force_x(i)
              end do       
          else if(Key_Dimension == 3)then      
          end if          
          close(202)  
    
      end select

      write(temp,'(I5)') isub
      c_File_name_2   =  trim(Full_Pathname)//'.'
     &                 //'fydf'//'_'//ADJUSTL(temp)        
      
      select case(Key_Data_Format)
      case(1:2) 
          open(302,file=c_File_name_2,status='unknown')
          if (Key_Dimension == 2) then      
              do i=1,all_Solid_Freedom/2
                  write(302, '(E20.12)') Force_y(i)
              end do       
          else if(Key_Dimension == 3)then      
          end if          
          close(302)      
      end select 
      
      RETURN
      END SUBROUTINE Save_Dof_Force
