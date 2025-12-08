 
  SUBROUTINE Save_Surface_Load_3D(isub,i_SL,F,all_Solid_Freedom)
  use Global_Float_Type      
  use Global_Common
  use Global_Model
  use Global_Filename
  use Global_Stress
  use Global_POST
  
  implicit none
  
  integer,intent(in)::isub,all_Solid_Freedom,i_SL
  real(kind=FT),intent(in)::F(all_Solid_Freedom)
  
  real(kind=FT) Force_x(all_Solid_Freedom/3)
  real(kind=FT) Force_y(all_Solid_Freedom/3)
  real(kind=FT) Force_z(all_Solid_Freedom/3)
  integer i
  character(200) c_File_name_1,c_File_name_2,c_File_name_3
  character(5) temp,temp2  
  
  if (Key_Save_Nothing==1) return 
  
  if(Key_Simple_Post==1)then
      return
  endif
  
  print *,'    Saving nodal force of pressure load ',i_SL
  

  do i=1,all_Solid_Freedom/3
      Force_x(i)=F(i*3-2)
      Force_y(i)=F(i*3-1)
      Force_z(i)=F(i*3)
  end do      

  write(temp,'(I5)') isub
  write(temp2,'(I5)') i_SL
  c_File_name_1   =  trim(Full_Pathname)//'.fxsl'//'_'//trim(ADJUSTL(temp2))//'_'//ADJUSTL(temp)          
  c_File_name_2   =  trim(Full_Pathname)//'.fysl'//'_'//trim(ADJUSTL(temp2))//'_'//ADJUSTL(temp)     
  c_File_name_3   =  trim(Full_Pathname)//'.fzsl'//'_'//trim(ADJUSTL(temp2))//'_'//ADJUSTL(temp)          
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

  RETURN
  END SUBROUTINE Save_Surface_Load_3D
