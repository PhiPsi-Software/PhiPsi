 
      SUBROUTINE Read_Disp_to(iter,Variable,n)
      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      
      implicit none
      
      integer,intent(in)::iter,n
      real(kind=FT),intent(out)::Variable(n)
      
      integer i
      character(200) c_File_name_1
      character(5) temp
      
      write(temp,'(I5)') iter      
                  
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'disp'//'_'//ADJUSTL(temp)      
     
      select case(Key_Data_Format)
      case(1)
          open(101,file=c_File_name_1,status='unknown')         
          read(101,'(E20.12)') (Variable(i),i=1,n)
          close(101) 
      case(2)
          open(101,file=c_File_name_1,status='unknown',
     &                                form='unformatted')         
          read(101) (Variable(i),i=1,n)
          close(101) 
      end select
      
      RETURN
      END SUBROUTINE Read_Disp_to
