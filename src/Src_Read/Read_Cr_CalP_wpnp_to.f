 
      SUBROUTINE Read_Cr_CalP_Wpnp_to(Result_File_Num,Variable,m,n)
      
      use Global_Float_Type      
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Filename
      
      implicit none
      
      integer,intent(in)::m,n,Result_File_Num
      real(kind=FT),intent(out)::Variable(m,n)
      
      integer i,j
      character(200) c_File_name_1
      character(5) temp
                     
      print *,'    Reading wpnp of crack points......'
      write(temp,'(I5)') Result_File_Num
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'wpnp'//'_'//ADJUSTL(temp)            
      
     
      select case(Key_Data_Format)
      case(1)
          open(101,file=c_File_name_1,status='unknown')         
          do i=1,num_Crack
              read(101,'(2000E20.12)') (Variable(i,
     &                                   j),j=1,Cracks_CalP_Num(i))
          end do
          close(101) 
      case(2)
          open(101,file=c_File_name_1,status='unknown',
     &                                form='unformatted')         
          do i=1,num_Crack
              read(101) (Variable(i,j),j=1,Cracks_CalP_Num(i))
          end do
          close(101) 
      end select
      
      RETURN
      END SUBROUTINE Read_Cr_CalP_Wpnp_to
