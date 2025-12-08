 
      SUBROUTINE Read_Cr_CalP_Aper_to(iter,Variable,m,n)
      use Global_Float_Type      
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Filename
      
      implicit none
      
      integer,intent(in)::iter,m,n
      real(kind=FT),intent(out)::Variable(m,n)
      
      integer i,j
      character(200) c_File_name_1
      character(5) temp
      
      
      Variable(1:m,1:n) = ZR
      
      write(temp,'(I5)') iter      
                     
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'cape'//'_'//ADJUSTL(temp)      
     
     
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
      END SUBROUTINE Read_Cr_CalP_Aper_to
