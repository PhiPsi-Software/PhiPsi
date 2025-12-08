 
      SUBROUTINE Save_Ifrac_Num_Iter(ifra_num_iter,Num_Frac)

      use Global_Float_Type      
      use Global_Common
      use Global_Filename
      use Global_POST
      
      implicit none
      
      integer,intent(in):: ifra_num_iter(Num_Frac)
      integer,intent(in):: Num_Frac
      integer i
      character(200) c_File_name_1
      
      if (Key_Save_Nothing==1) return
      
      print *,'    Saving num of iteration of each time step...'
      c_File_name_1   =  trim(Full_Pathname)//'.iite'   
      open(201,file=c_File_name_1,status='unknown')  
      write(201,*) '    Num of iteration of each time step' 
      do i=1,Num_Frac
          if(ifra_num_iter(i)>0)then
              write(201, '(I5)') ifra_num_iter(i)
          endif
      end do
      close(201)
 
      RETURN
      END SUBROUTINE Save_Ifrac_Num_Iter
