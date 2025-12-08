 
      SUBROUTINE Save_Ifrac_LS_CPU_time(ifra_CPU_time,Num_Frac)

      use Global_Float_Type      
      use Global_Common
      use Global_Filename
      use Global_POST
      
      implicit none
      
      real(kind=FT),intent(in):: ifra_CPU_time(Num_Frac)
      integer,intent(in):: Num_Frac
      integer i
      character(200) c_File_name_1
      
      if (Key_Save_Nothing==1) return 
      
      print *,'    Saving LS CPU time (s) of each time step...'
      c_File_name_1   =  trim(Full_Pathname)//'.icpt_LS'   
      open(201,file=c_File_name_1,status='unknown')   
      write(201,*) '    LS CPU time (s) of each fracturing step:' 
      do i=1,Num_Frac
          write(201, '(I10,E15.5)') i,ifra_CPU_time(i)
      end do
      close(201)    
 
      RETURN
      END SUBROUTINE Save_Ifrac_LS_CPU_time
