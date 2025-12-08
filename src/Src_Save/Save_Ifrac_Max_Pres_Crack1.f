 
      SUBROUTINE Save_Ifrac_Max_Pres_Crack1(ifra_Max_Press,Num_Frac)

      use Global_Float_Type      
      use Global_Common
      use Global_Filename
      use Global_POST
      
      implicit none
      
      real(kind=FT),intent(in):: ifra_Max_Press(Num_Frac)
      integer,intent(in):: Num_Frac
      integer i
      character(200) c_File_name_1
      
      if (Key_Save_Nothing==1) return 
      
      print *,'    Saving max pressure of crack 1 of each time step...'
      c_File_name_1   =  trim(Full_Pathname)//'.ipre'   
      open(201,file=c_File_name_1,status='unknown') 
      write(201,*) '    Max pressure (MPa )of crack 1'
     &                  //' of each time step' 
      do i=1,Num_Frac
          if(ifra_Max_Press(i)>ZR)then
              write(201, '(F12.5)') ifra_Max_Press(i)
          endif
      end do
      close(201)   
 
      RETURN
      END SUBROUTINE Save_Ifrac_Max_Pres_Crack1
