 
      SUBROUTINE Save_Crack_Tan_Relative_Disp(isub,Cracks_CalP_Tan_Aper)

      use Global_Float_Type      
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Filename
      use Global_POST
      
      implicit none
      real(kind=FT),intent(in)::Cracks_CalP_Tan_Aper(num_Crack,
     &                                                Max_Num_Cr_CalP)
      integer isub
      integer i,j
      character(200) c_File_name_1
      character(5) temp  
      real(kind=FT) tem_1,tem_2 
      
      
      if (Key_Save_Nothing==1) return
      
      tem_1 = 1000.0D0
      tem_2 = 1.0D6

      print *,'    Saving tangent apertures of calculation points...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'ctap'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown') 
      if(Key_Unit_System==1)then   
          do i=1,num_Crack
              write(101,'(2000E20.12)') (Cracks_CalP_Tan_Aper(i,
     &                                   j),j=1,Cracks_CalP_Num(i))
          end do
      elseif(Key_Unit_System==2)then
          do i=1,num_Crack
              write(101,'(2000E20.12)') (Cracks_CalP_Tan_Aper(i,
     &                             j)/tem_1,j=1,Cracks_CalP_Num(i))
          end do
      endif
      close(101) 
 
      RETURN
      END SUBROUTINE Save_Crack_Tan_Relative_Disp
