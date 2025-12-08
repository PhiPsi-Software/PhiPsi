 
      SUBROUTINE Save_Files_Cr_CalP_3D(isub)

      use Global_Float_Type      
      use Global_Common
      use Global_Crack_Common
      use Global_Crack_3D
      use Global_Model
      use Global_Filename
      use Global_POST
      
      implicit none
      
      integer isub
      integer i,j
      character(200) c_File_name_1
      character(5) temp  
      real(kind=FT) tem_1 ,tem_2
      
      if (Key_Save_Nothing==1) return
      tem_1 = 1000.0D0
      tem_2 = 1.0D6
      

      print *,'    Saving apertures of calculation points...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'cape'//'_'//ADJUSTL(temp)    
      open(101,file=c_File_name_1,status='unknown') 
      do i=1,num_Crack
          write(101,'(50000E20.12)') (Cracks_CalP_Aper_3D(i)%row(j),
     &                         j=1,Cracks_CalP_Num_3D(i))
      end do
      close(101) 


      print *,'    Saving pressure of calculation points...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'cpre'//'_'//ADJUSTL(temp)    
      open(101,file=c_File_name_1,status='unknown') 
      do i=1,num_Crack
          write(101,'(50000E20.12)') (Cracks_CalP_Pres_3D(i)%row(j),
     &                         j=1,Cracks_CalP_Num_3D(i)) 
      end do
      close(101) 
      

      print *,'    Saving apertures of fluid element...'
      write(temp,'(I5)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'capf'//'_'//ADJUSTL(temp)    
      open(102,file=c_File_name_1,status='unknown') 
      do i=1,num_Crack
          write(102,'(50000E20.12)') 
     &                  (Cracks_FluidEle_Aper_3D(i)%row(j),
     &                               j=1,Cracks_FluidEle_num_3D(i))      
      end do
      close(102) 

      
      RETURN
      END SUBROUTINE Save_Files_Cr_CalP_3D
