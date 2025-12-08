 
      SUBROUTINE Save_HF_Natural_Cracks_Coors(c_num_Na_Cr,
     &                                        c_Na_Crack_Coor)

      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Model
      use Global_Filename
      use Global_POST
      
      implicit none
      integer,intent(in)::c_num_Na_Cr              
      real(kind=FT),intent(in)::c_Na_Crack_Coor(c_num_Na_Cr,2,2)      
      integer i
      character(200) c_File_name_1,c_File_name_2
      real(kind=FT) tem_1 
      
      tem_1 = 1000.0D0
      
      if (Key_Save_Nothing==1) return 
      
      print *,'    Saving coordinates of natural cracks...'
      c_File_name_1   =  trim(Full_Pathname)//'.ncrx'
      c_File_name_2   =  trim(Full_Pathname)//'.ncry'      
      open(101,file=c_File_name_1,status='unknown')    
      open(102,file=c_File_name_2,status='unknown')  
      if(Key_Unit_System==1)then     
          do i=1,c_num_Na_Cr
            write(101,'(E20.12,A,E20.12)')
     &            c_Na_Crack_Coor(i,1,1),',',c_Na_Crack_Coor(i,2,1)
            write(102,'(E20.12,A,E20.12)')
     &            c_Na_Crack_Coor(i,1,2),',',c_Na_Crack_Coor(i,2,2)
          end do
      elseif(Key_Unit_System==2)then 
          do i=1,c_num_Na_Cr
            write(101,'(E20.12,A,E20.12)')
     &      c_Na_Crack_Coor(i,1,1)/tem_1,',',
     &      c_Na_Crack_Coor(i,2,1)/tem_1
            write(102,'(E20.12,A,E20.12)')
     &      c_Na_Crack_Coor(i,1,2)/tem_1,',',
     &      c_Na_Crack_Coor(i,2,2)/tem_1
          end do
      endif
      
      close(101)
      close(102)   
      RETURN
      END SUBROUTINE Save_HF_Natural_Cracks_Coors
