 
      SUBROUTINE Save_Killed_Elements(isub)
      use Global_Float_Type
      use Global_Common
      use Global_Filename
      use Global_Model
      use Global_POST
      
      implicit none
      
      integer isub
      integer i,j,i_E
      character(200) c_File_name_1
      character(5) temp   
      integer num_Ele_Killed
   
      if (Key_Save_Nothing==1) return
      
      write(temp,'(I5)') isub
      print *,'    Saving killed elements...'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'kiel'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown')         
      do i=1,isub
          num_Ele_Killed = count(Ele_Killed_Each_Load_Step(i,:)>0)
          if(num_Ele_Killed>=1)then
              do j=1,num_Ele_Killed
                write(101, '(I10)') Ele_Killed_Each_Load_Step(i,j)
              enddo
          endif        
      end do
      do i_E=1,Num_Elem
        if(Elem_Break(i_E))then
            write(101, '(I10)') i_E
        endif
      enddo        
      close(101)       
      

     

      RETURN
      END SUBROUTINE Save_Killed_Elements
