 
      SUBROUTINE Save_Files_Inclusions(isub)
      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_Model
      use Global_Filename
      use Global_Inclusion
      use Global_POST
      
      implicit none
      integer isub
      integer i,j
      character(200) c_File_name_1
      character(5) temp   

      
      if (Key_Save_Nothing==1) return 
      
      if(num_Circ_Incl/=0)then
          select case(Key_Data_Format)
          case(1:2)           
              print *,'    Saving coordinates of circle inclusions...'
              write(temp,'(I5)') isub
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                     //'jzcr'   
              open(101,file=c_File_name_1,status='unknown')    
              
              do i=1,num_Circ_Incl
                  write(101, '(3E20.12)') Circ_Inclu_Coor(i,1:3)
              end do
              close(101)     
              print *,'    Saving ennj file of inclusion...'
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'ennj'//'_'//ADJUSTL(temp)       
              open(103,file=c_File_name_1,status='unknown')         
              do i=1,Num_Node
                  write(103, '(200I10)') (Enriched_Node_Type_Incl(i,
     &                                    j),j=1,num_Circ_Incl)
              end do
              close(103)           
              print *,'    Saving eltj file of inclusion...'
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'eltj'//'_'//ADJUSTL(temp)       
              open(104,file=c_File_name_1,status='unknown')         
              do i=1,Num_Elem
                  write(104, '(200I10)') 
     &               (Elem_Type_Incl(i,j),j=1,num_Circ_Incl)
              end do
              close(104)   
              print *,'    Saving posj file of inclusion...'
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'posj'//'_'//ADJUSTL(temp)         
              open(110,file=c_File_name_1,status='unknown')      
              do i=1,Num_Node
                  write(110, '(200I10)')
     &                 (c_POS_Incl(i,j),j=1,num_Circ_Incl)
              end do
              close(110) 
          end select
      endif
      
      if(num_Poly_Incl/=0)then
          select case(Key_Data_Format)
          case(1:2)         
              print *,'    Saving x coordinates of poly inclusions...'
              write(temp,'(I5)') isub
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                     //'jzpx'   
              open(101,file=c_File_name_1,status='unknown')    
              do i=1,num_Poly_Incl
                  write(101, '(200E20.12)') (Poly_Incl_Coor_x(i,
     &                                 j),j=1,Poly_Inclu_Edges_Num(i))
              end do
              close(101)     
              print *,'    Saving y coordinates of poly inclusions...'
              write(temp,'(I5)') isub
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                     //'jzpy'   
              open(101,file=c_File_name_1,status='unknown')    
              do i=1,num_Poly_Incl
                  write(101, '(200E20.12)') (Poly_Incl_Coor_y(i,
     &                                 j),j=1,Poly_Inclu_Edges_Num(i))
              end do
              close(101)     
              print *,'    Saving ennj file of inclusion...'
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'ennj'//'_'//ADJUSTL(temp)       
              open(103,file=c_File_name_1,status='unknown')         
              do i=1,Num_Node
                  write(103, '(200I10)') (Enriched_Node_Type_Incl(i,
     &                                    j),j=1,num_Poly_Incl)
              end do
              close(103)           
              print *,'    Saving eltj file of inclusion...'
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'eltj'//'_'//ADJUSTL(temp)       
              open(104,file=c_File_name_1,status='unknown')         
              do i=1,Num_Elem
                  write(104, '(200I10)') 
     &               (Elem_Type_Incl(i,j),j=1,num_Poly_Incl)
              end do
              close(104)   
              print *,'    Saving posj file of inclusion...'
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'posj'//'_'//ADJUSTL(temp)         
              open(110,file=c_File_name_1,status='unknown')      
              do i=1,Num_Node
                  write(110, '(200I10)')
     &                 (c_POS_Incl(i,j),j=1,num_Poly_Incl)
              end do
              close(110) 
          
          end select
      endif      
      
      RETURN
      END SUBROUTINE Save_Files_Inclusions
