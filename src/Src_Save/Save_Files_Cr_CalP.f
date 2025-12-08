 
      SUBROUTINE Save_Files_Cr_CalP(isub)

      use Global_Float_Type      
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Filename
      use Global_POST
      
      implicit none
      integer isub
      integer i,j
      character(200) c_File_name_1,c_File_name_2
      character(5) temp   
      real(kind=FT) tem_1 ,tem_2 
      
      if (Key_Save_Nothing==1) return
      
      
      tem_1 = 1000.0D0
      tem_2 = 1.0D6
      select case(Key_Data_Format)
      case(1:2)      
          print *,'    Saving coordinates of calculation points...'
          write(temp,'(I5)') isub
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'apex'//'_'//ADJUSTL(temp)  
          c_File_name_2   =  trim(Full_Pathname)//'.'
     &                 //'apey'//'_'//ADJUSTL(temp)      
          open(101,file=c_File_name_1,status='unknown')    
          open(102,file=c_File_name_2,status='unknown')  
          if(Key_Unit_System==1)then    
              do i=1,num_Crack
                  write(101, '(2000E20.12)') (Cracks_CalP_Coors(i,
     &                                    j,1),j=1,Cracks_CalP_Num(i))
                  write(102, '(2000E20.12)') (Cracks_CalP_Coors(i,
     &                                    j,2),j=1,Cracks_CalP_Num(i))
              end do
          elseif(Key_Unit_System==2)then
              do i=1,num_Crack
                  write(101, '(2000E20.12)') (Cracks_CalP_Coors(i,
     &                                j,1)/tem_1,j=1,Cracks_CalP_Num(i))
                  write(102, '(2000E20.12)') (Cracks_CalP_Coors(i,
     &                                j,2)/tem_1,j=1,Cracks_CalP_Num(i))
              end do
          endif
          close(101)
          close(102)            
          print *,'    Saving orientations of calculation points...'
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'cori'//'_'//ADJUSTL(temp)      
          open(101,file=c_File_name_1,status='unknown')         
          do i=1,num_Crack
              write(101,'(2000E20.12)') (Cracks_CalP_Orient(i,
     &                                   j),j=1,Cracks_CalP_Num(i))
          end do
          close(101)    
          print *,'    Saving apertures of calculation points...'
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'cape'//'_'//ADJUSTL(temp)      
          open(101,file=c_File_name_1,status='unknown') 
          if(Key_Unit_System==1)then    
              do i=1,num_Crack
                  write(101,'(2000E20.12)') (Cracks_CalP_Aper(i,
     &                                       j),j=1,Cracks_CalP_Num(i))
              end do
          elseif(Key_Unit_System==2)then
              do i=1,num_Crack
                  write(101,'(2000E20.12)') (Cracks_CalP_Aper(i,
     &                                 j)/tem_1,j=1,Cracks_CalP_Num(i))
              end do
          endif
          close(101) 
          
          if ((Key_Analysis_Type==3).or.(Key_Analysis_Type==4) .or.
     &        (Key_Analysis_Type==5) ) then
              print *,'    Saving pressures of calculation points...'
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                     //'cpre'//'_'//ADJUSTL(temp)      
              open(101,file=c_File_name_1,status='unknown')         
              if(Key_Unit_System==1)then    
                do i=1,num_Crack
                  write(101,'(2000E20.12)') (Cracks_CalP_Pres(i,
     &                                       j),j=1,Cracks_CalP_Num(i))
                end do
              elseif(Key_Unit_System==2)then 
                do i=1,num_Crack
                  write(101,'(2000E20.12)') (Cracks_CalP_Pres(i,
     &                                  j)*tem_2,j=1,Cracks_CalP_Num(i))
                end do
              endif
              close(101)       
              print *,'    Saving flow velocity of '
     &                     //'calculation points...'
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                     //'cvel'//'_'//ADJUSTL(temp)      
              open(111,file=c_File_name_1,status='unknown')         
              
              if(Key_Unit_System==1)then    
                do i=1,num_Crack
                  write(111,'(2000E20.12)') (Cracks_CalP_Velo(i,
     &                                       j),j=1,Cracks_CalP_Num(i))
                end do
              elseif(Key_Unit_System==2)then 
                  do i=1,num_Crack
                      write(111,'(2000E20.12)') (Cracks_CalP_Velo(i,
     &                                  j)/tem_1,j=1,Cracks_CalP_Num(i))
                  end do
              endif
              
              close(111)   
              print *,'    Saving quantity of flow of '
     &                     //'calculation points...'
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                     //'cqua'//'_'//ADJUSTL(temp)      
              open(112,file=c_File_name_1,status='unknown')         
              if(Key_Unit_System==1)then     
                do i=1,num_Crack
                  write(112,'(2000E20.12)') (Cracks_CalP_Quan(i,
     &                                       j),j=1,Cracks_CalP_Num(i))
                end do
              elseif(Key_Unit_System==2)then 
                do i=1,num_Crack
                  write(112,'(2000E20.12)') (Cracks_CalP_Quan(i,
     &                                 j)/tem_2,j=1,Cracks_CalP_Num(i))
                end do
              endif
              close(112)   
          end if
          if (Key_Analysis_Type==1 .and. Key_TipEnrich ==4 ) then
              print *,'    Saving tractions of calculation'
     &                     //' points of cohesive crack...'
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                     //'cohx'//'_'//ADJUSTL(temp)      
              open(101,file=c_File_name_1,status='unknown')   
              do i=1,num_Crack
                  write(101,'(2000E20.12)') (Cracks_CalP_Tractions(i,
     &                                   j,1),j=1,Cracks_CalP_Num(i))
              enddo
              close(101)    
              c_File_name_2   =  trim(Full_Pathname)//'.'
     &                     //'cohy'//'_'//ADJUSTL(temp)      
              open(102,file=c_File_name_2,status='unknown')    
              do i=1,num_Crack
                  write(102,'(2000E20.12)') (Cracks_CalP_Tractions(i,
     &                                   j,2),j=1,Cracks_CalP_Num(i))
              enddo
              close(102)   
          endif
      end select
 
      RETURN
      END SUBROUTINE Save_Files_Cr_CalP
