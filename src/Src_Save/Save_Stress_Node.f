 
      SUBROUTINE Save_Stress_Node(isub,Key_CoorSys)

      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_Stress
      use Global_POST
      
      implicit none
      
      integer isub,Key_CoorSys
      integer i
      character(200) c_File_name_1
      character(5) temp  
      real(kind=FT) tem_2 
      
      
      if (Key_Save_Nothing==1) return 
      
      tem_2 = 1.0D6

      if (Key_CoorSys==1) then      
          print *,'    Saving stress of all nodes...'
          write(temp,'(I5)') isub
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                     //'strn'//'_'//ADJUSTL(temp)        
          
          select case(Key_Data_Format)
          case(1)  
              open(202,file=c_File_name_1,status='unknown')
              if (Key_Dimension == 2) then 
                  if(Key_Unit_System==1)then     
                      do i=1,num_Node
                          write(202, '(I8,4E20.12)') i,
     &                                Stress_xx_Node(i),
     &                                Stress_yy_Node(i),
     &                                Stress_xy_Node(i),    
     &                                Stress_vm_Node(i)
                      end do
                  elseif(Key_Unit_System==2)then
                      do i=1,num_Node
                          write(202, '(I8,4E20.12)') i,
     &                                Stress_xx_Node(i)*tem_2,
     &                                Stress_yy_Node(i)*tem_2,
     &                                Stress_xy_Node(i)*tem_2,    
     &                                Stress_vm_Node(i)*tem_2
                      end do
                  endif        
              else if(Key_Dimension == 3)then
                  do i=1,num_Node
                      write(202, '(6E20.12)') 
     &                            Stress_xx_Node(i),
     &                            Stress_yy_Node(i),
     &                            Stress_zz_Node(i),
     &                            Stress_xy_Node(i),
     &                            Stress_yz_Node(i),
     &                            Stress_xz_Node(i)
                  end do          
              end if          
              close(202)  
          case(2)  
              open(203,file=c_File_name_1,status='unknown',
     &                 form='unformatted',access='stream')     
              if (Key_Dimension == 2) then        
                  write(203) (Stress_xx_Node(i),
     &                        Stress_yy_Node(i),
     &                        Stress_xy_Node(i),    
     &                        Stress_vm_Node(i),i=1,num_Node)
              else if(Key_Dimension == 3)then
                  write(203)     (Stress_xx_Node(i),
     &                            Stress_yy_Node(i),
     &                            Stress_zz_Node(i),
     &                            Stress_xy_Node(i),
     &                            Stress_yz_Node(i),
     &                            Stress_xz_Node(i),       
     &                            Stress_vm_Node(i),i=1,num_Node)

              end if          
              close(203)      
          end select
          
          if(Key_Thermal_Stress==1)then
              print *,'    Saving thermal stress of all nodes...'
              write(temp,'(I5)') isub
              c_File_name_1   =  trim(Full_Pathname)//'.'
     &                         //'sttn'//'_'//ADJUSTL(temp)        
              
              select case(Key_Data_Format)
              case(1)  
                  open(202,file=c_File_name_1,status='unknown')
                  if (Key_Dimension == 2) then 
                      if(Key_Unit_System==1)then     
                          do i=1,num_Node
                              write(202, '(I8,4E20.12)') i,
     &                                    TStress_xx_Node(i),
     &                                    TStress_yy_Node(i),
     &                                    TStress_xy_Node(i),    
     &                                    TStress_vm_Node(i)
                          end do
                      elseif(Key_Unit_System==2)then 
                          do i=1,num_Node
                              write(202, '(I8,4E20.12)') i,
     &                                    TStress_xx_Node(i)*tem_2,
     &                                    TStress_yy_Node(i)*tem_2,
     &                                    TStress_xy_Node(i)*tem_2,    
     &                                    TStress_vm_Node(i)*tem_2
                          end do
                      endif       
                  else if(Key_Dimension == 3)then
        
                  end if          
                  close(202)  
              case(2)  
                  open(203,file=c_File_name_1,status='unknown',
     &                     form='unformatted',access='stream')     
                  if (Key_Dimension == 2) then        
                      write(203)     (TStress_xx_Node(i),
     &                                TStress_yy_Node(i),
     &                                TStress_vm_Node(i),i=1,num_Node)         
                  else if(Key_Dimension == 3)then
                  end if          
                  close(203)      
              end select
          endif

      elseif(Key_CoorSys==2) then
          print *,'    Saving nodes stress in cylindrical coordinate...'
          write(temp,'(I5)') isub
          c_File_name_1   =  trim(Full_Pathname)//'.'
     &                     //'stnc'//'_'//ADJUSTL(temp)        
          
          select case(Key_Data_Format)
          case(1)  
              open(202,file=c_File_name_1,status='unknown')
              if (Key_Dimension == 2) then 
        
              else if(Key_Dimension == 3)then
                  do i=1,num_Node
                      write(202, '(I8,7E20.12)') i,
     &                            Stress_Crr_Node(i),
     &                            Stress_Ctt_Node(i),
     &                            Stress_Czz_Node(i),
     &                            Stress_Crt_Node(i),
     &                            Stress_Ctz_Node(i),
     &                            Stress_Crz_Node(i),       
     &                            Stress_Cvm_Node(i)   
                  end do          
              end if          
              close(202)  
          case(2)  
              open(203,file=c_File_name_1,status='unknown',
     &                 form='unformatted',access='stream')     
              if (Key_Dimension == 2) then        
       
              else if(Key_Dimension == 3)then
                  write(203)     (Stress_Crr_Node(i),
     &                            Stress_Ctt_Node(i),
     &                            Stress_Czz_Node(i),
     &                            Stress_Crt_Node(i),
     &                            Stress_Ctz_Node(i),
     &                            Stress_Crz_Node(i),       
     &                            Stress_Cvm_Node(i),i=1,num_Node)

              end if          
              close(203)      
          end select
      endif
      RETURN
      END SUBROUTINE Save_Stress_Node
