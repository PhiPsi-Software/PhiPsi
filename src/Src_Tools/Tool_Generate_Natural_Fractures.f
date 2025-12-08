 
      subroutine Tool_Generate_Natural_Fractures(
     &               M_x_min,M_x_max,M_y_min,M_y_max,
     &               Length_Cr,Orien_Cr,Length_Cr_Delta,Orien_Cr_Delta,
     &               c_num_Na_Cr,c_R,c_Na_Cr)

      use Global_Float_Type
      use Global_Common
      use Global_Model
      
      implicit none
      integer,intent(in)::c_num_Na_Cr
      real(kind=FT),intent(in)::M_x_min,M_x_max,M_y_min,M_y_max
      real(kind=FT),intent(in)::Length_Cr,Orien_Cr
      real(kind=FT),intent(in)::Length_Cr_Delta,Orien_Cr_Delta
      real(kind=FT),intent(in)::c_R
      real(kind=FT),intent(out)::c_Na_Cr(c_num_Na_Cr,2,2)
      real(kind=FT) c_Na_Cr_Center(c_num_Na_Cr,2)
      integer num_Cr,i_Try,i_Check
      real(kind=FT) c_x,c_y,Ran_Num,check_x,check_y,tem_Dis
      real(kind=FT) c_Tip_1(2),c_Tip_2(2),Orien_Cr_Rad
      real(kind=FT) Orien_Cr_Delta_Rad
      logical Yes_In_Circle,Yes_Tip1_In_Model,Yes_Tip2_In_Model
      real(kind=FT) U_1,U_2,UUU1,UUU2,R_Length_Cr,R_Orien_Cr_Rad
      
      print *,'    Generating natural fractures......'
      
      c_Na_Cr(1:c_num_Na_Cr,1:2,1:2) = ZR
      c_Na_Cr_Center(1:c_num_Na_Cr,1:2) = ZR
      num_Cr = 0
      Orien_Cr_Rad = Orien_Cr*pi/180.0D0
      Orien_Cr_Delta_Rad = Orien_Cr_Delta*pi/180.0D0
      
      
      if (Length_Cr > ((M_x_max-M_x_min)+(M_y_max-M_y_min)))then
           print *,'    Error :: too long intial natural fracture!'
           call Warning_Message('S',Keywords_Blank) 
      endif
      
      
      if(key_Random==1)then
          call Init_Random_Seed()
      endif
      
      do i_Try = 1,30000000
          
          if(key_Random==1 .or. key_Random==0)then
              call random_number (Ran_Num)
          elseif(key_Random==2 .or. key_Random==3)then
              call Tool_Generate_Random_Number(Ran_Num)
          endif
          
          c_x =  M_x_min +(M_x_max-M_x_min)*Ran_Num
          
          if(key_Random==1 .or. key_Random==0)then
              call random_number (Ran_Num)
          elseif(key_Random==2 .or. key_Random==3)then
              call Tool_Generate_Random_Number(Ran_Num)
          endif
          
          c_y =  M_y_min +(M_y_max-M_y_min)*Ran_Num
 
          Yes_In_Circle = .False.
          do i_Check=1,num_Cr
              check_x = c_Na_Cr_Center(i_Check,1)
              check_y = c_Na_Cr_Center(i_Check,2)
              tem_Dis = sqrt((c_x-check_x)**2 + (c_y-check_y)**2)
              if(tem_Dis<=c_R)then
                  Yes_In_Circle = .True.
                  exit
              endif
          enddo
          if(Yes_In_Circle .eqv. .False.)then
              Yes_Tip1_In_Model =  .False.
              Yes_Tip2_In_Model =  .False.
              
              if(key_Random==1 .or. key_Random==0)then
                  call random_number (U_1) 
                  call random_number (U_2)
              elseif(key_Random==2 .or. key_Random==3)then
                  call Tool_Generate_Random_Number(U_1)
                  call Tool_Generate_Random_Number(U_2)
              endif
              
              UUU1 = sqrt(-TWO*log(U_1)) * cos(TWO*PI*U_2)
              UUU2 = sqrt(-TWO*log(U_1)) * sin(TWO*PI*U_2)
              if(UUU1 >THR)UUU1=THR
              if(UUU1 <-THR)UUU1=-THR
              if(UUU2 >THR)UUU2=THR
              if(UUU2 <-THR)UUU2=-THR
              R_Length_Cr = Length_Cr+Length_Cr_Delta*UUU1/THR
              R_Orien_Cr_Rad=Orien_Cr_Rad+Orien_Cr_Delta_Rad*UUU2/THR
              c_Tip_1(1) = c_x + R_Length_Cr*HLF*cos(R_Orien_Cr_Rad)
              c_Tip_1(2) = c_y + R_Length_Cr*HLF*sin(R_Orien_Cr_Rad)
              c_Tip_2(1) = c_x - R_Length_Cr*HLF*cos(R_Orien_Cr_Rad)
              c_Tip_2(2) = c_y - R_Length_Cr*HLF*sin(R_Orien_Cr_Rad)
              if(c_Tip_1(1) > M_x_min .and. c_Tip_1(1)< M_x_max .and.
     &           c_Tip_1(2) > M_y_min .and. c_Tip_1(2)< M_y_max )then
                  Yes_Tip1_In_Model = .True.
              endif
              if(c_Tip_2(1) > M_x_min .and. c_Tip_2(1)< M_x_max .and.
     &           c_Tip_2(2) > M_y_min .and. c_Tip_2(2)< M_y_max )then
                  Yes_Tip2_In_Model = .True.
              endif
              if(Yes_Tip1_In_Model .and. Yes_Tip2_In_Model)then
                  num_Cr = num_Cr +1
                  c_Na_Cr_Center(num_Cr,1) = c_x 
                  c_Na_Cr_Center(num_Cr,2) = c_y
                  c_Na_Cr(num_Cr,1,1:2) =c_Tip_1(1:2)
                  c_Na_Cr(num_Cr,2,1:2) =c_Tip_2(1:2)
              endif
          endif
          if(num_Cr>=c_num_Na_Cr)then
              exit
          endif
      enddo
      if(num_Cr < c_num_Na_Cr)then
          print *,'    Warning :: natural fractures generation failed!'
          print *,'               natural fractures are too many or'
     &                       //' too long!'
          call Warning_Message('S',Keywords_Blank) 
      endif
      
      return 
      end SUBROUTINE Tool_Generate_Natural_Fractures