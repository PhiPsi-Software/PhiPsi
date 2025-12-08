 
      subroutine Tool_Generate_Circles(
     &               M_x_min,M_x_max,M_y_min,M_y_max,
     &               c_Ave_R_Ball,c_Delta_R_Ball,
     &               c_num_Cir,check_R,Max_Dis_to_Edge,
     &               c_Center_Cir,c_R_Cir)

      use Global_Float_Type
      use Global_Common
      use Global_Model
      
      implicit none
      integer,intent(in)::c_num_Cir
      real(kind=FT),intent(in)::M_x_min,M_x_max,M_y_min,M_y_max
      real(kind=FT),intent(in)::c_Ave_R_Ball
      real(kind=FT),intent(in)::c_Delta_R_Ball
      real(kind=FT),intent(in)::check_R
      real(kind=FT),intent(out)::c_Center_Cir(c_num_Cir,2)
      real(kind=FT),intent(in)::Max_Dis_to_Edge
      real(kind=FT) c_R_Cir(c_num_Cir)
      
      integer num_B,i_Try,i_Check
      real(kind=FT) c_x,c_y,Ran_Num,check_x,check_y,tem_Dis
      logical Yes_In_Circle
      real(kind=FT) U_1,U_2,UUU1
      real(kind=FT)  try_R,try_x,try_y
      logical Yes_In_Model_x,Yes_In_Model_y
      
      print *,'    Generating random circles......'
      
      c_Center_Cir(1:c_num_Cir,1:2) = ZR
      c_R_Cir(1:c_num_Cir) = ZR
      num_B = 0
      
      if(key_Random==1)then
          call Init_Random_Seed()
      endif
      do i_Try = 1,5000000
          
          
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
          do i_Check=1,num_B
              check_x = c_Center_Cir(i_Check,1)
              check_y = c_Center_Cir(i_Check,2)
              tem_Dis = sqrt((c_x-check_x)**2 + (c_y-check_y)**2)
              if(tem_Dis<=check_R)then
                  Yes_In_Circle = .True.
                  exit
              endif
          enddo
          if(Yes_In_Circle .eqv. .False.)then
              
              if(key_Random==1 .or. key_Random==0)then
                  call random_number (U_1) 
                  call random_number (U_2)
              elseif(key_Random==2 .or. key_Random==3)then
                  call Tool_Generate_Random_Number(U_1)
                  call Tool_Generate_Random_Number(U_2)
              endif
              
              UUU1 = sqrt(-TWO*log(U_1)) * cos(TWO*PI*U_2)
              if(UUU1 >THR)UUU1=THR
              if(UUU1 <-THR)UUU1=-THR
              try_R = c_Ave_R_Ball + c_Delta_R_Ball*UUU1/THR
              Try_x = c_x ; Try_y = c_y 
              Yes_In_Model_x = .False.
              Yes_In_Model_y = .False.
              if (Try_x+try_R<M_x_max-Max_Dis_to_Edge .and.
     &            Try_x-try_R > M_x_min+Max_Dis_to_Edge) then
                  Yes_In_Model_x = .True.
              endif
              if (Try_y+try_R<M_y_max -Max_Dis_to_Edge .and.
     &            Try_y-try_R > M_y_min +Max_Dis_to_Edge) then
                  Yes_In_Model_y = .True.
              endif
              if(Yes_In_Model_x .and. Yes_In_Model_y)then
                  num_B = num_B +1
                  c_Center_Cir(num_B,1) = c_x 
                  c_Center_Cir(num_B,2) = c_y
                  c_R_Cir(num_B)        = try_R
              endif
          endif
          if(num_B>=c_num_Cir)then
              exit
          endif
      enddo
      if(num_B < c_num_Cir)then
          print *,'    Error :: circles generation failed!'
          call Warning_Message('S',Keywords_Blank)
      endif
      
      return 
      end SUBROUTINE Tool_Generate_Circles 
      

