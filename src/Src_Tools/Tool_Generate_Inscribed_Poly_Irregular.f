 
      subroutine Tool_Generate_Inscribed_Poly_Irregular(
     &                 M_x_min,M_x_max,M_y_min,M_y_max,
     &                 c_num_Rand_Poly_Incl,
     &                 c_num_Vert_Poly_Incl,
     &                 c_num_Rand_Poly_Incl_for_Each_Type,
     &                 c_Rand_Poly_Incl_R_Min_and_Max,
     &                 c_Irregular_R_Delta_Factor,
     &                 c_Irregular_Angle_Factor,
     &                 c_Irregular_Extension_Factor,
     &                 c_Irregular_Inclination,
     %                 check_R_Factor,
     &                 check_Edge_Dis_Factor,
     &                 Poly_Incl_Coor_x,
     &                 Poly_Incl_Coor_y)


      use Global_Float_Type
      use Global_Common
      use Global_Model
      
      
      implicit none
      integer,intent(in)::c_num_Rand_Poly_Incl
      integer,intent(in)::c_num_Vert_Poly_Incl
      real(kind=FT),intent(in)::M_x_min,M_x_max,M_y_min,M_y_max
      integer,intent(in):: c_num_Rand_Poly_Incl_for_Each_Type(10)
      real(kind=FT),intent(in)::c_Rand_Poly_Incl_R_Min_and_Max(10,2)    
      real(kind=FT),intent(in)::c_Irregular_R_Delta_Factor
      real(kind=FT),intent(in)::c_Irregular_Angle_Factor
      real(kind=FT),intent(in)::c_Irregular_Extension_Factor
      real(kind=FT),intent(in)::c_Irregular_Inclination
      real(kind=FT) c_Delta_R
      real(kind=FT),intent(in)::check_R_Factor
      real(kind=FT),intent(out)::Poly_Incl_Coor_x(c_num_Rand_Poly_Incl,
     &                                            c_num_Vert_Poly_Incl)   
      real(kind=FT),intent(out)::Poly_Incl_Coor_y(c_num_Rand_Poly_Incl,
     &                                            c_num_Vert_Poly_Incl)   
      real(kind=FT),intent(in)::check_Edge_Dis_Factor
      
      integer num_Cir,i_Try,i_Check
      real(kind=FT) c_x,c_y,Ran_Num,check_x,check_y,tem_Dis
      logical Yes_In_Circle
      real(kind=FT) U_1,U_2,UUU1
      real(kind=FT)  try_R,try_x,try_y
      logical Yes_In_Model_x,Yes_In_Model_y
      real(kind=FT) c_Center_Cir(c_num_Rand_Poly_Incl,2),
     &              c_R_Cir(c_num_Rand_Poly_Incl)
      integer i_point
      real(kind=FT) delta_Theta
      integer i_Type,num_Circle_this_Type,num_Poly_this_Type
      real(kind=FT) min_R_this_Type,max_R_this_Type,Ave_R_this_Type
      
      real(kind=FT) Max_R
      real(kind=FT) check_Edge_Dis 
      real(kind=FT) check_R_Dis
      real(kind=FT) check_R 
      integer i_Cir,c_Cir
      real(kind=FT) Ran_n_Num(c_num_Vert_Poly_Incl),
     &              Ran_n_R(c_num_Vert_Poly_Incl),
     &              Ran_n_Angle(c_num_Vert_Poly_Incl)
      real(kind=FT) sum_Angle
      integer n_Vertex
      real(kind=FT) c_Theta
      integer i
      real(kind=FT) Elongation,Rotation,tem1,tem2,tem3,tem4,c_angle
      
      print *,'    Generating random irregular polygon......'
      
      Poly_Incl_Coor_x(1:c_num_Rand_Poly_Incl,1:c_num_Vert_Poly_Incl)=ZR
      Poly_Incl_Coor_y(1:c_num_Rand_Poly_Incl,1:c_num_Vert_Poly_Incl)=ZR
      c_Center_Cir(1:c_num_Rand_Poly_Incl,1:2) = ZR
      c_R_Cir(1:c_num_Rand_Poly_Incl)           = ZR
      num_Cir = 0
      Max_R = maxval(c_Rand_Poly_Incl_R_Min_and_Max(1:10,1:2))
      Elongation = c_Irregular_Extension_Factor
      Rotation   = c_Irregular_Inclination*pi/180.0D0
      
      if(key_Random==1)then
          call Init_Random_Seed()
      endif
      
      do i_Type = 1,10
          num_Circle_this_Type = 0
          num_Poly_this_Type=c_num_Rand_Poly_Incl_for_Each_Type(i_Type)
          min_R_this_Type = c_Rand_Poly_Incl_R_Min_and_Max(i_Type,1)
          max_R_this_Type = c_Rand_Poly_Incl_R_Min_and_Max(i_Type,2)
          Ave_R_this_Type = (min_R_this_Type + max_R_this_Type)/TWO
          if (num_Poly_this_Type>=1)then
            do i_Try = 1,50000000
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
              do i_Check=1,num_Cir
                  check_x = c_Center_Cir(i_Check,1)
                  check_y = c_Center_Cir(i_Check,2)
                  check_R = c_R_Cir(i_Check)  
                  tem_Dis = sqrt((c_x-check_x)**2 + (c_y-check_y)**2)
                  check_R_Dis=check_R_Factor*(Ave_R_this_Type+check_R)      
                  if(tem_Dis<=check_R_Dis)then
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
                  c_Delta_R = (max_R_this_Type-min_R_this_Type)/TWO
                  try_R = Ave_R_this_Type + c_Delta_R*UUU1/THR
                  Try_x = c_x
                  Try_y = c_y 
                  Yes_In_Model_x = .False.
                  Yes_In_Model_y = .False.
                  check_Edge_Dis = check_Edge_Dis_Factor*Max_R
                  if (Try_x+try_R<M_x_max-check_Edge_Dis .and.
     &                Try_x-try_R > M_x_min+check_Edge_Dis) then
                      Yes_In_Model_x = .True.
                  endif
                  if (Try_y+try_R<M_y_max -check_Edge_Dis .and.
     &                Try_y-try_R > M_y_min +check_Edge_Dis) then
                      Yes_In_Model_y = .True.
                  endif
                  if(Yes_In_Model_x .and. Yes_In_Model_y)then
                      num_Cir = num_Cir +1
                      num_Circle_this_Type = num_Circle_this_Type +1
                      c_Center_Cir(num_Cir,1) = c_x 
                      c_Center_Cir(num_Cir,2) = c_y
                      c_R_Cir(num_Cir)        = try_R
                  endif
              endif
              if(num_Circle_this_Type >= num_Poly_this_Type)then
                  exit
              endif
            enddo
          endif
      enddo
      if(num_Cir < c_num_Rand_Poly_Incl)then
          print *,'    Error :: Inscribed polygons generation failed!'
          call Warning_Message('S',Keywords_Blank)
      endif
      
      
      delta_Theta = TWO*pi/c_num_Vert_Poly_Incl
      n_Vertex = c_num_Vert_Poly_Incl
      c_Cir = 0
      do i_Type = 1,10
          num_Circle_this_Type = 0
          num_Poly_this_Type=c_num_Rand_Poly_Incl_for_Each_Type(i_Type)
          min_R_this_Type = c_Rand_Poly_Incl_R_Min_and_Max(i_Type,1)
          max_R_this_Type = c_Rand_Poly_Incl_R_Min_and_Max(i_Type,2)
          Ave_R_this_Type = (min_R_this_Type + max_R_this_Type)/TWO
          if (num_Poly_this_Type>=1)then
              do i_Cir = 1,num_Poly_this_Type
                  c_Cir = c_Cir +1
                  if(key_Random==1 .or. key_Random==0)then
                      call random_number (Ran_n_Num(1:n_Vertex)) 
                  elseif(key_Random==2 .or. key_Random==3)then
                      do i=1,n_Vertex
                        call Tool_Generate_Random_Number(Ran_n_Num(i))
                      enddo
                  endif
     
                  Ran_n_R(1:n_Vertex)= Ave_R_this_Type +
     &                         (TWO*Ran_n_Num(1:n_Vertex)-ONE)*
     &                         (max_R_this_Type -min_R_this_Type)*
     &                         c_Irregular_R_Delta_Factor             
                  if(key_Random==1 .or. key_Random==0)then
                      call random_number (Ran_n_Num(1:n_Vertex)) 
                  elseif(key_Random==2 .or. key_Random==3)then
                      do i=1,n_Vertex
                        call Tool_Generate_Random_Number(Ran_n_Num(i)) 
                      enddo
                  endif
                  
                  
                  Ran_n_Angle(1:n_Vertex)= delta_Theta +
     &                         (TWO*Ran_n_Num(1:n_Vertex)-ONE)*
     &                         delta_Theta*c_Irregular_Angle_Factor
                  sum_Angle =sum(Ran_n_Angle(1:n_Vertex))
                  Ran_n_Angle(1:n_Vertex)=
     &                   Ran_n_Angle(1:n_Vertex)*TWO*pi/sum_Angle
                  c_Theta =  ZR
                  do i_Point =1,c_num_Vert_Poly_Incl
                      c_Theta = c_Theta + Ran_n_Angle(i_Point)
                      tem1=
     &                   Ran_n_R(i_Point)*cos(c_Theta)*sqrt(Elongation)
                      tem2=
     &                   Ran_n_R(i_Point)*sin(c_Theta)/sqrt(Elongation)    
                      c_angle = atan2(tem2,tem1)
                      tem3 = sqrt(tem1**2+tem2**2)*cos(c_angle+Rotation)
                      tem4 = sqrt(tem1**2+tem2**2)*sin(c_angle+Rotation)
                      Poly_Incl_Coor_x(c_Cir,i_Point) = 
     &                  c_Center_Cir(c_Cir,1) + tem3                                
                      Poly_Incl_Coor_y(c_Cir,i_Point) = 
     &                    c_Center_Cir(c_Cir,2) + tem4   
                  enddo
                  
              enddo
          endif
      enddo
      
      return 
      end SUBROUTINE Tool_Generate_Inscribed_Poly_Irregular 
      

