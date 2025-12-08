 
subroutine Tool_Generate_Natural_Fractures_3D





use Global_Float_Type
use Global_Common
use Global_Model
use Global_Crack_Common
use Global_Crack_3D

implicit none
real(kind=FT) Orien_Cr_Delta_Rad
real(kind=FT) U_1,U_2,UUU
real(kind=FT) V_1,V_2,VVV
real(kind=FT) W_1,W_2,WWW
real(kind=FT) X_1,X_2,XXX
real(kind=FT) Y_1,Y_2,YYY
real(kind=FT) c_Na_Cr(num_Rand_Na_Crack,4,3)
real(kind=FT) temp_size(3)
real(kind=FT) M_x_min,M_x_max,M_y_min,M_y_max,M_z_min,M_z_max 
real(kind=FT) c_x,c_y,c_z
real(kind=FT) Ran_Num
integer num_Cr,i_Try,i_Check
real(kind=FT) c_Na_Cr_Center(num_Rand_Na_Crack,3)
real(kind=FT) check_x,check_y,check_z,tem_Dis
real(kind=FT) c_Crack_Size,c_n_Vector(3),c_Crack_Center(3)
real(kind=FT) c_Crack3D_Coor(200,3)
real(kind=FT) c_Point(3)
integer i_P
logical Yes_Points_In_Model,Yes_In_Circle
real(kind=FT) c_Rotate_Degree
real(kind=FT) c_Rotate_Degree_for_Longside_Vector
real(kind=FT) Modified_Vector(3)
real(kind=FT) Rotate_Axis(3)
real(kind=FT) Rand_Point(3)
integer c_num_Crack
real(kind=FT) c_NaCr_3D_Rect_L,c_NaCr_3D_Rect_W
real(kind=FT) NaCr_3D_Rect_Longside_Vector_Modified(3)

if(num_Rand_Na_Crack<=0)then
  return
endif    

print *,'    Generating natural fractures......'

if(Key_NaCr_Type_3D==1)then
  allocate(Na_Crack3D_Coor(num_Rand_Na_Crack,4,3))   
  allocate(Each_NaCr3D_Poi_Num(num_Rand_Na_Crack))  
  Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack) =  4
  allocate(NaCr3D_Status(num_Rand_Na_Crack,5))
  NaCr3D_Status(1:num_Rand_Na_Crack,1:5) = 0
  
elseif(Key_NaCr_Type_3D==2)then
  allocate(Na_Crack3D_Cir_Coor(num_Rand_Na_Crack,7))   
  allocate(Na_Crack3D_Coor(num_Rand_Na_Crack,Circle_3D_Eqv_Polygon_Resolution,3))   
  allocate(Each_NaCr3D_Poi_Num(num_Rand_Na_Crack))  
  Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack) = Circle_3D_Eqv_Polygon_Resolution
  allocate(NaCr3D_Status(num_Rand_Na_Crack,5))
  NaCr3D_Status(1:num_Rand_Na_Crack,1:5) = 0
  
elseif(Key_NaCr_Type_3D==3)then
  allocate(Na_Crack3D_Coor(num_Rand_Na_Crack,Num_Poly_Edges_NaCr,3))   
  allocate(Each_NaCr3D_Poi_Num(num_Rand_Na_Crack))  
  Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack) =  Num_Poly_Edges_NaCr
  allocate(NaCr3D_Status(num_Rand_Na_Crack,5))
  NaCr3D_Status(1:num_Rand_Na_Crack,1:5) = 0
  
elseif(Key_NaCr_Type_3D==4)then
  allocate(Na_Crack3D_Coor(num_Rand_Na_Crack,4,3))   
  allocate(Each_NaCr3D_Poi_Num(num_Rand_Na_Crack))  
  Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack) =  4
  allocate(NaCr3D_Status(num_Rand_Na_Crack,5))
  NaCr3D_Status(1:num_Rand_Na_Crack,1:5) = 0
endif 


c_Na_Cr(1:num_Rand_Na_Crack,1:4,1:3) = ZR
num_Cr = 0
Orien_Cr_Delta_Rad = NaCr_3D_n_Vector_Delta*pi/180.0D0
if (Key_Fracture_Zone==0)then
  M_x_min = Min_X_Coor
  M_x_max = Max_X_Coor
  M_y_min = Min_Y_Coor
  M_y_max = Max_Y_Coor
  M_z_min = Min_Z_Coor
  M_z_max = Max_Z_Coor
elseif(Key_Fracture_Zone==1)then
  M_x_min = Frac_Zone_Minx
  M_x_max = Frac_Zone_Maxx
  M_y_min = Frac_Zone_Miny
  M_y_max = Frac_Zone_Maxy
  M_z_min = Frac_Zone_Minz
  M_z_max = Frac_Zone_Maxz     
endif

c_num_Crack = 0

temp_size(1) = Frac_Zone_Maxx - Frac_Zone_Minx
temp_size(2) = Frac_Zone_Maxy - Frac_Zone_Miny
temp_size(3) = Frac_Zone_Maxz - Frac_Zone_Minz
if (Key_Fracture_Zone==1)then
  if (NaCr_3D_Size > minval(temp_size(1:3)))then
       print *,'    Error :: too long intial natural fracture!'
       call Warning_Message('S',Keywords_Blank) 
  endif
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

  if(key_Random==1 .or. key_Random==0)then
      call random_number (Ran_Num)
  elseif(key_Random==2 .or. key_Random==3)then
      call Tool_Generate_Random_Number(Ran_Num)
  endif
  
  c_z =  M_z_min +(M_z_max-M_z_min)*Ran_Num
  
  c_Crack_Center(1:3) = [c_x,c_y,c_z]
  
  Yes_In_Circle = .False.
  if(Key_NaCr_Cross==0) then
      do i_Check=1,num_Cr
          check_x = c_Na_Cr_Center(i_Check,1)
          check_y = c_Na_Cr_Center(i_Check,2)
          check_z = c_Na_Cr_Center(i_Check,3)
          tem_Dis = sqrt((c_x-check_x)**2 + (c_y-check_y)**2 +(c_z-check_z)**2)
          if(tem_Dis<= NaCr_3D_Check_R)then
              Yes_In_Circle = .True.
              exit
          endif
      enddo
  elseif(Key_NaCr_Cross==1)then
      Yes_In_Circle = .False.
  endif
  
  if(Yes_In_Circle .eqv. .False.)then
      if(key_Random==1 .or. key_Random==0)then
          call random_number (U_1) 
          call random_number (U_2)
          call random_number (V_1) 
          call random_number (V_2)     
          if(Key_NaCr_Type_3D==4)then
              call random_number (W_1) 
              call random_number (W_2)  
              call random_number (X_1) 
              call random_number (X_2)  
          endif
      elseif(key_Random==2 .or. key_Random==3)then
          call Tool_Generate_Random_Number(U_1)
          call Tool_Generate_Random_Number(U_2)
          call Tool_Generate_Random_Number(V_1)
          call Tool_Generate_Random_Number(V_2)  
          if(Key_NaCr_Type_3D==4)then
              call Tool_Generate_Random_Number(W_1)
              call Tool_Generate_Random_Number(W_2)
              call Tool_Generate_Random_Number(X_1)
              call Tool_Generate_Random_Number(X_2) 
              call Tool_Generate_Random_Number(Y_1)
              call Tool_Generate_Random_Number(Y_2)   
          endif
      endif
      
      UUU = sqrt(-TWO*log(U_1)) * cos(TWO*PI*U_2)
      VVV = sqrt(-TWO*log(V_1)) * cos(TWO*PI*V_2)
      if(UUU > THR) UUU=  THR
      if(UUU <-THR) UUU= -THR
      if(VVV > THR) VVV=  THR
      if(VVV <-THR) VVV= -THR
      
      if(Key_NaCr_Type_3D==4)then
          WWW = sqrt(-TWO*log(W_1)) * cos(TWO*PI*W_2)
          XXX = sqrt(-TWO*log(X_1)) * cos(TWO*PI*X_2)
          YYY = sqrt(-TWO*log(Y_1)) * cos(TWO*PI*Y_2)
          if(WWW > THR) WWW=  THR
          if(WWW <-THR) WWW= -THR
          if(XXX > THR) XXX=  THR
          if(XXX <-THR) XXX= -THR
          if(YYY > THR) YYY=  THR
          if(YYY <-THR) YYY= -THR
      endif
      
      c_Crack_Size = NaCr_3D_Size +NaCr_3D_Sz_Delta*UUU/THR
      
      if(Key_NaCr_Type_3D==4)then
          c_NaCr_3D_Rect_L= NaCr_3D_Rect_L +NaCr_3D_Rect_L_Delta*WWW/THR
          c_NaCr_3D_Rect_W= NaCr_3D_Rect_W +NaCr_3D_Rect_W_Delta*XXX/THR
          if(abs(NaCr_3D_Rect_Longside_Vector_Delta)<=Tol_10)then
              NaCr_3D_Rect_Longside_Vector_Modified = NaCr_3D_Rect_Longside_Vector
          else
              c_Rotate_Degree_for_Longside_Vector = NaCr_3D_Rect_Longside_Vector_Delta*YYY/THR 
              c_Rotate_Degree_for_Longside_Vector = c_Rotate_Degree_for_Longside_Vector*pi/Con_180
              call Tool_Generate_Rand_Point_in_3D_Circle(c_Crack_Center(1:3), &
                     NaCr_3D_Rect_Longside_Vector(1:3),TWO,Rand_Point)
              Rotate_Axis(1:3)  = Rand_Point - c_Crack_Center(1:3)
              call Vector_Normalize(3,Rotate_Axis(1:3))   
              call Tool_Vector_a_Rotate_around_Vector_b(NaCr_3D_Rect_Longside_Vector(1:3), &
                      Rotate_Axis,c_Rotate_Degree_for_Longside_Vector,NaCr_3D_Rect_Longside_Vector_Modified)
          endif
      endif
      
      c_Rotate_Degree = NaCr_3D_n_Vector_Delta*VVV/THR
      c_Rotate_Degree = c_Rotate_Degree*pi/Con_180
      
      
      
      call Tool_Generate_Rand_Point_in_3D_Circle(c_Crack_Center(1:3), &
                 NaCr_3D_n_Vector(1:3),TWO,Rand_Point)
      
      Rotate_Axis(1:3)  = Rand_Point - c_Crack_Center(1:3)
      call Vector_Normalize(3,Rotate_Axis(1:3))   
      
      
      call Tool_Vector_a_Rotate_around_Vector_b(NaCr_3D_n_Vector(1:3), &
                  Rotate_Axis,c_Rotate_Degree,Modified_Vector)
      c_n_Vector(1:3) = Modified_Vector(1:3)
      
          
      if(Key_NaCr_Type_3D==1)then
          call Tool_Get_3D_Rect_by_n_Vector_and_Center(c_n_Vector(1:3), &
                          c_Crack_Center(1:3),c_Crack_Size,c_Crack3D_Coor(1:4,1:3))              
          Yes_Points_In_Model = .True.
          do  i_P = 1,4
              c_Point = c_Crack3D_Coor(i_P,1:3)
              if(c_Point(1) > M_x_min .and. c_Point(1) < M_x_max .and. &
                 c_Point(2) > M_y_min .and. c_Point(2) < M_y_max .and. &
                 c_Point(3) > M_z_min .and.c_Point(3) < M_z_max)then
              else
                  Yes_Points_In_Model = .False.
                  exit
              endif
          enddo
          if( Yes_Points_In_Model .eqv..True.)then
              num_Cr = num_Cr +1
              c_Na_Cr_Center(num_Cr,1:3) = c_Crack_Center(1:3)
              
             
             c_num_Crack = c_num_Crack + 1
             Na_Crack3D_Coor(c_num_Crack,1,1:3)=c_Crack3D_Coor(1,1:3)
             Na_Crack3D_Coor(c_num_Crack,2,1:3)=c_Crack3D_Coor(2,1:3)
             Na_Crack3D_Coor(c_num_Crack,3,1:3)=c_Crack3D_Coor(3,1:3)
             Na_Crack3D_Coor(c_num_Crack,4,1:3)=c_Crack3D_Coor(4,1:3)    
          endif
      elseif(Key_NaCr_Type_3D==2)then

          
          call Tool_Get_3D_Polygon_by_n_Vector_and_Center(c_n_Vector(1:3),  &
                           Circle_3D_Eqv_Polygon_Resolution, &
                           c_Crack_Center(1:3),c_Crack_Size, &
                           c_Crack3D_Coor(1:Circle_3D_Eqv_Polygon_Resolution,1:3))    

          Yes_Points_In_Model = .True.
          do  i_P = 1,Circle_3D_Eqv_Polygon_Resolution
              c_Point = c_Crack3D_Coor(i_P,1:3)
              if(c_Point(1) > M_x_min .and. c_Point(1) < M_x_max .and. &
                 c_Point(2) > M_y_min .and. c_Point(2) < M_y_max .and. &
                 c_Point(3) > M_z_min .and. c_Point(3) < M_z_max)then
              else
                  Yes_Points_In_Model = .False.
                  exit
              endif
          enddo
          if( Yes_Points_In_Model .eqv..True.)then
              num_Cr = num_Cr +1
              c_Na_Cr_Center(num_Cr,1:3) = c_Crack_Center(1:3)
              c_num_Crack = c_num_Crack + 1
              Na_Crack3D_Coor(c_num_Crack,1:Circle_3D_Eqv_Polygon_Resolution,1:3)= &
                             c_Crack3D_Coor(1:Circle_3D_Eqv_Polygon_Resolution,1:3)

          endif          
          
          
      elseif(Key_NaCr_Type_3D==3)then
          call Tool_Get_3D_Polygon_by_n_Vector_and_Center(c_n_Vector(1:3),Num_Poly_Edges_NaCr, &
                     c_Crack_Center(1:3),c_Crack_Size,c_Crack3D_Coor(1:Num_Poly_Edges_NaCr,1:3))    
          
          Yes_Points_In_Model = .True.
          do  i_P = 1,Num_Poly_Edges_NaCr
              c_Point = c_Crack3D_Coor(i_P,1:3)
              if(c_Point(1) > M_x_min .and. c_Point(1) < M_x_max .and. &
                 c_Point(2) > M_y_min .and. c_Point(2) < M_y_max .and. &
                 c_Point(3) > M_z_min .and. c_Point(3) < M_z_max)then
              else
                  Yes_Points_In_Model = .False.
                  exit
              endif
          enddo
          if( Yes_Points_In_Model .eqv..True.)then
              num_Cr = num_Cr +1
              c_Na_Cr_Center(num_Cr,1:3) = c_Crack_Center(1:3)
              
             
             c_num_Crack = c_num_Crack + 1
             Na_Crack3D_Coor(c_num_Crack,1:Num_Poly_Edges_NaCr,1:3)= c_Crack3D_Coor(1:Num_Poly_Edges_NaCr,1:3)     
          endif
      elseif(Key_NaCr_Type_3D==4)then    
          call Tool_Get_3D_Rect_by_n_Vector_Center_LVector_LW_Ratio(c_n_Vector(1:3),&
                 c_Crack_Center(1:3),c_NaCr_3D_Rect_L,c_NaCr_3D_Rect_W,&
                 NaCr_3D_Rect_Longside_Vector_Modified(1:3),c_Crack3D_Coor(1:4,1:3))
          Yes_Points_In_Model = .True.
          do  i_P = 1,4
              c_Point = c_Crack3D_Coor(i_P,1:3)
              if(c_Point(1) > M_x_min .and. c_Point(1) < M_x_max .and. &
                 c_Point(2) > M_y_min .and. c_Point(2) < M_y_max .and. &
                 c_Point(3) > M_z_min .and.c_Point(3) < M_z_max)then
              else
                  Yes_Points_In_Model = .False.
                  exit
              endif
          enddo
          if( Yes_Points_In_Model .eqv..True.)then
              num_Cr = num_Cr +1
              c_Na_Cr_Center(num_Cr,1:3) = c_Crack_Center(1:3)
              
             c_num_Crack = c_num_Crack + 1
             Na_Crack3D_Coor(c_num_Crack,1,1:3)=c_Crack3D_Coor(1,1:3)
             Na_Crack3D_Coor(c_num_Crack,2,1:3)=c_Crack3D_Coor(2,1:3)
             Na_Crack3D_Coor(c_num_Crack,3,1:3)=c_Crack3D_Coor(3,1:3)
             Na_Crack3D_Coor(c_num_Crack,4,1:3)=c_Crack3D_Coor(4,1:3)    
          endif
      endif
  endif
  if(num_Cr>=num_Rand_Na_Crack)then
      exit
  endif
enddo

if(num_Cr < num_Rand_Na_Crack)then
  print *,'    Warning :: natural fractures generation failed!'
  print *,'               natural fractures are too many or too long!'
  call Warning_Message('S',Keywords_Blank) 
endif

if(Key_NaCr_Type_3D==1 .or. Key_NaCr_Type_3D==4 .or. Key_NaCr_Type_3D==3)then
    call Save_Files_NF_3D(1,1)
endif

if(Key_NaCr_Active_Scheme_3D == 1) then
    if(Key_NaCr_Type_3D==1 .or. Key_NaCr_Type_3D==4)then
        Crack3D_Coor(num_Crack+1:num_Crack+num_Rand_Na_Crack,1:4,1:3)=Na_Crack3D_Coor(1:num_Rand_Na_Crack,1:4,1:3)
        Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,1) = 2
        Key_CS_Crack(num_Crack+1:num_Crack+num_Rand_Na_Crack) = 1
        if(Key_NaCr_Growth==0)then
            Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,3)=0
        elseif(Key_NaCr_Growth==1) then
            Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,3)=1
        endif
        Each_Cr_Poi_Num(num_Crack+1:num_Crack+num_Rand_Na_Crack)   = 4    
        num_Crack = num_Crack + num_Rand_Na_Crack
    elseif(Key_NaCr_Type_3D==2)then
        
        Crack3D_Coor(num_Crack+1:num_Crack+num_Rand_Na_Crack,1:Circle_3D_Eqv_Polygon_Resolution,1:3)=&
             Na_Crack3D_Coor(1:num_Rand_Na_Crack,1:Circle_3D_Eqv_Polygon_Resolution,1:3)
        Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,1) = 2
        Key_CS_Crack(num_Crack+1:num_Crack+num_Rand_Na_Crack) = 1
        if(Key_NaCr_Growth==0)then
            Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,3)=0
        elseif(Key_NaCr_Growth==1) then
            Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,3)=1
        endif
        Each_Cr_Poi_Num(num_Crack+1:num_Crack+num_Rand_Na_Crack)   = Circle_3D_Eqv_Polygon_Resolution   
        num_Crack = num_Crack + num_Rand_Na_Crack
        
    elseif(Key_NaCr_Type_3D==3)then
        Crack3D_Coor(num_Crack+1:num_Crack+num_Rand_Na_Crack,1:Num_Poly_Edges_NaCr,1:3)=&
             Na_Crack3D_Coor(1:num_Rand_Na_Crack,1:Num_Poly_Edges_NaCr,1:3)
        Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,1) = 2
        Key_CS_Crack(num_Crack+1:num_Crack+num_Rand_Na_Crack) = 1
        if(Key_NaCr_Growth==0)then
            Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,3)=0
        elseif(Key_NaCr_Growth==1) then
            Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,3)=1
        endif
        Each_Cr_Poi_Num(num_Crack+1:num_Crack+num_Rand_Na_Crack)   = Num_Poly_Edges_NaCr   
        num_Crack = num_Crack + num_Rand_Na_Crack
    endif 

    NaCr3D_Status(1:num_Rand_Na_Crack,1) = 1
elseif(Key_NaCr_Active_Scheme_3D == 2) then

elseif(Key_NaCr_Active_Scheme_3D == 3) then  
  
  allocate(Na_Crack3D_Ele_List(num_Rand_Na_Crack))
endif


return 
end SUBROUTINE Tool_Generate_Natural_Fractures_3D