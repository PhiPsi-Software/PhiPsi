 
subroutine Tool_Set_Natural_Fractures_by_kpp_3D






use Global_Float_Type
use Global_Common
use Global_Model
use Global_Crack_Common
use Global_Crack_3D
use Global_Read_kpp

implicit none

if(num_Rand_Na_Crack<=0)then
  return
endif    

print *,'    Set natural fractures accroding to kpp files......'

if(Key_NaCr_Type_3D==1)then
  allocate(Na_Crack3D_Coor(num_Rand_Na_Crack,4,3))   
  allocate(Each_NaCr3D_Poi_Num(num_Rand_Na_Crack))  
  Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack) =  4
  allocate(NaCr3D_Status(num_Rand_Na_Crack,5))
  NaCr3D_Status(1:num_Rand_Na_Crack,1:5) = 0
  Na_Crack3D_Coor(1:num_Rand_Na_Crack,1:4,1:3) = Read_kpp_Na_Crack3D_Coor(1:num_Rand_Na_Crack,1:4,1:3)
  Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack)     = Read_kpp_Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack)    
  
elseif(Key_NaCr_Type_3D==2)then
  allocate(Na_Crack3D_Cir_Coor(num_Rand_Na_Crack,7))   
  allocate(NaCr3D_Status(num_Rand_Na_Crack,5))
  NaCr3D_Status(1:num_Rand_Na_Crack,1:5) = 0
  Na_Crack3D_Cir_Coor(1:num_Rand_Na_Crack,1:7) = Read_kpp_Na_Crack3D_Cir_Coor(1:num_Rand_Na_Crack,1:7) 
  
elseif(Key_NaCr_Type_3D==3)then
  allocate(Na_Crack3D_Coor(num_Rand_Na_Crack,Num_Poly_Edges_NaCr,3))   
  allocate(Each_NaCr3D_Poi_Num(num_Rand_Na_Crack))  
  Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack) =  Num_Poly_Edges_NaCr
  allocate(NaCr3D_Status(num_Rand_Na_Crack,5))
  NaCr3D_Status(1:num_Rand_Na_Crack,1:5) = 0
  Na_Crack3D_Coor(1:num_Rand_Na_Crack,1:4,1:3) = Read_kpp_Na_Crack3D_Coor(1:num_Rand_Na_Crack,1:4,1:3)
  Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack)     = Read_kpp_Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack)   
  
elseif(Key_NaCr_Type_3D==4)then
  allocate(Na_Crack3D_Coor(num_Rand_Na_Crack,4,3))   
  allocate(Each_NaCr3D_Poi_Num(num_Rand_Na_Crack))  
  Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack) =  4
  allocate(NaCr3D_Status(num_Rand_Na_Crack,5))
  NaCr3D_Status(1:num_Rand_Na_Crack,1:5) = 0
  Na_Crack3D_Coor(1:num_Rand_Na_Crack,1:4,1:3) = Read_kpp_Na_Crack3D_Coor(1:num_Rand_Na_Crack,1:4,1:3)
  Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack)     = Read_kpp_Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack)   
endif 

if(Key_NaCr_Type_3D==1 .or. Key_NaCr_Type_3D==4 .or. Key_NaCr_Type_3D==3)then
    call Save_Files_NF_3D(1,1)
endif

if(Key_NaCr_Active_Scheme_3D == 1) then
    if(Key_NaCr_Type_3D==1 .or. Key_NaCr_Type_3D==4)then
        Crack3D_Coor(num_Crack+1:num_Crack+num_Rand_Na_Crack,1:4,1:3)= Na_Crack3D_Coor(1:num_Rand_Na_Crack,1:4,1:3)
        Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,1) = 2
        Key_CS_Crack(num_Crack+1:num_Crack+num_Rand_Na_Crack) = 1
        if(Key_NaCr_Growth==0)then
            Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,3)=0
        elseif(Key_NaCr_Growth==1) then
            Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,3)=1
        endif
        Each_Cr_Poi_Num(num_Crack+1:num_Crack+num_Rand_Na_Crack)   = Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack)    
        num_Crack = num_Crack + num_Rand_Na_Crack
    elseif(Key_NaCr_Type_3D==2)then
        Crack3D_Cir_Coor(num_Crack+1:num_Crack+num_Rand_Na_Crack,1:7) = Na_Crack3D_Cir_Coor(1:num_Rand_Na_Crack,1:7)
        Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,1) = 2
        Key_CS_Crack(num_Crack+1:num_Crack+num_Rand_Na_Crack) = 1
        if(Key_NaCr_Growth==0)then
            Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,3)=0
        elseif(Key_NaCr_Growth==1) then
            Crack_Type_Status_3D(num_Crack+1:num_Crack+num_Rand_Na_Crack,3)=1
        endif   
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
        Each_Cr_Poi_Num(num_Crack+1:num_Crack+num_Rand_Na_Crack)   = Each_NaCr3D_Poi_Num(1:num_Rand_Na_Crack)  
        num_Crack = num_Crack + num_Rand_Na_Crack
    endif 

    NaCr3D_Status(1:num_Rand_Na_Crack,1) = 1
    
elseif(Key_NaCr_Active_Scheme_3D == 2) then

elseif(Key_NaCr_Active_Scheme_3D == 3) then  
  
  allocate(Na_Crack3D_Ele_List(num_Rand_Na_Crack))
endif
  
return 
end SUBROUTINE Tool_Set_Natural_Fractures_by_kpp_3D