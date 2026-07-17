!-----------------------------------------------------------
! Brief: Initialize 3D natural fractures from KPP keyword data.
!
! Parameters:
!   (none)
!
! Notes:   Allocates the global natural-fracture arrays and
!   copies coordinates and status flags from the kpp read
!   buffers, dispatching on Key_NaCr_Type_3D to rectangular,
!   circular, polygonal or narrow-rectangular shapes.
!-----------------------------------------------------------

subroutine Tool_Set_Natural_Fractures_by_kpp_3D

! Set natural fractures based on the data in the KPP file. 2023-08-24.
!NEWFTU2023082401.
! Modified from Tool_Generate_Natural_Fractures_3D.f90.


! The global variable Crack_Type_Status_3D(i_C,10) is used to indicate the type and status of
! cracks.
! Column 1 (Fracture Type): =1, HF fracture; =2, natural fracture; =3, post-fracturing hydraulic
! fracture
! Note: Natural fractures and propped hydraulic fractures may potentially turn into HF fractures.
! Column 2 (Fracture Status): =1, HF fracturing not completed; =2, HF fracturing completed
! Column 3 (Can the crack continue to propagate): =1, yes; =0, no
! Column 4 (Whether the fracture has obtained a fluid node): =1, Yes; =0, No
! Column 5 (Did the crack propagate in the previous step?): =1, Yes; =0, No



!***************************
! Variable Type Declaration
!***************************
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
