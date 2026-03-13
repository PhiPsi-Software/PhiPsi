!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
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