 
subroutine Cal_Crack_Aperture_3D(isub,c_DISP)

use Global_Float_Type
use Global_Crack_Common
use Global_Crack_3D
use Global_Model
use Global_Elem_Area_Vol
use Global_Common
use Global_HF
use Global_Post
use Global_Filename
use omp_lib
use Global_Crack

use Global_Cal_Ele_Num_by_Coors_3D  

implicit none
real(kind=FT),intent(in)::c_DISP(Total_FD)
integer,intent(in)::isub
integer i_C,i_FluidEle
integer num_Flu_Nodes
real(kind=FT) Flu_Ele_Area,Flu_Ele_Vector(3),ori_n(3)
real(kind=FT) Flu_Ele_Centroid(3),delta_L
real(kind=FT) up_offset_P(3),lw_offset_P(3)
integer up_Elem_num,lw_Elem_num
real(kind=FT) up_Kesi,up_Yita,up_Zeta
real(kind=FT) lw_Kesi,lw_Yita,lw_Zeta
real(kind=FT) up_P_Disp(3),lw_P_Disp(3),Relative_Disp(3)
real(kind=FT) norm_Flu_Ele_Vector
real(kind=FT) c_Aperture
integer i_FluNode,c_FluNode
real(kind=FT) delta_L_Factor
real(kind=FT) c_El_Ave_Aperture
real(kind=FT) c_Node_Coor(3)
integer Fl_El_Nodes(15),num_All_Flu_Nodes
integer count_All_Flu_Nodes(Max_Max_N_CalP_3D) 
integer i_Node
logical Yes_Warning_1,Yes_Warning_2
integer c_node_1,c_node_2,c_node_3
real(kind=FT) c_Ape_1,c_Ape_2,c_Ape_3,c_Ave_Ape
real(kind=FT) offset_dis,x_new,y_new,z_new
integer Ele_Num_Cache
integer c_Max_N_FluEl_3D,c_Max_N_CalP_3D
character(200) c_File_name_1
logical alive

delta_L_Factor = 0.02D0
Cracks_Volume(1:Max_Num_Cr_3D) = ZR



!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,i_FluidEle,Flu_Ele_Area,  &
!$OMP                  num_Flu_Nodes,Flu_Ele_Centroid,                  &
!$OMP                  Flu_Ele_Vector,delta_L,up_offset_P,lw_offset_P,  &
!$OMP                  up_Elem_num,up_Kesi,up_Yita,up_Zeta,             &
!$OMP                  lw_Elem_num,lw_Kesi,lw_Yita,lw_Zeta,             &
!$OMP                  up_P_Disp,lw_P_Disp,Relative_Disp,c_Max_N_CalP_3D,  &
!$OMP                  norm_Flu_Ele_Vector,Ele_Num_Cache,c_Aperture,c_Max_N_FluEl_3D)    &         
!$OMP             SCHEDULE(static)     
do i_C = 1,num_Crack
  
  Ele_Num_Cache = 1
  c_Max_N_FluEl_3D = size(Cracks_FluidEle_CalP_3D(i_C)%row,1)
  Cracks_FluidEle_Aper_3D(i_C)%row(1:c_Max_N_FluEl_3D) = ZR
  
  c_Max_N_CalP_3D = size(Cracks_CalP_Coors_3D(i_C)%row,1)
  Cracks_CalP_Aper_3D(i_C)%row(1:c_Max_N_CalP_3D)       = ZR
  Cracks_CalP_LowDis_3D(i_C)%row(1:c_Max_N_CalP_3D,1:3) = ZR
  Cracks_CalP_LowDis_3D(i_C)%row(1:c_Max_N_CalP_3D,1:3) = ZR
  
  do i_FluidEle = 1,Cracks_FluidEle_num_3D(i_C)
      Flu_Ele_Area =Cracks_FluidEle_Area_3D(i_C)%row(i_FluidEle)
      num_Flu_Nodes=Cracks_FluidEle_num_CalP_3D(i_C)%row(i_FluidEle)
      Flu_Ele_Centroid(1:3) = Cracks_FluidEle_Centroid_3D(i_C)%row(i_FluidEle,1:3)
      Flu_Ele_Vector(1:3) = Cracks_FluidEle_Vector_3D(i_C)%row(i_FluidEle,1:3)
      
      if(Key_Crack_Aperture_Method==1) then
          call Cal_Crack_Point_Aperture_3D(c_DISP,i_C,Flu_Ele_Centroid,Relative_Disp,i_FluidEle,0,0,0)
          
      elseif(Key_Crack_Aperture_Method==2) then
          delta_L = delta_L_Factor*Ave_Elem_L_Enrich
          up_offset_P(1:3) = Flu_Ele_Centroid + delta_L*Flu_Ele_Vector(1:3)
          lw_offset_P(1:3)= Flu_Ele_Centroid - delta_L*Flu_Ele_Vector(1:3)

          call Cal_Ele_Num_by_Coors_3D(up_offset_P(1),up_offset_P(2),up_offset_P(3),Ele_Num_Cache,up_Elem_num)        
          call Cal_Ele_Num_by_Coors_3D(lw_offset_P(1),lw_offset_P(2),lw_offset_P(3),Ele_Num_Cache,lw_Elem_num)  
          
          call Cal_KesiYita_by_Coor_3D(up_offset_P,up_Elem_num,up_Kesi,up_Yita,up_Zeta)
          call Cal_KesiYita_by_Coor_3D(lw_offset_P,lw_Elem_num,lw_Kesi,lw_Yita,lw_Zeta)  
          if(up_Elem_num<=0)then
              print *,"    WARNING :: illegal up_Elem_num in Cal_Crack_Aperture_3D.f"
              Cracks_FluidEle_Aper_3D(i_C)%row(i_FluidEle) =ZR
              cycle
          endif
          if(lw_Elem_num<=0)then
              print *,"    WARNING :: illegal lw_Elem_num in Cal_Crack_Aperture_3D.f"
              Cracks_FluidEle_Aper_3D(i_C)%row(i_FluidEle) =ZR
              cycle
          endif
          call Cal_Any_Point_Disp_KesiYita_3D(up_Elem_num,up_Kesi,up_Yita,up_Zeta,c_DISP,up_P_Disp)
          call Cal_Any_Point_Disp_KesiYita_3D(lw_Elem_num,lw_Kesi,lw_Yita,lw_Zeta,c_DISP,lw_P_Disp)  
          Relative_Disp(1:3) = up_P_Disp- lw_P_Disp   
      endif
      
      norm_Flu_Ele_Vector = ONE
      c_Aperture = (Relative_Disp(1)* Flu_Ele_Vector(1) &
                  + Relative_Disp(2)* Flu_Ele_Vector(2) &
                  + Relative_Disp(3)* Flu_Ele_Vector(3))/norm_Flu_Ele_Vector
                  
      if(Key_Contact ==5 .or. Key_Non_Negtive_Aperture_3D  == 1) then
          if(c_Aperture<=ZR) c_Aperture=ZR
      endif
      
      Cracks_FluidEle_Aper_3D(i_C)%row(i_FluidEle) = c_Aperture    
  end do
end do  
!$omp end parallel do  

!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,i_FluidEle,c_Aperture,  &
!$OMP                 num_All_Flu_Nodes,count_All_Flu_Nodes, &
!$OMP                 Fl_El_Nodes,i_FluNode,num_Flu_Nodes)       &         
!$OMP            SCHEDULE(static)     
do i_C = 1,num_Crack
  
  num_All_Flu_Nodes = Cracks_CalP_Num_3D(i_C)
  Cracks_CalP_Aper_3D(i_C)%row(1:num_All_Flu_Nodes)  = ZR
  count_All_Flu_Nodes(1:num_All_Flu_Nodes) = 0
  do i_FluidEle = 1,Cracks_FluidEle_num_3D(i_C)
      c_Aperture =  Cracks_FluidEle_Aper_3D(i_C)%row(i_FluidEle)
      num_Flu_Nodes=Cracks_FluidEle_num_CalP_3D(i_C)%row(i_FluidEle)
      Fl_El_Nodes(1:num_Flu_Nodes) = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,1:num_Flu_Nodes)
      Cracks_CalP_Aper_3D(i_C)%row(Fl_El_Nodes(1:num_Flu_Nodes)) = &
              Cracks_CalP_Aper_3D(i_C)%row(Fl_El_Nodes(1:num_Flu_Nodes)) +c_Aperture
      count_All_Flu_Nodes(Fl_El_Nodes(1:num_Flu_Nodes)) = count_All_Flu_Nodes(Fl_El_Nodes(1:num_Flu_Nodes)) +1
  enddo
  
  
  do i_FluNode=1,num_All_Flu_Nodes
      Cracks_CalP_Aper_3D(i_C)%row(i_FluNode) = Cracks_CalP_Aper_3D(i_C)%row(i_FluNode)/count_All_Flu_Nodes(i_FluNode)
  enddo
  
end do  
!$omp end parallel do       

Cracks_Volume(1:num_Crack) = ZR
!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,i_FluidEle,Flu_Ele_Area,c_El_Ave_Aperture)SCHEDULE(static)       
do i_C = 1,num_Crack
  
  do i_FluidEle = 1,Cracks_FluidEle_num_3D(i_C)
    Flu_Ele_Area = Cracks_FluidEle_Area_3D(i_C)%row(i_FluidEle)
    c_El_Ave_Aperture = Cracks_FluidEle_Aper_3D(i_C)%row(i_FluidEle)
    Cracks_Volume(i_C) = Cracks_Volume(i_C)  +c_El_Ave_Aperture*Flu_Ele_Area
#ifndef Silverfrost
    if (isnan(Cracks_Volume(i_C))) then
        write(*,*) '    Warn-2022121601 :: Volume of crack ',i_C,'is NaN!'
        Cracks_Volume(i_C) = ZR
        exit
    endif
#endif    
  enddo
end do
!$omp end parallel do  

if (Key_Save_Nothing==1) goto 555

if (Key_Save_avol_file==0) goto 555
 
c_File_name_1   =  trim(Full_Pathname)//'.avol' 
if(isub==1)then
    inquire(file=c_File_name_1, exist=alive)
    if(alive.EQV..True.)then
          OPEN  (UNIT=105, FILE=c_File_name_1, STATUS='OLD') 
          CLOSE (UNIT=105, STATUS='DELETE')
    endif
endif
open(301,file=c_File_name_1,status='unknown',position='append',action='write') 
if(isub==1)then
  write(301,*) '    i_Step   |   All Volume (m3)  |   Only Positive (m3)'   
endif   
write(301, '(I10,F18.5,5X,F18.5)') isub,sum(Cracks_Volume(1:num_Crack)),&
                                        sum(Cracks_Volume(1:num_Crack),mask=(Cracks_Volume>ZR))
close(301)    
555 continue

!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,i_Node,ori_n,c_Node_Coor,&
!$OMP               delta_L,up_offset_P,lw_offset_P,up_Elem_num,       &
!$OMP               lw_Elem_num,c_Aperture,up_Kesi,up_Yita,up_Zeta,    &
!$OMP               lw_Kesi,lw_Yita,lw_Zeta,up_P_Disp,lw_P_Disp,       &
!$OMP               Relative_Disp,norm_Flu_Ele_Vector,Ele_Num_Cache)   &
!$OMP             SCHEDULE(static)   
do i_C = 1,num_Crack
  
  Ele_Num_Cache = 1
  do i_Node =1,Crack3D_Meshed_Node_num(i_C)      
      ori_n =Crack3D_Meshed_Node_Nor_Vector(i_C)%row(i_Node,1:3)
      c_Node_Coor = Crack3D_Meshed_Node(i_C)%row(i_Node,1:3)

      if(Key_Crack_Aperture_Method==1) then
          call Cal_Crack_Point_Aperture_3D(c_DISP,i_C,c_Node_Coor,Relative_Disp,0,0,0,i_Node)
      elseif(Key_Crack_Aperture_Method==2) then
          delta_L = delta_L_Factor*Ave_Elem_L_Enrich
          up_offset_P(1:3) = c_Node_Coor + delta_L*ori_n(1:3)
          lw_offset_P(1:3) = c_Node_Coor - delta_L*ori_n(1:3)
        
          call Cal_Ele_Num_by_Coors_3D(up_offset_P(1),up_offset_P(2),up_offset_P(3),Ele_Num_Cache,up_Elem_num)        
          call Cal_Ele_Num_by_Coors_3D(lw_offset_P(1),lw_offset_P(2),lw_offset_P(3),Ele_Num_Cache,lw_Elem_num)  

          if(up_Elem_num<=0)then
              c_Aperture = ZR
              Crack3D_Meshed_Node_Value(i_C)%row(i_Node,1)  = ZR
              cycle
          endif
          if(lw_Elem_num<=0)then
              c_Aperture = ZR
              Crack3D_Meshed_Node_Value(i_C)%row(i_Node,1)  = ZR 
              cycle
          endif
          call Cal_KesiYita_by_Coor_3D(up_offset_P,up_Elem_num,up_Kesi,up_Yita,up_Zeta)
          call Cal_KesiYita_by_Coor_3D(lw_offset_P,lw_Elem_num,lw_Kesi,lw_Yita,lw_Zeta)                
          call Cal_Any_Point_Disp_KesiYita_3D(up_Elem_num,up_Kesi,up_Yita,up_Zeta,c_DISP,up_P_Disp)
          call Cal_Any_Point_Disp_KesiYita_3D(lw_Elem_num,lw_Kesi,lw_Yita,lw_Zeta,c_DISP,lw_P_Disp)    

          Relative_Disp(1:3) = up_P_Disp- lw_P_Disp   

      endif
      
      
      
      norm_Flu_Ele_Vector = ONE
      c_Aperture=(Relative_Disp(1)*ori_n(1)+Relative_Disp(2)*ori_n(2)+Relative_Disp(3)*ori_n(3))/norm_Flu_Ele_Vector
      
      if(Key_Contact ==5 .or. Key_Non_Negtive_Aperture_3D  == 1) then
          if(c_Aperture<=ZR) c_Aperture=ZR
      endif
      
      Crack3D_Meshed_Node_Value(i_C)%row(i_Node,1)  = c_Aperture  
  enddo
end do   
!$omp end parallel do        

!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,i_FluidEle,c_node_1, &
!$OMP       c_node_2,c_node_3,c_Ape_1,c_Ape_2,c_Ape_3,c_Ave_Ape) SCHEDULE(static)        
do i_C = 1,num_Crack
  do i_FluidEle =1,Cracks_FluidEle_num_3D(i_C)
      c_node_1 = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,1)
      c_node_2 = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,2)
      c_node_3 = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,3)
      c_Ape_1  = Cracks_CalP_Aper_3D(i_C)%row(c_node_1) 
      c_Ape_2  = Cracks_CalP_Aper_3D(i_C)%row(c_node_2) 
      c_Ape_3  = Cracks_CalP_Aper_3D(i_C)%row(c_node_3) 
      c_Ave_Ape= (c_Ape_1+c_Ape_2+c_Ape_3)/THR
      Cracks_FluidEle_Aper_3D(i_C)%row(i_FluidEle) = c_Ave_Ape 
  enddo
enddo
!$omp end parallel do      

do i_C = 1,num_Crack
    Crack_Max_Min_Aperture(i_C,1) = maxval(Cracks_FluidEle_Aper_3D(i_C)%row(1:Cracks_FluidEle_num_3D(i_C))) 
    Crack_Max_Min_Aperture(i_C,2) = minval(Cracks_FluidEle_Aper_3D(i_C)%row(1:Cracks_FluidEle_num_3D(i_C))) 
    Crack_Max_Min_Aperture(i_C,3) = sum(Cracks_FluidEle_Aper_3D(i_C)%row(1:Cracks_FluidEle_num_3D(i_C)))/num_Crack
enddo
   
return 
end SUBROUTINE Cal_Crack_Aperture_3D               
