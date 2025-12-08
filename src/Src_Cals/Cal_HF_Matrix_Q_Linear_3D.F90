 
subroutine Cal_HF_Matrix_Q_Linear_3D(Counter,c_Total_FD,c_num_Tol_CalP_Water)

use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D

implicit none
integer,intent(in)::Counter,c_Total_FD,c_num_Tol_CalP_Water
integer i_N,i_C,i_HF_Elem,i_HF_GAUSS    
real(kind=FT) Kesi,Yita,Zeta
real(kind=FT) Orient1,Orient2,x1,y1,x2,y2
integer c_NN(8),Num_Div_Points 
real(kind=FT) N(3,24),tem_n(8)


real(kind=FT) HF_GAUSS_x,HF_GAUSS_y,HF_GAUSS_z
integer HF_GAUSS_Elem

integer,ALLOCATABLE::  Local_P(:)
real(kind=FT),ALLOCATABLE::   N_HF(:)
real(kind=FT),ALLOCATABLE::  tem_Q(:,:),tem(:,:)


integer Local_W(60),num_Local_W                           
real(kind=FT) N_W(3,60)
real(kind=FT) ori_n(3)
real(kind=FT) c_Sign,ori_HF_elem,ele_Sign
integer Num_Enriched_Node,j_E,c_Adj_Ele,c_Junc_Ele
integer i_FluidEle,i_CalP
integer c_CalP,num_Flu_Nodes
real(kind=FT) Flu_Ele_Area
real(kind=FT) Flu_Ele_Centroid(3)
real(kind=FT) c_Center(3)
real(kind=FT) BaseLine_Mid(3)
real(kind=FT) Gauss_Coor_Local(3)
integer ref_elem
real(kind=FT) Tip_T_Matrx(3,3),r_Gauss,theta_Gauss,T(3,3)   
integer Solid_El
real(kind=FT)BaseLine_A(3),BaseLine_B(3)
real(kind=FT)BaseLine_x_Vec(3),BaseLine_y_Vec(3),BaseLine_z_Vec(3)    
integer c_Cr_Location
integer i_Dof
integer num_HF_Cracks,HF_Cracks_List(num_Crack)
integer c_Crack

#ifndef Silverfrost
print *,'    Calculating coupled matrix Q......'

HF_Cracks_List = 0
num_HF_Cracks  = 0
do i_C=1,num_Crack
    if (Crack_Type_Status_3D(i_C,1) == 1) then  
        num_HF_Cracks = num_HF_Cracks +1
        HF_Cracks_List(num_HF_Cracks) = i_C
    endif
enddo



!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,c_Crack,i_FluidEle,tem,tem_Q,num_Flu_Nodes, & 
!$OMP                  Local_P,Flu_Ele_Centroid,c_Center,ori_n,T,Solid_El,Flu_Ele_Area,   & 
!$OMP                  HF_GAUSS_Elem,Kesi,Yita,Zeta,c_NN,N,tem_n,N_HF,num_Local_W,N_W,    & 
!$OMP                  i_N,Local_W,ref_elem,c_Cr_Location,BaseLine_A,BaseLine_B,          &
!$OMP                  BaseLine_Mid,BaseLine_x_Vec,BaseLine_y_Vec,BaseLine_z_Vec,         &
!$OMP                  Tip_T_Matrx,Gauss_Coor_Local,r_Gauss,theta_Gauss,i_Dof)            &          
!$OMP            SCHEDULE(static)  

do i_C=1,num_HF_Cracks  
    c_Crack = HF_Cracks_List(i_C)

    
    ALLOCATE(Local_P(c_num_Tol_CalP_Water))  
    ALLOCATE(N_HF(c_num_Tol_CalP_Water))
    ALLOCATE(tem(60,c_num_Tol_CalP_Water)) 
    Local_P(1:c_num_Tol_CalP_Water) =0
    N_HF(1:c_num_Tol_CalP_Water) = ZR
    tem(1:60,1:c_num_Tol_CalP_Water)  = ZR
    
    do i_FluidEle = 1,Cracks_FluidEle_num_3D(c_Crack)
        tem(1:60,1:c_num_Tol_CalP_Water)  = ZR
        num_Flu_Nodes = Cracks_FluidEle_num_CalP_3D(c_Crack)%row(i_FluidEle)
        
        Local_P(1:num_Flu_Nodes) = Cracks_FluidEle_Glo_CalP_3D(c_Crack)%row(i_FluidEle,1:num_Flu_Nodes) 


        Flu_Ele_Centroid(1:3) = Cracks_FluidEle_Centroid_3D(c_Crack)%row(i_FluidEle,1:3)    
        c_Center = Cracks_FluidEle_Centroid_3D(c_Crack)%row(i_FluidEle,1:3)
        ori_n    = Cracks_FluidEle_Vector_3D(c_Crack)%row(i_FluidEle,1:3)
        T = Cracks_FluidEle_LCS_T_3D(c_Crack)%row(i_FluidEle,1:3,1:3)
        Solid_El = Cracks_FluidEle_EleNum_3D(c_Crack)%row(i_FluidEle)     
        Flu_Ele_Area = Cracks_FluidEle_Area_3D(c_Crack)%row(i_FluidEle)    

        
        
        if (Flu_Ele_Area<=0.05D0*Ave_Elem_Area_Enrich) then
            cycle
        endif
        
        HF_GAUSS_Elem = Cracks_FluidEle_EleNum_3D(c_Crack)%row(i_FluidEle)
        call Cal_KesiYita_by_Coor_3D(Flu_Ele_Centroid,HF_GAUSS_Elem,Kesi,Yita,Zeta)
        c_NN(1:8)    = G_NN(1:8,HF_GAUSS_Elem)
        call Cal_N_3D(Kesi,Yita,Zeta,N)     
        tem_n(1) = N(1,1);    tem_n(2) = N(1,4)
        tem_n(3) = N(1,7);    tem_n(4) = N(1,10)      
        tem_n(5) = N(1,13);   tem_n(6) = N(1,16)
        tem_n(7) = N(1,19);   tem_n(8) = N(1,22)    
        N_HF(1:num_Flu_Nodes) = ONE/num_Flu_Nodes   
        num_Local_W = 0
        N_W(1:3,1:60) =ZR
        
        do i_N = 1,8
          if (Enriched_Node_Type_3D(c_NN(i_N),c_Crack) ==2) then
              num_Local_W = num_Local_W+1
              Local_W(3*num_Local_W-2)=3*c_POS_3D(c_NN(i_N),c_Crack)-2
              Local_W(3*num_Local_W-1)=3*c_POS_3D(c_NN(i_N),c_Crack)-1     
              Local_W(3*num_Local_W)  =3*c_POS_3D(c_NN(i_N),c_Crack)
              N_W(1,3*num_Local_W-2)   = TWO*tem_n(i_N)
              N_W(2,3*num_Local_W-1)   = TWO*tem_n(i_N)
              N_W(3,3*num_Local_W)     = TWO*tem_n(i_N)           
          elseif(Enriched_Node_Type_3D(c_NN(i_N),c_Crack) ==3)then
              num_Local_W = num_Local_W+1
              Local_W(3*num_Local_W-2)=3*c_POS_3D(c_NN(i_N),c_Crack)-2
              Local_W(3*num_Local_W-1)=3*c_POS_3D(c_NN(i_N),c_Crack)-1     
              Local_W(3*num_Local_W)  =3*c_POS_3D(c_NN(i_N),c_Crack)
              N_W(1,3*num_Local_W-2)   = TWO*tem_n(i_N)*1.0D0
              N_W(2,3*num_Local_W-1)   = TWO*tem_n(i_N)*1.0D0
              N_W(3,3*num_Local_W)     = TWO*tem_n(i_N)*1.0D0 
          elseif(Enriched_Node_Type_3D(c_NN(i_N),c_Crack)==1) then
              if ((Elem_Type_3D(Solid_El,c_Crack).eq.1 )) then
                  ref_elem=Solid_El
              else
                  ref_elem=Ele_Num_Tip_Enriched_Node_3D(c_NN(i_N))%row(c_Crack)
              endif
              call Vector_Location_Int_v2(Solid_El_Max_num_Crs,Solid_El_Crs(ref_elem,1:Solid_El_Max_num_Crs), &
                      c_Crack,c_Cr_Location)
              BaseLine_A = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location,1,1:3)
              BaseLine_B = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location,2,1:3)
              BaseLine_Mid = (BaseLine_A+BaseLine_B)/TWO
              BaseLine_x_Vec=Solid_El_Tip_BaseLine_x_Vec(ref_elem)%row(c_Cr_Location,1:3)
              BaseLine_y_Vec=Solid_El_Tip_BaseLine_y_Vec(ref_elem)%row(c_Cr_Location,1:3)
              BaseLine_z_Vec=Solid_El_Tip_BaseLine_z_Vec(ref_elem)%row(c_Cr_Location,1:3)
              Tip_T_Matrx(1:3,1:3)=Solid_El_Tip_BaseLine_T_Matrix(ref_elem)%row(c_Cr_Location,1:3,1:3)

              Gauss_Coor_Local= MATMUL(Tip_T_Matrx,c_Center-BaseLine_Mid)
              r_Gauss = sqrt(Gauss_Coor_Local(1)**2 + Gauss_Coor_Local(2)**2)
              theta_Gauss = atan2(Gauss_Coor_Local(2),Gauss_Coor_Local(1))    

              num_Local_W = num_Local_W+1
              Local_W(3*num_Local_W-2)=3*c_POS_3D(c_NN(i_N),c_Crack)-2
              Local_W(3*num_Local_W-1)=3*c_POS_3D(c_NN(i_N),c_Crack)-1     
              Local_W(3*num_Local_W)  =3*c_POS_3D(c_NN(i_N),c_Crack) 

              N_W(1,3*num_Local_W-2)= TWO*sqrt(r_Gauss)*tem_n(i_N)
              N_W(2,3*num_Local_W-1)= TWO*sqrt(r_Gauss)*tem_n(i_N)
              N_W(3,3*num_Local_W)  = TWO*sqrt(r_Gauss)*tem_n(i_N)     
          endif
        end do
        
        
        call Vectors_Multi(MATMUL(transpose(N_W(1:3,1:3*num_Local_W)),ori_n), &
                3*num_Local_W,N_HF(1:num_Flu_Nodes),num_Flu_Nodes,            &
                tem(1:3*num_Local_W,1:num_Flu_Nodes)) 


        

        
        do i_Dof = 1,3*num_Local_W
          if(.not. allocated(Coupled_Q_3D(Local_W(i_Dof))%row))then  
            allocate(Coupled_Q_3D(Local_W(i_Dof))%row(3))
            Coupled_Q_3D(Local_W(i_Dof))%row(1:3) = ZR
            
            allocate(Coupled_Q_3D_Index(Local_W(i_Dof))%row(3))
            Coupled_Q_3D_Index(Local_W(i_Dof))%row(1:3) = 0
            



            
            
            Coupled_Q_3D(Local_W(i_Dof))%row(1:3) = tem(i_Dof,1:3)*Flu_Ele_Area 
            Coupled_Q_3D_Index(Local_W(i_Dof))%row(1:3) = Local_P(1:3)
          else
            
            Coupled_Q_3D(Local_W(i_Dof))%row(1:3) = &
                         Coupled_Q_3D(Local_W(i_Dof))%row(1:3) +&
                         tem(i_Dof,1:3)*Flu_Ele_Area 
          endif
        enddo
              
              
    enddo
    DEALLOCATE(Local_P)  
    DEALLOCATE(N_HF)
    DEALLOCATE(tem) 
end do
!$omp end parallel do   

#endif

return 
end SUBROUTINE Cal_HF_Matrix_Q_Linear_3D           
