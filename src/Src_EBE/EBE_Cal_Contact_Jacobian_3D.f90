 
subroutine EBE_Cal_Contact_Jacobian_3D(isub,i_NR_P,num_freeD,freeDOF,size_local_0, &
       all_local_0,diag_precon_no_invert_0,diag_precon_no_invert,                   &
       Kn,Kn_Gauss,Kt1_Gauss,Kt2_Gauss,CT_State_Gauss)


use Global_Float_Type
use Global_Common
use Global_Model
use Global_XFEM_Elements
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Global_Material
use OMP_LIB

implicit none
integer,intent(in)::isub,i_NR_P
integer,intent(in)::num_FreeD
integer,intent(in)::freeDOF(num_FreeD)
real(kind=FT),intent(in)::Kn
real(kind=FT),intent(in)::diag_precon_no_invert_0(0:num_FreeD)
integer,intent(in)::size_local_0(Num_Elem)
integer,intent(in)::all_local_0(MDOF_3D,Num_Elem)      
real(kind=FT),intent(out)::diag_precon_no_invert(0:num_FreeD)
real(kind=FT),intent(in):: Kt1_Gauss(num_Crack,Max_Max_N_FluEl_3D),  &
      Kt2_Gauss(num_Crack,Max_Max_N_FluEl_3D),         &
      Kn_Gauss(num_Crack,Max_Max_N_FluEl_3D)
integer,intent(in)::CT_State_Gauss(num_Crack,Max_Max_N_FluEl_3D)
integer i_N,i_C
integer c_NN(8)
integer i_CT_Elem
integer Local_W(MDOF_3D),num_Local_W
real(kind=FT) tem_n(8),N_W(3,MDOF_3D),N(3,24)
real(kind=FT) D_ep(3,3),D_ep_o(3,3),T(3,3)
real(kind=FT) c_kn
real(kind=FT) c_Center(3),ori_n(3),Flu_Ele_Area,c_kt1,c_kt2
integer Solid_El
real(kind=FT) Kesi,Yita,Zeta,N_HF(3)
real(kind=FT)BaseLine_A(3),BaseLine_B(3)
real(kind=FT)BaseLine_x_Vec(3),BaseLine_y_Vec(3),BaseLine_z_Vec(3)
real(kind=FT) BaseLine_Mid(3),Gauss_Coor_Local(3)
integer ref_elem
real(kind=FT) Tip_T_Matrx(3,3),r_Gauss,theta_Gauss
integer k,c_loca,c_loca_W,local_FreeD(MDOF_3D),j,True_local_FreeD(MDOF_3D)
integer local(MDOF_3D)
integer num_Loc_ESM,c_Ele_Location,jj
real(kind=FT),POINTER::diag_precon_no_invert_thread(:,:)
integer i_Thread,c_thread,max_threads
integer c_Cr_Location
integer i_E, c_Elem,cEle_Loc,c_size
integer cnt
integer c_LocXFEM(3)
integer c_DOFs(3)
real(kind=FT) c_Beta(3) 
integer i_beta,j_beta
real(kind=FT) c_Penalty_CS
integer i_f
real(kind=FT) c_Aperture
real(kind=FT),ALLOCATABLE::localK_Temp(:,:)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_Elem,cEle_Loc,c_size) 
do i_E = 1,num_XFEM_Elem
    c_Elem = XFEM_Elem_List(i_E)
    cEle_Loc = Elem_Location(c_Elem,1) 
    c_size = size(storK_XFEM(cEle_Loc)%row(:,:),1)
    if(allocated(storK_XFEM_Updated(cEle_Loc)%row)) deallocate(storK_XFEM_Updated(cEle_Loc)%row)
    allocate(storK_XFEM_Updated(cEle_Loc)%row(c_size,c_size))
    storK_XFEM_Updated(cEle_Loc)%row(:,:) = storK_XFEM(cEle_Loc)%row(:,:)
enddo 
!$OMP END PARALLEL DO

diag_precon_no_invert(0:num_FreeD) = diag_precon_no_invert_0(0:num_FreeD)

max_threads = omp_get_max_threads()
if(Key_EBE_Precondition == 1)then
  ALLOCATE(diag_precon_no_invert_thread(0:num_FreeD,1:max_threads))
  diag_precon_no_invert_thread(0:num_FreeD,1:max_threads)= ZR 
endif


if(Key_Contact==7) then
    goto 100
endif
      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_C,i_CT_Elem,c_Center,ori_n,T,Solid_El,Flu_Ele_Area,c_kt1,c_kt2,c_kn,  &
!$OMP             D_ep_o,Kesi,Yita,Zeta,c_NN,num_Loc_ESM,N_W,i_N,Local_W,ref_elem,BaseLine_A,BaseLine_B,BaseLine_Mid,   &
!$OMP             BaseLine_x_Vec,BaseLine_y_Vec,BaseLine_z_Vec, Tip_T_Matrx,Gauss_Coor_Local,r_Gauss,theta_Gauss,       &
!$OMP             k,c_loca_W,local_FreeD,True_local_FreeD,localK_Temp,c_Ele_Location,num_Local_W,                            &
!$OMP             j,c_loca,local,N,tem_n,N_HF,D_ep,c_Cr_Location) 
c_thread = omp_get_thread_num()+1
!$OMP do SCHEDULE(static) 
do i_C=1,num_Crack
  do i_CT_Elem = 1,Cracks_FluidEle_num_3D(i_C)
      c_Center = Cracks_FluidEle_Centroid_3D(i_C)%row(i_CT_Elem,1:3)
      ori_n    = Cracks_FluidEle_Vector_3D(i_C)%row(i_CT_Elem,1:3)
      T=Cracks_FluidEle_LCS_T_3D(i_C)%row(i_CT_Elem,1:3,1:3)
      Solid_El = Cracks_FluidEle_EleNum_3D(i_C)%row(i_CT_Elem)     
      Flu_Ele_Area = Cracks_FluidEle_Area_3D(i_C)%row(i_CT_Elem)  
      
      
      
      if(CT_State_Gauss(i_C,i_CT_Elem) /= 0)then
          c_kt1 = ZR
          c_kt2 = ZR
          c_kn = Kn_Gauss(i_C,i_CT_Elem)
          D_ep_o(1:3,1:3) =  ZR
          D_ep_o(1,1:3) = [c_kt1,    ZR,       ZR]
          D_ep_o(2,1:3) = [ZR,    c_kt2,       ZR]
          D_ep_o(3,1:3) = [ZR,       ZR,     c_kn]
          
          call Cal_KesiYita_by_Coor_3D(c_Center,Solid_El,Kesi,Yita,Zeta)
          c_NN(1:8)    = G_NN(1:8,Solid_El)
          
          num_Loc_ESM = size_local_0(Solid_El)
          local(1:num_Loc_ESM)= all_local_0(1:num_Loc_ESM,Solid_El) 

          call Cal_N_3D(Kesi,Yita,Zeta,N)     
          tem_n(1) = N(1,1);    tem_n(2) = N(1,4)
          tem_n(3) = N(1,7);    tem_n(4) = N(1,10)      
          tem_n(5) = N(1,13);   tem_n(6) = N(1,16)
          tem_n(7) = N(1,19);   tem_n(8) = N(1,22)  
          
          N_HF(1:3) = ONE/THR      
          
          D_ep = MATMUL(MATMUL(transpose(T),D_ep_o),T)
          
          num_Local_W = 0
          N_W(1:3,1:MDOF_3D) =ZR
          do i_N = 1,8
              if(Enriched_Node_Type_3D(c_NN(i_N),i_C) ==2)then
                  num_Local_W = num_Local_W+1
                  Local_W(3*num_Local_W-2)=3*c_POS_3D(c_NN(i_N),i_C)-2
                  Local_W(3*num_Local_W-1)=3*c_POS_3D(c_NN(i_N),i_C)-1     
                  Local_W(3*num_Local_W)  =3*c_POS_3D(c_NN(i_N),i_C)
                  N_W(1,3*num_Local_W-2)   = TWO*tem_n(i_N)
                  N_W(2,3*num_Local_W-1)   = TWO*tem_n(i_N)
                  N_W(3,3*num_Local_W)     = TWO*tem_n(i_N)     
              elseif(Enriched_Node_Type_3D(c_NN(i_N),i_C) ==3)then

              elseif(Enriched_Node_Type_3D(c_NN(i_N),i_C)==1) then
                  if ((Elem_Type_3D(Solid_El,i_C).eq.1 )) then
                      ref_elem=Solid_El
                  else
                      ref_elem=Ele_Num_Tip_Enriched_Node_3D(c_NN(i_N))%row(i_C)
                  endif
  
                  call Vector_Location_Int_v2(Solid_El_Max_num_Crs,Solid_El_Crs(ref_elem,1:Solid_El_Max_num_Crs),i_C,c_Cr_Location)
     
                  BaseLine_A = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location,1,1:3)
                  BaseLine_B = Solid_El_Tip_BaseLine(ref_elem)%row(c_Cr_Location,2,1:3)
                  BaseLine_Mid = (BaseLine_A+BaseLine_B)/TWO
                  BaseLine_x_Vec = Solid_El_Tip_BaseLine_x_Vec(ref_elem)%row(c_Cr_Location,1:3)
                  BaseLine_y_Vec = Solid_El_Tip_BaseLine_y_Vec(ref_elem)%row(c_Cr_Location,1:3)
                  BaseLine_z_Vec = Solid_El_Tip_BaseLine_z_Vec(ref_elem)%row(c_Cr_Location,1:3)
                  Tip_T_Matrx(1:3,1:3) = Solid_El_Tip_BaseLine_T_Matrix(ref_elem)%row(c_Cr_Location,1:3,1:3)

                  Gauss_Coor_Local= MATMUL(Tip_T_Matrx,c_Center-BaseLine_Mid)
                  r_Gauss = sqrt(Gauss_Coor_Local(1)**2 + Gauss_Coor_Local(2)**2)
                  theta_Gauss = atan2(Gauss_Coor_Local(2),Gauss_Coor_Local(1))    
                  num_Local_W = num_Local_W+1
                  Local_W(3*num_Local_W-2)=3*c_POS_3D(c_NN(i_N),i_C)-2
                  Local_W(3*num_Local_W-1)=3*c_POS_3D(c_NN(i_N),i_C)-1     
                  Local_W(3*num_Local_W)  =3*c_POS_3D(c_NN(i_N),i_C) 
 
                  N_W(1,3*num_Local_W-2)= TWO*sqrt(r_Gauss)*tem_n(i_N)
                  N_W(2,3*num_Local_W-1)= TWO*sqrt(r_Gauss)*tem_n(i_N)
                  N_W(3,3*num_Local_W)  = TWO*sqrt(r_Gauss)*tem_n(i_N)     
                  
              endif
          enddo
          DO k=1,3*num_Local_W
              
              c_loca_W =  Location_FreeDOF(Local_W(k))
              
              local_FreeD(k) = c_loca_W
          enddo
          
          DO k=1,3*num_Local_W
              c_loca_W = minloc(local(1:num_Loc_ESM),1,MASK = (local(1:num_Loc_ESM).eq.local_FreeD(k)))
              True_local_FreeD(k) = c_loca_W
          enddo 
          
          allocate(localK_Temp(MDOF_3D,MDOF_3D))
          localK_Temp(1:MDOF_3D,1:MDOF_3D) = ZR
          localK_Temp(1:3*num_Local_W,1:3*num_Local_W) = &
              MATMUL(MATMUL(transpose(N_W(1:3,1:3*num_Local_W)),D_ep),N_W(1:3,1:3*num_Local_W))*Flu_Ele_Area 

          if (Elem_XFEM_Flag(Solid_El) ==1)then  
              c_Ele_Location = Elem_Location(Solid_El,1) 
              !$OMP CRITICAL
              storK_XFEM_Updated(c_Ele_Location)%row(True_local_FreeD(1:3*num_Local_W),True_local_FreeD(1:3*num_Local_W))=        &
                  storK_XFEM_Updated(c_Ele_Location)%row(True_local_FreeD(1:3*num_Local_W),True_local_FreeD(1:3*num_Local_W)) + & 
                    localK_Temp(1:3*num_Local_W,1:3*num_Local_W)
              !$OMP END CRITICAL
          elseif (Elem_XFEM_Flag(Solid_El) ==0)then
              print *,'Error-2023082701 :: Contact should not be related to FEM elements, in EBE_Cal_Contact_Jacobian_3D.f90.'
              call Warning_Message('S',Keywords_Blank)  
          endif

          if (Key_EBE_Precondition==1)then
              DO j=1,3*num_Local_W
                  c_loca=local_FreeD(j)        
                  diag_precon_no_invert_thread(c_loca,c_thread)=diag_precon_no_invert_thread(c_loca,c_thread) +localK_Temp(j,j)
              END DO 
          endif
          
          deallocate(localK_Temp)
      endif
  enddo
enddo
!$omp end do
!$omp end parallel     

100 continue

if(Key_Contact==7)then
    do i_C=1,num_Crack
    do i_CT_Elem = 1,Cracks_FluidEle_num_3D(i_C)
          c_Center = Cracks_FluidEle_Centroid_3D(i_C)%row(i_CT_Elem,1:3)
          ori_n    = Cracks_FluidEle_Vector_3D(i_C)%row(i_CT_Elem,1:3)
          T=Cracks_FluidEle_LCS_T_3D(i_C)%row(i_CT_Elem,1:3,1:3)
          Solid_El = Cracks_FluidEle_EleNum_3D(i_C)%row(i_CT_Elem)     
          Flu_Ele_Area = Cracks_FluidEle_Area_3D(i_C)%row(i_CT_Elem)  
          
          c_Aperture  = Cracks_FluidEle_Aper_3D(i_C)%row(i_CT_Elem)  
          
          
          c_Penalty_CS = Penalty_CS_Natural_Crack

          if(CT_State_Gauss(i_C,i_CT_Elem) /= 0)then
              c_kt1 = ZR
              c_kt2 = ZR
              c_kn = Kn_Gauss(i_C,i_CT_Elem)
              D_ep_o(1:3,1:3) =  ZR
              D_ep_o(1,1:3) = [c_kt1,    ZR,       ZR]
              D_ep_o(2,1:3) = [ZR,    c_kt2,       ZR]
              D_ep_o(3,1:3) = [ZR,       ZR,     c_kn]
              
              call Cal_KesiYita_by_Coor_3D(c_Center,Solid_El,Kesi,Yita,Zeta)
              c_NN(1:8)    = G_NN(1:8,Solid_El)
              
              num_Loc_ESM = size_local_0(Solid_El)
              local(1:num_Loc_ESM)= all_local_0(1:num_Loc_ESM,Solid_El) 

              call Cal_N_3D(Kesi,Yita,Zeta,N)     
              tem_n(1) = N(1,1);    tem_n(2) = N(1,4)
              tem_n(3) = N(1,7);    tem_n(4) = N(1,10)      
              tem_n(5) = N(1,13);   tem_n(6) = N(1,16)
              tem_n(7) = N(1,19);   tem_n(8) = N(1,22)    
              N_HF(1:3) = ONE/THR      
              D_ep = MATMUL(MATMUL(transpose(T),D_ep_o),T)
              
              num_Local_W = 0
              N_W(1:3,1:MDOF_3D) =ZR
              cnt = 0
              c_LocXFEM(1:3) = 0
              cEle_Loc = Elem_Location(Solid_El,1)      
              do i_N = 1,8
                  if(Enriched_Node_Type_3D(c_NN(i_N),i_C) ==2)then
                      num_Local_W = num_Local_W+1
                      Local_W(3*num_Local_W-2)=3*c_POS_3D(c_NN(i_N),i_C)-2
                      Local_W(3*num_Local_W-1)=3*c_POS_3D(c_NN(i_N),i_C)-1     
                      Local_W(3*num_Local_W)  =3*c_POS_3D(c_NN(i_N),i_C)
                      N_W(1,3*num_Local_W-2)   = TWO*tem_n(i_N)
                      N_W(2,3*num_Local_W-1)   = TWO*tem_n(i_N)
                      N_W(3,3*num_Local_W)     = TWO*tem_n(i_N)    
                      
                      cnt = cnt + 1
                      c_DOFs(1) = 3*c_POS_3D(c_NN(i_N),i_C) - 2
                      c_DOFs(2) = 3*c_POS_3D(c_NN(i_N),i_C) - 1
                      c_DOFs(3) = 3*c_POS_3D(c_NN(i_N),i_C)  
                    
                      c_LocXFEM(1) = 24 + 3*cnt -2
                      c_LocXFEM(2) = 24 + 3*cnt -1
                      c_LocXFEM(3) = 24 + 3*cnt                    
                      
                      c_Beta(1:3) = Enriched_Node_Crack_n_Vector_3D(c_NN(i_N))%row(i_C,1:3)
                      if(Key_Penalty_CS_Method==1) then
                          do j_beta=1,3
                              do i_beta=1,3       
                                  storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i_beta),c_LocXFEM(j_beta)) = &
                                    storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i_beta),c_LocXFEM(j_beta)) &
                                       + c_Penalty_CS*c_Beta(i_beta)*c_Beta(j_beta)
                              enddo
                          enddo
                      elseif(Key_Penalty_CS_Method==2) then
                          do j_beta=3,3
                              do i_beta=3,3           
                                  storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i_beta),c_LocXFEM(j_beta)) = &
                                  storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i_beta),c_LocXFEM(j_beta)) &
                                       + c_Penalty_CS*c_Beta(i_beta)*c_Beta(j_beta)
                              enddo
                          enddo
                      endif
                  elseif(Enriched_Node_Type_3D(c_NN(i_N),i_C) ==3)then

                  elseif(Enriched_Node_Type_3D(c_NN(i_N),i_C)==1) then
                      if ((Elem_Type_3D(Solid_El,i_C).eq.1)) then
                          ref_elem=Solid_El
                      else
                          ref_elem=Ele_Num_Tip_Enriched_Node_3D(c_NN(i_N))%row(i_C)
                      endif
                      
                      do i_f=1,Num_F_Functions
                        cnt = cnt + 1
                        c_DOFs(1) = 3*(c_POS_3D(c_NN(i_N),i_C)+i_f-1) - 2
                        c_DOFs(2) = 3*(c_POS_3D(c_NN(i_N),i_C)+i_f-1) - 1
                        c_DOFs(3) = 3*(c_POS_3D(c_NN(i_N),i_C)+i_f-1)   
                        
                        c_LocXFEM(1) = 24 + 3*cnt -2
                        c_LocXFEM(2) = 24 + 3*cnt -1
                        c_LocXFEM(3) = 24 + 3*cnt                    
                        
                          c_Beta(1:3) = Enriched_Node_Crack_n_Vector_3D(c_NN(i_N))%row(i_C,1:3)
                          if(Key_Penalty_CS_Method==1) then
                              do j_beta=1,3
                                  do i_beta=1,3       
                                      storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i_beta),c_LocXFEM(j_beta)) = &
                                        storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i_beta),c_LocXFEM(j_beta)) &
                                           + c_Penalty_CS*c_Beta(i_beta)*c_Beta(j_beta)
                                  enddo
                              enddo
                          elseif(Key_Penalty_CS_Method==2) then
                              do j_beta=3,3
                                  do i_beta=3,3           
                                      storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i_beta),c_LocXFEM(j_beta)) = &
                                      storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i_beta),c_LocXFEM(j_beta)) &
                                           + c_Penalty_CS*c_Beta(i_beta)*c_Beta(j_beta)
                                  enddo
                              enddo
                          endif             
                          
                      end do  
                  endif
              enddo
             

              if (Elem_XFEM_Flag(Solid_El) ==1)then  

              elseif (Elem_XFEM_Flag(Solid_El) ==0)then
                  print *,'Error-2023082701 :: Contact should not be related to FEM elements, in EBE_Cal_Contact_Jacobian_3D.f90.'
                  call Warning_Message('S',Keywords_Blank)  
              endif
          endif
    enddo
    enddo
endif

if(Key_EBE_Precondition == 1)then
  DO i_Thread = 1,omp_get_max_threads()
     diag_precon_no_invert(0:num_FreeD)  = diag_precon_no_invert(0:num_FreeD)  + diag_precon_no_invert_thread(0:num_FreeD,i_Thread)
  ENDDO      
  deALLOCATE(diag_precon_no_invert_thread)
  
  forall(jj=1:num_FreeD,abs(diag_precon_no_invert(jj))<=Tol_10)
      diag_precon_no_invert(jj)=Tol_10
  end forall             
endif


return 
end SUBROUTINE EBE_Cal_Contact_Jacobian_3D       
