 
subroutine EBE_XFEM_PCG_3D_Modify_K_CS(isub,num_freeD,freeDOF,size_local_0, &
       all_local_0,diag_precon_no_invert)

use Global_Float_Type
use Global_Common
use Global_Model
use Global_XFEM_Elements
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Global_Material
use Global_Contact
use OMP_LIB

implicit none
integer,intent(in)::isub
integer,intent(in)::num_FreeD
integer,intent(in)::freeDOF(num_FreeD)
integer,intent(in)::size_local_0(Num_Elem)
integer,intent(in)::all_local_0(MDOF_3D,Num_Elem)      
real(kind=FT),intent(out)::diag_precon_no_invert(0:num_FreeD)
integer i_N,i_C
integer c_NN(8)
integer k,c_loca,c_loca_W,j
real(kind=FT) localK(MDOF_3D,MDOF_3D)
integer local(MDOF_3D)
integer num_Loc_ESM,jj
real(kind=FT),POINTER::diag_precon_no_invert_thread(:,:)
integer i_Thread,c_thread,max_threads
integer i_E, c_Elem,cEle_Loc,c_size
integer Yes_DOF_Dealed(Total_FD)
real(kind=FT) c_Penalty_CS
real(kind=FT) c_Beta(3) 
integer c_LocXFEM(3),c_DOFs(3)
integer cnt,i_f
integer i
integer kk
real(kind=FT) Diag_localK(MDOF_3D)

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


max_threads = omp_get_max_threads()
if(Key_EBE_Precondition == 1)then
  ALLOCATE(diag_precon_no_invert_thread(0:num_FreeD,1:max_threads))
  diag_precon_no_invert_thread(0:num_FreeD,1:max_threads)= ZR 
endif
      
print *,"    PCG-EBE: get c_Penalty_CS for Penalty function treatment..."   
Yes_DOF_Dealed(1:Total_FD) = 0

c_Penalty_CS =  Penalty_CS_Natural_Crack


print *,"    PCG-EBE: Penalty function treatment for XFEM elements..."   


do i_C=1,num_Crack
  do i_E = 1,num_XFEM_Elem
      c_Elem = XFEM_Elem_List(i_E)
      if(Elem_Conta_Sta(c_Elem,i_C)==1) then
          cEle_Loc = Elem_Location(c_Elem,1)           
          
          if(Elem_Type_3D(c_Elem,i_C)==0) then
              cycle
          endif
      
          c_NN(1:8)    = G_NN(1:8,c_Elem)

          cnt = 0
          c_LocXFEM(1:3) = 0
          do i_N = 1,8
            if (Enriched_Node_Type_3D(c_NN(i_N),i_C) .eq. 2)then  
                cnt = cnt + 1
                c_DOFs(1) = 3*c_POS_3D(c_NN(i_N),i_C) - 2
                c_DOFs(2) = 3*c_POS_3D(c_NN(i_N),i_C) - 1
                c_DOFs(3) = 3*c_POS_3D(c_NN(i_N),i_C)  
                
                c_LocXFEM(1) = 24 + 3*cnt -2
                c_LocXFEM(2) = 24 + 3*cnt -1
                c_LocXFEM(3) = 24 + 3*cnt                    
                
                  Yes_DOF_Dealed(c_DOFs(1))=1
                  c_Beta(1:3) = Enriched_Node_Crack_n_Vector_3D(c_NN(i_N))%row(i_C,1:3)
                  if(Key_Penalty_CS_Method==1) then
                      do j=1,3
                          do i=1,3       
                              storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) = &
                                storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) &
                                   + c_Penalty_CS*c_Beta(i)*c_Beta(j)
                          enddo
                      enddo
                  elseif(Key_Penalty_CS_Method==2) then
                      do j=3,3
                          do i=3,3           
                              storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) = &
                              storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) &
                                   + c_Penalty_CS*c_Beta(i)*c_Beta(j)
                          enddo
                      enddo
                  endif
            elseif(Enriched_Node_Type_3D(c_NN(i_N),i_C).eq.1)then 
                do i_f=1,Num_F_Functions
                    cnt = cnt + 1
                    c_DOFs(1) = 3*(c_POS_3D(c_NN(i_N),i_C)+i_f-1) - 2
                    c_DOFs(2) = 3*(c_POS_3D(c_NN(i_N),i_C)+i_f-1) - 1
                    c_DOFs(3) = 3*(c_POS_3D(c_NN(i_N),i_C)+i_f-1)   
                    
                    c_LocXFEM(1) = 24 + 3*cnt -2
                    c_LocXFEM(2) = 24 + 3*cnt -1
                    c_LocXFEM(3) = 24 + 3*cnt                    
                    
                    if (Yes_DOF_Dealed(c_DOFs(1))==0) then
                      Yes_DOF_Dealed(c_DOFs(1))=1
                      c_Beta(1:3) = Enriched_Node_Crack_n_Vector_3D(c_NN(i_N))%row(i_C,1:3)
                      if(Key_Penalty_CS_Method==1) then
                          do j=1,3
                              do i=1,3       
                                  storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) = &
                                    storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) &
                                       + c_Penalty_CS*c_Beta(i)*c_Beta(j)
                              enddo
                          enddo
                      elseif(Key_Penalty_CS_Method==2) then
                          do j=3,3
                              do i=3,3           
                                  storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) = &
                                  storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) &
                                       + c_Penalty_CS*c_Beta(i)*c_Beta(j)
                              enddo
                          enddo
                      endif    
                    endif
                end do  
              
            elseif(Enriched_Node_Type_3D(c_NN(i_N),i_C).eq.3)then 
                cnt = cnt + 1
                c_DOFs(1) = 3*c_POS_3D(c_NN(i_N),i_C) - 2
                c_DOFs(2) = 3*c_POS_3D(c_NN(i_N),i_C) - 1
                c_DOFs(3) = 3*c_POS_3D(c_NN(i_N),i_C)        
                
                c_LocXFEM(1) = 24 + 3*cnt -2
                c_LocXFEM(2) = 24 + 3*cnt -1
                c_LocXFEM(3) = 24 + 3*cnt                    
                
                  Yes_DOF_Dealed(c_DOFs(1))=1
                  c_Beta(1:3) = Enriched_Node_Crack_n_Vector_3D(c_NN(i_N))%row(i_C,1:3)
                  if(Key_Penalty_CS_Method==1) then
                      do j=1,3
                          do i=1,3       
                              storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) = &
                                storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) &
                                   + c_Penalty_CS*c_Beta(i)*c_Beta(j)
                          enddo
                      enddo
                  elseif(Key_Penalty_CS_Method==2) then
                      do j=3,3
                          do i=3,3           
                              storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) = &
                              storK_XFEM_Updated(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) &
                                   + c_Penalty_CS*c_Beta(i)*c_Beta(j)
                          enddo
                      enddo
                  endif              
            else
                cycle
            end if
          end do
      endif
  enddo
enddo



if(Key_EBE_Precondition == 1)then
  
  diag_precon_no_invert(0:num_FreeD)= ZR      
  print *,"    PCG-EBE: geting preconditioner for FEM elements..."          
  diag_precon_no_invert_thread(0:num_FreeD,1:max_threads)= ZR 
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,c_Elem,i_E,num_Loc_ESM,local,localK,kk,Diag_localK) 
  c_thread = omp_get_thread_num()+1
  !$OMP do  
  do i_E=1,num_FEM_Elem
      c_Elem = FEM_Elem_List(i_E)
      num_Loc_ESM = size_local_0(c_Elem)
      local(1:num_Loc_ESM)=all_local_0(1:num_Loc_ESM,c_Elem) 
      
      if(Key_EBE_Sym_Storage_K==0)then
          DO kk=1,num_Loc_ESM     
              Diag_localK(kk) =storK_FEM(kk,kk,Elem_Location(c_Elem,2)) 
          enddo       
      elseif(Key_EBE_Sym_Storage_K==1)then
          DO kk=1,num_Loc_ESM     
              Diag_localK(kk) =storK_FEM_Sym((kk-1)*24 -(kk-1)*kk/2 +kk,Elem_Location(c_Elem,2)) 
          enddo   
      endif
      
      diag_precon_no_invert_thread(local(1:num_Loc_ESM),c_thread)  =   &
                 diag_precon_no_invert_thread(local(1:num_Loc_ESM),c_thread) + Diag_localK(1:num_Loc_ESM)       
  enddo
  !$omp end do
  !$omp end parallel    
  
  DO i_Thread = 1,omp_get_max_threads()
     diag_precon_no_invert(0:num_FreeD)  =  diag_precon_no_invert(0:num_FreeD)  + diag_precon_no_invert_thread(0:num_FreeD,i_Thread)
  ENDDO
  
  print *,"    PCG-EBE: geting preconditioner for XFEM elements..."          
  diag_precon_no_invert_thread(0:num_FreeD,1:max_threads)= ZR 
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,c_Elem,i_E,num_Loc_ESM,local,localK,kk,Diag_localK) 
  c_thread = omp_get_thread_num()+1
  !$OMP do
  do i_E=1,num_XFEM_Elem
      c_Elem = XFEM_Elem_List(i_E)
      num_Loc_ESM = size_local_0(c_Elem)
      local(1:num_Loc_ESM)=all_local_0(1:num_Loc_ESM,c_Elem) 
      DO kk=1,num_Loc_ESM     
          Diag_localK(kk) =storK_XFEM_Updated(Elem_Location(c_Elem,1))%row(kk,kk) 
      enddo         
      diag_precon_no_invert_thread(local(1:num_Loc_ESM),c_thread)  =   &
                 diag_precon_no_invert_thread(local(1:num_Loc_ESM),c_thread) + Diag_localK(1:num_Loc_ESM)       
  enddo
  !$omp end do
  !$omp end parallel    
  
  DO i_Thread = 1,omp_get_max_threads()
     diag_precon_no_invert(0:num_FreeD)  =  diag_precon_no_invert(0:num_FreeD)  + diag_precon_no_invert_thread(0:num_FreeD,i_Thread)
  ENDDO  
  
  forall(jj=1:num_FreeD,abs(diag_precon_no_invert(jj))<=Tol_10)
      diag_precon_no_invert(jj)=Tol_10
  end forall
endif

return 
end SUBROUTINE EBE_XFEM_PCG_3D_Modify_K_CS     
