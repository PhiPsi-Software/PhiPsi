 
SUBROUTINE EBE_Determine_Contact_State_by_Iteration_3D(isub,c_cg_tol,max_num_PCG,num_FreeD,freeDOF, &
               globalF,DISP,diag_precon_no_invert,Num_Indent)
                                  
use Global_Float_Type
use Global_Common
use Global_Filename
use Global_Model
use Global_XFEM_Elements
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Global_Contact
use Global_Inter_Cal_Contact_Red_Resid
use omp_lib
use Global_POST
      
use Global_EBE_Cal_Contact_Jacobian_3D
use Global_Cal_Contact_Resid_3D   
use Global_Cal_Contact_Contact_State_Gauss_3D
use Global_Cal_Contact_PN_and_PT_3D
use Global_Cal_Contact_Conve_Factor

implicit none
integer,intent(in)::isub,max_num_PCG
integer,intent(in)::num_FreeD
integer,intent(in)::freeDOF(num_FreeD)
real(kind=FT),intent(inout)::DISP(Total_FD)
real(kind=FT),intent(in)::diag_precon_no_invert(0:num_FreeD)
real(kind=FT),intent(in)::c_cg_tol
real(kind=FT),intent(in)::globalF(Total_FD)
integer Num_Indent
real(kind=FT)::diag_precon_no_invert_Updated(0:num_FreeD)
real(kind=FT)::diag_precon_no_invert_0(0:num_FreeD)
real(kind=FT)::diag_precon(0:num_FreeD)
integer:: size_local_0(Num_Elem)
integer:: all_local_0(MDOF_3D,Num_Elem)
      
integer i_C,i_E,i_N,c_Edge_Elem
logical Yes_Contact
integer i_NR_P
real(kind=FT) R_PSI(Total_FD)
real(kind=FT) Last_R_PSI(Total_FD)
real(kind=FT) NR_DISP(Total_FD)
real(kind=FT) Conve_Tolerance,Conve_Factor
real(kind=FT) fric_mu
logical Yes_Conve
real(kind=FT) Contact_DISP_0(Total_FD)
real(kind=FT) B(3,80),tem_B(3,80),detJ,tem_Val
integer num_Loc_ESM,num_Loc_ESM_C_Crack      
real(kind=FT),ALLOCATABLE::p(:),loads(:),x(:),xnew(:),u(:), pcg_d(:)

real(kind=FT) Saved_Conv_Factor(Max_Contact_Iter)
logical Yes_Oscill
integer c_Node,cg_iters,i_PCG,i_Node,i,c_DOF
real(kind=FT) localK(MDOF_3D,MDOF_3D)
integer local(MDOF_3D)
real(kind=FT) alpha,beta,up,Tol
real(kind=FT) delta_disp(Total_FD)
real(kind=FT) PC_Gauss_x(num_Crack,Max_Max_N_FluEl_3D)
real(kind=FT) PC_Gauss_y(num_Crack,Max_Max_N_FluEl_3D)
real(kind=FT) PC_Gauss_z(num_Crack,Max_Max_N_FluEl_3D)
integer CT_State_Gauss(num_Crack,Max_Max_N_FluEl_3D)
real(kind=FT) Kn,Kt
real(kind=FT) Kn_Gauss(num_Crack,Max_Max_N_FluEl_3D)
real(kind=FT) Kt1_Gauss(num_Crack,Max_Max_N_FluEl_3D)
real(kind=FT) Kt2_Gauss(num_Crack,Max_Max_N_FluEl_3D)
character(200) c_File_name_1,c_File_name_2
character(200) c_File_name_3
character(5) temp
integer j
integer c_Ele_Location
real(kind=FT),allocatable::u_thread(:,:)
integer i_Thread,c_Thread,max_threads
integer c_Elem    
integer c_i,c_j
  
4031 FORMAT(8X,'---------------------------------------') 
4032 FORMAT(8X,'Contact iteration done after ',I3,' tries!') 
4033 FORMAT(8X,'Number of contact elements is ',I5) 
4132 FORMAT(5X,'Contact iteration done after ',I3,' tries!') 
4133 FORMAT(5X,'Number of contact elements is ',I5) 
4131 FORMAT(8X,'---------------------------------------') 
2003 FORMAT(8X,'  WARNING :: Oscillation detected!') 
2013 FORMAT(8X,'  WARNING :: Oscillation detected!') 
4001 FORMAT(8X,'>>  Contact NR iteration ',I3,' of ',I3,' started:') 
4002 FORMAT(5X,'>>  Contact NR iteration ',I3,' of ',I3,' started:')  
4022 FORMAT(12X,'Convergence factor is', E12.4 ,' | ',E12.4,'.')   
4023 FORMAT(9X,'Convergence factor is', E12.4 ,' | ',E12.4,'.')   
5001 FORMAT(12X,'Sticking solid elements:  ',I5,'| Sliding solid elements:  ',I5) 
5002 FORMAT(12X,'Sticking contact elements:',I5,'| Sliding contact elements:',I5) 
5101 FORMAT(9X,'Sticking solid elements:  ',I5,'| Sliding solid elements:  ',I5) 
5102 FORMAT(9X,'Sticking contact elements:',I5,'| Sliding contact elements:',I5) 
101 FORMAT(12X,'PCG-EBE: number of iterations to convergence is ',I5)  
102 FORMAT(12X,'PCG-EBE: convergence factor is ',E12.5,' / ',E12.5) 
103 FORMAT(17X,'CG_Iters:',I5, ' -> CG_Factor:', E12.5,' / ',E12.5)   
111 FORMAT(9X,'PCG-EBE: number of iterations to convergence is ',I5)  
112 FORMAT(9X,'PCG-EBE: convergence factor is ',E12.5,' / ',E12.5) 
113 FORMAT(14X,'CG_Iters:',I5, ' -> CG_Factor:', E12.5,' / ',E12.5)   
      
      
Contact_DISP_0(1:Total_FD) = ZR
Contact_DISP_0 = DISP
Saved_Conv_Factor(1:Max_Contact_Iter)  = ZR
diag_precon_no_invert_0(0:num_FreeD) =diag_precon_no_invert(0:num_FreeD)
size_local_0  = Size_Local_3D
all_local_0 =All_Local_3D

IF(ALLOCATED(Elem_Conta_Sta)) DEALLOCATE(Elem_Conta_Sta)  
ALLOCATE(Elem_Conta_Sta(Num_Elem,num_Crack))
ALLOCATE(p(0:num_FreeD),loads(0:num_FreeD),x(0:num_FreeD),xnew(0:num_FreeD),u(0:num_FreeD), pcg_d(0:num_FreeD))   
if(.not. allocated(storK_XFEM_Updated)) allocate(storK_XFEM_Updated(num_XFEM_Elem))   

max_threads = omp_get_max_threads()
if(.not. allocated(u_thread)) ALLOCATE(u_thread(0:num_FreeD,max_threads))
     
kn = kn_Cont_Penalty
kt = kt_Cont_Penalty
fric_mu = fric_mu_Cont
Conve_Tolerance = Conve_Tol_Penalty

select case (Key_Contact)
  case(1)
      do i_NR_P = 1,Max_Contact_Iter
          if (Num_Indent==8) then
              write(*,4001) i_NR_P,Max_Contact_Iter
          elseif(Num_Indent==5)then
              write(*,4002) i_NR_P,Max_Contact_Iter
          endif
          if(i_NR_P==1)then 
              NR_DISP(1:Total_FD)    = ZR
              delta_disp(1:Total_FD) = ZR                
              Last_R_PSI(1:Total_FD) = ZR
              Elem_Conta_Sta(1:Num_Elem,1:num_Crack)  = 0  
              Kn_Gauss(1:num_Crack,1:Max_Max_N_FluEl_3D)  = kn
              Kt1_Gauss(1:num_Crack,1:Max_Max_N_FluEl_3D) = kt
              Kt2_Gauss(1:num_Crack,1:Max_Max_N_FluEl_3D) = kt
        
              call Cal_Contact_Contact_State_Gauss_3D(isub,1,1,i_NR_P,Contact_DISP_0,Yes_Contact,Elem_Conta_Sta,   &
                    CT_State_Gauss)   
              if(Yes_Contact .eqv. .False.)then
                  if (Num_Indent==8) then
                      print *, '           No penetration detected, leaving contact iteration!'
                  elseif(Num_Indent==5)then
                      print *, '        No penetration detected, leaving contact iteration!'
                  endif
                  goto 9999
              endif
              
          endif              
          
          if (Num_Indent==8) then
              print *,'           Get contact force and update Kt.'
          elseif(Num_Indent==5)then
              print *,'        Get contact force and update Kt.'
          endif
          call Cal_Contact_PN_and_PT_3D(isub,1,1,i_NR_P,Total_FD,num_freeD,           &
                    Kn,Kn_Gauss,Kt1_Gauss,Kt2_Gauss,fric_mu,NR_DISP,delta_disp,       &
                    CT_State_Gauss,Elem_Conta_Sta,PC_Gauss_x,PC_Gauss_y,PC_Gauss_z,Num_Indent)
                    
          if (Num_Indent==8) then
              write(*,5001) count(Elem_Conta_Sta==1),count(Elem_Conta_Sta==2)
              write(*,5002) count(CT_State_Gauss==1),count(CT_State_Gauss==2)
          elseif(Num_Indent==5)then
              write(*,5101) count(Elem_Conta_Sta==1),count(Elem_Conta_Sta==2)
              write(*,5102) count(CT_State_Gauss==1),count(CT_State_Gauss==2)
          endif   

          
          if (Num_Indent==8) then
              print *,'           Assemble Jacobian matrix.'
          elseif(Num_Indent==5)then
              print *,'        Assemble Jacobian matrix.'
          endif
          
          call EBE_Cal_Contact_Jacobian_3D(isub,i_NR_P,num_freeD,freeDOF,size_local_0,              &
                    all_local_0,diag_precon_no_invert_0,diag_precon_no_invert_Updated,              &
                    Kn,Kn_Gauss,Kt1_Gauss,Kt2_Gauss,CT_State_Gauss)  
          
          
          if (Num_Indent==8) then
              print *,'           Get residual.'
          elseif(Num_Indent==5)then
              print *,'        Get residual.'
          endif

          call Cal_Contact_Resid_3D(isub,1,1,i_NR_P,Total_FD,num_freeD,globalF,     &
                           NR_DISP,freeDOF,PC_Gauss_x,PC_Gauss_y,PC_Gauss_z,R_PSI)     

          if (Num_Indent==8) then
              print *, "           Get the preconditioner and get starting loads..."
          elseif(Num_Indent==5)then
              print *, "        Get the preconditioner and get starting loads..."
          endif
          p    = ZR
          p(0) = ZR
          loads= ZR 
          loads(1:num_FreeD) = -R_PSI(freeDOF(1:num_FreeD))
          loads(0)= ZR 
              
          if (Key_EBE_Precondition==0) then
              diag_precon=ONE
              diag_precon(0) =ZR 
          elseif(Key_EBE_Precondition==1) then
              diag_precon(1:)=ONE/diag_precon_no_invert_Updated(1:) 
              diag_precon(0) =ZR 
          endif
          
          pcg_d=diag_precon*loads 
          p    =pcg_d
          x    =ZR
          cg_iters=0
          
          if (Num_Indent==8) then
              print *, '           PCG equation solution...'
          elseif(Num_Indent==5)then
              print *, '        PCG equation solution...'
          endif
          do i_PCG =1,max_num_PCG
              cg_iters=cg_iters+1 
              u=ZR     
              
              
              select case(Key_EBE_Sym_Storage_K)
              
              case(0)
                  u_thread(0:num_FreeD,1:max_threads) = ZR   
                  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,c_Ele_Location)
                  c_thread = omp_get_thread_num()+1
                  !$OMP DO SCHEDULE(static)  
                  do i_E=1,num_FEM_Elem
                      c_Elem = FEM_Elem_List(i_E)
                      num_Loc_ESM = Size_Local_3D(c_Elem)
                      local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
                      c_Ele_Location = Elem_Location(c_Elem,2) 
                      u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+ &
                                                   MATMUL(storK_FEM(1:num_Loc_ESM,1:num_Loc_ESM,c_Ele_Location),&
                                                          p(local(1:num_Loc_ESM)))
                  enddo
                  !$omp end do 
                  !$omp end parallel
              
              case(1)
                  u_thread(0:num_FreeD,1:max_threads) = ZR   
                  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,localK,c_Ele_Location,c_i,c_j)
                  c_thread = omp_get_thread_num()+1
                  !$OMP DO SCHEDULE(static)   
                  do i_E=1,num_FEM_Elem
                      c_Elem = FEM_Elem_List(i_E)
                      num_Loc_ESM = Size_Local_3D(c_Elem)
                      local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
                      c_Ele_Location = Elem_Location(c_Elem,2) 
                      
                      do c_i=1,num_Loc_ESM
                          do c_j=c_i,num_Loc_ESM
                              localK(c_i,c_j) = storK_FEM_Sym((c_i-1)*24 -(c_i-1)*c_i/2 +c_j,c_Ele_Location)  
                          enddo
                      enddo
                      do c_i=1,num_Loc_ESM
                          do c_j=1,c_i-1
                              localK(c_i,c_j) =  localK(c_j,c_i) 
                          enddo
                      enddo
                      u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+ &
                                                   MATMUL(localK(1:num_Loc_ESM,1:num_Loc_ESM),p(local(1:num_Loc_ESM))) 
                  enddo
                  !$omp end do 
                  !$omp end parallel
              end select
              
              DO i_Thread = 1,omp_get_max_threads()
                 u(0:num_FreeD)  =  u(0:num_FreeD)  + u_thread(0:num_FreeD,i_Thread)
              ENDDO  

              
              u_thread(0:num_FreeD,1:max_threads) = ZR   
              !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,c_Ele_Location)
              c_thread = omp_get_thread_num()+1
              !$OMP DO SCHEDULE(static) 
              do i_E=1,num_XFEM_Elem
                  c_Elem = XFEM_Elem_List(i_E)
                  num_Loc_ESM = Size_Local_3D(c_Elem)
                  local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
                  c_Ele_Location = Elem_Location(c_Elem,1)               
                  u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+ &
                                               MATMUL(storK_XFEM_Updated(c_Ele_Location)%row(1:num_Loc_ESM,1:num_Loc_ESM),&
                                                      p(local(1:num_Loc_ESM)))
              enddo
              !$omp end do 
              !$omp end parallel
              DO i_Thread = 1,omp_get_max_threads()
                 u(0:num_FreeD)  =  u(0:num_FreeD)  + u_thread(0:num_FreeD,i_Thread)
              ENDDO  
              
              up=DOT_PRODUCT(loads,pcg_d)
            

              
              alpha=up/DOT_PRODUCT(p,u) 

              
              xnew=x+p*alpha 
              loads=loads-u*alpha
              pcg_d=diag_precon*loads 
              
              beta=DOT_PRODUCT(loads,pcg_d)/up 

              
              p=pcg_d+p*beta

              Tol = MAXVAL(ABS(x(0:num_FreeD)-xnew(0:num_FreeD)))/MAXVAL(ABS(x(0:num_FreeD)))
              
              x=xnew
              
              if(Tol<c_cg_tol)then  
                  exit
              endif
              
              if(cg_iters==max_num_PCG)then
                  if (Num_Indent==8) then
                      print *, '           -------------------------------------'
                      print *, '               Warning :: PCG-EBE failed!       '  
                      print *, '           -------------------------------------'
                  elseif(Num_Indent==5)then
                      print *, '        -------------------------------------'
                      print *, '            Warning :: PCG-EBE failed!       '  
                      print *, '        -------------------------------------'
                  endif
                  exit
              endif 
      
          enddo
          
 
          if (Num_Indent==8) then
              write(*,101)  cg_iters      
              write(*,102)  Tol,c_cg_tol     
          elseif(Num_Indent==5)then
              write(*,111)  cg_iters      
              write(*,112)  Tol,c_cg_tol     
          endif
                  
          loads=xnew
          delta_disp = ZR
          if (Num_Indent==8) then
              print *, "           Update nodal displacement..."
          elseif(Num_Indent==5)then
              print *, "        Update nodal displacement..."
          endif  
          delta_disp(freeDOF(1:num_FreeD)) =loads(1:num_FreeD)
          NR_DISP(freeDOF(1:num_FreeD)) = NR_DISP(freeDOF(1:num_FreeD)) + delta_disp(freeDOF(1:num_FreeD))
        
          if (Num_Indent==8) then
              print *,'           Check convergence.'
          elseif(Num_Indent==5)then
              print *,'        Check convergence.'
          endif  
          call Cal_Contact_Conve_Factor(isub,1,1,i_NR_P,Conve_Tolerance,Total_FD,freeDOF,num_freeD, &
                 globalF,R_PSI,Last_R_PSI,delta_disp,NR_DISP,Contact_DISP_0,Yes_Conve,Conve_Factor)
 
          if (Num_Indent==8) then
              write(*,4022) Conve_Factor,Conve_Tolerance
          elseif(Num_Indent==5)then
              write(*,4023) Conve_Factor,Conve_Tolerance
          endif  
          Saved_Conv_Factor(i_NR_P) =Conve_Factor
          if(Yes_Conve)then
              if (Num_Indent==8) then
                  write(*,4031) 
                  write(*,4032) i_NR_P
                  write(*,4033) count(Elem_Conta_Sta/=0)
                  write(*,4031)
              elseif(Num_Indent==5)then
                  write(*,4131)
                  write(*,4132) i_NR_P
                  write(*,4133) count(Elem_Conta_Sta/=0)
                  write(*,4131)
              endif  
              exit
          endif
          if (i_NR_P>=6)then
              if (Num_Indent==8) then
                  print *,'           Check oscillation.'
              elseif(Num_Indent==5)then
                  print *,'        Check oscillation.'
              endif  
              call Tool_Check_Oscillation_by_6_Variables(Saved_Conv_Factor(i_NR_P-5:i_NR_P),Yes_Oscill)
              if (Yes_Oscill) then
                  if (Num_Indent==8) then
                      write(*,2003)
                  elseif(Num_Indent==5)then
                      write(*,2013)
                  endif 
                  exit
              endif
          endif
          Last_R_PSI = R_PSI
      enddo
          
      if(.not. Yes_Conve)then
          if (Num_Indent==8) then
              print *, '       Contact iteration loop done!'
          elseif(Num_Indent==5)then
              print *, '    Contact iteration loop done!'
          endif  
      endif
      
      DISP = NR_DISP
      DEALLOCATE(p,loads,x,xnew,u,pcg_d)         
      if(allocated(u_thread)) deallocate(u_thread)
  case(2)
  
  case(6)
      call Cal_Contact_Contact_State_Gauss_3D(isub,1,1,i_NR_P,Contact_DISP_0,Yes_Contact,Elem_Conta_Sta,   &
            CT_State_Gauss)   

      
      if(Yes_Contact .eqv. .False.)then
          if (Num_Indent==8) then
              print *, '           No penetration detected, leaving contact iteration!'
          elseif(Num_Indent==5)then
              print *, '        No penetration detected, leaving contact iteration!'
          endif  
          
          goto 9999
      endif     
      
      if (Num_Indent==8) then
          print *,"           Update storK_XFEM..."  
      elseif(Num_Indent==5)then
          print *,"        Update storK_XFEM..."  
      endif  
       
      call EBE_XFEM_PCG_3D_Modify_K_CS(isub,num_freeD,freeDOF,size_local_0,              &
                    all_local_0,diag_precon_no_invert_Updated)  
                    
      if (Num_Indent==8) then
          print *, "           Get the preconditioner and get starting loads..."
      elseif(Num_Indent==5)then
          print *, "        Get the preconditioner and get starting loads..."
      endif  
      p    = ZR
      loads=ZR 
      loads(1:num_FreeD) =  globalF(freeDOF(1:num_FreeD))  


      if (Key_EBE_Precondition==0) then
          diag_precon=ONE
          diag_precon(0) =ZR 
      elseif(Key_EBE_Precondition==1) then
          diag_precon(1:)=ONE/diag_precon_no_invert_Updated(1:) 
          diag_precon(0) =ZR 
      endif
      
      pcg_d=diag_precon*loads 
      p=pcg_d
      x=ZR
      cg_iters=0
      
      if (Num_Indent==8) then
          print *, "           PCG equation solution..."
      elseif(Num_Indent==5)then
          print *, "        PCG equation solution..."
      endif  
      
      do i_PCG =1,max_num_PCG
          cg_iters=cg_iters+1 
          u=ZR     

          u_thread(0:num_FreeD,1:max_threads) = ZR   
          c_thread = omp_get_thread_num()+1 

          select case(Key_EBE_Sym_Storage_K)
          
          case(0)
              do i_E=1,num_FEM_Elem
                  c_Elem = FEM_Elem_List(i_E)
                  num_Loc_ESM = Size_Local_3D(c_Elem)
                  local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
                  c_Ele_Location = Elem_Location(c_Elem,2) 
                  u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+ &
                                               MATMUL(storK_FEM(1:num_Loc_ESM,1:num_Loc_ESM,c_Ele_Location),&
                                                      p(local(1:num_Loc_ESM)))
              enddo
          case(1)
              do i_E=1,num_FEM_Elem
                  c_Elem = FEM_Elem_List(i_E)
                  num_Loc_ESM = Size_Local_3D(c_Elem)
                  local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
                  c_Ele_Location = Elem_Location(c_Elem,2) 
                  do c_i=1,num_Loc_ESM
                      do c_j=c_i,num_Loc_ESM
                          localK(c_i,c_j) = storK_FEM_Sym((c_i-1)*24 -(c_i-1)*c_i/2 +c_j,c_Ele_Location)  
                      enddo
                  enddo
                  do c_i=1,num_Loc_ESM
                      do c_j=1,c_i-1
                          localK(c_i,c_j) =  localK(c_j,c_i) 
                      enddo
                  enddo
                  u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+ &
                                               MATMUL(localK(1:num_Loc_ESM,1:num_Loc_ESM),p(local(1:num_Loc_ESM))) 
              enddo
          
          end select
              
          DO i_Thread = 1,omp_get_max_threads()
             u(0:num_FreeD)  =  u(0:num_FreeD)  + u_thread(0:num_FreeD,i_Thread)
          ENDDO  
          
          u_thread(0:num_FreeD,1:max_threads) = ZR   
          c_thread = omp_get_thread_num()+1
          do i_E=1,num_XFEM_Elem
              c_Elem = XFEM_Elem_List(i_E)
              num_Loc_ESM = Size_Local_3D(c_Elem)
              local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
              c_Ele_Location = Elem_Location(c_Elem,1) 
              u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+ &
                                           MATMUL(storK_XFEM_Updated(c_Ele_Location)%row(1:num_Loc_ESM,1:num_Loc_ESM),&
                                                  p(local(1:num_Loc_ESM)))
          enddo
          DO i_Thread = 1,omp_get_max_threads()
             u(0:num_FreeD)  =  u(0:num_FreeD)  + u_thread(0:num_FreeD,i_Thread)
          ENDDO  
          
          
          
          up=DOT_PRODUCT(loads,pcg_d)
          

          
          alpha=up/DOT_PRODUCT(p,u) 

          
          xnew=x+p*alpha 
          loads=loads-u*alpha
          pcg_d=diag_precon*loads 
          
          beta=DOT_PRODUCT(loads,pcg_d)/up 

          
          p=pcg_d+p*beta

          Tol = MAXVAL(ABS(x(0:num_FreeD)-xnew(0:num_FreeD)))/MAXVAL(ABS(x(0:num_FreeD)))
          
          x=xnew
     
          
          if(Tol<c_cg_tol.OR.cg_iters==max_num_PCG)then
              exit
          endif
      enddo
      if (Num_Indent==8) then
          print *, "           Number of CG iterations to convergence was",cg_iters
      elseif(Num_Indent==5)then
          print *, "        Number of CG iterations to convergence was",cg_iters
      endif  
      
      
      loads=xnew
      if (Num_Indent==8) then
          print *, "           Update nodal displacement..." 
      elseif(Num_Indent==5)then
          print *, "        Update nodal displacement..." 
      endif  
          
      NR_DISP(freeDOF(1:num_FreeD)) = loads(1:num_FreeD)
      DISP(freeDOF(1:num_FreeD)) = NR_DISP(freeDOF(1:num_FreeD))
      DEALLOCATE(p,loads,x,xnew,u,pcg_d)         
      if(allocated(u_thread)) deallocate(u_thread)
      
  case(7)
      do i_NR_P = 1,Max_Contact_Iter
          write(*,4001) i_NR_P,Max_Contact_Iter
          if(i_NR_P==1)then 
              NR_DISP(1:Total_FD)    = ZR
              delta_disp(1:Total_FD)  = ZR                
              Last_R_PSI(1:Total_FD) = ZR
              Elem_Conta_Sta(1:Num_Elem,1:num_Crack)  = 0  
              Kn_Gauss(1:num_Crack,1:Max_Max_N_FluEl_3D)  = kn
              Kt1_Gauss(1:num_Crack,1:Max_Max_N_FluEl_3D) = kt
              Kt2_Gauss(1:num_Crack,1:Max_Max_N_FluEl_3D) = kt
        
              call Cal_Contact_Contact_State_Gauss_3D(isub,1,1,i_NR_P,Contact_DISP_0,Yes_Contact,Elem_Conta_Sta,   &
                    CT_State_Gauss)   
              
              if(Yes_Contact .eqv. .False.)then
                  if (Num_Indent==8) then
                      print *, '           No penetration detected, leaving contact iteration!'
                  elseif(Num_Indent==5)then
                      print *, '        No penetration detected, leaving contact iteration!'
                  endif  
                  
                  goto 9999
              endif
          endif              
          
          if (Num_Indent==8) then
              print *,'           Get contact force and update Kt.'
          elseif(Num_Indent==5)then
              print *,'        Get contact force and update Kt.'
          endif  
          
          call Cal_Contact_PN_and_PT_3D(isub,1,1,i_NR_P,Total_FD,num_freeD,           &
                    Kn,Kn_Gauss,Kt1_Gauss,Kt2_Gauss,fric_mu,NR_DISP,delta_disp,       &
                    CT_State_Gauss,Elem_Conta_Sta,PC_Gauss_x,PC_Gauss_y,PC_Gauss_z,Num_Indent)
                    
          if (Num_Indent==8) then
              write(*,5001) count(Elem_Conta_Sta==1),count(Elem_Conta_Sta==2)
          elseif(Num_Indent==5)then
              write(*,5101) count(Elem_Conta_Sta==1),count(Elem_Conta_Sta==2)
          endif  
          
          if (Num_Indent==8) then
              print *,'           Assemble Jacobian matrix.'
          elseif(Num_Indent==5)then
              print *,'        Assemble Jacobian matrix.'
          endif  
          
          call EBE_Cal_Contact_Jacobian_3D(isub,i_NR_P,num_freeD,freeDOF,size_local_0,              &
                    all_local_0,diag_precon_no_invert_0,diag_precon_no_invert_Updated,              &
                    Kn,Kn_Gauss,Kt1_Gauss,Kt2_Gauss,CT_State_Gauss)  
          
          
          if (Num_Indent==8) then
              print *,'           Get residual.'
          elseif(Num_Indent==5)then
              print *,'        Get residual.'
          endif  
          
          call Cal_Contact_Resid_3D(isub,1,1,i_NR_P,Total_FD,num_freeD,globalF,     &
                           NR_DISP,freeDOF,PC_Gauss_x,PC_Gauss_y,PC_Gauss_z,R_PSI)     

          if (Num_Indent==8) then
              print *, "           Get the preconditioner and get starting loads..."
          elseif(Num_Indent==5)then
              print *, "        Get the preconditioner and get starting loads..."
          endif  
          
          p    = ZR
          p(0) = ZR
          loads= ZR 
          loads(1:num_FreeD) = -R_PSI(freeDOF(1:num_FreeD))
          loads(0)= ZR 

              
          if (Key_EBE_Precondition==0) then
              diag_precon=ONE
              diag_precon(0) =ZR 
          elseif(Key_EBE_Precondition==1) then
              diag_precon(1:)=ONE/diag_precon_no_invert_Updated(1:) 
              diag_precon(0) =ZR 
          endif
          
          pcg_d=diag_precon*loads 
          p    =pcg_d
          x    =ZR
          cg_iters=0
          
          if (Num_Indent==8) then
              print *, '           PCG equation solution...'
          elseif(Num_Indent==5)then
              print *, '        PCG equation solution...'
          endif  
          
          do i_PCG =1,max_num_PCG
              cg_iters=cg_iters+1 
              u=ZR     
              
              select case(Key_EBE_Sym_Storage_K)
              
              case(0)
                  u_thread(0:num_FreeD,1:max_threads) = ZR   
                  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,c_Ele_Location)
                  c_thread = omp_get_thread_num()+1
                  !$OMP DO SCHEDULE(static)  
                  do i_E=1,num_FEM_Elem
                      c_Elem = FEM_Elem_List(i_E)
                      num_Loc_ESM = Size_Local_3D(c_Elem)
                      local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
                      c_Ele_Location = Elem_Location(c_Elem,2)  
                      u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+ &
                                                   MATMUL(storK_FEM(1:num_Loc_ESM,1:num_Loc_ESM,c_Ele_Location),&
                                                          p(local(1:num_Loc_ESM)))
                  enddo
                  !$omp end do 
                  !$omp end parallel
              
              case(1)
                  u_thread(0:num_FreeD,1:max_threads) = ZR   
                  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,localK,c_Ele_Location,c_i,c_j)
                  c_thread = omp_get_thread_num()+1
                  !$OMP DO SCHEDULE(static)   
                  do i_E=1,num_FEM_Elem
                      c_Elem = FEM_Elem_List(i_E)
                      num_Loc_ESM = Size_Local_3D(c_Elem)
                      local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
                      c_Ele_Location = Elem_Location(c_Elem,2) 
                      do c_i=1,num_Loc_ESM
                          do c_j=c_i,num_Loc_ESM
                              localK(c_i,c_j) = storK_FEM_Sym((c_i-1)*24 -(c_i-1)*c_i/2 +c_j,c_Ele_Location)  
                          enddo
                      enddo
                      do c_i=1,num_Loc_ESM
                          do c_j=1,c_i-1
                              localK(c_i,c_j) =  localK(c_j,c_i) 
                          enddo
                      enddo
                      u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+ &
                                                   MATMUL(localK(1:num_Loc_ESM,1:num_Loc_ESM),p(local(1:num_Loc_ESM))) 
                  enddo
                  !$omp end do 
                  !$omp end parallel
              end select
              
              DO i_Thread = 1,omp_get_max_threads()
                 u(0:num_FreeD)  =  u(0:num_FreeD)  + u_thread(0:num_FreeD,i_Thread)
              ENDDO  

              u_thread(0:num_FreeD,1:max_threads) = ZR   
              !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,c_Ele_Location)
              c_thread = omp_get_thread_num()+1
              !$OMP DO SCHEDULE(static) 
              do i_E=1,num_XFEM_Elem
                  c_Elem = XFEM_Elem_List(i_E)
                  num_Loc_ESM = Size_Local_3D(c_Elem)
                  local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
                  c_Ele_Location = Elem_Location(c_Elem,1) 
                  u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+ &
                                               MATMUL(storK_XFEM_Updated(c_Ele_Location)%row(1:num_Loc_ESM,1:num_Loc_ESM),&
                                                      p(local(1:num_Loc_ESM))) 
              enddo
              !$omp end do 
              !$omp end parallel
              DO i_Thread = 1,omp_get_max_threads()
                 u(0:num_FreeD)  =  u(0:num_FreeD)  + u_thread(0:num_FreeD,i_Thread)
              ENDDO  
              
              up=DOT_PRODUCT(loads,pcg_d)

              alpha=up/DOT_PRODUCT(p,u) 
              xnew=x+p*alpha 
              loads=loads-u*alpha
              pcg_d=diag_precon*loads 
              
              beta=DOT_PRODUCT(loads,pcg_d)/up 
              
              p=pcg_d+p*beta

              Tol = MAXVAL(ABS(x(0:num_FreeD)-xnew(0:num_FreeD)))/MAXVAL(ABS(x(0:num_FreeD)))
              
              x=xnew
              
              if(cg_iters<=300)then
                if(mod(cg_iters,50) == 0) then
                    if (Num_Indent==8) then
                        write(*,103)  cg_iters,Tol,c_cg_tol 
                    elseif(Num_Indent==5)then
                        write(*,113)  cg_iters,Tol,c_cg_tol 
                    endif  
                endif
              else
                if(mod(cg_iters,100) == 0) then
                    if (Num_Indent==8) then
                        write(*,103)  cg_iters,Tol,c_cg_tol 
                    elseif(Num_Indent==5)then
                        write(*,113)  cg_iters,Tol,c_cg_tol 
                    endif  
                endif
              endif   
              
              if(Tol<c_cg_tol)then  
                  exit
              endif
              
              if(cg_iters==max_num_PCG)then
                  if (Num_Indent==8) then
                      print *, '           -------------------------------------'
                      print *, '               Warning :: PCG-EBE failed!       '  
                      print *, '           -------------------------------------'
                  elseif(Num_Indent==5)then
                      print *, '        -------------------------------------'
                      print *, '            Warning :: PCG-EBE failed!       '  
                      print *, '        -------------------------------------'
                  endif
                  exit
              endif 
      
          enddo
          
          if (Num_Indent==8) then
              write(*,101)  cg_iters      
              write(*,102)  Tol,c_cg_tol 
          elseif(Num_Indent==5)then
              write(*,111)  cg_iters      
              write(*,112)  Tol,c_cg_tol 
          endif
          
          loads=xnew
          delta_disp = ZR
          if (Num_Indent==8) then
              print *, "           Update nodal displacement..." 
          elseif(Num_Indent==5)then
              print *, "        Update nodal displacement..." 
          endif  
              
          delta_disp(freeDOF(1:num_FreeD)) =loads(1:num_FreeD)
          NR_DISP(freeDOF(1:num_FreeD)) = NR_DISP(freeDOF(1:num_FreeD)) + delta_disp(freeDOF(1:num_FreeD))
        
          if (Num_Indent==8) then
              print *,'           Check convergence.'
          elseif(Num_Indent==5)then
              print *,'        Check convergence.'
          endif  
          
          call Cal_Contact_Conve_Factor(isub,1,1,i_NR_P,Conve_Tolerance,Total_FD,freeDOF,num_freeD, &
                 globalF,R_PSI,Last_R_PSI,delta_disp,NR_DISP,Contact_DISP_0,Yes_Conve,Conve_Factor)
 
          
          if (Num_Indent==8) then
              write(*,4022) Conve_Factor,Conve_Tolerance
          elseif(Num_Indent==5)then
              write(*,4023) Conve_Factor,Conve_Tolerance
          endif 
          
          Saved_Conv_Factor(i_NR_P) =Conve_Factor
          if(Yes_Conve)then 
              if (Num_Indent==8) then
                  write(*,4031) 
                  write(*,4032) i_NR_P
                  write(*,4033) count(Elem_Conta_Sta/=0)
                  write(*,4031)
              elseif(Num_Indent==5)then
                  write(*,4131)
                  write(*,4132) i_NR_P
                  write(*,4133) count(Elem_Conta_Sta/=0)
                  write(*,4131)
              endif 
              exit
          endif
          if (i_NR_P>=6)then
              if (Num_Indent==8) then
                  print *,'           Check oscillation.'
              elseif(Num_Indent==5) then
                  print *,'        Check oscillation.'
              endif  
              
              call Tool_Check_Oscillation_by_6_Variables(Saved_Conv_Factor(i_NR_P-5:i_NR_P),Yes_Oscill)
              if (Yes_Oscill) then
                  if (Num_Indent==8) then
                      write(*,2003)
                  elseif(Num_Indent==5)then
                      write(*,2013)
                  endif 
                  exit
              endif
          endif
          Last_R_PSI = R_PSI
      enddo
          
      if(.not. Yes_Conve)then
          if (Num_Indent==8) then
              print *, '       Contact iteration loop done!'
          elseif(Num_Indent==5)then
              print *, '    Contact iteration loop done!'
          endif  
      endif
      
      DISP = NR_DISP
      DEALLOCATE(p,loads,x,xnew,u,pcg_d)         
      if(allocated(u_thread)) deallocate(u_thread)
end select

9999 continue


if(Key_Save_Nothing ==0 .and. Key_Simple_Post==0)then
    if(Key_Contact==1) then
        write(temp,'(I5)') isub
        c_File_name_1=trim(Full_Pathname)//'.cfrx_'//ADJUSTL(temp)     
        c_File_name_2=trim(Full_Pathname)//'.cfry_'//ADJUSTL(temp)    
        c_File_name_3=trim(Full_Pathname)//'.cfrz_'//ADJUSTL(temp)  
        open(801,file=c_File_name_1,status='unknown') 
        open(802,file=c_File_name_2,status='unknown') 
        open(803,file=c_File_name_3,status='unknown') 
        do i_C=1,num_Crack
              write(801, '(5000E20.12)') (PC_Gauss_x(i_C,j) ,j=1,Cracks_FluidEle_num_3D(i_C))       
              write(802, '(5000E20.12)') (PC_Gauss_y(i_C,j) ,j=1,Cracks_FluidEle_num_3D(i_C))      
              write(803, '(5000E20.12)') (PC_Gauss_z(i_C,j) ,j=1,Cracks_FluidEle_num_3D(i_C))         
        end do
        close(801)          
        close(802)    
        close(803) 
        c_File_name_1=trim(Full_Pathname)//'.csce_'//ADJUSTL(temp)     
        open(901,file=c_File_name_1,status='unknown') 
        do i_C=1,num_Crack
            write(901, '(5000I3)') (CT_State_Gauss(i_C,j),j=1,Cracks_FluidEle_num_3D(i_C))    
        end do
        close(901)    
    endif
endif

if (allocated(storK_XFEM_Updated)) deallocate(storK_XFEM_Updated)

RETURN
END SUBROUTINE EBE_Determine_Contact_State_by_Iteration_3D
