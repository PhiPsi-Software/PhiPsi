 
SUBROUTINE EBE_Determine_Contact_State_by_Iteration(isub, &
                  c_cg_tol,max_num_PCG, &
                  num_FreeD,freeDOF,globalF,DISP,      &
                  storK,size_local,all_local, &
                  diag_precon)


                          
use Global_Float_Type
use Global_Common
use Global_Filename
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack
use Global_Crack_Common
use Global_HF
use Global_Contact
use Global_Inter_Cal_Contact_Red_Resid
use omp_lib


implicit none
integer,intent(in)::isub,max_num_PCG
integer,intent(in)::num_FreeD
integer,intent(in)::freeDOF(num_FreeD)
real(kind=FT),intent(inout)::DISP(Total_FD)
real(kind=FT),intent(inout)::storK(MDOF_2D,MDOF_2D,Num_Elem)
real(kind=FT),intent(inout)::diag_precon(0:num_FreeD)
integer,intent(inout):: size_local(Num_Elem)
integer,intent(inout):: all_local(MDOF_2D,Num_Elem)
real(kind=FT),intent(in)::c_cg_tol
real(kind=FT),intent(in)::globalF(Total_FD)
real(kind=FT)::storK_0(MDOF_2D,MDOF_2D,Num_Elem)
real(kind=FT)::diag_precon_0(0:num_FreeD)
integer:: size_local_0(Num_Elem)
integer:: all_local_0(MDOF_2D,Num_Elem)
integer i_Thread,c_Thread,max_threads
real(kind=FT),ALLOCATABLE::u_thread(:,:) 
integer i_E
logical Yes_Contact
integer i_NR_P
real(kind=FT) R_PSI(Total_FD)
real(kind=FT) Last_R_PSI(Total_FD)
real(kind=FT) NR_DISP(Total_FD)
real(kind=FT) PC_Gauss_x(num_Crack,Max_Num_Cr_CalP-1,2)
real(kind=FT) PC_Gauss_y(num_Crack,Max_Num_Cr_CalP-1,2)
integer CT_State_Gauss(num_Crack,Max_Num_Cr_CalP-1,2)
real(kind=FT) Kn,Kt
real(kind=FT) Kn_Gauss(num_Crack,Max_Num_Cr_CalP-1,2)
real(kind=FT) Kt_Gauss(num_Crack,Max_Num_Cr_CalP-1,2)
real(kind=FT) Conve_Tolerance,Conve_Factor
real(kind=FT) fric_mu
logical Yes_Conve
real(kind=FT) Contact_DISP_0(Total_FD)
integer num_Loc_ESM
real(kind=FT),ALLOCATABLE::p(:),loads(:),x(:),xnew(:),u(:),pcg_d(:)
real(kind=FT) Saved_Conv_Factor(Max_Contact_Iter)
logical Yes_Oscill
integer cg_iters,i_PCG
integer local(MDOF_2D)
real(kind=FT) alpha,beta,up,Tol
real(kind=FT) delta_disp(Total_FD)
real(kind=FT) ddot
real(kind=FT) tem_value_1

1997 FORMAT(7X,'\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!')  
1998 FORMAT(7X,'////////////////////////////////////////!')  
2002 FORMAT(7X,'  PhiPsi failed to converge for contact iteration!')   
2003 FORMAT(7X,'  WARNING :: Oscillation detected!')   
4001 FORMAT(7X,'+++  Contact NR iteration ',I3,' of ',I3,' started:') 
4012 FORMAT(7X,'No penetration detected, leaving contact iteration!') 
4022 FORMAT(12X,'Convergence factor is', E16.7 ,'.')  
4032 FORMAT(7X,'Contact iteration done after ',I3,' tries!') 
4033 FORMAT(7X,'Number of contact elements is ',I5) 
5001 FORMAT(12X,'Sticking elements:',I5,'| Sliding elements:',I5) 
print *,'    Determine contact states by iteration......' 

Contact_DISP_0(1:Total_FD) = ZR
Contact_DISP_0 = DISP
Saved_Conv_Factor(1:Max_Contact_Iter)  = ZR
storK_0 = storK
diag_precon_0(0:num_FreeD) =diag_precon(0:num_FreeD)
size_local_0  = size_local
all_local_0   = all_local
ALLOCATE(p(0:num_FreeD),loads(0:num_FreeD),  &
          x(0:num_FreeD),xnew(0:num_FreeD),u(0:num_FreeD), &
          pcg_d(0:num_FreeD))   
max_threads = omp_get_max_threads()
if (allocated(u_thread)) deallocate(u_thread)
ALLOCATE(u_thread(0:num_FreeD,max_threads))

kn = kn_Cont_Penalty
kt = kt_Cont_Penalty
fric_mu = fric_mu_Cont
Conve_Tolerance = Conve_Tol_Penalty
select case (Key_Contact)
case(1)
  do i_NR_P = 1,Max_Contact_Iter
      write(*,4001) i_NR_P,Max_Contact_Iter
      if(i_NR_P==1)then 
          NR_DISP(1:Total_FD)    = ZR
          delta_disp(1:Total_FD)  = ZR
          
          
          Last_R_PSI(1:Total_FD) = ZR
          Elem_Conta_Sta(1:Num_Elem,1:Max_Num_Cr) = 0 
          Kt_Gauss(1:num_Crack,1:Max_Num_Cr_CalP-1,1:2) = kt
          Kn_Gauss(1:num_Crack,1:Max_Num_Cr_CalP-1,1:2) = kn
          call Cal_Contact_Contact_State_Gauss(isub,1,1,i_NR_P,   &
                    Contact_DISP_0,Yes_Contact,Elem_Conta_Sta,CT_State_Gauss)    
          if(Yes_Contact.eqv..False.)then
              write(*,4012)
              goto 9999
          endif
      endif
      
      print *,'           Get contact force and update Kt.'
      call Cal_Contact_PN_and_PT(isub,1,1,  &
                  i_NR_P,Total_FD,num_freeD, &
                  Kn,Kn_Gauss,Kt_Gauss,fric_mu,NR_DISP,delta_disp,&
                  CT_State_Gauss,Elem_Conta_Sta,PC_Gauss_x,PC_Gauss_y)

      
      write(*,5001) count(Elem_Conta_Sta==1),count(Elem_Conta_Sta==2)

      print *,'           Assemble Jacobian matrix.'
      call EBE_Cal_Contact_Jacobian(     &
                  isub,i_NR_P, &
                  num_freeD,freeDOF,&
                  storK_0,size_local_0,all_local_0,diag_precon_0,&
                  storK,diag_precon,     &
                  Kn,Kn_Gauss,Kt_Gauss,&
                  CT_State_Gauss)
      
      
      print *,'           Get residual.'
      call Cal_Contact_Resid(isub,1,1,i_NR_P,  &
           Total_FD,num_freeD,globalF, &
           NR_DISP,freeDOF(1:num_FreeD),PC_Gauss_x,PC_Gauss_y,R_PSI)
     
      

      print *, "           Get the preconditioner and get starting loads..."
      p = ZR
      loads=ZR 
      loads(1:num_FreeD) = -R_PSI(freeDOF(1:num_FreeD))
      
      diag_precon(1:)=ONE/diag_precon(1:) 
      diag_precon(0) =ZR 


      
      pcg_d=diag_precon*loads 
      p=pcg_d
      x=ZR
      cg_iters=0
      
      print *, "           PCG equation solution..."
      do i_PCG =1,max_num_PCG
          cg_iters=cg_iters+1 
          u=ZR  
          u_thread(0:num_FreeD,1:max_threads)= ZR  
          !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,num_Loc_ESM,local) 
          c_thread = omp_get_thread_num()+1
          !$OMP DO         
          do i_E=1,Num_Elem
                 num_Loc_ESM = size_local(i_E)
                 local(1:num_Loc_ESM)=all_local(1:num_Loc_ESM,i_E) 
                 u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+ &
                       MATMUL(storK(1:num_Loc_ESM,1:num_Loc_ESM,i_E) ,p(local(1:num_Loc_ESM)))  
          enddo  
          !$omp end do
          !$omp end parallel       
            
          DO i_Thread = 1,omp_get_max_threads()
              u  =  u  + u_thread(:,i_Thread)
          ENDDO   
  
#ifdef Silverfrost
          up=DOT_PRODUCT(loads,pcg_d)
#endif
#ifndef Silverfrost
          up=ddot(num_FreeD,loads(1:num_FreeD),1,pcg_d(1:num_FreeD),1)
#endif
          
#ifdef Silverfrost
          alpha=up/DOT_PRODUCT(p,u) 
#endif
#ifndef Silverfrost
          alpha=up/ddot(num_FreeD,p(1:num_FreeD),1,u(1:num_FreeD),1)
#endif          
          xnew=x+p*alpha 
          loads=loads-u*alpha
          pcg_d=diag_precon*loads 
          
#ifdef Silverfrost
          beta=DOT_PRODUCT(loads,pcg_d)/up 
#endif
#ifndef Silverfrost
          beta=ddot(num_FreeD,loads(1:num_FreeD),1,pcg_d(1:num_FreeD),1)/up
#endif          
          p=pcg_d+p*beta

          
          
          
          tem_value_1 = MAXVAL(ABS(x(0:num_FreeD)))
          if(tem_value_1==ZR) tem_value_1=Tol_20
          Tol = MAXVAL(ABS(x(0:num_FreeD)-xnew(0:num_FreeD)))/tem_value_1
          
          
          x=xnew
     
          
          if(Tol<c_cg_tol.OR.cg_iters==max_num_PCG)then
              exit
          endif
      enddo
      print *, "           Number of CG iterations to convergence was",cg_iters
      loads=xnew
      delta_disp = ZR
      print *, "           Update nodal displacement..."     
      delta_disp(freeDOF(1:num_FreeD)) =loads(1:num_FreeD)
      NR_DISP = NR_DISP + delta_disp

      print *,'           Check convergence.'
      call Cal_Contact_Conve_Factor( &
                isub,1,1,i_NR_P,Conve_Tolerance, &
                Total_FD,freeDOF(1:num_freeD),num_freeD, &
                globalF,R_PSI,Last_R_PSI, &
                delta_disp,NR_DISP,Contact_DISP_0, &
                Yes_Conve,Conve_Factor)

      write(*,4022) Conve_Factor
      Saved_Conv_Factor(i_NR_P) =Conve_Factor
      if(Yes_Conve)then
          write(*,1997)
          write(*,4032) i_NR_P
          write(*,4033) count(Elem_Conta_Sta/=0)
          write(*,1997)
          exit
      endif
      if (i_NR_P>=6)then
          print *,'           Check oscillation.'
          call Tool_Check_Oscillation_by_6_Variables(Saved_Conv_Factor(i_NR_P-5:i_NR_P),Yes_Oscill)
          if (Yes_Oscill) then
              write(*,1998)
              write(*,2003)
              write(*,1998)
              exit
          endif
      endif
      Last_R_PSI = R_PSI
  enddo
  
  if(.not.Yes_Conve)then
      write(*,2002)
  endif
  
  DISP = NR_DISP
  DEALLOCATE(p,loads,x,xnew,u,pcg_d)            
case(2)

end select

9999 continue


RETURN
END SUBROUTINE EBE_Determine_Contact_State_by_Iteration
