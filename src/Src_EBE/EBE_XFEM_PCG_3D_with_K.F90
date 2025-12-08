 
SUBROUTINE EBE_XFEM_PCG_3D_with_K(isub,Lambda,c_cg_tol,max_num_PCG,num_FreeD,freeDOF,F,disp,  &
                                  diag_precon_no_invert)

use Global_Float_Type
use Global_Model
use Global_Filename
use Global_Common
use Global_Material
use Global_Crack_3D
use Function_MVMUL_LP
use Global_XFEM_Elements
use omp_lib
      
implicit none
integer,intent(in)::isub,max_num_PCG,num_FreeD
real(kind=FT),intent(in)::Lambda,c_cg_tol,F(num_FreeD)
integer,intent(in)::freeDOF(num_FreeD)
real(kind=FT),intent(in)::diag_precon_no_invert(0:num_FreeD)
real(kind=FT),intent(out)::disp(Total_FD)

integer jj,kk
integer i_E,i_G,c_NN(8)
real(kind=FT) localK(MDOF_3D,MDOF_3D)
real(kind=FT) Diag_localK(MDOF_3D)
integer local(MDOF_3D)
real(kind=FT),ALLOCATABLE::p(:),loads(:),x(:),xnew(:),u(:),diag_precon(:),pcg_d(:)

integer cg_iters,i_PCG
real(kind=FT) alpha,beta,up,Tol
integer i_DOF
real(kind=FT) DOT_PRODUCT_pu,DOT_PRODUCT_lp
integer num_Loc_ESM
integer cEle_Loc


real(kind=FT),ALLOCATABLE::u_thread(:,:)
integer i_Thread,c_thread,max_threads
integer c_Elem


real(kind=FT) ddot
EXTERNAL ddot
real(kind=FT) dot_p_u,dot_l_p
real(kind=FT),ALLOCATABLE:: Vector_x_xnew(:)

integer c_i,c_j


ALLOCATE(p(0:num_FreeD),loads(0:num_FreeD),x(0:num_FreeD),xnew(0:num_FreeD),u(0:num_FreeD), &
          diag_precon(0:num_FreeD),pcg_d(0:num_FreeD))   

ALLOCATE(Vector_x_xnew(0:num_FreeD))

max_threads = omp_get_max_threads()
ALLOCATE(u_thread(0:num_FreeD,1:max_threads))
      


loads=ZR 
loads(1:num_FreeD) = F(1:num_FreeD)



if (Key_EBE_Precondition==0) then
  diag_precon=ONE
  diag_precon(0) =ZR 
elseif(Key_EBE_Precondition==1) then
  diag_precon(1:)=ONE/diag_precon_no_invert(1:)
  diag_precon(0) =ZR 
endif

pcg_d=diag_precon*loads 
p=pcg_d
x=ZR
cg_iters=0
      
      
      
if (Key_BLAS==1) then
    goto 100
endif 
do i_PCG =1,max_num_PCG
    cg_iters=cg_iters+1 
    u=ZR

    select case(Key_EBE_Sym_Storage_K)
  
    case(0)
    u_thread(0:num_FreeD,1:max_threads)= ZR   
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,localK,cEle_Loc) 
    c_thread = omp_get_thread_num()+1
    !$OMP DO 
    do i_E=1,num_FEM_Elem
          c_Elem = FEM_Elem_List(i_E)
          num_Loc_ESM = Size_Local_3D(c_Elem)
          local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
          cEle_Loc = Elem_Location(c_Elem,2) 
          localK(1:num_Loc_ESM,1:num_Loc_ESM)= storK_FEM(1:num_Loc_ESM,1:num_Loc_ESM,cEle_Loc)      
          u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+&
                           MATMUL(localK(1:num_Loc_ESM,1:num_Loc_ESM),p(local(1:num_Loc_ESM)))
    enddo
    !$omp end do
    !$omp end parallel   
    
    case(1)
    u_thread(0:num_FreeD,1:max_threads)= ZR   
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,localK,cEle_Loc,c_i,c_j) 
    c_thread = omp_get_thread_num()+1
    !$OMP DO 
    do i_E=1,num_FEM_Elem
          c_Elem = FEM_Elem_List(i_E)
          num_Loc_ESM = Size_Local_3D(c_Elem)
          local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
          cEle_Loc = Elem_Location(c_Elem,2) 
          do c_i=1,num_Loc_ESM
              do c_j=c_i,num_Loc_ESM
                  localK(c_i,c_j) = storK_FEM_Sym((c_i-1)*24 -(c_i-1)*c_i/2 +c_j,cEle_Loc)  
              enddo
          enddo
          do c_i=1,num_Loc_ESM
              do c_j=1,c_i-1
                  localK(c_i,c_j) =  localK(c_j,c_i) 
              enddo
          enddo
          u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+&
                           MATMUL(localK(1:num_Loc_ESM,1:num_Loc_ESM),p(local(1:num_Loc_ESM)))
    enddo
    !$omp end do
    !$omp end parallel 
    
    endselect
    
    DO i_Thread = 1,omp_get_max_threads()
         u(0:num_FreeD)  =  u(0:num_FreeD)  + u_thread(0:num_FreeD,i_Thread)
    ENDDO

    u_thread(0:num_FreeD,1:max_threads)= ZR   
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,localK,cEle_Loc) 
    c_thread = omp_get_thread_num()+1
    !$OMP DO
    do i_E=1,num_XFEM_Elem
          c_Elem = XFEM_Elem_List(i_E)
          num_Loc_ESM = Size_Local_3D(c_Elem)
          local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
          cEle_Loc = Elem_Location(c_Elem,1) 
          
          localK(1:num_Loc_ESM,1:num_Loc_ESM)= storK_XFEM(cEle_Loc)%row(1:num_Loc_ESM,1:num_Loc_ESM)  
          
          u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+&
                           MATMUL(localK(1:num_Loc_ESM,1:num_Loc_ESM),p(local(1:num_Loc_ESM)))
    enddo
    !$omp end do
    !$omp end parallel   
    DO i_Thread = 1,omp_get_max_threads()
         u(0:num_FreeD)  =  u(0:num_FreeD)  + u_thread(0:num_FreeD,i_Thread)
    ENDDO

       
      up=DOT_PRODUCT(loads,pcg_d)
      alpha=up/DOT_PRODUCT(p,u) 
      
      xnew = x + p*alpha 

      
      loads=loads-u*alpha
      pcg_d=diag_precon*loads 

      
      beta=DOT_PRODUCT(loads,pcg_d)/up 
      
       p=pcg_d+p*beta

      if(i_PCG==1) then
          Tol = 100.0D0
      else
          Tol = MAXVAL(ABS(x(0:num_FreeD)-xnew(0:num_FreeD)))/MAXVAL(ABS(x(0:num_FreeD)))
      endif
      
      x=xnew

      if(cg_iters<=300)then
        if(mod(cg_iters,50) == 0) then
            write(*,103)  cg_iters,Tol,c_cg_tol 
        endif
      else
        if(mod(cg_iters,100) == 0) then
            write(*,103)  cg_iters,Tol,c_cg_tol 
        endif
      endif   
      
      if(Tol<c_cg_tol)then  
          exit
      endif
      
      if(cg_iters==max_num_PCG)then
          print *, '    -------------------------------------'
          print *, '        Warning :: PCG-EBE failed!       '  
          print *, '    -------------------------------------'
          goto 200
      endif          
enddo







100 continue
do i_PCG =1,max_num_PCG
    cg_iters=cg_iters+1 
    u=ZR
    select case(Key_EBE_Sym_Storage_K)
  
    case(0)
    u_thread(0:num_FreeD,1:max_threads)= ZR   
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,localK,cEle_Loc) 
    c_thread = omp_get_thread_num()+1
    !$OMP DO 
    do i_E=1,num_FEM_Elem
          c_Elem = FEM_Elem_List(i_E)
          num_Loc_ESM = Size_Local_3D(c_Elem)
          local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
          cEle_Loc = Elem_Location(c_Elem,2) 
          localK(1:num_Loc_ESM,1:num_Loc_ESM)= storK_FEM(1:num_Loc_ESM,1:num_Loc_ESM,cEle_Loc)      
          u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+&
                           MATMUL(localK(1:num_Loc_ESM,1:num_Loc_ESM),p(local(1:num_Loc_ESM)))
    enddo
    !$omp end do
    !$omp end parallel   
    
    case(1)
    u_thread(0:num_FreeD,1:max_threads)= ZR   
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,localK,cEle_Loc,c_i,c_j) 
    c_thread = omp_get_thread_num()+1
    !$OMP DO 
    do i_E=1,num_FEM_Elem
          c_Elem = FEM_Elem_List(i_E)
          num_Loc_ESM = Size_Local_3D(c_Elem)
          local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
          cEle_Loc = Elem_Location(c_Elem,2) 
          do c_i=1,num_Loc_ESM
              do c_j=c_i,num_Loc_ESM
                  localK(c_i,c_j) = storK_FEM_Sym((c_i-1)*24 -(c_i-1)*c_i/2 +c_j,cEle_Loc)  
              enddo
          enddo
          do c_i=1,num_Loc_ESM
              do c_j=1,c_i-1
                  localK(c_i,c_j) =  localK(c_j,c_i) 
              enddo
          enddo
          u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+&
                           MATMUL(localK(1:num_Loc_ESM,1:num_Loc_ESM),p(local(1:num_Loc_ESM)))
    enddo
    !$omp end do
    !$omp end parallel       
    endselect
    
    DO i_Thread = 1,omp_get_max_threads()
         
#ifdef Silverfrost
         u(0:num_FreeD)  =  u(0:num_FreeD)  + u_thread(0:num_FreeD,i_Thread)
#endif
#ifndef Silverfrost
         call DAXPY(num_FreeD,ONE,u_thread(1:num_FreeD,i_Thread),1,u(1:num_FreeD),1) 
         u(0) = u(0) + u_thread(0,i_Thread)
#endif
    ENDDO

    u_thread(0:num_FreeD,1:max_threads)= ZR   
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,localK,cEle_Loc) 
    c_thread = omp_get_thread_num()+1
    !$OMP DO
    do i_E=1,num_XFEM_Elem
          c_Elem = XFEM_Elem_List(i_E)
          num_Loc_ESM = Size_Local_3D(c_Elem)
          local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
          cEle_Loc = Elem_Location(c_Elem,1) 
          
          localK(1:num_Loc_ESM,1:num_Loc_ESM)= storK_XFEM(cEle_Loc)%row(1:num_Loc_ESM,1:num_Loc_ESM)
          
          u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+&
                           MATMUL(localK(1:num_Loc_ESM,1:num_Loc_ESM),p(local(1:num_Loc_ESM)))
    enddo
    !$omp end do
    !$omp end parallel   
    DO i_Thread = 1,omp_get_max_threads()
         
#ifdef Silverfrost  
         u(0:num_FreeD)  =  u(0:num_FreeD)  + u_thread(0:num_FreeD,i_Thread)
#endif
#ifndef Silverfrost  
         call DAXPY(num_FreeD,ONE,u_thread(1:num_FreeD,i_Thread),1,u(1:num_FreeD),1) 
         u(0) = u(0) + u_thread(0,i_Thread)    
#endif
    ENDDO


       
#ifdef Silverfrost
      up=DOT_PRODUCT(loads,pcg_d)
      alpha=up/DOT_PRODUCT(p,u)   
#endif
#ifndef Silverfrost
      up      = ddot(num_FreeD,loads(1:num_FreeD),1,pcg_d(1:num_FreeD),1)      
      up      = up + loads(0)*pcg_d(0)
      dot_p_u = ddot(num_FreeD,p(1:num_FreeD),1,u(1:num_FreeD),1) 
      dot_p_u = dot_p_u + p(0)*u(0)
      alpha   = up/dot_p_u  
#endif
      
#ifdef Silverfrost
      xnew = x + p*alpha 
#endif
#ifndef Silverfrost
      call DCOPY(num_FreeD,x(1:num_FreeD),1,xnew(1:num_FreeD),1)
      xnew(0) = x(0)
      call DAXPY(num_FreeD,alpha,p(1:num_FreeD),1,xnew(1:num_FreeD),1)
      xnew(0) = xnew(0) + alpha*p(0)
#endif
#ifdef Silverfrost
      loads=loads-u*alpha
      pcg_d=diag_precon*loads  
#endif
#ifndef Silverfrost
      call DAXPY(num_FreeD,-alpha,u(1:num_FreeD),1,loads(1:num_FreeD),1)  
      loads(0) = loads(0) - alpha*u(0)
      pcg_d=diag_precon*loads 
#endif
       beta=DOT_PRODUCT(loads,pcg_d)/up    
      
      
#ifdef Silverfrost
      p=pcg_d+p*beta     
#endif      
#ifndef Silverfrost
      call DSCAL(num_FreeD,beta,p(1:num_FreeD),1)
      p(0) = p(0)*beta
      call DAXPY(num_FreeD,ONE,pcg_d(1:num_FreeD),1,p(1:num_FreeD),1)
      p(0) = p(0)+pcg_d(0)
#endif
      
#ifdef Silverfrost
      Tol = MAXVAL(ABS(x(0:num_FreeD)-xnew(0:num_FreeD)))/MAXVAL(ABS(x(0:num_FreeD)))
#endif
#ifndef Silverfrost
      call DCOPY(num_FreeD,x(1:num_FreeD),1,Vector_x_xnew(1:num_FreeD),1)
      Vector_x_xnew(0) = x(0)
      call DAXPY(num_FreeD,-ONE,xnew(1:num_FreeD),1,Vector_x_xnew(1:num_FreeD),1)   
      Vector_x_xnew(0) = Vector_x_xnew(0) - xnew(0)
      Tol = MAXVAL(ABS(Vector_x_xnew(0:num_FreeD)))/MAXVAL(ABS(x(0:num_FreeD)))
#endif
      
      
#ifdef Silverfrost
      x=xnew
#endif
#ifndef Silverfrost
      call DCOPY(num_FreeD,xnew(1:num_FreeD),1,x(1:num_FreeD),1)
      x(0) = xnew(0)
#endif

      if(cg_iters<=300)then
        if(mod(cg_iters,50) == 0) then
            write(*,103)  cg_iters,Tol,c_cg_tol 
        endif
      else
        if(mod(cg_iters,100) == 0) then
            write(*,103)  cg_iters,Tol,c_cg_tol 
        endif
      endif   
      
  
      if(Tol<c_cg_tol)then  
          exit
      endif
      
      if(cg_iters==max_num_PCG)then
          print *, "    -------------------------------------"
          print *, "    Warning :: PCG-EBE failed!           "     
          print *, "    -------------------------------------"
          exit
      endif          
enddo

200 continue

101 FORMAT(8X,'PCG-EBE: number of iterations to convergence is ',I5)  
102 FORMAT(8X,'PCG-EBE: convergence factor is ',E12.5,' / ',E12.5) 
103 FORMAT(17X,'CG_Iters:',I5, ' -> CG_Factor:', E12.5,' / ',E12.5)   
write(*,101)  cg_iters      
write(*,102)  Tol,c_cg_tol       

loads=xnew
disp = ZR
disp(freeDOF(1:num_FreeD)) =loads(1:num_FreeD)      
DEALLOCATE(p,loads,x,xnew,u,diag_precon,pcg_d)  
if(allocated(u_thread)) deallocate(u_thread)

RETURN
END SUBROUTINE EBE_XFEM_PCG_3D_with_K
