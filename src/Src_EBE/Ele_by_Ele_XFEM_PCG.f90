 
SUBROUTINE Ele_by_Ele_XFEM_PCG(isub,Lambda,c_cg_tol,max_num_PCG,  &
                             num_FreeD,freeDOF,F,disp,  &
                             Total_Num_G_P,  &
                             storK,size_local,all_local,  &
                             diag_precon_no_invert)

use Global_Float_Type
use Global_Crack
use Global_Crack_Common
use Global_Model
use Global_Filename
use Global_Common
use Global_Material
use Global_HF
use Global_Contact
use Global_Inclusion
use Global_Cross
use omp_lib

implicit none
integer,intent(in)::isub,max_num_PCG,num_FreeD
real(kind=FT),intent(in)::Lambda,c_cg_tol,F(num_FreeD)
integer,intent(in)::freeDOF(num_FreeD)
integer,intent(out)::Total_Num_G_P
real(kind=FT),intent(out)::disp(Total_FD)
real(kind=FT),intent(out)::storK(MDOF_2D,MDOF_2D,Num_Elem)
real(kind=FT),intent(out)::diag_precon_no_invert(0:num_FreeD)
integer,intent(out):: size_local(Num_Elem)
integer,intent(out):: all_local(MDOF_2D,Num_Elem)
real(kind=FT)::diag_precon(0:num_FreeD)
integer j
integer i_E,i_G,c_NN(4),mat_num,k
real(kind=FT) c_D(3,3)
real(kind=FT) c_X_NODES(4),c_Y_NODES(4)
integer local(MDOF_2D)

real(kind=FT),ALLOCATABLE::p(:),loads(:),x(:),xnew(:),u(:),pcg_d(:)
integer cg_iters,i_PCG
real(kind=FT) alpha,beta,up,Tol
real(kind=FT) kesi_Enr(Num_Gauss_Points),yita_Enr(Num_Gauss_Points),weight_Enr(Num_Gauss_Points)
real(kind=FT) kesi_N_Enr(Num_Gauss_P_FEM),yita_N_Enr(Num_Gauss_P_FEM), weight_N_Enr(Num_Gauss_P_FEM)
integer c_Num_Gauss_Point,i_C
real(kind=FT) kesi(900),yita(900),weight(900)
integer:: Location_ESM(MDOF_2D)
integer::Location_ESM_C_Crack(MDOF_2D),Location_ESM_C_Cr_NoFEM(MDOF_2D)
integer num_Loc_ESM_C_Cr_NoFEM
real(kind=FT) B(3,80),tem_B(3,80),detJ
integer num_B,num_tem_B
integer num_Loc_ESM,num_Loc_ESM_C_Crack

integer c_loca
integer j_C,i_Incl,i_H,i_Cross
real(kind=FT) c_N(2,8),c_G_x,c_G_y,c_thick
logical Yes_Gauss_in_Incl
integer c_Incl_Num,c_MatNum
integer i_Thread,c_Thread,max_threads
real(kind=FT),ALLOCATABLE::u_thread(:,:)

real(kind=FT),ALLOCATABLE::localK_FEM_or_XFEM(:,:)
real(kind=FT) tem_value

Total_Num_G_P = 0
print *, "    >>>> Start of element by element PCG solver <<<<"
print *, "    Step 1: prepare data..."

ALLOCATE(p(0:num_FreeD),loads(0:num_FreeD), &
     x(0:num_FreeD),xnew(0:num_FreeD),u(0:num_FreeD),pcg_d(0:num_FreeD))
diag_precon= ZR

max_threads = omp_get_max_threads()
if (allocated(u_thread)) deallocate(u_thread)
ALLOCATE(u_thread(0:num_FreeD,max_threads))

print *, "    Step 2: get and store element stiffness matrix..."
if (Key_Integral_Sol  == 2)then
  call Cal_Gauss_Points_QUAD(Num_Gauss_Points,kesi_Enr,yita_Enr,weight_Enr)
  call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,kesi_N_Enr,yita_N_Enr,weight_N_Enr)
elseif (Key_Integral_Sol  == 3)then
  call Cal_Gauss_Points_QUAD_for_SUBQUAD(Num_Sub_Quads,kesi_Enr,yita_Enr,weight_Enr)
  call Cal_Gauss_Points_QUAD(Num_Gauss_P_FEM,kesi_N_Enr,yita_N_Enr,weight_N_Enr)
endif

EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Gauss_Points)= .false.
storK(1:MDOF_2D,1:MDOF_2D,1:Num_Elem) =ZR
all_local(1:MDOF_2D,1:Num_Elem) = 0

do i_E=1,Num_Elem
  Ele_GP_Start_Num(i_E) = Total_Num_G_P + 1
  c_NN    = G_NN(1:4,i_E)
  if (num_Crack/= 0 .and. (maxval(Enriched_Node_Type(c_NN,1:num_Crack)).ne.0)) then
      c_Num_Gauss_Point = Num_Gauss_Points
      Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
  elseif(num_Hole/= 0 .and. (maxval(Enriched_Node_Type_Hl(c_NN,1:num_Hole)).ne.0))then
      c_Num_Gauss_Point = Num_Gauss_Points
      Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
  elseif(num_Cross/= 0 .and. (maxval(Enriched_Node_Type_Cross(c_NN,1:num_Cross)).ne.0))then
      c_Num_Gauss_Point = Num_Gauss_Points
      Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
  elseif(num_Inclusion/= 0 .and.(maxval(Enriched_Node_Type_Incl(c_NN,1:num_Inclusion)).ne.0))then
      c_Num_Gauss_Point = Num_Gauss_Points
      Total_Num_G_P     = Total_Num_G_P + c_Num_Gauss_Point
  else
      c_Num_Gauss_Point = Num_Gauss_P_FEM
      Total_Num_G_P     = Total_Num_G_P +c_Num_Gauss_Point
  end if
  num_GP_Elem(i_E) = c_Num_Gauss_Point
enddo
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,mat_num,c_thick,c_D,  &
!$OMP            c_X_NODES,c_Y_NODES,c_NN,kesi,yita,weight,  &
!$OMP            Location_ESM_C_Crack,num_Loc_ESM_C_Crack,  &
!$OMP            Location_ESM_C_Cr_NoFEM,c_G_x,c_G_y,c_MatNum,  &
!$OMP            c_Num_Gauss_Point,i_C,j_C,Location_ESM,num_Loc_ESM,B,  &
!$OMP            num_Loc_ESM_C_Cr_NoFEM,i_G,i_H,i_Cross,i_Incl,  &
!$OMP            k,c_loca,localK_FEM_or_XFEM,num_B,tem_B,num_tem_B,detJ,j,local,  &
!$OMP            Yes_Gauss_in_Incl,c_Incl_Num)
do i_E=1,Num_Elem
  mat_num = Elem_Mat(i_E)
  c_thick = thick(Elem_Mat(i_E))

  c_D(1:3,1:3)     = D(Elem_Mat(i_E),1:3,1:3)

  if(Flag_Weibull_E)then
      if (Key_Weibull_E(Elem_Mat(i_E)) ==1)then
          c_D = Weibull_Elements_D_Matrix(i_E,1:3,1:3)
      endif
  endif

  c_NN    = G_NN(1:4,i_E)
  c_X_NODES = G_X_NODES(1:4,i_E)
  c_Y_NODES = G_Y_NODES(1:4,i_E)
  if (num_Crack/= 0 .and. (maxval(Enriched_Node_Type(c_NN,1:num_Crack)).ne.0))then
      kesi(1:Num_Gauss_Points)    = kesi_Enr
      yita(1:Num_Gauss_Points)    = yita_Enr
      weight(1:Num_Gauss_Points)  = weight_Enr
      c_Num_Gauss_Point = Num_Gauss_Points
  elseif(num_Hole/= 0 .and.(maxval(Enriched_Node_Type_Hl(c_NN,1:num_Hole)).ne.0))then
      kesi(1:Num_Gauss_Points)    = kesi_Enr
      yita(1:Num_Gauss_Points)    = yita_Enr
      weight(1:Num_Gauss_Points)  = weight_Enr
      c_Num_Gauss_Point = Num_Gauss_Points
  elseif(num_Cross/= 0 .and.(maxval(Enriched_Node_Type_Cross(c_NN,1:num_Cross)).ne.0))then
      kesi(1:Num_Gauss_Points)    = kesi_Enr
      yita(1:Num_Gauss_Points)    = yita_Enr
      weight(1:Num_Gauss_Points)  = weight_Enr
      c_Num_Gauss_Point = Num_Gauss_Points
  elseif(num_Inclusion/= 0 .and. (maxval(Enriched_Node_Type_Incl(c_NN,1:num_Inclusion)).ne.0))then
      kesi(1:Num_Gauss_Points)    = kesi_Enr
      yita(1:Num_Gauss_Points)    = yita_Enr
      weight(1:Num_Gauss_Points)  = weight_Enr
      c_Num_Gauss_Point = Num_Gauss_Points
  else
      kesi(1:Num_Gauss_P_FEM)    = kesi_N_Enr
      yita(1:Num_Gauss_P_FEM)    = yita_N_Enr
      weight(1:Num_Gauss_P_FEM)  = weight_N_Enr
      c_Num_Gauss_Point = Num_Gauss_P_FEM
  end if
  Location_ESM(1:MDOF_2D)   = 0
  num_Loc_ESM           = 0
  if(num_Crack/=0)then
    do i_C =1,num_Crack
      call Location_Element_Stiff_Matrix(i_E,i_C,c_POS(1:Num_Node,i_C),&
                         Location_ESM_C_Crack, num_Loc_ESM_C_Crack,&
                          Location_ESM_C_Cr_NoFEM, num_Loc_ESM_C_Cr_NoFEM)
      Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = &
        Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
      num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
    end do
  endif
  if(num_Hole/=0)then
    do i_H =1,num_Hole
      call Location_Element_Stiff_Matrix_Hl(i_E,i_H,c_POS_Hl(1:Num_Node,i_H),Location_ESM_C_Crack, &
                       num_Loc_ESM_C_Crack,Location_ESM_C_Cr_NoFEM,num_Loc_ESM_C_Cr_NoFEM)
      Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
      num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
    end do
  endif
  if(num_Cross/=0)then
    do i_Cross =1,num_Cross
      call Location_Element_Stiff_Matrix_Cross(i_E,i_Cross,c_POS_Cross(1:Num_Node,i_Cross),Location_ESM_C_Crack,  &
                         num_Loc_ESM_C_Crack,Location_ESM_C_Cr_NoFEM,num_Loc_ESM_C_Cr_NoFEM)
      Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
      num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
    end do
  endif
  if(num_Inclusion/=0)then
    do i_Incl =1,num_Inclusion
      call Location_Element_Stiff_Matrix_Incl(i_E,i_Incl,c_POS_Incl(1:Num_Node,i_Incl), &
                            Location_ESM_C_Crack,num_Loc_ESM_C_Crack,Location_ESM_C_Cr_NoFEM, &
                            num_Loc_ESM_C_Cr_NoFEM)
      Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
      num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
    end do
  endif
  size_local(i_E) = num_Loc_ESM
  DO k=1,num_Loc_ESM
      if (any(freeDOF(1:num_FreeD)== Location_ESM(k)))then
          c_loca = minloc(freeDOF(1:num_FreeD),1, MASK = (freeDOF(1:num_FreeD).eq.Location_ESM(k)))
      else
          c_loca =0
      endif

      local(k) = c_loca
  enddo
  all_local(1:num_Loc_ESM,i_E) = local(1:num_Loc_ESM)
  allocate(localK_FEM_or_XFEM(MDOF_2D,MDOF_2D))
  localK_FEM_or_XFEM(1:MDOF_2D,1:MDOF_2D) = ZR

  do i_G = 1,c_Num_Gauss_Point
      B(1:3,1:80) = ZR
      num_B = 0
      call Cal_N(kesi(i_G),yita(i_G),c_N)
      call Cal_detJ(kesi(i_G),yita(i_G),c_X_NODES,c_Y_NODES,detJ)
      c_G_x = DOT_PRODUCT(c_N(1,1:7:2),c_X_NODES(1:4))
      c_G_y = DOT_PRODUCT(c_N(1,1:7:2),c_Y_NODES(1:4))
      if(num_Crack/=0)then
        do i_C =1,num_Crack
          call Cal_B_Matrix_Crack(kesi(i_G),yita(i_G),i_C,i_E,i_G,c_NN,c_X_NODES,c_Y_NODES,  &
                             tem_B,num_tem_B)
          B(1:3,num_B+1:num_B+num_tem_B) =  tem_B(1:3,1:num_tem_B)
          num_B = num_B + num_tem_B
        end do
      endif
      if(num_Hole/=0)then
          do i_H =1,num_Hole
              call Cal_B_Matrix_Hl(kesi(i_G),yita(i_G),i_H,i_E,i_G, &
                                 c_NN,c_X_NODES,c_Y_NODES,tem_B,num_tem_B)
              B(1:3,num_B+1:num_B+num_tem_B) =  tem_B(1:3,1:num_tem_B)
              num_B = num_B + num_tem_B
          end do
      endif
      if(num_Cross/=0)then
          do i_Cross =1,num_Cross
              call Cal_B_Matrix_Cross(kesi(i_G),yita(i_G),i_Cross,i_E,i_G, &
                                 c_NN,c_X_NODES,c_Y_NODES,tem_B,num_tem_B)
              B(1:3,num_B+1:num_B+num_tem_B) =  tem_B(1:3,1:num_tem_B)
              num_B = num_B + num_tem_B
          end do
      endif
      if(num_Circ_Incl/=0)then
          do i_Incl =1,num_Circ_Incl
              call Cal_B_Matrix_Circ_Incl(kesi(i_G),yita(i_G),i_Incl,i_E,i_G, &
                                 c_NN,c_X_NODES,c_Y_NODES,tem_B,num_tem_B)
              B(1:3,num_B+1:num_B+num_tem_B) =  tem_B(1:3,1:num_tem_B)
              num_B = num_B + num_tem_B
          end do
          call Tool_Yes_Point_in_Inclusions(c_G_x,c_G_y,Yes_Gauss_in_Incl,c_Incl_Num)
          if(Yes_Gauss_in_Incl)then
              c_MatNum = Circ_Inclu_Mat_Num(c_Incl_Num)
              c_D     = D(c_MatNum,:,:)
          endif
      endif
      if(num_Poly_Incl/=0)then
          do i_Incl =1,num_Poly_Incl
              call Cal_B_Matrix_Poly_Incl(kesi(i_G),  &
                                yita(i_G),i_Incl,i_E,i_G,c_NN,c_X_NODES,c_Y_NODES,tem_B,num_tem_B)
              B(1:3,num_B+1:num_B+num_tem_B) =  tem_B(1:3,1:num_tem_B)
              num_B = num_B + num_tem_B
          end do
          call Tool_Yes_Point_in_Inclusions(c_G_x,c_G_y,Yes_Gauss_in_Incl,c_Incl_Num)
          if(Yes_Gauss_in_Incl)then
              c_MatNum = Poly_Inclu_Mat_Num(c_Incl_Num)
              c_D     = D(c_MatNum,:,:)
          endif
      endif
      localK_FEM_or_XFEM(1:num_B,1:num_B) = localK_FEM_or_XFEM(1:num_B,1:num_B) +  &
                    c_thick*weight(i_G)*detJ*MATMUL(MATMUL(transpose  &
                        (B(1:3,1:num_B)),c_D), B(1:3,1:num_B))
  enddo
  storK(1:num_B,1:num_B,i_E)=localK_FEM_or_XFEM(1:num_B,1:num_B)
  DO j=1,num_B
      c_loca=local(j)
      if(abs(localK_FEM_or_XFEM(j,j))<=Tol_10)then
          localK_FEM_or_XFEM(j,j)  = 1.0D0
      endif
      diag_precon(c_loca)=diag_precon(c_loca) + localK_FEM_or_XFEM(j,j)
  END DO
  deallocate(localK_FEM_or_XFEM)
enddo
!$omp end parallel do

diag_precon_no_invert(1:num_FreeD)= diag_precon(1:num_FreeD)
diag_precon_no_invert(0) = ZR


print *, "    Step 3: invert the preconditioner and get starting loads..."
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

print *, "    Step 4: pcg equation solution..."
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
      u_thread(local(1:num_Loc_ESM),c_thread) = u_thread(local(1:num_Loc_ESM),c_thread)+ &
                  MATMUL(storK(1:num_Loc_ESM,1:num_Loc_ESM,i_E),p(local(1:num_Loc_ESM)))
  enddo
  !$omp end do
  !$omp end parallel

  DO i_Thread = 1,omp_get_max_threads()
      u  =  u  + u_thread(:,i_Thread)
  ENDDO

  up=DOT_PRODUCT(loads,pcg_d)
  alpha=up/DOT_PRODUCT(p,u)
  xnew=x+p*alpha
  loads=loads-u*alpha
  pcg_d=diag_precon*loads
  beta=DOT_PRODUCT(loads,pcg_d)/up
  p=pcg_d+p*beta


  tem_value = MAXVAL(ABS(x(0:num_FreeD)))
  if(tem_value==ZR) tem_value = Tol_20
  Tol = MAXVAL(ABS(x(0:num_FreeD)-xnew(0:num_FreeD)))/tem_value

  x=xnew
  if(Tol<c_cg_tol.OR.cg_iters==max_num_PCG)then
      exit
  endif
enddo
print *, "    Number of CG iterations to convergence was",cg_iters

loads=xnew
disp = ZR
print *, "    Step 5: retrive nodal displacement..."
disp(freeDOF(1:num_FreeD)) =loads(1:num_FreeD)
DEALLOCATE(p,loads,x,xnew,u,pcg_d)
print *, "    >>>>  End of element by element PCG solver  <<<<"
RETURN
END SUBROUTINE Ele_by_Ele_XFEM_PCG
