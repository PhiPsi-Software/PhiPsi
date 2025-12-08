 
SUBROUTINE Ele_by_Ele_XFEM_PCG_3D(isub,Lambda,c_cg_tol,max_num_PCG,num_FreeD,freeDOF,F,disp, &
                                  diag_precon_no_invert)

use Global_Float_Type
use Global_Model
use Global_Filename
use Global_Common
use Global_Material
use Global_Crack_Common
use Global_Crack_3D
use Function_MATMUL_LP 
use Global_XFEM_Elements
use Global_POST
use omp_lib
use Global_Cal_Ele_Stiffness_Matrix_3D_8nodes
use Global_Inter_Location_Element_Stiff_Matrix_3D
implicit none
integer,intent(in)::isub,max_num_PCG,num_FreeD
real(kind=FT),intent(in)::Lambda,c_cg_tol,F(num_FreeD)
integer,intent(in)::freeDOF(num_FreeD)
real(kind=FT),intent(out)::disp(Total_FD)
real(kind=FT),intent(out)::diag_precon_no_invert(0:num_FreeD)
integer jj,kk
integer i_E,i_G,c_NN(8),mat_num,k
real(kind=FT) c_D(6,6)
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
real(kind=FT) localK(MDOF_3D,MDOF_3D)
real(kind=FT) Diag_localK(MDOF_3D)
integer local(MDOF_3D)
real(kind=FT),ALLOCATABLE::p(:),loads(:),x(:),xnew(:),u(:), diag_precon(:),pcg_d(:)
integer cg_iters,i_PCG
real(kind=FT) alpha,beta,up,Tol
real(kind=FT) Rot_c_D_Comp(6,6),c_D_Comp(6,6),Volume_Ratio
real(kind=FT) T_Matrix(6,6),TT_Matrix(6,6)
real(kind=FT) kesi_Enr(Num_Gau_Points_3D),yita_Enr(Num_Gau_Points_3D),zeta_Enr(Num_Gau_Points_3D),weight_Enr(Num_Gau_Points_3D)
real(kind=FT) kesi_Enr_MC(Num_Gau_Points_3D_MC),yita_Enr_MC(Num_Gau_Points_3D_MC), &
              zeta_Enr_MC(Num_Gau_Points_3D_MC),weight_Enr_MC(Num_Gau_Points_3D_MC)
real(kind=FT) kesi_No_Enr(Num_Gauss_P_FEM_3D),yita_No_Enr(Num_Gauss_P_FEM_3D),&
              zeta_No_Enr(Num_Gauss_P_FEM_3D),weight_No_Enr(Num_Gauss_P_FEM_3D)

real(kind=FT),ALLOCATABLE::kesi_Enr_Cubes(:),yita_Enr_Cubes(:),zeta_Enr_Cubes(:),weight_Enr_Cubes(:)
integer c_Num_Gauss_Point,i_C
real(kind=FT),ALLOCATABLE::kesi(:),yita(:),zeta(:),weight(:)
integer:: Location_ESM(MDOF_3D)
integer::Location_ESM_C_Crack(MDOF_3D),Location_ESM_C_Cr_NoFEM(MDOF_3D)
integer num_Loc_ESM_C_Cr_NoFEM
real(kind=FT) B(6,MDOF_3D),tem_B(6,MDOF_3D),detJ
integer num_B,num_tem_B
integer num_Loc_ESM,num_Loc_ESM_C_Crack
integer c_loca,j_C,i_Node
integer cEle_Loc,i_N
integer Yes_Node_Dealed(num_Node)
integer Yes_DOF_Dealed(Total_FD)
integer num_t,num_h,num_j,cnt,tem,i_f
integer c_LocXFEM(3),c_DOFs(3)
integer i,j
real(kind=FT) c_Beta(3) 
real(kind=FT) c_Penalty_CS
logical Need_XFEM_Assem_Flag
real(kind=FT) c_Vector(MDOF_3D*MDOF_3D)
integer c_Rec
real(kind=FT) DOT_PRODUCT_pu,DOT_PRODUCT_lp
integer c_Elem
integer ndof,i_ndof
real(kind=FT) dot_p_u,dot_l_p
real(kind=FT),ALLOCATABLE:: Vector_x_xnew(:)
real(kind=FT),ALLOCATABLE::u_thread(:,:)
real(kind=FT),ALLOCATABLE::diag_precon_thread(:,:)
integer i_Thread,c_Thread,max_threads
integer c_POS_3D_c_Ele(8)
integer c_i,c_j
real(kind=FT),ALLOCATABLE::localK_FEM(:,:)
real(kind=FT),ALLOCATABLE::localK_XFEM(:,:)
integer date_time(8)
integer(LIT) c_S_Time,F_time
character*10  current_data
5003 FORMAT('     Elapsed CPU time of EBE-PCG solution - ',I8,' s, about ',F10.4,' mins')

print *, "    >>>> Start of element by element PCG solver <<<<"
print *,'    PCG-EBE: prepare data...' 

ALLOCATE(p(0:num_FreeD),loads(0:num_FreeD),x(0:num_FreeD),xnew(0:num_FreeD),u(0:num_FreeD), &
           diag_precon(0:num_FreeD),pcg_d(0:num_FreeD))   

ALLOCATE(Vector_x_xnew(0:num_FreeD))

max_threads = omp_get_max_threads()
if (allocated(u_thread)) deallocate(u_thread)
ALLOCATE(u_thread(0:num_FreeD,max_threads))
if(Key_EBE_Precondition == 1)then
    ALLOCATE(diag_precon_thread(0:num_FreeD,max_threads))
endif


print *, "    PCG-EBE: get and store element stiffness matrix..."
call Cal_Gauss_Points_3D_8nodes(Num_Gau_Points_3D,kesi_Enr,yita_Enr,zeta_Enr,weight_Enr)      
call Cal_Gauss_Points_3D_8nodes(Num_Gau_Points_3D_MC,kesi_Enr_MC,yita_Enr_MC,zeta_Enr_MC,weight_Enr_MC)
call Cal_Gauss_Points_3D_8nodes(Num_Gauss_P_FEM_3D,kesi_No_Enr,yita_No_Enr,zeta_No_Enr,weight_No_Enr)
 
if (Key_Integral_Sol  == 3) then     
    allocate(kesi_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    allocate(yita_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    allocate(zeta_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    allocate(weight_Enr_Cubes(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
    call Cal_Gauss_Points_3D_for_SUBCUBES(Num_Sub_3D_Cubes,Num_Gau_Points_3D_Cube,&
                                          kesi_Enr_Cubes,yita_Enr_Cubes,zeta_Enr_Cubes,weight_Enr_Cubes)                                
endif
      
      
EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Max_3D_gauss)  = .false.
if (isub==1) then
    if(Key_EBE_Sym_Storage_K==0)then
        storK_FEM(1:24,1:24,1:num_FEM_Elem) =ZR
    elseif(Key_EBE_Sym_Storage_K==1) then
        storK_FEM_Sym(1:300,1:num_FEM_Elem)     =ZR
    endif 
    
    All_Local_3D(1:MDOF_3D,1:Num_Elem) = 0
    print *, '    PCG-EBE: First step may take more CPU time...'
endif
      
if (isub==1 .or. First_XFEM_Step==1) then
    if (First_XFEM_Step ==1) then
        First_XFEM_Step = -1
    endif
    
    print *,'    PCG-EBE: Assembling K for FEM elements...'    
    
    if(allocated(kesi)) deallocate(kesi)  
    if(allocated(yita)) deallocate(yita)
    if(allocated(zeta)) deallocate(zeta)
    if(allocated(weight)) deallocate(weight)
    allocate(kesi(Num_Gauss_P_FEM_3D))
    allocate(yita(Num_Gauss_P_FEM_3D))
    allocate(zeta(Num_Gauss_P_FEM_3D))
    allocate(weight(Num_Gauss_P_FEM_3D))              
    kesi(1:Num_Gauss_P_FEM_3D)   = kesi_No_Enr
    yita(1:Num_Gauss_P_FEM_3D)   = yita_No_Enr
    zeta(1:Num_Gauss_P_FEM_3D)   = zeta_No_Enr
    weight(1:Num_Gauss_P_FEM_3D) = weight_No_Enr
    c_Num_Gauss_Point            = Num_Gauss_P_FEM_3D 
    
    ndof = 3*8
    allocate(localK_FEM(ndof,ndof))
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_Elem,mat_num,c_D,Volume_Ratio,c_X_NODES,c_D_Comp,T_Matrix,TT_Matrix,       &
    !$OMP            Rot_c_D_Comp,c_Y_NODES,c_Z_NODES,c_NN,Location_ESM,k,c_loca,local,localK_FEM,num_Loc_ESM,c_i,c_j)    &
    !$OMP            SCHEDULE(static)                                     
    do i_E = 1,num_FEM_Elem
          c_Elem = FEM_Elem_List(i_E)
          
          if(Key_XA/=2) then
              mat_num = Elem_Mat(c_Elem)
              c_D(1:6,1:6)     = D(Elem_Mat(c_Elem),1:6,1:6)        
              if (Material_Type(mat_num)==5)then
                  Volume_Ratio = Material_Para_Added(mat_num,10)
                  c_D_Comp = D_Comp(mat_num,1:6,1:6)
                  T_Matrix = Ele_ComMat_RotMatrix(c_Elem,1:6,1:6)
                  TT_Matrix= TRANSPOSE(T_Matrix)
                  Rot_c_D_Comp = MATMUL(TT_Matrix,c_D_Comp)
                  Rot_c_D_Comp = MATMUL(Rot_c_D_Comp,T_Matrix)
                  c_D =(ONE-Volume_Ratio)*c_D + Volume_Ratio*Rot_c_D_Comp
              endif
          else
              c_D       = Elem_D_XA(c_Elem,1:6,1:6)
          endif
          
          c_NN    = G_NN(1:8,c_Elem)
          c_X_NODES = G_X_NODES(1:8,c_Elem)
          c_Y_NODES = G_Y_NODES(1:8,c_Elem)    
          c_Z_NODES = G_Z_NODES(1:8,c_Elem)  
          if(Key_Post_Elements_Gauss_Num==1) then
              Elements_Gauss_Num(c_Elem)      =  c_Num_Gauss_Point
          endif
          
          num_Loc_ESM = ndof
          
          Location_ESM(1:num_Loc_ESM)=[c_NN(1)*3-2,c_NN(1)*3-1,c_NN(1)*3,c_NN(2)*3-2,c_NN(2)*3-1,c_NN(2)*3,  &
                                       c_NN(3)*3-2,c_NN(3)*3-1,c_NN(3)*3,c_NN(4)*3-2,c_NN(4)*3-1,c_NN(4)*3,  &
                                       c_NN(5)*3-2,c_NN(5)*3-1,c_NN(5)*3,c_NN(6)*3-2,c_NN(6)*3-1,c_NN(6)*3,  &
                                       c_NN(7)*3-2,c_NN(7)*3-1,c_NN(7)*3,c_NN(8)*3-2,c_NN(8)*3-1,c_NN(8)*3] 
                                
          Size_Local_3D(c_Elem)  = num_Loc_ESM
          
          DO k=1,num_Loc_ESM 
              if (Flag_FreeDOF(Location_ESM(k)) == 1 )then
                  
                  c_loca = Location_FreeDOF(Location_ESM(k))
              else
                  c_loca =0
              endif               
              local(k) = c_loca
          enddo
          All_Local_3D(1:num_Loc_ESM,c_Elem) = local(1:num_Loc_ESM)
          localK_FEM = ZR
          
          call Cal_Ele_Stiffness_Matrix_3D_8nodes(c_Elem,Num_Gauss_P_FEM_3D,c_X_NODES,c_Y_NODES,c_Z_NODES, &
                                c_D,kesi,yita,zeta,weight,localK_FEM(1:num_Loc_ESM,1:num_Loc_ESM))  

          
          if(Key_EBE_Sym_Storage_K==0)then
              storK_FEM(1:ndof,1:ndof,Elem_Location(c_Elem,2))=localK_FEM(1:ndof,1:ndof)
          elseif(Key_EBE_Sym_Storage_K==1)then
              do c_i=1,ndof
                  do c_j=c_i,ndof
                      storK_FEM_Sym((c_i-1)*24 -(c_i-1)*c_i/2 +c_j,Elem_Location(c_Elem,2))  = localK_FEM(c_i,c_j)
                  enddo
              enddo
          endif
          
    enddo
    !$omp end parallel do  
    deallocate(kesi)  
    deallocate(yita)
    deallocate(zeta)
    deallocate(weight)      
    deallocate(localK_FEM)
endif


if (isub>=2 .and. Num_Rollbacked_FEM_Elements >0) then      
    print *,'    PCG-EBE: Assembling K for rollbacked FEM elements...'
    ndof = 3*8
    allocate(localK_FEM(ndof,ndof))
    allocate(kesi(Num_Gauss_P_FEM_3D))
    allocate(yita(Num_Gauss_P_FEM_3D))
    allocate(zeta(Num_Gauss_P_FEM_3D))
    allocate(weight(Num_Gauss_P_FEM_3D))              
    kesi(1:Num_Gauss_P_FEM_3D)   = kesi_No_Enr
    yita(1:Num_Gauss_P_FEM_3D)   = yita_No_Enr
    zeta(1:Num_Gauss_P_FEM_3D)   = zeta_No_Enr
    weight(1:Num_Gauss_P_FEM_3D) = weight_No_Enr
    c_Num_Gauss_Point            = Num_Gauss_P_FEM_3D 
    

    
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_Elem,c_D,Volume_Ratio,c_D_Comp,T_Matrix,TT_Matrix,  &
    !$OMP            Rot_c_D_Comp,Location_ESM,k,c_loca,local,localK_fem,c_i,c_j)    &
    !$OMP            SCHEDULE(static)                                     
    do i_E = 1,Num_Rollbacked_FEM_Elements
          c_Elem = Rollbacked_FEM_Elements(i_E)
          
          if(Key_XA/=2) then
              c_D(1:6,1:6)     = D(Elem_Mat(c_Elem),1:6,1:6)        
              if (Material_Type(Elem_Mat(c_Elem))==5)then
                  Volume_Ratio = Material_Para_Added(Elem_Mat(c_Elem),10)
                  c_D_Comp = D_Comp(Elem_Mat(c_Elem),1:6,1:6)
                  T_Matrix = Ele_ComMat_RotMatrix(c_Elem,1:6,1:6)
                  TT_Matrix= TRANSPOSE(T_Matrix)
                  Rot_c_D_Comp = MATMUL(TT_Matrix,c_D_Comp)
                  Rot_c_D_Comp = MATMUL(Rot_c_D_Comp,T_Matrix)
                  c_D =(ONE-Volume_Ratio)*c_D + Volume_Ratio*Rot_c_D_Comp
              endif   
          else
              c_D       = Elem_D_XA(c_Elem,1:6,1:6)
          endif
          
          Location_ESM(1:ndof)=&
            [G_NN(1,c_Elem)*3-2,G_NN(1,c_Elem)*3-1,G_NN(1,c_Elem)*3,G_NN(2,c_Elem)*3-2,G_NN(2,c_Elem)*3-1,G_NN(2,c_Elem)*3,  &
             G_NN(3,c_Elem)*3-2,G_NN(3,c_Elem)*3-1,G_NN(3,c_Elem)*3,G_NN(4,c_Elem)*3-2,G_NN(4,c_Elem)*3-1,G_NN(4,c_Elem)*3,  &
             G_NN(5,c_Elem)*3-2,G_NN(5,c_Elem)*3-1,G_NN(5,c_Elem)*3,G_NN(6,c_Elem)*3-2,G_NN(6,c_Elem)*3-1,G_NN(6,c_Elem)*3,  &
             G_NN(7,c_Elem)*3-2,G_NN(7,c_Elem)*3-1,G_NN(7,c_Elem)*3,G_NN(8,c_Elem)*3-2,G_NN(8,c_Elem)*3-1,G_NN(8,c_Elem)*3] 
                                
          Size_Local_3D(c_Elem)  = ndof
          
          DO k=1,ndof 
              if (Flag_FreeDOF(Location_ESM(k)) == 1 )then
                  
                  c_loca = Location_FreeDOF(Location_ESM(k))
              else
                  c_loca =0
              endif               
              local(k) = c_loca
          enddo
          All_Local_3D(1:ndof,c_Elem) = local(1:ndof)
          call Cal_Ele_Stiffness_Matrix_3D_8nodes(c_Elem,Num_Gauss_P_FEM_3D,&
                                G_X_NODES(1:8,c_Elem),G_Y_NODES(1:8,c_Elem),G_Z_NODES(1:8,c_Elem), &
                                c_D,kesi,yita,zeta,weight,localK_fem(1:ndof,1:ndof))
                                
          if(Key_EBE_Sym_Storage_K==0)then
              storK_FEM(1:ndof,1:ndof,Elem_Location(c_Elem,2))=localK_fem(1:ndof,1:ndof)
          elseif(Key_EBE_Sym_Storage_K==1)then
              do c_i=1,ndof
                  do c_j=c_i,ndof
                      storK_FEM_Sym((c_i-1)*24 -(c_i-1)*c_i/2 +c_j,Elem_Location(c_Elem,2))  = localK_fem(c_i,c_j)
                  enddo
              enddo
          endif
    enddo
    !$omp end parallel do  
    
    if(Key_Save_Nothing==0) then
        if(Key_Post_Elements_Gauss_Num==1) then
            Elements_Gauss_Num(1:num_FEM_Elem)      =  c_Num_Gauss_Point
        endif
    endif
    
    deallocate(kesi)  
    deallocate(yita)
    deallocate(zeta)
    deallocate(weight)    
    deallocate(localK_FEM)
endif


print *,'    PCG-EBE: Assembling K for XFEM elements...'



!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(c_Elem,                                                               &
!$OMP            i_E,mat_num,c_D,Volume_Ratio,c_X_NODES,c_D_Comp,T_Matrix,TT_Matrix,Rot_c_D_Comp,               &
!$OMP            c_Y_NODES,c_Z_NODES,c_NN,kesi,yita,zeta,weight,Location_ESM_C_Crack,num_Loc_ESM_C_Crack,       &
!$OMP            Location_ESM_C_Cr_NoFEM,c_Num_Gauss_Point,i_C,j_C,Location_ESM,num_Loc_ESM,B,                  &
!$OMP            num_Loc_ESM_C_Cr_NoFEM,i_G,k,c_loca,localK_XFEM,num_B,tem_B,num_tem_B,detJ,local,              &
!$OMP            cEle_Loc,Need_XFEM_Assem_Flag,c_Rec,c_Vector,c_POS_3D_c_Ele)                                   &
!$OMP            SCHEDULE(static)                                    

do i_E = 1,num_XFEM_Elem
      allocate(localK_XFEM(MDOF_3D,MDOF_3D))
      c_Elem = XFEM_Elem_List(i_E)
      if(isub==1) then
          Need_XFEM_Assem_Flag = .True.
      endif
      if(isub>=2) then
          Need_XFEM_Assem_Flag = .False.
          if(Elem_New_XFEM_Flag(c_Elem) == 1)then  
              Need_XFEM_Assem_Flag = .True.
          endif
      endif
      

      if(Key_XA/=2) then
          mat_num = Elem_Mat(c_Elem)
          c_D(1:6,1:6)     = D(Elem_Mat(c_Elem),1:6,1:6)        
          if (Material_Type(mat_num)==5)then
              Volume_Ratio = Material_Para_Added(mat_num,10)
              c_D_Comp = D_Comp(mat_num,1:6,1:6)
              T_Matrix = Ele_ComMat_RotMatrix(c_Elem,1:6,1:6)
              TT_Matrix= TRANSPOSE(T_Matrix)
              Rot_c_D_Comp = MATMUL(TT_Matrix,c_D_Comp)
              Rot_c_D_Comp = MATMUL(Rot_c_D_Comp,T_Matrix)
              c_D =(ONE-Volume_Ratio)*c_D + Volume_Ratio*Rot_c_D_Comp
          endif 
      else
          c_D       = Elem_D_XA(c_Elem,1:6,1:6)
      endif
      
      c_NN    = G_NN(1:8,c_Elem)
      c_X_NODES = G_X_NODES(1:8,c_Elem)
      c_Y_NODES = G_Y_NODES(1:8,c_Elem)    
      c_Z_NODES = G_Z_NODES(1:8,c_Elem)  
      
      select case(Key_Integral_Sol)
      case(2)
        if (num_Crack/= 0 .and. (maxval(Enriched_Node_Type_3D(c_NN,1:num_Crack)).ne.0))then
              if(allocated(kesi)) deallocate(kesi)
              if(allocated(yita)) deallocate(yita)
              if(allocated(zeta)) deallocate(zeta)
              if(allocated(weight)) deallocate(weight)              
              allocate(kesi(Num_Gau_Points_3D))
              allocate(yita(Num_Gau_Points_3D))
              allocate(zeta(Num_Gau_Points_3D))
              allocate(weight(Num_Gau_Points_3D))
              kesi(1:Num_Gau_Points_3D)    = kesi_Enr
              yita(1:Num_Gau_Points_3D)    = yita_Enr
              zeta(1:Num_Gau_Points_3D)    = zeta_Enr
              weight(1:Num_Gau_Points_3D)  = weight_Enr
              c_Num_Gauss_Point            = Num_Gau_Points_3D
        else 
              if(allocated(kesi)) deallocate(kesi)
              if(allocated(yita)) deallocate(yita)
              if(allocated(zeta)) deallocate(zeta)
              if(allocated(weight)) deallocate(weight)
              allocate(kesi(Num_Gauss_P_FEM_3D))
              allocate(yita(Num_Gauss_P_FEM_3D))
              allocate(zeta(Num_Gauss_P_FEM_3D))
              allocate(weight(Num_Gauss_P_FEM_3D))              
              kesi(1:Num_Gauss_P_FEM_3D)   = kesi_No_Enr
              yita(1:Num_Gauss_P_FEM_3D)   = yita_No_Enr
              zeta(1:Num_Gauss_P_FEM_3D)   = zeta_No_Enr
              weight(1:Num_Gauss_P_FEM_3D) = weight_No_Enr
              c_Num_Gauss_Point            = Num_Gauss_P_FEM_3D
        end if
                    
        if(Elem_num_Related_Cracks(c_Elem)>=2) then  
              if(allocated(kesi)) deallocate(kesi)
              if(allocated(yita)) deallocate(yita)
              if(allocated(zeta)) deallocate(zeta)
              if(allocated(weight)) deallocate(weight)
              allocate(kesi(Num_Gau_Points_3D_MC))
              allocate(yita(Num_Gau_Points_3D_MC))
              allocate(zeta(Num_Gau_Points_3D_MC))
              allocate(weight(Num_Gau_Points_3D_MC))    
              kesi(1:Num_Gau_Points_3D_MC)    = kesi_Enr_MC
              yita(1:Num_Gau_Points_3D_MC)    = yita_Enr_MC
              zeta(1:Num_Gau_Points_3D_MC)    = zeta_Enr_MC
              weight(1:Num_Gau_Points_3D_MC)  = weight_Enr_MC
              c_Num_Gauss_Point               = Num_Gau_Points_3D_MC
        endif
      case(3)
        if (num_Crack/= 0 .and. (maxval(Enriched_Node_Type_3D(c_NN,1:num_Crack)).ne.0))then
              if(allocated(kesi)) deallocate(kesi) 
              if(allocated(yita)) deallocate(yita)
              if(allocated(zeta)) deallocate(zeta)
              if(allocated(weight)) deallocate(weight)              
              allocate(kesi(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
              allocate(yita(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
              allocate(zeta(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
              allocate(weight(Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube))
              kesi(1:Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube)    = kesi_Enr_Cubes
              yita(1:Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube)    = yita_Enr_Cubes
              zeta(1:Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube)    = zeta_Enr_Cubes
              weight(1:Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube)  = weight_Enr_Cubes
              c_Num_Gauss_Point            = Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube
        else 
              if(allocated(kesi)) deallocate(kesi)  
              if(allocated(yita)) deallocate(yita)
              if(allocated(zeta)) deallocate(zeta)
              if(allocated(weight)) deallocate(weight)
              allocate(kesi(Num_Gauss_P_FEM_3D))
              allocate(yita(Num_Gauss_P_FEM_3D))
              allocate(zeta(Num_Gauss_P_FEM_3D))
              allocate(weight(Num_Gauss_P_FEM_3D))              
              kesi(1:Num_Gauss_P_FEM_3D)   = kesi_No_Enr
              yita(1:Num_Gauss_P_FEM_3D)   = yita_No_Enr
              zeta(1:Num_Gauss_P_FEM_3D)   = zeta_No_Enr
              weight(1:Num_Gauss_P_FEM_3D) = weight_No_Enr
              c_Num_Gauss_Point            = Num_Gauss_P_FEM_3D
        end if
        
      case(4)
      end select

      if(Key_Post_Elements_Gauss_Num==1) then
          Elements_Gauss_Num(c_Elem)      =  c_Num_Gauss_Point
      endif
      
      Location_ESM(1:MDOF_3D)   = 0
      num_Loc_ESM               = 0
      
      do i_C =1,num_Crack
          c_POS_3D_c_Ele(1:8) = c_POS_3D(c_NN,i_C)
          
          if(i_C >1 .and. sum(c_POS_3D_c_Ele(1:8))==0) cycle
          
          call Location_Element_Stiff_Matrix_3D(c_Elem,i_C,c_POS_3D_c_Ele(1:8),Location_ESM_C_Crack,    &
                                     num_Loc_ESM_C_Crack,Location_ESM_C_Cr_NoFEM,num_Loc_ESM_C_Cr_NoFEM)
          Location_ESM(num_Loc_ESM+1:num_Loc_ESM+num_Loc_ESM_C_Crack) = Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
          num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack
      end do
          
       
          
          
      Size_Local_3D(c_Elem)  = num_Loc_ESM
      
      
      if(isub>=2)then
          if (allocated(size_local_Old))then
              if(Size_Local_3D(c_Elem) /= size_local_Old(c_Elem)) then
                  Need_XFEM_Assem_Flag =.True.
              endif    
          endif
      endif
      
      local(1:MDOF_3D) = 0
      
      DO k=1,num_Loc_ESM
          if (Flag_FreeDOF(Location_ESM(k)) == 1 )then
              
              c_loca = Location_FreeDOF(Location_ESM(k))
          else
              c_loca =0
          endif 
          local(k) = c_loca              
      enddo
                  
      All_Local_3D(1:num_Loc_ESM,c_Elem) = local(1:num_Loc_ESM)
      localK_XFEM(1:MDOF_3D,1:MDOF_3D) = ZR
      
      if (Need_XFEM_Assem_Flag .eqv. .True.) then  

        select case(Key_Integral_Sol)
        case(2)
            do i_G = 1,c_Num_Gauss_Point
              B(1:6,1:MDOF_3D) = ZR
              num_B = 0
              call Cal_detJ_3D(kesi(i_G),yita(i_G),zeta(i_G),c_X_NODES,c_Y_NODES,c_Z_NODES,detJ)  
              if(num_Crack/=0)then
                  do j_C =1,num_Crack 
                      c_POS_3D_c_Ele(1:8) = c_POS_3D(c_NN,j_C)
                      if(j_C >1 .and. sum(c_POS_3D_c_Ele(1:8))==0) cycle
                      
                      
                      tem_B(1:6,1:MDOF_3D) = ZR
                      call Cal_B_Matrix_Crack_3D(kesi(i_G),yita(i_G),zeta(i_G),j_C,c_Elem,i_G,c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES, &
                                   tem_B,num_tem_B)     
                      B(1:6,num_B+1:num_B+num_tem_B)=tem_B(1:6,1:num_tem_B)
                      num_B = num_B + num_tem_B                  
                  end do
              endif

              localK_XFEM(1:num_B,1:num_B) = localK_XFEM(1:num_B,1:num_B) + weight(i_G)*detJ* &
                    MATMUL(MATMUL(transpose(B(1:6,1:num_B)),c_D),B(1:6,1:num_B))     
            end do
        case(3)
            do i_G = 1,c_Num_Gauss_Point
              B(1:6,1:MDOF_3D) = ZR
              num_B = 0
              call Cal_detJ_3D(kesi(i_G),yita(i_G),zeta(i_G),c_X_NODES,c_Y_NODES,c_Z_NODES,detJ)  
              if(num_Crack/=0)then
                  do j_C =1,num_Crack 
                      c_POS_3D_c_Ele(1:8) = c_POS_3D(c_NN,j_C)
                      if(j_C >1 .and. sum(c_POS_3D_c_Ele(1:8))==0) cycle
                      
                      tem_B(1:6,1:MDOF_3D) = ZR
                      call Cal_B_Matrix_Crack_3D(kesi(i_G),yita(i_G),zeta(i_G),j_C,c_Elem,i_G,c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES, &
                                   tem_B,num_tem_B)     
                      B(1:6,num_B+1:num_B+num_tem_B)=tem_B(1:6,1:num_tem_B)
                      num_B = num_B + num_tem_B                  
                  end do
              endif

              localK_XFEM(1:num_B,1:num_B) = localK_XFEM(1:num_B,1:num_B) + weight(i_G)*detJ* &
                    MATMUL(MATMUL(transpose(B(1:6,1:num_B)),c_D),B(1:6,1:num_B))     
            end do
        case(4)

        end select

      elseif (Need_XFEM_Assem_Flag .eqv. .False.) then
          c_Rec = Elem_Location_old(c_Elem,1)
          
          num_B = size_local_Old(c_Elem)
          localK_XFEM(1:num_B,1:num_B) = storK_XFEM_Old(c_Rec)%row(1:num_B,1:num_B)
          
          
      endif
   
      cEle_Loc = Elem_Location(c_Elem,1) 
      if(.not. allocated(storK_XFEM(cEle_Loc)%row)) then
          allocate(storK_XFEM(cEle_Loc)%row(num_B,num_B))
      endif
      storK_XFEM(cEle_Loc)%row(1:num_B,1:num_B)=localK_XFEM(1:num_B,1:num_B)
      deallocate(localK_XFEM)
      
enddo
!$omp end parallel do   

if (allocated(size_local_Old))DEALLOCATE(size_local_Old)
allocate(size_local_Old(num_Elem))
size_local_Old = 0
size_local_Old = Size_Local_3D
      
115 FORMAT(8X,'Element ',I6,' condition number:',E12.5)  
if (Key_EBE_Condition_Number == 1)then
  print *,'            Calculating condition number of XFEM elements...'
  if (allocated(EBE_Condition_Number))then
      DEALLOCATE(EBE_Condition_Number)
  endif
  allocate(EBE_Condition_Number(num_Elem))
  EBE_Condition_Number(1:Num_Elem) = ZR
  do i_E=1,Num_Elem
      if (Elem_XFEM_Flag(i_E) ==1)then  
          cEle_Loc = Elem_Location(i_E,1) 
          num_Loc_ESM = Size_Local_3D(i_E)
          call Matrix_Condition_Number_v2(2,num_Loc_ESM,          &
            storK_XFEM(cEle_Loc)%row(1:num_Loc_ESM,1:num_Loc_ESM),EBE_Condition_Number(i_E)) 
          write(*,115) i_E,EBE_Condition_Number(i_E)
      endif
  enddo           
endif 


if ( .not. any(Key_CS_Crack(:) == 1))then  
    goto 1000
endif

if (Key_CS_Natural_Crack/=1) then  
    goto 1000
endif

Yes_Node_Dealed(1:num_Node) =0
Yes_DOF_Dealed(1:Total_FD) = 0
print *,"    PCG-EBE: get c_Penalty_CS for Penalty function treatment..." 


c_Penalty_CS =  Penalty_CS_Natural_Crack  

print *,"    PCG-EBE: Penalty function treatment for XFEM elements..."   

!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,i_E,c_Elem,cEle_Loc,c_NN,cnt,c_LocXFEM, &
!$OMP         i_N,c_DOFs,i_f,c_Beta,i,j) &
!$OMP SCHEDULE(static) 
do i_C=1,num_Crack
  if(Key_CS_Crack(i_C)==1)then
    do i_E = 1,num_XFEM_Elem
      c_Elem = XFEM_Elem_List(i_E)
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
              if (Yes_DOF_Dealed(c_DOFs(1))==0) then
                  Yes_DOF_Dealed(c_DOFs(1))=1
                  c_Beta(1:3)=Enriched_Node_Crack_n_Vector_3D(c_NN(i_N))%row(i_C,1:3)
                  if(Key_Penalty_CS_Method==1) then
                      do j=1,3
                        do i=1,3
                            storK_XFEM(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) = &
                                  storK_XFEM(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) &
                                    + c_Penalty_CS*c_Beta(i)*c_Beta(j)
                        enddo
                      enddo
                  elseif(Key_Penalty_CS_Method==2) then
                      do j=3,3
                        do i=3,3
                            storK_XFEM(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) = &
                                  storK_XFEM(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) &
                                    + c_Penalty_CS*c_Beta(i)*c_Beta(j)
                        enddo
                      enddo
                  endif             
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
                  c_Beta(1:3)=Enriched_Node_Crack_n_Vector_3D(c_NN(i_N))%row(i_C,1:3)
                  if(Key_Penalty_CS_Method==1) then
                      do j=1,3
                        do i=1,3
                            storK_XFEM(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) = &
                                  storK_XFEM(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) &
                                    + c_Penalty_CS*c_Beta(i)*c_Beta(j)
                        enddo
                      enddo
                  elseif(Key_Penalty_CS_Method==2) then
                      do j=3,3
                        do i=3,3
                            storK_XFEM(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) = &
                                  storK_XFEM(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) &
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
              if (Yes_DOF_Dealed(c_DOFs(1))==0) then
                  Yes_DOF_Dealed(c_DOFs(1))=1
                  c_Beta(1:3)=Enriched_Node_Crack_n_Vector_3D(c_NN(i_N))%row(i_C,1:3)
                  if(Key_Penalty_CS_Method==1) then
                      do j=1,3
                        do i=1,3
                            storK_XFEM(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) = &
                                  storK_XFEM(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) &
                                    + c_Penalty_CS*c_Beta(i)*c_Beta(j)
                        enddo
                      enddo
                  elseif(Key_Penalty_CS_Method==2) then
                      do j=3,3
                        do i=3,3
                            storK_XFEM(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) = &
                                  storK_XFEM(cEle_Loc)%row(c_LocXFEM(i),c_LocXFEM(j)) &
                                    + c_Penalty_CS*c_Beta(i)*c_Beta(j)
                        enddo
                      enddo
                  endif             
              endif
              
            else
                cycle
            end if
            
          end do
    enddo
  endif
enddo
!$OMP END PARALLEL do 

1000 continue
      
if(Key_EBE_Precondition == 1)then
  diag_precon(0:num_FreeD)= ZR      
  print *,"    PCG-EBE: geting preconditioner for FEM elements..."          
  diag_precon_thread(0:num_FreeD,1:max_threads)= ZR 
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,c_Elem,i_E,num_Loc_ESM,local,localK,kk,Diag_localK) 
  c_thread = omp_get_thread_num()+1
  !$OMP DO SCHEDULE(static) 
  do i_E =1,num_FEM_Elem
      c_Elem = FEM_Elem_List(i_E)
      num_Loc_ESM = Size_Local_3D(c_Elem)
      if(num_Loc_ESM>24) then 
          print *,'    ERROR-2023021501 :: num_Loc_ESM > 24 for FEM element!'
          print *,'                        In EBE_XFEM_PCG_3D_Get_K.f90!'
          print *,'                        num_Loc_ESM:',num_Loc_ESM
          print *,'                        c_Elem:',c_Elem
          call Warning_Message('S',Keywords_Blank)
      endif
      
      local(1:num_Loc_ESM) = All_Local_3D(1:num_Loc_ESM,c_Elem) 
      
      if(Key_EBE_Sym_Storage_K==0)then
          DO kk=1,num_Loc_ESM     
              Diag_localK(kk) =storK_FEM(kk,kk,Elem_Location(c_Elem,2)) 
          enddo       
      elseif(Key_EBE_Sym_Storage_K==1)then
          DO kk=1,num_Loc_ESM     
              Diag_localK(kk) =storK_FEM_Sym((kk-1)*24 -(kk-1)*kk/2 +kk,Elem_Location(c_Elem,2)) 
          enddo   
      endif
      
      diag_precon_thread(local(1:num_Loc_ESM),c_thread)  =   &
                 diag_precon_thread(local(1:num_Loc_ESM),c_thread) + Diag_localK(1:num_Loc_ESM)       
  enddo
  !$omp end do
  !$omp end parallel    
  
  DO i_Thread = 1,omp_get_max_threads()
     diag_precon(0:num_FreeD)  =  diag_precon(0:num_FreeD)  + diag_precon_thread(0:num_FreeD,i_Thread)
  ENDDO
  
  
  
  print *,"    PCG-EBE: geting preconditioner for XFEM elements..."          
  diag_precon_thread(0:num_FreeD,1:max_threads)= ZR 
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,c_Elem,i_E,num_Loc_ESM,local,localK,kk,Diag_localK) 
  c_thread = omp_get_thread_num()+1
  !$OMP DO SCHEDULE(static) 
  do i_E=1,num_XFEM_Elem
      c_Elem = XFEM_Elem_List(i_E)
      num_Loc_ESM = Size_Local_3D(c_Elem)
      local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
      DO kk=1,num_Loc_ESM     
          Diag_localK(kk) =storK_XFEM(Elem_Location(c_Elem,1))%row(kk,kk) 
          
      enddo         
      diag_precon_thread(local(1:num_Loc_ESM),c_thread)  =   &
                 diag_precon_thread(local(1:num_Loc_ESM),c_thread) + Diag_localK(1:num_Loc_ESM)       
  enddo
  !$omp end do
  !$omp end parallel    

  DO i_Thread = 1,omp_get_max_threads()
     diag_precon(0:num_FreeD)  =  diag_precon(0:num_FreeD)  + diag_precon_thread(0:num_FreeD,i_Thread)
  ENDDO  
  
  if(allocated(diag_precon_thread)) deallocate(diag_precon_thread)  
      
  forall(jj=1:num_FreeD,abs(diag_precon(jj))<=Tol_10)
      diag_precon(jj)=Tol_10
  endforall
endif
      

diag_precon_no_invert(1:num_FreeD)= diag_precon(1:num_FreeD) 
diag_precon_no_invert(0) = ZR

      
loads=ZR 
loads(1:num_FreeD) = F(1:num_FreeD)



if (Key_EBE_Precondition==0) then
      diag_precon=ONE
      diag_precon(0) =ZR 
elseif(Key_EBE_Precondition==1) then
      print *, "    PCG-EBE: invert preconditioner and get starting loads..."
      diag_precon(1:)=ONE/diag_precon(1:) 
      diag_precon(0) =ZR 
endif
    
pcg_d=diag_precon*loads 
p=pcg_d
x=ZR
cg_iters=0
      
      
if(Key_Print_EBEPCG_Solution_Time==1)then
  call Tool_Get_Current_Time(current_data,date_time,c_S_Time)
endif
  
  
print *, '    PCG-EBE: PCG equation solution...'
do i_PCG =1,max_num_PCG
  cg_iters=cg_iters+1 
  u=ZR

    select case(Key_EBE_Sym_Storage_K)
  
    case(0)
    
        u_thread(0:num_FreeD,1:max_threads)= ZR  
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,cEle_Loc) 
        c_thread = omp_get_thread_num()+1
        !$OMP DO SCHEDULE(static) 
        do i_E=1,num_FEM_Elem
              c_Elem = FEM_Elem_List(i_E)
              num_Loc_ESM = Size_Local_3D(c_Elem)
              local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
              cEle_Loc = Elem_Location(c_Elem,2) 
              u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+&
                               MATMUL(storK_FEM(1:num_Loc_ESM,1:num_Loc_ESM,cEle_Loc),p(local(1:num_Loc_ESM)))
        enddo
        !$omp end do
        !$omp end parallel   
    
    case(1)
    
        u_thread(0:num_FreeD,1:max_threads)= ZR  
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,localK,cEle_Loc,c_i,c_j) 
        c_thread = omp_get_thread_num()+1
        !$OMP DO SCHEDULE(static)  
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
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_E,c_Elem,num_Loc_ESM,local,cEle_Loc) 
    c_thread = omp_get_thread_num()+1
    !$OMP DO SCHEDULE(static) 
    do i_E=1,num_XFEM_Elem
          c_Elem = XFEM_Elem_List(i_E)
          num_Loc_ESM = Size_Local_3D(c_Elem)
          local(1:num_Loc_ESM)=All_Local_3D(1:num_Loc_ESM,c_Elem) 
          cEle_Loc = Elem_Location(c_Elem,1) 
          u_thread(local(1:num_Loc_ESM),c_thread)=u_thread(local(1:num_Loc_ESM),c_thread)+&
                           MATMUL(storK_XFEM(cEle_Loc)%row(1:num_Loc_ESM,1:num_Loc_ESM),&
                                  p(local(1:num_Loc_ESM)))
                            
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
      exit
  endif          
enddo

101 FORMAT(14X,'Number of iterations to convergence is ',I5)  
102 FORMAT(14X,'Convergence factor is ',E12.5,' / ',E12.5)  
103 FORMAT(14X,'CG_Iters:',I5, ' -> CG_Factor:', E12.5,' / ',E12.5)  
write(*,101)  cg_iters      
write(*,102)  Tol,c_cg_tol   
      
if(Key_Print_EBEPCG_Solution_Time==1)then
    call Tool_Get_Current_Time(current_data,date_time,F_time)
    WRITE(*,5003) F_time-c_S_Time,(dble(F_time)-dble(c_S_Time))/Con_60
endif

loads=xnew
disp = ZR
print *, "    PCG-EBE: retrive nodal displacement..."           
disp(freeDOF(1:num_FreeD)) =loads(1:num_FreeD)      
DEALLOCATE(p,loads,x,xnew,u,diag_precon,pcg_d)  
print *, "    >>>>  End of element by element PCG solver  <<<<"


if (allocated(storK_XFEM_Old)) deallocate(storK_XFEM_Old)
allocate(storK_XFEM_Old(num_XFEM_Elem))
!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_E,c_Elem,cEle_Loc) 
do i_E = 1,num_XFEM_Elem
    c_Elem = XFEM_Elem_List(i_E)
    cEle_Loc = Elem_Location(c_Elem,1) 
    allocate(storK_XFEM_Old(cEle_Loc)%row(Size_Local_3D(c_Elem),Size_Local_3D(c_Elem)))
    storK_XFEM_Old(cEle_Loc)%row(:,:) = storK_XFEM(cEle_Loc)%row(:,:)
enddo 
!$OMP END PARALLEL do
RETURN
END SUBROUTINE Ele_by_Ele_XFEM_PCG_3D
