 
subroutine Cal_Crack_Permeability_3D(isub)

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

implicit none
integer,intent(in)::isub
integer i_C,num_All_Flu_Nodes,i_FluNode,i_FluidEle
real(kind=FT) c_Aperture,c_n_vector(3)
real(kind=FT) c_k
real(kind=FT) c_T_Matrx(3,3),c_T_Matrx_T(3,3)
real(kind=FT) c_ThetaX,c_ThetaY,c_ThetaZ
integer c_node_1,c_node_2,c_node_3
integer c_fluid_node
real(kind=FT) k_Matrix_Crack(3,3),k_Matrix_Global(3,3)
real(kind=FT) k_xx,k_yy,k_zz,k_xy,k_yz,k_xz
integer c_Max_N_Node_3D
character(200) c_File_name_1,c_File_name_2
character(5) temp
integer i,j
real(kind=FT),ALLOCATABLE:: tem_Cracks_CalP_k_3D(:,:,:)
integer max_fluid_node
integer c_Ele
real(kind=FT) c_flu_k_6(6)
integer i_Ele
integer i_Thread,c_Thread,max_threads
real(kind=FT),ALLOCATABLE::Ele_Permeability_3D_Thread(:,:,:)
real(kind=FT) c_Flu_Area,c_Ele_Area,c_Flu_w,c_Ele_L,c_Ele_V
real(kind=FT),ALLOCATABLE::Ele_VolumeRatio_3D_Thread(:,:)
real(kind=FT) c_VolumeRatio
if(Key_Get_Permeability==0) return   


max_fluid_node = maxval(Cracks_CalP_Num_3D(1:num_Crack))
allocate(tem_Cracks_CalP_k_3D(num_Crack,max_fluid_node,6))
tem_Cracks_CalP_k_3D = ZR
!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,i_FluidEle,c_T_Matrx,c_T_Matrx_T,  &
!$OMP      i_FluNode,c_fluid_node,c_Aperture,c_k,k_Matrix_Crack,k_Matrix_Global, &
!$OMP      k_xx,k_yy,k_zz,k_xy,k_yz,k_xz)                                        &
!$OMP SCHEDULE(static)    
do i_C = 1,num_Crack
    
    do i_FluidEle =1,Cracks_FluidEle_num_3D(i_C) 
        c_T_Matrx   = Cracks_FluidEle_LCS_T_3D(i_C)%row(i_FluidEle,1:3,1:3)
        c_T_Matrx_T = transpose(c_T_Matrx)
        
        do i_FluNode=1,3
            c_fluid_node = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,i_FluNode)
            c_Aperture = Cracks_CalP_Aper_3D(i_C)%row(c_fluid_node)
            if(c_Aperture<=ZR) then
                c_Aperture=ZR
            endif
            c_k        = c_Aperture**2/12.0D0
            k_Matrix_Crack(1:3,1:3) = ZR
            k_Matrix_Crack(1,1) = c_k
            k_Matrix_Crack(2,2) = c_k
            k_Matrix_Crack(3,3) = ZR
            k_Matrix_Global(1:3,1:3) = ZR
            
            k_Matrix_Global= MATMUL(c_T_Matrx_T,k_Matrix_Crack)
            k_Matrix_Global= MATMUL(k_Matrix_Global,c_T_Matrx)
            
            k_xx = k_Matrix_Global(1,1)
            k_yy = k_Matrix_Global(2,2)
            k_zz = k_Matrix_Global(3,3)
            k_xy = k_Matrix_Global(1,2)
            k_yz = k_Matrix_Global(2,3)
            k_xz = k_Matrix_Global(1,3)
             
             tem_Cracks_CalP_k_3D(i_C,c_fluid_node,1)  = k_xx
             tem_Cracks_CalP_k_3D(i_C,c_fluid_node,2)  = k_yy
             tem_Cracks_CalP_k_3D(i_C,c_fluid_node,3)  = k_zz
             tem_Cracks_CalP_k_3D(i_C,c_fluid_node,4)  = k_xy
             tem_Cracks_CalP_k_3D(i_C,c_fluid_node,5)  = k_yz
             tem_Cracks_CalP_k_3D(i_C,c_fluid_node,6)  = k_xz
        enddo
    enddo
enddo   
!$omp end parallel do  

if(Key_Save_Nothing /= 1)then
    write(temp,'(I5)') isub
    c_File_name_1   =  trim(Full_Pathname)//'.ckxx'//'_'//ADJUSTL(temp)    
    open(101,file=c_File_name_1,status='unknown') 
    do i_C=1,num_Crack
    write(101,'(50000E20.12)') (tem_Cracks_CalP_k_3D(i_C,j,1),j=1,Cracks_CalP_Num_3D(i_C))     
    end do
    close(101) 
    c_File_name_1   =  trim(Full_Pathname)//'.ckyy'//'_'//ADJUSTL(temp)    
    open(102,file=c_File_name_1,status='unknown') 
    do i_C=1,num_Crack
    write(102,'(50000E20.12)') (tem_Cracks_CalP_k_3D(i_C,j,2),j=1,Cracks_CalP_Num_3D(i_C))     
    end do
    close(102) 
    c_File_name_1   =  trim(Full_Pathname)//'.ckzz'//'_'//ADJUSTL(temp)    
    open(103,file=c_File_name_1,status='unknown') 
    do i_C=1,num_Crack
    write(103,'(50000E20.12)') (tem_Cracks_CalP_k_3D(i_C,j,3),j=1,Cracks_CalP_Num_3D(i_C))     
    end do
    close(103) 
    c_File_name_1   =  trim(Full_Pathname)//'.ckxy'//'_'//ADJUSTL(temp)    
    open(104,file=c_File_name_1,status='unknown') 
    do i_C=1,num_Crack
    write(104,'(50000E20.12)') (tem_Cracks_CalP_k_3D(i_C,j,4),j=1,Cracks_CalP_Num_3D(i_C))     
    end do
    close(104) 
    c_File_name_1   =  trim(Full_Pathname)//'.ckyz'//'_'//ADJUSTL(temp)    
    open(105,file=c_File_name_1,status='unknown') 
    do i_C=1,num_Crack
    write(105,'(50000E20.12)') (tem_Cracks_CalP_k_3D(i_C,j,5),j=1,Cracks_CalP_Num_3D(i_C))     
    end do
    close(105) 
    c_File_name_1   =  trim(Full_Pathname)//'.ckxz'//'_'//ADJUSTL(temp)    
    open(106,file=c_File_name_1,status='unknown') 
    do i_C=1,num_Crack
    write(106,'(50000E20.12)') (tem_Cracks_CalP_k_3D(i_C,j,6),j=1,Cracks_CalP_Num_3D(i_C))     
    end do
    close(106)        
endif

if (allocated(Ele_Permeability_3D)) deallocate(Ele_Permeability_3D)
allocate(Ele_Permeability_3D(Num_Elem,6))
Ele_Permeability_3D(1:Num_Elem,1:6) = ZR
max_threads = omp_get_max_threads()
ALLOCATE(Ele_Permeability_3D_Thread(Num_Elem,6,max_threads))  
Ele_Permeability_3D_Thread(1:Num_Elem,1:6,1:max_threads) = ZR
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_C,i_FluidEle,c_Ele,c_flu_k_6,i_FluNode,c_fluid_node,c_Flu_Area,&
!$OMP                                  c_Flu_w,c_Ele_V) 
c_thread = omp_get_thread_num()+1 
!$OMP do  
do i_C = 1,num_Crack
    do i_FluidEle =1,Cracks_FluidEle_num_3D(i_C) 
        c_Ele = Cracks_FluidEle_EleNum_3D(i_C)%row(i_FluidEle)
        c_Flu_Area= Cracks_FluidEle_Area_3D(i_C)%row(i_FluidEle)
        c_Flu_w   = Cracks_FluidEle_Aper_3D(i_C)%row(i_FluidEle)
        
        if (c_Flu_w<=ZR) c_Flu_w = ZR
        
        c_Ele_V   = Elem_Vol(c_Ele)
        c_flu_k_6(1:6) = ZR
        do i_FluNode=1,3
            c_fluid_node   = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,i_FluNode)
            c_flu_k_6(1:6) = c_flu_k_6(1:6) + tem_Cracks_CalP_k_3D(i_C,c_fluid_node,1:6)
        enddo
        c_flu_k_6(1:6) = c_flu_k_6(1:6)/THR
        
        
        c_flu_k_6(1:6) = c_flu_k_6(1:6)*c_Flu_Area*c_Flu_w/c_Ele_V
        
        c_flu_k_6(1:6) = c_flu_k_6(1:6)*1.01325D15
        

        Ele_Permeability_3D_Thread(c_Ele,1:6,c_thread) = Ele_Permeability_3D_Thread(c_Ele,1:6,c_thread) + c_flu_k_6(1:6) 
    enddo
enddo   
!$omp end do
!$omp end parallel    
DO i_Thread = 1,omp_get_max_threads()
    Ele_Permeability_3D  =  Ele_Permeability_3D  + Ele_Permeability_3D_Thread(:,:,i_Thread)
ENDDO
deallocate(Ele_Permeability_3D_Thread) 
  
if(Key_Save_Nothing /= 1)then
    c_File_name_1   =  trim(Full_Pathname)//'.elek'//'_'//ADJUSTL(temp)    
    select case(Key_Data_Format)
    case(1)
        open(107,file=c_File_name_1,status='unknown') 
        do i_Ele =1,Num_Elem
            write(107,'(6E20.12)') Ele_Permeability_3D(i_Ele,1:6)     
        end do
        close(107) 
    case(2)
        open(107,file=c_File_name_1,status='unknown',form='unformatted',access='stream')    
        do i_Ele =1,Num_Elem
            write(107) Ele_Permeability_3D(i_Ele,1:6)   
        end do
        close(107) 
    endselect
endif


if (Key_Cpp_Call_Fortran_Lib==1) then
    if (allocated(Ele_VolumeRatio_3D)) deallocate(Ele_VolumeRatio_3D)
    allocate(Ele_VolumeRatio_3D(Num_Elem))           
    Ele_VolumeRatio_3D(1:Num_Elem) = ZR
    max_threads = omp_get_max_threads()
    ALLOCATE(Ele_VolumeRatio_3D_Thread(Num_Elem,max_threads))  
    Ele_VolumeRatio_3D_Thread(1:Num_Elem,1:max_threads) = ZR
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c_thread,i_C,i_FluidEle,c_Ele,c_VolumeRatio,c_Flu_Area,&
    !$OMP                                  c_Flu_w,c_Ele_V) 
    c_thread = omp_get_thread_num()+1 
    !$OMP do  
    do i_C = 1,num_Crack
        do i_FluidEle =1,Cracks_FluidEle_num_3D(i_C) 
            c_Ele = Cracks_FluidEle_EleNum_3D(i_C)%row(i_FluidEle)
            c_Flu_Area= Cracks_FluidEle_Area_3D(i_C)%row(i_FluidEle)
            c_Flu_w   = Cracks_FluidEle_Aper_3D(i_C)%row(i_FluidEle)
            
            if (c_Flu_w<=ZR) c_Flu_w = ZR
            
            c_Ele_V   = Elem_Vol(c_Ele)
            
            c_VolumeRatio = c_Flu_Area*c_Flu_w/c_Ele_V
            
            Ele_VolumeRatio_3D_Thread(c_Ele,c_thread) = Ele_VolumeRatio_3D_Thread(c_Ele,c_thread) + c_VolumeRatio
        enddo
    enddo   
    !$omp end do
    !$omp end parallel    
    DO i_Thread = 1,omp_get_max_threads()
        Ele_VolumeRatio_3D  =  Ele_VolumeRatio_3D  + Ele_VolumeRatio_3D_Thread(:,i_Thread)
    ENDDO
    deallocate(Ele_VolumeRatio_3D_Thread) 

endif



deallocate(tem_Cracks_CalP_k_3D)



   
return 
end SUBROUTINE Cal_Crack_Permeability_3D              
