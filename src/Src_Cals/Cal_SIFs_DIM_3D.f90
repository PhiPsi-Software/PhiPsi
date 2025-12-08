 
subroutine Cal_SIFs_DIM_3D(iter,c_DISP,Tab_Num)


use Global_Float_Type
use Global_Crack_Common
use Global_Crack_3D
use Global_DISP
use Global_Elem_Area_Vol
use Global_Model
use Global_HF
use Global_Common
use Global_Material
use Global_Filename
use Global_Ragged_Array_Real_Classs
use Global_POST
use Global_Cal_Ele_Num_by_Coors_3D  
use Global_Crack   
use Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh
   
implicit none

integer,intent(in)::iter,Tab_Num
real(kind=FT),intent(in)::c_DISP(Total_FD)
integer i_C
integer mat_num
real(kind=FT) offset_delta
real(kind=FT) Line_AB(2,3),new_Line_AB(2,3),new_Point(3)
integer up_Elem_num_1,dn_Elem_num_1
integer up_Elem_num_2,dn_Elem_num_2
real(kind=FT) up_Kesi_1,up_Yita_1,up_Zeta_1,dn_Kesi_1,dn_Yita_1,dn_Zeta_1
real(kind=FT) up_Kesi_2,up_Yita_2,up_Zeta_2,dn_Kesi_2,dn_Yita_2,dn_Zeta_2
real(kind=FT) up_Disp_1(3),dn_Disp_1(3)
real(kind=FT) up_Disp_2(3),dn_Disp_2(3)
real(kind=FT) Relative_Disp_1(3),Relative_Disp_2(3)
real(kind=FT) c_E,c_v,c_G,k
real(kind=FT) T_Matrx(3,3)
real(kind=FT) Relative_Disp_Cr_Sys_1(3),Relative_Disp_Cr_Sys_2(3)
real(kind=FT) r_1,r_2
real(kind=FT) cc_1,cc_2,A
character(5) temp
character(200) c_File_name_test
integer num_Cr_Edges,Edge_P1,Edge_P2,i_V,Meshed_Cr_Ele_Num
real(kind=FT) Edge_P1_Point(3),Edge_P2_Point(3),Egde_Last_P1_Point(3)
real(kind=FT) Vector_x(3),Vector_y(3),Vector_z(3)
integer Edge_P
real(kind=FT) Edge_P_Point(3)      
real(kind=FT) Offseted_Point1(3),Offseted_Point2(3)
real(kind=FT) Crack_Center(3)
real(kind=FT) Offseted_Point1_Up(3),Offseted_Point1_Dn(3)
real(kind=FT) Offseted_Point2_Up(3),Offseted_Point2_Dn(3)
integer c_OUT_Elem,num_Ver_in_Model
integer n_Sigma
integer Num_CrMesh_Outlines,i_Out_Node,c_Mesh_Node
integer in_Ele_Num
real(kind=FT) c_x,c_y,c_z
real(kind=FT) Tool_Function_2Point_Dis_3D
integer c_OUT_Elem_Old(5000)
real(kind=FT) Final_Point(5000,3) 
integer Final_Point_Prop(5000) 
real(kind=FT) c_x_old,c_y_old,c_z_old
real(kind=FT),ALLOCATABLE::Sai_Cr(:),Ld_Cr(:),Ks_Cr(:),D_Cr(:) 
integer num_of_check,i_Check_Theta 
real(kind=FT) c_KI_3D,c_KII_3D,c_KIII_3D
real(kind=FT),ALLOCATABLE::Check_Theta(:)
real(kind=FT),ALLOCATABLE::Spe_Pri_Stress(:)      
real(kind=FT) tem_part1,tem_part2,Stress_Theta,Tao_Theta
real(kind=FT) tem_root_2pir
real(kind=FT) c_Theta 
integer i_Theta_min
real(kind=FT) r_for_c_Theta_finding
real(kind=FT) cc_Part1,cc_Part2,cc_Part3
real(kind=FT) c_KIc 
real(kind=FT),ALLOCATABLE::KI_eq(:),Theta_All(:)
integer c_Mesh_Node_Next
real(kind=FT) Schollm_Max_Theta_in_pi
integer Ele_Num_Cache
real(kind=FT) Signed_Dis,Signed_Dis_v2,Check_Ball_R
logical c_Yes_Node_PER_in_FS,Yes_Found_Min_Signed_Dis
real(kind=FT) c_PER_Node_to_FS(3),c_n_Vector(3)
integer Num_of_SIFs_R_Points   
real(kind=FT) Relative_Disp_Cr_Sys(30,3)    
real(kind=FT) Offseted_Points(30,3)  
real(kind=FT) Relative_Disp(30,3) 
integer i_R_Point
real(kind=FT) offset_r,offset_r_value_x(30)
real(kind=FT) KI_before_fitting(30),KII_before_fitting(30),KIII_before_fitting(30)
real(kind=FT) c_fit_line_k(3),c_fit_line_b(3)
real(kind=FT) Offseted_Points_Up(30,3),Offseted_Points_Dn(30,3)    
integer up_Elem_num(30),dn_Elem_num(30)
real(kind=FT) up_Kesi(30),up_Yita(30),up_Zeta(30),dn_Kesi(30),dn_Yita(30),dn_Zeta(30)
real(kind=FT) up_Disp(30,3),dn_Disp(30,3)
real(kind=FT) r_k


      
offset_delta = SIFs_DIM_3D_Offset_Delta_Factor*Ave_Elem_L_Enrich  
r_1          = SIFs_DIM_3D_r_1_Factor*Ave_Elem_L_Enrich
r_2          = SIFs_DIM_3D_r_2_Factor*Ave_Elem_L_Enrich 
r_k          = SIFs_DIM_3D_r_k_Factor*Ave_Elem_L_Enrich 
Num_of_SIFs_R_Points = 10
Relative_Disp_Cr_Sys(1:30,1:3) = ZR

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,num_Cr_Edges,Crack_Center,  &
!$OMP      num_Ver_in_Model,i_V,Vector_x,Vector_y,Vector_z,T_Matrx, &
!$OMP      Edge_P,Edge_P_Point,c_OUT_Elem,Line_AB,new_Line_AB,&
!$OMP      new_Point,&
!$OMP      cc_1,cc_2,A,Ele_Num_Cache,&
!$OMP      c_E,c_v,c_G,k,Signed_Dis,Signed_Dis_v2,Check_Ball_R,&
!$OMP      c_Yes_Node_PER_in_FS,Yes_Found_Min_Signed_Dis,c_n_Vector,c_PER_Node_to_FS,&
!$OMP      KI_before_fitting,KII_before_fitting,KIII_before_fitting,&
!$OMP      c_fit_line_k,c_fit_line_b,&
!$OMP      Relative_Disp_Cr_Sys,Offseted_Points,Relative_Disp,i_R_Point,offset_r_value_x,offset_r,&
!$OMP      Offseted_Points_Up,Offseted_Points_Dn,up_Elem_num,dn_Elem_num,&
!$OMP      up_Kesi,up_Yita,up_Zeta,dn_Kesi,dn_Yita,dn_Zeta,&
!$OMP      up_Disp,dn_Disp)
      

do i_C =1,num_Crack
    Ele_Num_Cache = 1
    num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
    
    if(allocated(KI_3D(i_C)%row))    deallocate(KI_3D(i_C)%row)
    if(allocated(KII_3D(i_C)%row))   deallocate(KII_3D(i_C)%row)
    if(allocated(KIII_3D(i_C)%row))  deallocate(KIII_3D(i_C)%row)
    if(allocated(KI_eq_3D(i_C)%row)) deallocate(KI_eq_3D(i_C)%row)
    allocate(KI_3D(i_C)%row(num_Cr_Edges));    KI_3D(i_C)%row(:) = ZR
    allocate(KII_3D(i_C)%row(num_Cr_Edges));   KII_3D(i_C)%row(:) = ZR
    allocate(KIII_3D(i_C)%row(num_Cr_Edges));  KIII_3D(i_C)%row(:)=ZR 
    allocate(KI_eq_3D(i_C)%row(num_Cr_Edges)); KI_eq_3D(i_C)%row(:)=ZR
    
    Crack_Center = Crack3D_Centroid(i_C,1:3)
    num_Ver_in_Model = 0
    
    do i_V =1,num_Cr_Edges
      Vector_x = Crack3D_Meshed_Vertex_x_Vector(i_C)%row(i_V,1:3)
      Vector_y = Crack3D_Meshed_Vertex_y_Vector(i_C)%row(i_V,1:3)
      Vector_z = Crack3D_Meshed_Vertex_z_Vector(i_C)%row(i_V,1:3)

      Edge_P = Crack3D_Meshed_Outline(i_C)%row(i_V,1)
      Edge_P_Point = Crack3D_Meshed_Node(i_C)%row(Edge_P,1:3)
      call Cal_Ele_Num_by_Coors_3D(Edge_P_Point(1),Edge_P_Point(2),Edge_P_Point(3),Ele_Num_Cache,c_OUT_Elem)
      
      if(c_OUT_Elem ==0) then
          cycle
      else
          num_Ver_in_Model = num_Ver_in_Model + 1
      endif
      
      T_Matrx(1:3,1:3) = Crack3D_Meshed_Vertex_T_Matrx(i_C)%row(i_V,1:3,1:3)
      
      
      Line_AB(1,1:3) = Edge_P_Point- TWO*Vector_x
      Line_AB(2,1:3) = Edge_P_Point     
      offset_r_value_x = ZR
      do i_R_Point = 1,Num_of_SIFs_R_Points
          offset_r = -r_1 - (r_2-r_1)/(Num_of_SIFs_R_Points-1)*(i_R_Point-1)
          offset_r_value_x(i_R_Point) = abs(offset_r)
          call Tool_Shorten_or_Extend_Line_3D(Line_AB,offset_r,'B',new_Line_AB,new_Point)
          Offseted_Points(i_R_Point,1:3) = new_Point
      enddo
      
      if(Key_Crack_Aperture_Method==1) then
          Check_Ball_R  = 3.0D0*Node_Max_L(Elem_Node(c_OUT_Elem,1))  
          do i_R_Point = 1,Num_of_SIFs_R_Points
              call D3_Get_Signed_Dis_to_Crack_Mesh(Offseted_Points(i_R_Point,1:3),i_C,Check_Ball_R,&
                             Signed_Dis,Signed_Dis_v2,c_Yes_Node_PER_in_FS, &
                             c_PER_Node_to_FS,Yes_Found_Min_Signed_Dis,c_n_Vector)
              if(c_Yes_Node_PER_in_FS) then
                  Offseted_Points(i_R_Point,1:3) =c_PER_Node_to_FS
              else
                  print *,"    ERROR-2023081301 :: point on the crack surface not valid, in Cal_SIFs_DIM_3D.f90!"
                  print *,'                        Crack_Number:',i_C  
                  print *,'                        Offseted_Point:',Offseted_Points(i_R_Point,1:3)  
                  call Warning_Message('S',Keywords_Blank)  
              endif
          enddo
      endif
      
      if(Key_Crack_Aperture_Method==2) then    
          do i_R_Point = 1,Num_of_SIFs_R_Points
              call Tool_Point_Given_Point_Normal_and_Distance_3D(Offseted_Points(i_R_Point,1:3),&
                 Vector_y,offset_delta,Offseted_Points_Up(i_R_Point,1:3))
              call Tool_Point_Given_Point_Normal_and_Distance_3D(Offseted_Points(i_R_Point,1:3),&
                 -Vector_y,offset_delta,Offseted_Points_Dn(i_R_Point,1:3))   

              call Cal_Ele_Num_by_Coors_3D(Offseted_Points_Up(i_R_Point,1),Offseted_Points_Up(i_R_Point,2),&
                                           Offseted_Points_Up(i_R_Point,3),&
                                           Ele_Num_Cache,up_Elem_num(i_R_Point))
              call Cal_Ele_Num_by_Coors_3D(Offseted_Points_Dn(i_R_Point,1),Offseted_Points_Dn(i_R_Point,2),&
                                           Offseted_Points_Dn(i_R_Point,3),&
                                           Ele_Num_Cache,dn_Elem_num(i_R_Point))
              if(up_Elem_num(i_R_Point)==0)then
                  if (Key_Stop_Outside_Crack ==1) then
                      print *,'    WARN-2022110704 :: up_Elem_num(i_R_Point)=0 in Cal_SIFs_DIM_3D.f90!'
                      print *,'                       The coordinates of point:',Offseted_Points_Up(i_R_Point,1:3)
                      print *,'                       Crack ',i_C, ' will stop propagation!'
                      Crack_Type_Status_3D(i_C,3) = 0
                      goto 1037
                  else
                      print *,'    ERROR-2022110704 :: up_Elem_num(i_R_Point)=0 in Cal_SIFs_DIM_3D.f90!'
                      print *,'                        The coordinates of point:',Offseted_Points_Up(i_R_Point,1:3)       
                      call Warning_Message('S',Keywords_Blank)  
                  endif
              endif
              if(dn_Elem_num(i_R_Point)==0)then
                  if (Key_Stop_Outside_Crack ==1) then
                      print *,'    WARN-2022110701 :: dn_Elem_num(i_R_Point)=0 in Cal_SIFs_DIM_3D.f90!'
                      print *,'                       The coordinates of point:',Offseted_Points_Dn(i_R_Point,1:3) 
                      print *,'                       Crack ',i_C, ' will stop propagation!'
                      Crack_Type_Status_3D(i_C,3) = 0
                      goto 1037
                  else
                      print *,'    ERROR-2022110701 :: dn_Elem_num(i_R_Point)=0 in Cal_SIFs_DIM_3D.f90!'
                      print *,'             The coordinates of point:',Offseted_Points_Dn(i_R_Point,1:3)        
                      call Warning_Message('S',Keywords_Blank)  
                  endif
              endif
          enddo
      endif
      if(Key_Crack_Aperture_Method==1) then
          do i_R_Point = 1,Num_of_SIFs_R_Points
              call Cal_Crack_Point_Aperture_3D(c_DISP,i_C,Offseted_Points(i_R_Point,1:3),Relative_Disp(i_R_Point,1:3),0,0,0,0)
          enddo
      elseif(Key_Crack_Aperture_Method==2) then
          do i_R_Point = 1,Num_of_SIFs_R_Points
              call Cal_KesiYita_by_Coor_3D(Offseted_Points_Up(i_R_Point,1:3),up_Elem_num(i_R_Point),&
                                           up_Kesi(i_R_Point),up_Yita(i_R_Point),up_Zeta(i_R_Point))
              call Cal_KesiYita_by_Coor_3D(Offseted_Points_Dn(i_R_Point,1:3),dn_Elem_num(i_R_Point),&
                                           dn_Kesi(i_R_Point),dn_Yita(i_R_Point),dn_Zeta(i_R_Point))
              call Cal_Any_Point_Disp_KesiYita_3D(up_Elem_num(i_R_Point),&
                   up_Kesi(i_R_Point),up_Yita(i_R_Point),up_Zeta(i_R_Point),c_DISP,up_Disp(i_R_Point,1:3))
              call Cal_Any_Point_Disp_KesiYita_3D(dn_Elem_num(i_R_Point),&
                   dn_Kesi(i_R_Point),dn_Yita(i_R_Point),dn_Zeta(i_R_Point),c_DISP,dn_Disp(i_R_Point,1:3))
              Relative_Disp(i_R_Point,1:3) = up_Disp(i_R_Point,1:3) - dn_Disp(i_R_Point,1:3)
          enddo

      endif 
      
      do i_R_Point = 1,Num_of_SIFs_R_Points
          Relative_Disp_Cr_Sys(i_R_Point,1:3)=MATMUL(T_Matrx,Relative_Disp(i_R_Point,1:3))  
      enddo
      
      if(Key_XA/=2) then
        c_E  = Material_Para(Elem_Mat(c_OUT_Elem),1)
        c_v  = Material_Para(Elem_Mat(c_OUT_Elem),2)
      else
        c_E  = Elem_E_XA(c_OUT_Elem)
        c_v  = Elem_Mu_XA(c_OUT_Elem)
      endif
      c_G  = c_E/TWO/(ONE+c_v)
      if(Key_SIFs_DIM_Method == 1) then
          k    = (THR - c_v)/(1+c_v)
      elseif (Key_SIFs_DIM_Method ==2)then
          k    = THR - FOU*c_v
      endif      

      do i_R_Point = 1,Num_of_SIFs_R_Points
          cc_1 = Relative_Disp_Cr_Sys(i_R_Point,2)/sqrt(r_1)
          KI_before_fitting(i_R_Point) = sqrt(TWO*pi)*c_G/(1+k)*cc_1 
          cc_1 = Relative_Disp_Cr_Sys(i_R_Point,1)/sqrt(r_1)  
          KII_before_fitting(i_R_Point) = sqrt(TWO*pi)*c_G/(1+k)*cc_1 
          cc_1 = Relative_Disp_Cr_Sys(i_R_Point,3)/sqrt(r_1)
          KIII_before_fitting(i_R_Point) = sqrt(TWO*pi)*c_G/(1+k)*cc_1
      enddo
      
      call Tool_Fit_Least_Square_Straight_Line(offset_r_value_x(1:Num_of_SIFs_R_Points),&
         KI_before_fitting(1:Num_of_SIFs_R_Points), Num_of_SIFs_R_Points,c_fit_line_k(1),c_fit_line_b(1))
      call Tool_Fit_Least_Square_Straight_Line(offset_r_value_x(1:Num_of_SIFs_R_Points),&
         KII_before_fitting(1:Num_of_SIFs_R_Points), Num_of_SIFs_R_Points,c_fit_line_k(2),c_fit_line_b(2))   
      call Tool_Fit_Least_Square_Straight_Line(offset_r_value_x(1:Num_of_SIFs_R_Points),&
         KIII_before_fitting(1:Num_of_SIFs_R_Points), Num_of_SIFs_R_Points,c_fit_line_k(3),c_fit_line_b(3))   
      KI_3D(i_C)%row(i_V)   = c_fit_line_k(1)*r_k + c_fit_line_b(1)
      KII_3D(i_C)%row(i_V)  = c_fit_line_k(2)*r_k + c_fit_line_b(2)
      KIII_3D(i_C)%row(i_V) = c_fit_line_k(3)*r_k + c_fit_line_b(3)
    end do
    1037 continue
enddo
!$OMP END PARALLEL DO

      
if(Key_Save_Nothing  == 0) then
    if(Tab_Num==5)then
          print *,'    Saving ckr1,ckr2,ckr3 files...'
    elseif(Tab_Num==8)then
          print *,'       Saving ckr1,ckr2,ckr3 files...'
    endif

    write(temp,'(I5)') iter
    c_File_name_test   = trim(Full_Pathname)//'.ckr1_'//ADJUSTL(temp)  
    open(101,file=c_File_name_test,status='unknown')  
    c_File_name_test   = trim(Full_Pathname)//'.ckr2_'//ADJUSTL(temp)  
    open(102,file=c_File_name_test,status='unknown')  
    c_File_name_test   = trim(Full_Pathname)//'.ckr3_'//ADJUSTL(temp)  
    open(103,file=c_File_name_test,status='unknown')  
    do i_C = 1,num_Crack   
          num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
          write(101, '(50000E20.12)') (KI_3D(i_C)%row(i_V),i_v=1,num_Cr_Edges)
          write(102, '(50000E20.12)') (KII_3D(i_C)%row(i_V),i_v=1,num_Cr_Edges)
          write(103, '(50000E20.12)') (KIII_3D(i_C)%row(i_V),i_v=1,num_Cr_Edges)     
    enddo
    close(101)
    close(102)
    close(103)
endif

if (Key_Denoise_Vertex_Value>=1)then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,num_Cr_Edges,n_Sigma)  
  do i_C = 1,num_Crack   
      num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
      if (i_C==1) then
          if(Tab_Num==5)then
              print *,'    Denoising vertex values...'   
          elseif(Tab_Num==8)then
              print *,'       Denoising vertex values...'   
          endif                  
      endif
      n_Sigma =2 
      call Tool_Denoise_Data(KI_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,Key_Denoise_Vertex_Value,n_Sigma,.True.) 
      call Tool_Denoise_Data(KII_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,Key_Denoise_Vertex_Value,n_Sigma,.True.)       
      call Tool_Denoise_Data(KIII_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,Key_Denoise_Vertex_Value,n_Sigma,.True.)      
  enddo
!$OMP END PARALLEL DO
endif
      
if (Key_Smooth_Vertex_Value>=1)then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,num_Cr_Edges,n_Sigma)        
      do i_C = 1,num_Crack     
          num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
          if (i_C==1) then
              if(Tab_Num==5)then
                  print *,'    Smoothing vertex values...' 
              elseif(Tab_Num==8)then
                  print *,'       Smoothing vertex values...'   
              endif                
          endif
          call Tool_Smooth_Data(KI_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,Key_Smooth_Vertex_Value,Smooth_Vertex_n,.True.)
          call Tool_Smooth_Data(KII_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,Key_Smooth_Vertex_Value,Smooth_Vertex_n,.True.)   
          call Tool_Smooth_Data(KIII_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,Key_Smooth_Vertex_Value,Smooth_Vertex_n,.True.)      
      enddo
!$OMP END PARALLEL DO 
             
      if (Key_Smooth_Vertex_Value2>=1)then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,num_Cr_Edges,n_Sigma)           
          do i_C = 1,num_Crack     
              num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
              if (i_C==1) then
                  if(Tab_Num==5)then
                      print *,'    Smoothing vertex values...'
                  elseif(Tab_Num==8)then
                      print *,'       Smoothing vertex values...'   
                  endif 
              endif
              call Tool_Smooth_Data(KI_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,Key_Smooth_Vertex_Value2,Smooth_Vertex_n2,.True.) 
              call Tool_Smooth_Data(KII_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges,Key_Smooth_Vertex_Value2,Smooth_Vertex_n2,.True.) 
              call Tool_Smooth_Data(KIII_3D(i_C)%row(1:num_Cr_Edges),num_Cr_Edges, Key_Smooth_Vertex_Value2,Smooth_Vertex_n2,.True.) 
          enddo
!$OMP END PARALLEL DO               
      endif           
endif  

 
      
write(temp,'(I5)') iter
if(Key_Save_Nothing /= 1) then 
    if(Tab_Num==5)then
          print *,'    Saving cvk1,cvk2,cvk3 files...'
    elseif(Tab_Num==8)then
          print *,'       Saving cvk1,cvk2,cvk3 files...'
    endif 
    c_File_name_test   = trim(Full_Pathname)//'.cvk1_'//ADJUSTL(temp)  
    open(101,file=c_File_name_test,status='unknown')  
    c_File_name_test   = trim(Full_Pathname)//'.cvk2_'//ADJUSTL(temp)  
    open(102,file=c_File_name_test,status='unknown')  
    c_File_name_test   = trim(Full_Pathname)//'.cvk3_'//ADJUSTL(temp)  
    open(103,file=c_File_name_test,status='unknown')  
    do i_C = 1,num_Crack   
          num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
          write(101, '(50000E20.12)') (KI_3D(i_C)%row(i_V),i_v=1,num_Cr_Edges)
          write(102, '(50000E20.12)') (KII_3D(i_C)%row(i_V),i_v=1,num_Cr_Edges)
          write(103, '(50000E20.12)') (KIII_3D(i_C)%row(i_V),i_v=1,num_Cr_Edges)     
    enddo
    close(101)
    close(102)
    close(103)
endif      
      
      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,Num_CrMesh_Outlines,    &
!$OMP      c_OUT_Elem_Old,Final_Point,Final_Point_Prop,&
!$OMP      i_Out_Node,c_Mesh_Node,c_x_old,c_y_old,c_z_old,&
!$OMP      num_of_check,Check_Theta,Spe_Pri_Stress,Sai_Cr,Ld_Cr,&
!$OMP      Ks_Cr,D_Cr,KI_eq,Theta_All,r_for_c_Theta_finding,&
!$OMP      tem_root_2pir,c_Mesh_Node_Next,in_Ele_Num,mat_num,c_KIc,&
!$OMP      c_x,c_y,c_z,c_OUT_Elem,c_KI_3D,c_KII_3D,c_KIII_3D,&
!$OMP      i_Check_Theta,Schollm_Max_Theta_in_pi,tem_part1,tem_part2,&
!$OMP      cc_Part1,cc_Part2,cc_Part3,Stress_Theta,Tao_Theta,&
!$OMP      i_Theta_min,c_Theta,Ele_Num_Cache)    
do i_C =1,num_Crack  
    Ele_Num_Cache = 1
    Num_CrMesh_Outlines = Crack3D_Meshed_Outline_num(i_C)
    c_OUT_Elem_Old(1:5000)    = 1
    Final_Point(1:5000,1:3)   = ZR 
    Final_Point_Prop(1:5000)  = 0
    do i_Out_Node = 1,Num_CrMesh_Outlines
      c_Mesh_Node = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
      c_x_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1) 
      c_y_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,2) 
      c_z_old  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,3)
      Final_Point(i_Out_Node,1)=c_x_old
      Final_Point(i_Out_Node,2)=c_y_old
      Final_Point(i_Out_Node,3)=c_z_old
    enddo
    num_of_check = 720
    if(allocated(Check_Theta)) deallocate(Check_Theta)
    if(allocated(Spe_Pri_Stress)) deallocate(Spe_Pri_Stress)
    if(allocated(Sai_Cr)) deallocate(Sai_Cr)
    if(allocated(Ld_Cr)) deallocate(Ld_Cr)
    if(allocated(Ks_Cr)) deallocate(Ks_Cr)
    if(allocated(D_Cr)) deallocate(D_Cr)
    if(allocated(KI_eq)) deallocate(KI_eq)
    if(allocated(Theta_All)) deallocate(Theta_All)
    allocate(Check_Theta(num_of_check+1))
    allocate(Spe_Pri_Stress(num_of_check+1))
    allocate(Sai_Cr(Num_CrMesh_Outlines))
    allocate(Ld_Cr(Num_CrMesh_Outlines))
    allocate(Ks_Cr(Num_CrMesh_Outlines))  
    allocate(D_Cr(Num_CrMesh_Outlines))  
    allocate(KI_eq(Num_CrMesh_Outlines))
    allocate(Theta_All(Num_CrMesh_Outlines))
    
    r_for_c_Theta_finding = 1.0*Ave_Elem_L
    tem_root_2pir = sqrt(TWO*pi*r_for_c_Theta_finding)
    
    
    do i_Out_Node = 1,Num_CrMesh_Outlines
      c_Mesh_Node = Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
      c_Mesh_Node_Next =Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,2)
      in_Ele_Num  = Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(c_Mesh_Node)
      if(in_Ele_Num ==0 ) then
          cycle
      endif
  
      if(Key_XA/=2) then
          mat_num = Elem_Mat(in_Ele_Num)
          c_KIc  = Material_Para(mat_num,6)
      else
          c_KIc  = Elem_KIc_XA(in_Ele_Num)
      endif
          
      c_x  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1) 
      c_y  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,2) 
      c_z  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,3) 
      
      call Cal_Ele_Num_by_Coors_3D(c_x,c_y,c_z,Ele_Num_Cache,c_OUT_Elem)
      if(c_OUT_Elem ==0) then
          cycle
      endif
  
      Ld_Cr(i_Out_Node) = Tool_Function_2Point_Dis_3D(Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1:3),&
                                                      Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node_Next,1:3))

      Ks_Cr(i_Out_Node) = ONE/Ld_Cr(i_Out_Node)
      D_Cr(i_Out_Node) = Ld_Cr(i_Out_Node)**3
      
      c_KI_3D   = KI_3D(i_C)%row(i_Out_Node)
      c_KII_3D  = KII_3D(i_C)%row(i_Out_Node)
      c_KIII_3D = KIII_3D(i_C)%row(i_Out_Node)
      do i_Check_Theta =1,num_of_check+1
        Schollm_Max_Theta_in_pi = Schollm_Max_Theta*pi/Con_180
        Check_Theta(i_Check_Theta)=-Schollm_Max_Theta_in_pi+TWO*Schollm_Max_Theta_in_pi/num_of_check*(i_Check_Theta-1)   
 
        tem_part1 = THR*cos(Check_Theta(i_Check_Theta)/TWO) +cos(THR*Check_Theta(i_Check_Theta)/TWO)
        tem_part2 = THR*sin(Check_Theta(i_Check_Theta)/TWO) +THR*sin(THR*Check_Theta(i_Check_Theta)/TWO)     
        Stress_Theta = c_KI_3D/FOU/tem_root_2pir*tem_part1 -c_KII_3D/FOU/tem_root_2pir*tem_part2
        Tao_Theta=c_KIII_3D*cos(Check_Theta(i_Check_Theta)/TWO)/tem_root_2pir
        Spe_Pri_Stress(i_Check_Theta) = Stress_Theta/TWO + ZP5*sqrt(Stress_Theta**2 + FOU*Tao_Theta**2)
      enddo 
      i_Theta_min = maxloc(Spe_Pri_Stress(1:num_of_check+1),1)
      c_Theta = Check_Theta(i_Theta_min)
      Theta_All(i_Out_Node) = c_Theta
      
      cc_Part1 = c_KI_3D*(cos(c_Theta/TWO))**2
      cc_Part2 = THR/TWO*c_KII_3D*sin(c_Theta)
      cc_Part3 =sqrt((cc_Part1 - cc_Part2)**2 + FOU*c_KIII_3D**2)
      KI_eq(i_Out_Node)  =ZP5*cos(c_Theta/TWO)*(cc_Part1-cc_Part2+cc_Part3)
      KI_eq_3D(i_C)%row(i_Out_Node) = KI_eq(i_Out_Node)
    end do
enddo      
!$OMP END PARALLEL DO   
      
if(Key_Save_Nothing  == 0) then
    if(Tab_Num==5)then
        print *,'    Saving cvke file for cracks...'
    elseif(Tab_Num==8)then
        print *,'       Saving cvke file for cracks...'
    endif

    write(temp,'(I5)') iter
    c_File_name_test = trim(Full_Pathname)//'.cvke_'//ADJUSTL(temp)  
    open(104,file=c_File_name_test,status='unknown')  
    do i_C = 1,num_Crack   
          num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
          write(104, '(50000E20.12)') (KI_eq_3D(i_C)%row(i_v) ,i_v=1,num_Cr_Edges)   
    enddo
    close(104)
endif    
      
RETURN
END SUBROUTINE Cal_SIFs_DIM_3D