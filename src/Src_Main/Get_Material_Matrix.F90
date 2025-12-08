!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
SUBROUTINE Get_Material_Matrix
! Calculate the material matrix D.

!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Dynamic
use Global_Material
use Global_HF
use Global_POST
use Global_Filename

!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
integer i_mat
real(kind=FT) D1,D2,D3,c_E,c_v,K_m,E_1,v_1,KIc_1,K_m_Part1,K_m_Part2,Inject_Q,c_v1,c_v2
integer i_E,c_mat_type
real(kind=FT) basic_E,shape_parameters,c0,c1
real(kind=FT) Weibull_Data(Num_Elem)
real(kind=FT) choosed_E
character(200) c_File_name_1
      
print *,' '
print *,'    Calculating material matrix D...'

ALLOCATE(D(num_of_Material,3,3))
ALLOCATE(D4(num_of_Material,4,4))
ALLOCATE(S(num_of_Material,3,3))
ALLOCATE(St(num_of_Material,2))
ALLOCATE(Sc(num_of_Material,2))
ALLOCATE(T_Alpha(num_of_Material))
ALLOCATE(KIc(num_of_Material,2))
ALLOCATE(E(num_of_Material,3))
ALLOCATE(v(num_of_Material,3))
ALLOCATE(thick(num_of_Material))
ALLOCATE(density(num_of_Material))
ALLOCATE(Lame_lambda(num_of_Material))
ALLOCATE(Lame_mu(num_of_Material))
ALLOCATE(MC_dt(num_of_Material))
ALLOCATE(MC_phi_deg(num_of_Material))
ALLOCATE(MC_phi_rad(num_of_Material))
ALLOCATE(MC_psi_deg(num_of_Material))
ALLOCATE(MC_psi_rad(num_of_Material))
ALLOCATE(MC_c(num_of_Material))
!-------------------------------------
! Recycling of all types of materials
!-------------------------------------
do i_mat =1,num_of_Material
  if (Material_Type(i_mat) .eq. 1 )then
      c_E  = Material_Para(i_mat,1)
      c_v  = Material_Para(i_mat,2)
      Lame_lambda(i_mat)  =c_E*c_v/(ONE+c_v)/(ONE-TWO*c_v)
      Lame_mu(i_mat)      =c_E/TWO/(ONE+c_v)
      density(i_mat)   = Material_Para(i_mat,3)
      St(i_mat,1)      = Material_Para(i_mat,5)
      Sc(i_mat,1)      = Material_Para(i_mat,7)
      T_Alpha(i_mat)   = Material_Para(i_mat,8)
      KIc(i_mat,1)     = Material_Para(i_mat,6)
      if(Key_Type_2D .eq. 1  )then
          D1 = c_E/(ONE-c_v**2)
          D2 = c_E*c_v/(ONE-c_v**2)
          D3 = c_E/TWO/(ONE+c_v)
          D(i_mat,1,:)  = [D1,   D2,   ZR]
          D(i_mat,2,:)  = [D2,   D1,   ZR]
          D(i_mat,3,:)  = [ZR,ZR,   D3]
          thick(i_mat)  = Material_Para(i_mat,4)
      elseif (Key_Type_2D .eq.2 )then
          D1 = c_E*(ONE-c_v)/(ONE+c_v)/(ONE-TWO*c_v)
          D2 = c_E*c_v/(ONE+c_v)/(ONE-TWO*c_v)
          D3 = c_E/TWO/(ONE+c_v)
          D(i_mat,1,:)  = [D1,   D2,   ZR]
          D(i_mat,2,:)  = [D2,   D1,   ZR]
          D(i_mat,3,:)  = [ZR,   ZR,   D3]
          thick(i_mat)  = ONE
          ! D4, 4x4 material matrix, Ref: MATLAB FEM Code - From Elasticity to Plasticity, Eq. 4.5
          D4(i_mat,1,1:4) =  [D1,   D2,   ZR,   D2]
          D4(i_mat,2,1:4) =  [D2,   D1,   ZR,   D2]
          D4(i_mat,3,1:4) =  [ZR,   ZR,   D3,   ZR]
          D4(i_mat,4,1:4) =  [D2,   D2,   ZR,   D1]
      end if
      ! S matrix
      !call Matrix_Inverse(D(i_mat,:,:),S(i_mat,:,:),3)
      call Matrix_Inverse_3x3(D(i_mat,1:3,1:3),S(i_mat,1:3,1:3))
      E(i_mat,1) = c_E
      v(i_mat,1) = c_v
  ! Materials other than elastic
else
      c_E  = Material_Para(i_mat,1)
      c_v  = Material_Para(i_mat,2)
      Lame_lambda(i_mat)  =c_E*c_v/(ONE+c_v)/(ONE-TWO*c_v)
      Lame_mu(i_mat)      =c_E/TWO/(ONE+c_v)
      density(i_mat)   = Material_Para(i_mat,3)
      St(i_mat,1)      = Material_Para(i_mat,5)
      Sc(i_mat,1)      = Material_Para(i_mat,7)
      T_Alpha(i_mat)   = Material_Para(i_mat,8)
      KIc(i_mat,1)     = Material_Para(i_mat,6)
      if(Key_Type_2D .eq. 1  )then
          D1 = c_E/(ONE-c_v**2)
          D2 = c_E*c_v/(ONE-c_v**2)
          D3 = c_E/TWO/(ONE+c_v)
          D(i_mat,1,:)  = [D1,   D2,   ZR]
          D(i_mat,2,:)  = [D2,   D1,   ZR]
          D(i_mat,3,:)  = [ZR,ZR,   D3]
          thick(i_mat)  = Material_Para(i_mat,4)
      elseif (Key_Type_2D .eq.2 )then
          D1 = c_E*(ONE-c_v)/(ONE+c_v)/(ONE-TWO*c_v)
          D2 = c_E*c_v/(ONE+c_v)/(ONE-TWO*c_v)
          D3 = c_E/TWO/(ONE+c_v)
          D(i_mat,1,:)  = [D1,   D2,   ZR]
          D(i_mat,2,:)  = [D2,   D1,   ZR]
          D(i_mat,3,:)  = [ZR,ZR,   D3]
          thick(i_mat)  = ONE
          ! D4, 4x4 material matrix, Ref: MATLAB FEM Code - From Elasticity to Plasticity, Eq. 4.5
          D4(i_mat,1,1:4) =  [D1,   D2,   ZR,   D2]
          D4(i_mat,2,1:4) =  [D2,   D1,   ZR,   D2]
          D4(i_mat,3,1:4) =  [ZR,   ZR,   D3,   ZR]
          D4(i_mat,4,1:4) =  [D2,   D2,   ZR,   D1]              
      end if
      ! S matrix
      !call Matrix_Inverse(D(i_mat,:,:),S(i_mat,:,:),3)
      call Matrix_Inverse_3x3(D(i_mat,:,:),S(i_mat,:,:))
      E(i_mat,1) = c_E
      v(i_mat,1) = c_v
      ! Mohr-Coulomb criterion related parameters, 2021-07-07
      MC_phi_deg(i_mat) =  Material_Para(i_mat,16)
      MC_phi_rad(i_mat) =  Material_Para(i_mat,16)*pi/Con_180
      MC_psi_deg(i_mat) =  Material_Para(i_mat,21)
      MC_psi_rad(i_mat) =  Material_Para(i_mat,21)*pi/Con_180
      MC_c(i_mat)       =  Material_Para(i_mat,15)
      ! Mohr-Coulomb criterion viscoplastic method pseudo time step, Ref: MATLAB FEM Code - From
      ! Elasticity to Plasticity, Eq.6.10
      MC_dt(i_mat)      = FOU*(ONE+c_v)*(ONE-TWO*c_v)/c_E/(ONE-TWO*c_v+(sin(MC_phi_rad(i_mat)))**2)
  end if
end do

!----------------------------------------------------------------------------------------------------
! If it is a hydraulic fracturing analysis, then K_m also needs to be calculated (using the
! parameters of material No. 1 for all calculations).
! and determine the type of hydraulic fracturing, see the formula in Paper_Bao_2014_A Coupled Finite
! Element method for the numerical simulation of hydraulic fracturing with a condensation
! technique_Eq.29.
!----------------------------------------------------------------------------------------------------
if (Key_Analysis_Type == 3) then
  ! Calculate K_m
  print *,'    Calculating K_m......'
  E_1 = Material_Para(1,1)
  v_1 = Material_Para(1,2)
  KIc_1 = Material_Para(1,6)
  Inject_Q = Inject_Q_Val(1)
  K_m_Part1 = FOU*sqrt(TWO/pi)*((KIc_1*(ONE-v_1**2)))/E_1
  K_m_Part2 = (E_1/(12.0D0*Viscosity*Inject_Q*(ONE-v_1**2)))**0.25D0
  K_m = K_m_Part1*K_m_Part2
  Global_K_m = K_m
  
  
  print *,'    Value of K_m:       ',K_m
  ! Determine Type
  if (K_m < ONE) then
      Type_of_HF = 1
      print *,'    Type of hydraulic fracturing: viscosity-dominated'
  elseif (K_m > FOU) then
      Type_of_HF = 2
      print *,'    Type of hydraulic fracturing: toughness-dominated'
  else
      Type_of_HF = 0
      print *,'    Type of hydraulic fracturing: not typical'
  end if
end if

!---------------------------------------------------
! Weibull processing. 2024-06-25. NEWFTU2024062402.
!---------------------------------------------------
#ifndef github
if(any(Key_Weibull_E == 1)) then
    Flag_Weibull_E = .True.
    allocate(Weibull_Elements_E(Num_Elem))
    Weibull_Elements_E(1:Num_Elem) = ZR
    allocate(Weibull_Elements_D_Matrix(Num_Elem,3,3))
    Weibull_Elements_D_Matrix(1:Num_Elem,1:3,1:3) = ZR
endif
if (Flag_Weibull_E) then
    do i_mat =1,num_of_Material
        if (Key_Weibull_E(i_mat) ==1)then
            ! Generate Weibull parameters.
            basic_E = Material_Para(i_mat,1) 
            c_v  = Material_Para(i_mat,2)
            shape_parameters = Weibull_Parameters_E(i_mat,1) 
            c0 = Weibull_Parameters_E(i_mat,2)
            c1 = Weibull_Parameters_E(i_mat,3)
            !Call Weibull Tool.
            call Tool_Generate_Weibull_Distribution_Data(Num_Elem,basic_E,shape_parameters,c0,c1,Weibull_Data(1:Num_Elem))
            
            do i_E=1,Num_Elem
                c_mat_type = Elem_Mat(i_E)  
                if(c_mat_type == i_mat)then
                    choosed_E = Weibull_Data(i_E)
                    Weibull_Elements_E(i_E) = choosed_E
                    if(Key_Type_2D .eq. 1  )then
                      D1 = choosed_E/(ONE-c_v**2)
                      D2 = choosed_E*c_v/(ONE-c_v**2)
                      D3 = choosed_E/TWO/(ONE+c_v)
                      Weibull_Elements_D_Matrix(i_E,1,1:3)  = [D1,   D2,   ZR]
                      Weibull_Elements_D_Matrix(i_E,2,1:3)  = [D2,   D1,   ZR]
                      Weibull_Elements_D_Matrix(i_E,3,1:3)  = [ZR,   ZR,   D3]
                    elseif (Key_Type_2D .eq.2 )then
                      D1 = choosed_E*(ONE-c_v)/(ONE+c_v)/(ONE-TWO*c_v)
                      D2 = choosed_E*c_v/(ONE+c_v)/(ONE-TWO*c_v)
                      D3 = choosed_E/TWO/(ONE+c_v)
                      Weibull_Elements_D_Matrix(i_E,1,1:3)   = [D1,   D2,   ZR]
                      Weibull_Elements_D_Matrix(i_E,2,1:3)   = [D2,   D1,   ZR]
                      Weibull_Elements_D_Matrix(i_E,3,1:3)   = [ZR,   ZR,   D3]        
                    end if
                endif
            enddo
        endif
    enddo
    
    ! Save the elee file, which contains the elastic modulus value for each element.
    if (Key_Save_Nothing/=1) then
        c_File_name_1 =trim(Full_Pathname)//'.elee'
        open(101,file=c_File_name_1,status='unknown') 
            do i_E=1,Num_Elem
                write(101, '(I8,E20.12)') i_E,Weibull_Elements_E(i_E)
            enddo
        close(101)
    endif
endif
#endif

RETURN
END SUBROUTINE Get_Material_Matrix
