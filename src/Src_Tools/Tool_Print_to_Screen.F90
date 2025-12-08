 
subroutine Tool_Print_to_Screen(print_type,Tab_number)
use Global_Float_Type
use Global_Common
use Global_Crack_Common
use Global_Crack_3D
use Global_DISP
use Global_Model
use Global_Stress
use Global_Strain

implicit none
integer,intent(in)::print_type,Tab_number
integer i_C
integer num_Cr_Edges
real(kind=FT) c_Max_KI,c_Min_KI,c_Ave_KI
real(kind=FT) c_Max_KII,c_Min_KII,c_Ave_KII
real(kind=FT) c_Max_KIII,c_Min_KIII,c_Ave_KIII
real(kind=FT) Max_Disp_x,Max_Disp_y,Min_Disp_x,Min_Disp_y,Max_Disp_z,Min_Disp_z
real(kind=FT) Max_Disp_x_Cyl,Max_Disp_y_Cyl,Min_Disp_x_Cyl,Min_Disp_y_Cyl,Max_Disp_z_Cyl,Min_Disp_z_Cyl
real(kind=FT) Max_Node_Stress_xx,Max_Node_Stress_yy,Max_Node_Stress_zz, &
              Min_Node_Stress_xx,Min_Node_Stress_yy,Min_Node_Stress_zz, &
              Max_Node_Stress_xy,Max_Node_Stress_yz,Max_Node_Stress_xz, &
              Min_Node_Stress_xy,Min_Node_Stress_yz,Min_Node_Stress_xz      
real(kind=FT) Max_Node_Strain_xx,Max_Node_Strain_yy,Max_Node_Strain_zz, &
              Min_Node_Strain_xx,Min_Node_Strain_yy,Min_Node_Strain_zz, & 
              Max_Node_Strain_xy,Max_Node_Strain_yz,Max_Node_Strain_xz, &
              Min_Node_Strain_xy,Min_Node_Strain_yz,Min_Node_Strain_xz           
real(kind=FT) Max_Node_Strain_Crr,Min_Node_Strain_Crr,                      &
              Max_Node_Strain_Ctt,Min_Node_Strain_Ctt,Max_Node_Strain_Czz,  &
              Min_Node_Strain_Czz,Max_Node_Strain_Crt,Min_Node_Strain_Crt,  &
              Max_Node_Strain_Ctz,Min_Node_Strain_Ctz,Max_Node_Strain_Crz,  &
              Min_Node_Strain_Crz
real(kind=FT) Max_Node_Stress_Crr,Min_Node_Stress_Crr,                        &
              Max_Node_Stress_Ctt,Min_Node_Stress_Ctt,Max_Node_Stress_Czz,    &
              Min_Node_Stress_Czz,Max_Node_Stress_Crt,Min_Node_Stress_Crt,    &
              Max_Node_Stress_Ctz,Min_Node_Stress_Ctz,Max_Node_Stress_Crz,    &
              Min_Node_Stress_Crz   
              
real(kind=FT), allocatable:: tem_vector(:)
   
5001 FORMAT(5X,'Ave_KI of_Crack   ',I4,': ',E13.6,' MPa*m^(1/2)')    
5002 FORMAT(5X,'Min_KI of_Crack   ',I4,': ',E13.6,' MPa*m^(1/2)')  
5003 FORMAT(5X,'Max_KI of crack   ',I4,': ',E13.6,' MPa*m^(1/2)')   
5011 FORMAT(5X,'Ave_KII of crack  ',I4,': ',E13.6,' MPa*m^(1/2)')        
5012 FORMAT(5X,'Min_KII of crack  ',I4,': ',E13.6,' MPa*m^(1/2)')  
5013 FORMAT(5X,'Max_KII of crack  ',I4,': ',E13.6,' MPa*m^(1/2)')   
5021 FORMAT(5X,'Ave_KIII of crack ',I4,': ',E13.6,' MPa*m^(1/2)')        
5022 FORMAT(5X,'Min_KIII of crack ',I4,': ',E13.6,' MPa*m^(1/2)')  
5023 FORMAT(5X,'Max_KIII of crack ',I4,': ',E13.6,' MPa*m^(1/2)')   

1021 FORMAT(5X,'Range of displacement x:  ',E16.9,' m to ',E16.9,' m')  
1022 FORMAT(5X,'Range of displacement y:  ',E16.9,' m to ',E16.9,' m')  
1023 FORMAT(5X,'Range of displacement z:  ',E16.9,' m to ',E16.9,' m')     

1121 FORMAT(5X,'Range of displacement rr:',E14.5,' m to ',E14.5,' m')  
1122 FORMAT(5X,'Range of displacement θθ:',E14.5,' m to ',E14.5,' m')  
1123 FORMAT(5X,'Range of displacement zz:',E14.5,' m to ',E14.5,' m')

1721 FORMAT(5X,'Range of stress xx:  ',E13.5,' MPa to ',E13.6,' MPa')  
1722 FORMAT(5X,'Range of stress yy:  ',E13.5,' MPa to ',E13.6,' MPa') 
1723 FORMAT(5X,'Range of stress zz:  ',E13.5,' MPa to ',E13.6,' MPa')  
1724 FORMAT(5X,'Range of stress xy:  ',E13.5,' MPa to ',E13.6,' MPa')  
1725 FORMAT(5X,'Range of stress yz:  ',E13.5,' MPa to ',E13.6,' MPa') 
1726 FORMAT(5X,'Range of stress xz:  ',E13.5,' MPa to ',E13.6,' MPa')  

1821 FORMAT(5X,'Range of strain xx:  ',F9.2,' με to ',F9.2,' με')  
1822 FORMAT(5X,'Range of strain yy:  ',F9.2,' με to ',F9.2,' με') 
1823 FORMAT(5X,'Range of strain zz:  ',F9.2,' με to ',F9.2,' με')  
1824 FORMAT(5X,'Range of strain xy:  ',F9.2,' με to ',F9.2,' με')  
1825 FORMAT(5X,'Range of strain yz:  ',F9.2,' με to ',F9.2,' με') 
1826 FORMAT(5X,'Range of strain xz:  ',F9.2,' με to ',F9.2,' με')   

1921 FORMAT(5X,'Range of strain rr:  ',F9.2,' με to ',F9.2,' με') 
1922 FORMAT(5X,'Range of strain θθ:  ',F9.2,' με to ',F9.2,' με') 
1923 FORMAT(5X,'Range of strain zz:  ',F9.2,' με to ',F9.2,' με')  
1924 FORMAT(5X,'Range of strain rθ:  ',F9.2,' με to ',F9.2,' με') 
1925 FORMAT(5X,'Range of strain θz:  ',F9.2,' με to ',F9.2,' με')    
1926 FORMAT(5X,'Range of strain rz:  ',F9.2,' με to ',F9.2,' με')     

1931 FORMAT(5X,'Range of stress rr:  ',E13.6,' MPa to ',E13.6,' MPa')  
1932 FORMAT(5X,'Range of stress θθ:  ',E13.6,' MPa to ',E13.6,' MPa')  
1933 FORMAT(5X,'Range of stress zz:  ',E13.6,' MPa to ',E13.6,' MPa')  
1934 FORMAT(5X,'Range of stress rθ:  ',E13.6,' MPa to ',E13.6,' MPa')  
1935 FORMAT(5X,'Range of stress θz:  ',E13.6,' MPa to ',E13.6,' MPa')  
1936 FORMAT(5X,'Range of stress rz:  ',E13.6,' MPa to ',E13.6,' MPa')    



select case(print_type)

case(1)
#ifndef Silverfrost
    if(Key_Print_SIFs_to_Screen==1) then
        if(Key_Dimension==3) then
            do i_C = 1,num_Crack
                num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
                c_Ave_KI = sum(KI_3D(i_C)%row(1:num_Cr_Edges))/num_Cr_Edges
                c_Min_KI = minval(KI_3D(i_C)%row(1:num_Cr_Edges))  
                c_Max_KI = maxval(KI_3D(i_C)%row(1:num_Cr_Edges))      
                c_Ave_KII = sum(KII_3D(i_C)%row(1:num_Cr_Edges))/num_Cr_Edges
                c_Min_KII = minval(KII_3D(i_C)%row(1:num_Cr_Edges))  
                c_Max_KII = maxval(KII_3D(i_C)%row(1:num_Cr_Edges))   
                c_Ave_KIII = sum(KIII_3D(i_C)%row(1:num_Cr_Edges))/num_Cr_Edges
                c_Min_KIII = minval(KIII_3D(i_C)%row(1:num_Cr_Edges))  
                c_Max_KIII = maxval(KIII_3D(i_C)%row(1:num_Cr_Edges))   
                if(Tab_number==5)then
                    write(*,5001) i_C,c_Ave_KI/Cn_M
                    write(*,5002) i_C,c_Min_KI/Cn_M
                    write(*,5003) i_C,c_Max_KI/Cn_M
                    write(*,5011) i_C,c_Ave_KII/Cn_M
                    write(*,5012) i_C,c_Min_KII/Cn_M
                    write(*,5013) i_C,c_Max_KII/Cn_M
                    write(*,5021) i_C,c_Ave_KIII/Cn_M
                    write(*,5022) i_C,c_Min_KIII/Cn_M
                    write(*,5023) i_C,c_Max_KIII/Cn_M
                elseif(Tab_number==7)then
                    
                else
                    print *,"     Nothing to be printed in Tool_Print_to_Screen.f90!"
                endif
            enddo
        endif
    endif
#endif

case(2)
    if(Key_Dimension==3) then
#ifndef Silverfrost
        Max_Disp_x = maxval(DISP(1:3*Num_Node:3))
        Min_Disp_x = minval(DISP(1:3*Num_Node:3))
        Max_Disp_y = maxval(DISP(2:3*Num_Node:3))
        Min_Disp_y = minval(DISP(2:3*Num_Node:3))
        Max_Disp_z = maxval(DISP(3:3*Num_Node:3))
        Min_Disp_z = minval(DISP(3:3*Num_Node:3))
#endif   
     
#ifdef Silverfrost
        allocate(tem_vector(Num_Node))
        
        tem_vector = DISP(1:3*Num_Node:3)
        Max_Disp_x = maxval(tem_vector)
        tem_vector = DISP(1:3*Num_Node:3)
        Min_Disp_x = minval(tem_vector)
        
        tem_vector = DISP(2:3*Num_Node:3)
        Max_Disp_y = maxval(tem_vector)
        tem_vector = DISP(2:3*Num_Node:3)
        Min_Disp_y = minval(tem_vector)
        
        tem_vector = DISP(3:3*Num_Node:3)
        Max_Disp_z = maxval(tem_vector)
        tem_vector = DISP(3:3*Num_Node:3)
        Min_Disp_z = minval(tem_vector)
        
        deallocate(tem_vector)
#endif
        
        if(Tab_number==5)then
            WRITE(*,1021) Min_Disp_x,Max_Disp_x
            WRITE(*,1022) Min_Disp_y,Max_Disp_y
            WRITE(*,1023) Min_Disp_z,Max_Disp_z   
        elseif(Tab_number==7)then
            
        else
            print *,"     Nothing to be printed in Tool_Print_to_Screen.f90!"
        endif
    endif

case(3)
    if(Key_Dimension==3) then
#ifndef Silverfrost
        Max_Disp_x_Cyl = maxval(DISP_Cylinder(1:Total_FD:3))
        Min_Disp_x_Cyl = minval(DISP_Cylinder(1:Total_FD:3))
        Max_Disp_y_Cyl = maxval(DISP_Cylinder(2:Total_FD:3))
        Min_Disp_y_Cyl = minval(DISP_Cylinder(2:Total_FD:3))
        Max_Disp_z_Cyl = maxval(DISP_Cylinder(3:Total_FD:3))
        Min_Disp_z_Cyl = minval(DISP_Cylinder(3:Total_FD:3))
#endif
        if(Tab_number==5)then
            WRITE(*,1121) Min_Disp_x_Cyl,Max_Disp_x_Cyl
            WRITE(*,1122) Min_Disp_y_Cyl,Max_Disp_y_Cyl
            WRITE(*,1123) Min_Disp_z_Cyl,Max_Disp_z_Cyl    
        elseif(Tab_number==7)then
            
        else
            print *,"     Nothing to be printed in Tool_Print_to_Screen.f90!"
        endif
    endif

case(4)
    if(Key_Dimension==3) then
        Max_Node_Stress_xx = maxval(Stress_xx_Node)
        Min_Node_Stress_xx = minval(Stress_xx_Node)
        Max_Node_Stress_yy = maxval(Stress_yy_Node)
        Min_Node_Stress_yy = minval(Stress_yy_Node)
        Max_Node_Stress_zz = maxval(Stress_zz_Node)
        Min_Node_Stress_zz = minval(Stress_zz_Node)
        Max_Node_Stress_xy = maxval(Stress_xy_Node)
        Min_Node_Stress_xy = minval(Stress_xy_Node)
        Max_Node_Stress_yz = maxval(Stress_yz_Node)
        Min_Node_Stress_yz = minval(Stress_yz_Node)
        Max_Node_Stress_xz = maxval(Stress_xz_Node)
        Min_Node_Stress_xz = minval(Stress_xz_Node)     
        if(Tab_number==5)then
            WRITE(*,1721)Min_Node_Stress_xx/Cn_M,Max_Node_Stress_xx/Cn_M
            WRITE(*,1722)Min_Node_Stress_yy/Cn_M,Max_Node_Stress_yy/Cn_M
            WRITE(*,1723)Min_Node_Stress_zz/Cn_M,Max_Node_Stress_zz/Cn_M    
            WRITE(*,1724)Min_Node_Stress_xy/Cn_M,Max_Node_Stress_xy/Cn_M
            WRITE(*,1725)Min_Node_Stress_yz/Cn_M,Max_Node_Stress_yz/Cn_M
            WRITE(*,1726)Min_Node_Stress_xz/Cn_M,Max_Node_Stress_xz/Cn_M     
        elseif(Tab_number==7)then
            
        else
            print *,"     Nothing to be printed in Tool_Print_to_Screen.f90!"
        endif
    endif

case(5)
    if(Key_Dimension==3) then
        Max_Node_Strain_xx = maxval(Strain_xx_Node)
        Min_Node_Strain_xx = minval(Strain_xx_Node)
        Max_Node_Strain_yy = maxval(Strain_yy_Node)
        Min_Node_Strain_yy = minval(Strain_yy_Node)
        Max_Node_Strain_zz = maxval(Strain_zz_Node)
        Min_Node_Strain_zz = minval(Strain_zz_Node)
        Max_Node_Strain_xy = maxval(Strain_xy_Node)
        Min_Node_Strain_xy = minval(Strain_xy_Node)
        Max_Node_Strain_yz = maxval(Strain_yz_Node)
        Min_Node_Strain_yz = minval(Strain_yz_Node)
        Max_Node_Strain_xz = maxval(Strain_xz_Node)
        Min_Node_Strain_xz = minval(Strain_xz_Node)             
        if(Tab_number==5)then
            WRITE(*,1821)Min_Node_Strain_xx*Cn_M,Max_Node_Strain_xx*Cn_M
            WRITE(*,1822)Min_Node_Strain_yy*Cn_M,Max_Node_Strain_yy*Cn_M
            WRITE(*,1823)Min_Node_Strain_zz*Cn_M,Max_Node_Strain_zz*Cn_M    
            WRITE(*,1824)Min_Node_Strain_xy*Cn_M,Max_Node_Strain_xy*Cn_M
            WRITE(*,1825)Min_Node_Strain_yz*Cn_M,Max_Node_Strain_yz*Cn_M
            WRITE(*,1826)Min_Node_Strain_xz*Cn_M,Max_Node_Strain_xz*Cn_M     
        elseif(Tab_number==7)then
            
        else
            print *,"     Nothing to be printed in Tool_Print_to_Screen.f90!"
        endif
    endif

case(6)
    if(Key_Dimension==3) then
        Max_Node_Stress_Crr = maxval(Stress_Crr_Node)
        Min_Node_Stress_Crr = minval(Stress_Crr_Node)
        Max_Node_Stress_Ctt = maxval(Stress_Ctt_Node)
        Min_Node_Stress_Ctt = minval(Stress_Ctt_Node)
        Max_Node_Stress_Czz = maxval(Stress_Czz_Node)
        Min_Node_Stress_Czz = minval(Stress_Czz_Node)
        Max_Node_Stress_Crt = maxval(Stress_Crt_Node)
        Min_Node_Stress_Crt = minval(Stress_Crt_Node)
        Max_Node_Stress_Ctz = maxval(Stress_Ctz_Node)
        Min_Node_Stress_Ctz = minval(Stress_Ctz_Node)
        Max_Node_Stress_Crz = maxval(Stress_Crz_Node)
        Min_Node_Stress_Crz = minval(Stress_Crz_Node)            
        if(Tab_number==5)then
            WRITE(*,1931)Min_Node_Stress_Crr/Cn_M,Max_Node_Stress_Crr/Cn_M
            WRITE(*,1932)Min_Node_Stress_Ctt/Cn_M,Max_Node_Stress_Ctt/Cn_M
            WRITE(*,1933)Min_Node_Stress_Czz/Cn_M,Max_Node_Stress_Czz/Cn_M    
            WRITE(*,1934)Min_Node_Stress_Crt/Cn_M,Max_Node_Stress_Crt/Cn_M
            WRITE(*,1935)Min_Node_Stress_Ctz/Cn_M,Max_Node_Stress_Ctz/Cn_M
            WRITE(*,1936)Min_Node_Stress_Crz/Cn_M,Max_Node_Stress_Crz/Cn_M      
        elseif(Tab_number==7)then
            
        else
            print *,"     Nothing to be printed in Tool_Print_to_Screen.f90!"
        endif
    endif

case(7)
    if(Key_Dimension==3) then
        Max_Node_Strain_Crr = maxval(Strain_Crr_Node)
        Min_Node_Strain_Crr = minval(Strain_Crr_Node)
        Max_Node_Strain_Ctt = maxval(Strain_Ctt_Node)
        Min_Node_Strain_Ctt = minval(Strain_Ctt_Node)
        Max_Node_Strain_Czz = maxval(Strain_Czz_Node)
        Min_Node_Strain_Czz = minval(Strain_Czz_Node)
        Max_Node_Strain_Crt = maxval(Strain_Crt_Node)
        Min_Node_Strain_Crt = minval(Strain_Crt_Node)
        Max_Node_Strain_Ctz = maxval(Strain_Ctz_Node)
        Min_Node_Strain_Ctz = minval(Strain_Ctz_Node)
        Max_Node_Strain_Crz = maxval(Strain_Crz_Node)
        Min_Node_Strain_Crz = minval(Strain_Crz_Node)             
        if(Tab_number==5)then
            WRITE(*,1921)Min_Node_Strain_Crr*Cn_M,Max_Node_Strain_Crr*Cn_M
            WRITE(*,1922)Min_Node_Strain_Ctt*Cn_M,Max_Node_Strain_Ctt*Cn_M
            WRITE(*,1923)Min_Node_Strain_Czz*Cn_M,Max_Node_Strain_Czz*Cn_M    
            WRITE(*,1924)Min_Node_Strain_Crt*Cn_M,Max_Node_Strain_Crt*Cn_M
            WRITE(*,1925)Min_Node_Strain_Ctz*Cn_M,Max_Node_Strain_Ctz*Cn_M
            WRITE(*,1926)Min_Node_Strain_Crz*Cn_M,Max_Node_Strain_Crz*Cn_M      
        elseif(Tab_number==7)then
            
        else
            print *,"     Nothing to be printed in Tool_Print_to_Screen.f90!"
        endif
    endif
end select
RETURN
end subroutine Tool_Print_to_Screen