!-----------------------------------------------------------
! Brief: Extract global max/min/avg equivalent SIF across HF cracks.
!
! Parameters:
!   Output: Max_KI_eq_3D    - global maximum equivalent KI
!           Min_KI_eq_3D    - global minimum equivalent KI
!           Ave_KI_eq_3D    - global average equivalent KI
!           c_Crack_Max_Num - crack index holding the max KIeq
!           c_Vertex_Max_Num- outline vertex index holding the max KIeq
!
! Notes:   Branches on Key_3D_Slip_HF_Keep_Pressure to choose between
!          standard HF KIeq averaging and pressure-preserving variant.
!-----------------------------------------------------------

subroutine D3_Get_Max_Min_Ave_KIeq_of_HF_Cracks( Max_KI_eq_3D,Min_KI_eq_3D,Ave_KI_eq_3D, &
c_Crack_Max_Num,c_Vertex_Max_Num)
! Extract the maximum, minimum, and average equivalent stress intensity factors of the HF fracture.
! Called by the main program PhiPsi3D_Static_HF_SlipWater.f.
!2022-07-05.   
!Updated 2023-05-16.

!============================
!Read public variable module
!============================
use Global_Float_Type                                                                
use Global_Crack_Common
use Global_Crack_3D
use Global_Ragged_Array_Real_Classs

!=====================
!Variable Declaration
!=====================
implicit none
real(kind=FT),intent(out)::Max_KI_eq_3D,Min_KI_eq_3D,Ave_KI_eq_3D
integer,intent(out)::c_Crack_Max_Num,c_Vertex_Max_Num
integer i_C
integer num_Cr_Edges
real(kind=FT) c_Max,c_Min, Sum_KIeq
integer Ver_Count
integer c_Max_location

!===============
!Fracture cycle
!===============
Max_KI_eq_3D = -Con_Big_15   
Min_KI_eq_3D =  Con_Big_15
Ver_Count = 0
Sum_KIeq  = ZR

do i_C = 1,num_Crack
    num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)

    !  ! If it is an HF crack
    !  if(Crack_Type_Status_3D(i_C,1)==1) then 
    !      c_Max = maxval(KI_eq_3D(i_C)%row(1:num_Cr_Edges))         !2022-09-02.
    !      c_Max_location=MAXLOC(KI_eq_3D(i_C)%row(1:num_Cr_Edges),1)!2023-05-16.
    !      c_Min = minval(KI_eq_3D(i_C)%row(1:num_Cr_Edges))         !2022-09-02.
    !      
    !      if(c_Max > Max_KI_eq_3D) then 
    !          Max_KI_eq_3D = c_Max
    !          c_Crack_Max_Num = i_C            !2023-05-16.
    !          c_Vertex_Max_Num= c_Max_location !2023-05-16.
    !      endif
    !      if(c_Min < Min_KI_eq_3D) then 
    !          Min_KI_eq_3D = c_Min
    !      endif
    !      Ver_Count = Ver_Count +num_Cr_Edges
    !      Sum_KIeq  = Sum_KIeq + sum(KI_eq_3D(i_C)%row(1:num_Cr_Edges)) !2022-09-02.
    !  endif

    !NEWFTU-2026022701. 
    if (Key_3D_Slip_HF_Keep_Pressure==0) then
        ! If it is an HF crack.
        if(Crack_Type_Status_3D(i_C,1)==1) then 
            c_Max = maxval(KI_eq_3D(i_C)%row(1:num_Cr_Edges))
            c_Max_location=MAXLOC(KI_eq_3D(i_C)%row(1:num_Cr_Edges),1)
            c_Min = minval(KI_eq_3D(i_C)%row(1:num_Cr_Edges))

            if(c_Max > Max_KI_eq_3D) then 
                Max_KI_eq_3D = c_Max
                c_Crack_Max_Num = i_C
                c_Vertex_Max_Num= c_Max_location
            endif
            if(c_Min < Min_KI_eq_3D) then 
                Min_KI_eq_3D = c_Min
            endif
            Ver_Count = Ver_Count +num_Cr_Edges
            Sum_KIeq  = Sum_KIeq + sum(KI_eq_3D(i_C)%row(1:num_Cr_Edges))
        endif 
        !If keep the fluid pressure.
    elseif (Key_3D_Slip_HF_Keep_Pressure==1) then
        ! If it is an HF crack
        if(Crack_Type_Status_3D(i_C,1)==1) then 
            ! If it is not a post HF crack.
            if (Crack_Type_Status_3D(i_C,2) /= 2) then 
                c_Max = maxval(KI_eq_3D(i_C)%row(1:num_Cr_Edges))
                c_Max_location=MAXLOC(KI_eq_3D(i_C)%row(1:num_Cr_Edges),1)
                c_Min = minval(KI_eq_3D(i_C)%row(1:num_Cr_Edges))

                if(c_Max > Max_KI_eq_3D) then 
                    Max_KI_eq_3D = c_Max
                    c_Crack_Max_Num = i_C
                    c_Vertex_Max_Num= c_Max_location
                endif
                if(c_Min < Min_KI_eq_3D) then 
                    Min_KI_eq_3D = c_Min
                endif
                Ver_Count = Ver_Count +num_Cr_Edges
                Sum_KIeq  = Sum_KIeq + sum(KI_eq_3D(i_C)%row(1:num_Cr_Edges))
                ! If it is a post HF crack.  
            elseif (Crack_Type_Status_3D(i_C,2) == 2) then 
                !do nothing.
            endif
        endif
    endif

enddo

Ave_KI_eq_3D = Sum_KIeq/Ver_Count

return 
end SUBROUTINE D3_Get_Max_Min_Ave_KIeq_of_HF_Cracks
