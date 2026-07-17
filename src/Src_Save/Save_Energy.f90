!------------------------------------------------------------------------
! Brief: Append the current energy quantities to the ener file.
!
! Parameters:
!   Input: iter                     - current time step index
!   Input: c_Elastic_Strain_Energy  - Elastic_Strain_Energy
!   Input: c_Fracture_Energy        - Fracture_Energy
!   Input: c_External_Work          - External_Work
!   Input: c_Residual_Energy        - Residual_Energy
!   Input: c_Normalized_Residual    - Normalized_Residual_Energy
!
! Notes:   Writes a header on iter=1 and appends one row per step.
!------------------------------------------------------------------------

SUBROUTINE Save_Energy(iter,c_Elastic_Strain_Energy,c_Fracture_Energy,c_External_Work, &
                       c_Residual_Energy,c_Normalized_Residual_Energy)
use Global_Float_Type
use Global_Filename
use Global_Common
use Global_POST
use Global_Crack_Common
use Global_Crack

implicit none

integer, intent(in) :: iter
real(kind=FT), intent(in) :: c_Elastic_Strain_Energy,c_Fracture_Energy,c_External_Work, &
                             c_Residual_Energy,c_Normalized_Residual_Energy

character(200) :: Filename_1
logical :: file_exists
integer i_C,n_pts,i_pt
real(kind=FT) Crack_Propagation_Length

if (Key_Save_Nothing==1) return

Filename_1 = trim(Full_Pathname)//'.ener'

do i_C = 1, num_Crack
    n_pts = Each_Cr_Poi_Num(i_C)
    if (n_pts < 2) cycle
    
    Crack_Propagation_Length = ZR

    do i_pt = 1, n_pts - 1
        Crack_Propagation_Length = Crack_Propagation_Length + &
            sqrt((Crack_Coor(i_C, i_pt+1, 1) - Crack_Coor(i_C, i_pt, 1))**2 + &
                 (Crack_Coor(i_C, i_pt+1, 2) - Crack_Coor(i_C, i_pt, 2))**2)
    end do

    if(Crack_Propagation_Length>=(Cracks_Growth_Length_Last_Step(i_C)-Tol_5)) then
        Crack_Propagation_Length = Crack_Propagation_Length - Cracks_Growth_Length_Last_Step(i_C)
    endif
    
    if(Crack_Propagation_Length >= (Initial_Cracks_Length(i_C)-Tol_5)) then
        Crack_Propagation_Length = Crack_Propagation_Length - Initial_Cracks_Length(i_C)
    endif
end do

inquire(file=Filename_1, exist=file_exists)

if (iter==1) then
    open(101,file=Filename_1,status='replace')
    write(101,'(A)') '# iter        Elastic_Strain_Energy Fracture_Energy       External_Work         Residual_Energy      '// &
         'Normalized_Residual    Crack_Propagation_Length'
    write(101,'(I10,6(2X,E20.12))') iter,c_Elastic_Strain_Energy,c_Fracture_Energy,c_External_Work, &
                       c_Residual_Energy,c_Normalized_Residual_Energy,Crack_Propagation_Length
    close(101)
else
    open(101,file=Filename_1,status='unknown',position='append')
    write(101,'(I10,6(2X,E20.12))') iter,c_Elastic_Strain_Energy,c_Fracture_Energy,c_External_Work, &
                       c_Residual_Energy,c_Normalized_Residual_Energy,Crack_Propagation_Length
    close(101)
endif

END SUBROUTINE Save_Energy
