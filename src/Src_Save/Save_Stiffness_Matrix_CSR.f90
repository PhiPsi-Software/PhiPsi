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
 
subroutine Save_Stiffness_Matrix_CSR(isub, K_CSR_NNZ, num_FreeD, XFEM_Start_DOF,&
                                     K_CSR_aa, K_CSR_ja, K_CSR_ia)
    use Global_Float_Type                                 
    use Global_Filename
    implicit none
    
    integer, intent(in) :: isub
    integer, intent(in) :: K_CSR_NNZ
    integer, intent(in) :: XFEM_Start_DOF
    integer, intent(in) :: num_FreeD
    real(kind=FT), intent(in) :: K_CSR_aa(K_CSR_NNZ)
    integer, intent(in) :: K_CSR_ja(K_CSR_NNZ)
    integer, intent(in) :: K_CSR_ia(num_FreeD+1)
    
    integer :: i, ios
    character(200) :: Filename_1, Filename_2, Filename_3, Filename_4
    character(5) :: temp
    
    print *, '    Saving CSR stiffness matrix...'
    
    write(temp,'(I5)') isub
    
    Filename_1 = trim(Full_Pathname)//'.csrn_'//ADJUSTL(temp)
    Filename_2 = trim(Full_Pathname)//'.csra_'//ADJUSTL(temp)
    Filename_3 = trim(Full_Pathname)//'.csrj_'//ADJUSTL(temp)
    Filename_4 = trim(Full_Pathname)//'.csri_'//ADJUSTL(temp)
    
    open(101, file=Filename_1, status='unknown', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Cannot create file ', trim(Filename_1)
        return
    endif
    write(101, '(I0)') K_CSR_NNZ
    write(101, '(I0)') num_FreeD
    write(101, '(I0)') XFEM_Start_DOF
    close(101)
    
    open(102, file=Filename_2, status='unknown', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Cannot create file ', trim(Filename_2)
        return
    endif
    do i = 1, K_CSR_NNZ
        write(102, '(E25.16E3)') K_CSR_aa(i)
    enddo
    close(102)
    
    open(103, file=Filename_3, status='unknown', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Cannot create file ', trim(Filename_3)
        return
    endif
    do i = 1, K_CSR_NNZ
        write(103, '(I0)') K_CSR_ja(i)
    enddo
    close(103)
    
    open(104, file=Filename_4, status='unknown', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Cannot create file ', trim(Filename_4)
        return
    endif
    do i = 1, num_FreeD + 1
        write(104, '(I0)') K_CSR_ia(i)
    enddo
    close(104)
    
    
    
    print *, '    CSR stiffness matrix saved successfully.'
    
end subroutine Save_Stiffness_Matrix_CSR
