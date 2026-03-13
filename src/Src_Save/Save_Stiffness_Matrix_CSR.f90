 
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
