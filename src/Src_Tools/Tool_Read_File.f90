!-----------------------------------------------------------
! Brief: Read a structured real-valued file into a 2D array.
!
! Parameters:
!   Input:  Full_Name  - Path of the file to read
!   Input:  Case_Type  - Tag indicating the file category
!   Input:  m,n        - Dimensions of the data block
!   Output: Temp_DATA  - Array receiving the file contents
!   Output: Flag_Blank - True if the file is empty
!
! Notes:   Logs a short status line based on the case tag, then
!   probes the file to detect blank content before reading data.
!-----------------------------------------------------------

subroutine Tool_Read_File(Full_Name,Case_Type,m,n,Temp_DATA,Flag_Blank)
!     This function reads the contents of a file
use Global_Float_Type
implicit none

character*200 Full_Name
character*4  Case_Type
integer i,j,m,n,istat
logical Flag_Blank
real(kind=FT) Temp_DATA(m,n)


open(112,file=Full_Name,status='old')
Read(112,*,IOSTAT=istat)
close(112)

if (istat /= 0) then
  Flag_Blank = .True.
else
  Flag_Blank = .False.
  open(122,file=Full_Name,status='old')
  read(122,*)((Temp_DATA(i,j),j=1,n),i=1,m)
  close(122)
end if

RETURN
END SUBROUTINE Tool_Read_File
