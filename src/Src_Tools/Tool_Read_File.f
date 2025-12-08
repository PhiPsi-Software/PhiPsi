 
      subroutine Tool_Read_File(Full_Name,Case_Type,m,n,Temp_DATA,
     &                          Flag_Blank)
      use Global_Float_Type
      implicit none

      character*200 Full_Name
      character*4  Case_Type
      integer i,j,m,n,istat
      logical Flag_Blank
      real(kind=FT) Temp_DATA(m,n)
      
      
      select case(Case_Type(1:4))
      case('node')    
          print *, "    Reading nodal files...."
      case('elem')    
          print *, "    Reading element files...." 
      case('boux')    
          print *, "    Reading boux files...." 
      case('bouy')    
          print *, "    Reading bouy files...." 
      case('bouz')    
          print *, "    Reading bouz files...." 
      case('focx')    
          print *, "    Reading focx files...." 
      case('focy')    
          print *, "    Reading focy files...." 
      case('focz')    
          print *, "    Reading focy files...." 
      case('blnk')    
      end select
      
      open(112,file=Full_Name,status='old')
      Read(112,*,IOSTAT=istat)
      close(112)
      
      if ( istat /= 0 ) then
          Flag_Blank = .True.
      else
          Flag_Blank = .False.
          open(122,file=Full_Name,status='old')
          read(122,*)((Temp_DATA(i,j),j=1,n),i=1,m)
          close(122)
      end if
      
      RETURN
      END SUBROUTINE Tool_Read_File