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
 
      subroutine Tool_Read_File(Full_Name,Case_Type,m,n,Temp_DATA,
     &                          Flag_Blank)
C     This function reads the contents of a file
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