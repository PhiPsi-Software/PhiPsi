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
 
      SUBROUTINE Read_EQ_Accel_data
      ! Read the earthquake acceleration and store it in the public variable EQ_Accel_data(10000)
      
      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_Dynamic 
      
      implicit none
      
      LOGICAL alive
      character*200 temp_name
      integer Tool_Count_Lines
      real(kind=FT),ALLOCATABLE::Temp_DATA(:,:)
      logical Flag_Blank
      
      print *, "    Trying to read eqac files...." 
      
      temp_name = trim(trim(Full_Pathname)//'.eqac')
      
      inquire(file=temp_name, exist=alive)  
      
      if(alive.EQV..FALSE.)then
          print *, "    ERROR :: Can not find eqac file," 
          print *, "             which is necessary when *Key_EQ=1!" 
          call Warning_Message('S',Keywords_Blank)
      else
          num_EQ_Accel = Tool_Count_Lines(temp_name) 
          ALLOCATE( EQ_Accel_data(num_EQ_Accel))
          ALLOCATE( Temp_DATA(num_EQ_Accel,1))
          Call Tool_Read_File(temp_name,"eqac",num_EQ_Accel,1,Temp_DATA,
     &                        Flag_Blank)
          EQ_Accel_data  = Temp_DATA(:,1)
          DEALLOCATE(Temp_DATA)
      endif  
      
      RETURN
      END SUBROUTINE Read_EQ_Accel_data
