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
 
      SUBROUTINE Save_HF_Injection_Pressure(ifra,c_Time)
c     Save water pressure & time at injection points (per step) to injp file.

      use Global_Float_Type      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      use Global_POST
      
      implicit none
      
      integer,intent(in)::ifra
      real(kind=FT),intent(in)::c_Time
      character(200) c_File_name_1
      logical alive
      real(kind=FT) c_Inj_Press
      integer i_C
      
      if (Key_Save_Nothing==1) return 
      
      print *,'    Saving injection pressure...'
      c_File_name_1   =  trim(Full_Pathname)//'.injp'   
      
      if(ifra==1)then
          inquire(file=c_File_name_1, exist=alive)  
          if(alive.EQV..True.)then
              OPEN  (UNIT=105, FILE=c_File_name_1, STATUS='OLD') 
              CLOSE (UNIT=105, STATUS='DELETE')
          endif
      endif
      

      c_Inj_Press = ZR
      do i_C=1,num_Crack
          if (i_C == Inject_Crack_Num) then                             
              c_Inj_Press = Cracks_CalP_Pres(i_C,CalP_num_InjP_Local)   
              Last_Inj_Pres = c_Inj_Press 
              c_Inj_Press = c_Inj_Press +
     &                    Cracks_CalP_Remo_Strs(i_C,CalP_num_InjP_Local) 
              exit
          endif
      enddo

      
      open(301,file=c_File_name_1,status='unknown',
     &             position='append',action='write') 
      if(ifra==1)then
          write(301,*) '    ifra   |   time   | injection pressure'   
      endif   
      write(301, '(I10,2E20.12)') ifra,c_Time,c_Inj_Press
      close(301)   
 
      RETURN
      END SUBROUTINE Save_HF_Injection_Pressure
