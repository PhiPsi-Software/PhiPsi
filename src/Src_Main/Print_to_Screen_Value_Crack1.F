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
 
      subroutine Print_to_Screen_Value_Crack1(num_Space,c_String)
c     Print content on the screen (information about Crack No. 1)
      use Global_Float_Type
      use Global_Common 
      use Global_Crack
      implicit none
      integer,intent(in):: num_Space
      character*5  c_String
      real(kind=FT),allocatable::tem_vector_real(:)
      
#ifndef Silverfrost
      if(c_String=='WIDTH')then
          if(Key_Unit_System==1)then
              print *,'    Max aperture of crack 1(mm):',
     &         maxval(Cracks_CalP_Aper(1,1:Cracks_CalP_Num(1)))
     &         *1000.0D0
          elseif(Key_Unit_System==2)then
              print *,'    Max aperture of crack 1(mm):',
     &         maxval(Cracks_CalP_Aper(1,1:Cracks_CalP_Num(1)))
          endif
      elseif(c_String=='PRESS')then    
          if(Key_Unit_System==1)then
              print *,'    Max pressure of crack 1(MPa):',
     &             maxval(Cracks_CalP_Pres(1,1:Cracks_CalP_Num(1))) 
     &                /1.0D6        
          elseif(Key_Unit_System==2)then
              print *,'    Max pressure of crack 1(MPa):',
     &             maxval(Cracks_CalP_Pres(1,1:Cracks_CalP_Num(1)))   
          endif
      elseif(c_String=='FLUXX')then  
          if(Key_Unit_System==1)then
              print *,'    Max flux of crack 1(m*2/s):',
     &            maxval(Cracks_CalP_Quan(1,1:Cracks_CalP_Num(1))) 
          elseif(Key_Unit_System==2)then
              print *,'    Max pressure of crack 1(m*2/s):',
     &            maxval(Cracks_CalP_Quan(1,1:Cracks_CalP_Num(1)))/1.0D6
          endif
      else
          print *, '    Error :: wrong input par for '
     &                             //'Print_to_Screen_Value_Crack1'
          call Warning_Message('S',Keywords_Blank) 
      endif
#endif
      
#ifdef Silverfrost
      if(c_String=='WIDTH')then
          if(Key_Unit_System==1)then
              allocate(tem_vector_real(Cracks_CalP_Num(1)))
              tem_vector_real=Cracks_CalP_Aper(1,1:Cracks_CalP_Num(1))
              print *,'    Max aperture of crack 1(mm):',
     &         maxval(tem_vector_real)*1000.0D0
              deallocate(tem_vector_real)
          elseif(Key_Unit_System==2)then
              print *,'    Max aperture of crack 1(mm):',
     &         maxval(Cracks_CalP_Aper(1,1:Cracks_CalP_Num(1)))
          endif
      elseif(c_String=='PRESS')then    
          if(Key_Unit_System==1)then
              allocate(tem_vector_real(Cracks_CalP_Num(1)))
              tem_vector_real=Cracks_CalP_Pres(1,1:Cracks_CalP_Num(1))

              print *,'    Max pressure of crack 1(MPa):',
     &             maxval(tem_vector_real)/1.0D6  
              deallocate(tem_vector_real)
          elseif(Key_Unit_System==2)then
              print *,'    Max pressure of crack 1(MPa):',
     &             maxval(Cracks_CalP_Pres(1,1:Cracks_CalP_Num(1)))   
          endif
      elseif(c_String=='FLUXX')then  
          if(Key_Unit_System==1)then
              allocate(tem_vector_real(Cracks_CalP_Num(1)))
              tem_vector_real=Cracks_CalP_Pres(1,1:Cracks_CalP_Num(1))
              
              print *,'    Max flux of crack 1(m*2/s):',
     &            maxval(tem_vector_real) 
              deallocate(tem_vector_real)
          elseif(Key_Unit_System==2)then
              print *,'    Max pressure of crack 1(m*2/s):',
     &            maxval(Cracks_CalP_Quan(1,1:Cracks_CalP_Num(1)))  
     &                 /1.0D6
          endif
      else
          print *, '    Error :: wrong input par for '
     &                             //'Print_to_Screen_Value_Crack1'
          call Warning_Message('S',Keywords_Blank) 
      endif
#endif

      return 
      end SUBROUTINE Print_to_Screen_Value_Crack1                
