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
 
      SUBROUTINE Stat_Crack_Connection(iter)
c     Determine the connectivity of the cracks.
c
c     Variable Description (see Module_Global for details):
c     integer Cracks_Cone_Num(Max_Num_Cr)
C     integer Cracks_Cone_Cr(Max_Num_Cr, Max_Num_Cone_Cr)
C     integer Cracks_Cone_NumTipCr(Max_Num_Cr)
C     integer Cracks_Cone_TipCrNum(Max_Num_Cr,2)
C     Cracks_Cone_TipType(i_C,1) -- The connected crack number of crack tip 1 of i_C
C     Cracks_Cone_TipType(i_C,2) -- The connected crack number at crack tip 2 of i_C
c     integer Cracks_Cone_TipJuEle(Max_Num_Cr,2)
c     real(kind=FT) Cracks_Cone_TipJuCor(Max_Num_Cr,2,2)
C     integer Cracks_Cone_NumMidCr(Max_Num_Cr)
C     integer Cracks_Cone_MidCrNum(Max_Num_Cr, Max_Num_Cone_Cr)
C     integer Cracks_Cone_MidCrTip(Max_Num_Cr, Max_Num_Cone_Cr)
c     integer Cracks_Cone_MidJuEle(Max_Num_Cr, Max_Num_Cone_Cr)
c     real(kind=FT) Cracks_Cone_MidJuCor(Max_Num_Cr, Max_Num_Cone_Cr, 2)

c     ****************************
c     Read public variable module
c     ****************************
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Common
      
c     *********************
c     Variable Declaration
c     *********************
      implicit none
      integer,intent(in)::iter
      integer i_C,Jun_Cr,Num_Cone_Cr
      
      integer Cr_num_Pr
      
      print *,"    Counting cracks' connections......"
      
C     ************************
C     Variable Initialization
C     ************************
      Cracks_Cone_Num(1:Num_Crack) =0
      Cracks_Cone_Cr(1:Num_Crack,1:Max_Num_Cone_Cr)       = 0
      !-----------------------------------------------------
      Cracks_Cone_NumTipCr(1:Max_Num_Cr)               = 0  
      Cracks_Cone_TipCrNum(1:Max_Num_Cr,1:2)           = 0
      Cracks_Cone_TipJuEle(1:Max_Num_Cr,1:2)           = 0
      Cracks_Cone_TipJuCor(Max_Num_Cr,1:2,1:2)         = ZR
      !-----------------------------------------------------
      Cracks_Cone_NumMidCr(1:Max_Num_Cr)               = 0
      Cracks_Cone_MidCrNum(1:Max_Num_Cr,1:Max_Num_Cone_Cr) = 0
      Cracks_Cone_MidCrTip(1:Max_Num_Cr,1:Max_Num_Cone_Cr) = 0
      Cracks_Cone_MidJuEle(1:Max_Num_Cr,1:Max_Num_Cone_Cr) = 0
      Cracks_Cone_MidJuCor(1:Max_Num_Cr,1:Max_Num_Cone_Cr,1:2)=ZR

C     *********************************************
C     Crack cycle
C     Algorithm: Count the junctions of each crack
c     *********************************************
      do i_C=1,num_Crack
          !-----------------------------------------------------------------------------------
          ! Check whether the number of crack associations exceeds the program's scale range.
          ! If it exceeds the limit, terminate the program.
          !-----------------------------------------------------------------------------------
          if (int(maxval(Cracks_Cone_Num(1:Num_Crack))*0.9) > 
     &                   Max_Num_Cone_Cr) then
              WRITE(*,1001)  
              print *, '    Error :: Fracture network is too'
     &                       //' complicated!'
              call Warning_Message('S',Keywords_Blank)           
          end if
          !------------------------------
          ! If Tip 1 is the Junction tip
          !------------------------------
          if (Crack_Tip_Type(i_C,1) ==1) then
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Get the crack numbers associated with the current crack
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              Jun_Cr = Crack_Jun_CrNum(i_C,1)
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! From the crack (connected at the crack tip)
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              Cracks_Cone_Num(i_C) = Cracks_Cone_Num(i_C) + 1
              Cracks_Cone_Cr(i_C,Cracks_Cone_Num(i_C)) =Jun_Cr
              Cracks_Cone_NumTipCr(i_C)  = Cracks_Cone_NumTipCr(i_C)+1
              Cracks_Cone_TipCrNum(i_C,1)   = Jun_Cr
              ! element number and coordinates of the junction point
              Cracks_Cone_TipJuEle(i_C,1)   = Crack_Jun_Elem(i_C,1)
              Cracks_Cone_TipJuCor(i_C,1,1:2)=Crack_Tip_Coor(i_C,1,1:2)
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Main Crack (Central Connection)
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              Cracks_Cone_Num(Jun_Cr)=Cracks_Cone_Num(Jun_Cr)+1
              Cracks_Cone_Cr(Jun_Cr,Cracks_Cone_Num(Jun_Cr))=i_C
              !Cracks_Cone_JunType(Jun_Cr,Cracks_Cone_Num(Jun_Cr))=3
              Cracks_Cone_NumMidCr(Jun_Cr) = 
     &                                  Cracks_Cone_NumMidCr(Jun_Cr) +1
              Cracks_Cone_MidCrNum(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr)) 
     &                                     = i_C
              Cracks_Cone_MidCrTip(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr)) 
     &                                     = 1
              ! element number and coordinates of the junction point
              Cracks_Cone_MidJuEle(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr)) 
     &                                     = Crack_Jun_Elem(i_C,1)
              Cracks_Cone_MidJuCor(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr),
     &                             1:2) =
     &                                 Crack_Tip_Coor(i_C,1,:)
          end if
          !------------------------------
          ! If Tip 2 is the Junction tip
          !------------------------------
          if (Crack_Tip_Type(i_C,2) ==1) then
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Get the crack numbers associated with the current crack
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              Jun_Cr = Crack_Jun_CrNum(i_C,2)
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! From the crack (connected at the crack tip)
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              Cracks_Cone_Num(i_C) = Cracks_Cone_Num(i_C) + 1
              Cracks_Cone_Cr(i_C,Cracks_Cone_Num(i_C)) =Jun_Cr
              Cracks_Cone_NumTipCr(i_C)  = Cracks_Cone_NumTipCr(i_C)+1
              Cracks_Cone_TipCrNum(i_C,2)   = Jun_Cr
              ! element number and coordinates of the junction point
              Cracks_Cone_TipJuEle(i_C,2)   = Crack_Jun_Elem(i_C,2)
              Cracks_Cone_TipJuCor(i_C,2,1:2)=Crack_Tip_Coor(i_C,2,1:2)
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Main Crack (Central Connection)
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              Cracks_Cone_Num(Jun_Cr)=Cracks_Cone_Num(Jun_Cr)+1
              Cracks_Cone_Cr(Jun_Cr,Cracks_Cone_Num(Jun_Cr))=i_C
              !Cracks_Cone_JunType(Jun_Cr,Cracks_Cone_Num(Jun_Cr))=3
              Cracks_Cone_NumMidCr(Jun_Cr) = 
     &                                  Cracks_Cone_NumMidCr(Jun_Cr) +1
              Cracks_Cone_MidCrNum(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr)) 
     &                                     = i_C
              Cracks_Cone_MidCrTip(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr)) 
     &                                     = 2
              ! element number and coordinates of the junction point
              Cracks_Cone_MidJuEle(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr)) 
     &                                     = Crack_Jun_Elem(i_C,2)
              Cracks_Cone_MidJuCor(Jun_Cr,Cracks_Cone_NumMidCr(Jun_Cr),
     &                             1:2) =
     &                                 Crack_Tip_Coor(i_C,2,:)
          end if
          
      end do
 1001 FORMAT('     !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!') 
 
C     *************************************
C     Check if it exceeds the program size
c     *************************************
 1013 FORMAT(5X,'Error :: crack ',I4,' has too many connected cracks!')
 1014 FORMAT(5X,'         The avaliable maximum number is ',I4,'!')
      ! Each crack can connect to at most Max_Num_Cone_Cr cracks.
      do i_C =1,num_Crack
          if(Cracks_Cone_NumMidCr(i_C)  > Max_Num_Cone_Cr)then
              write(*,1013) i_C
              write(*,1014) Max_Num_Cone_Cr
              call Warning_Message('S',Keywords_Blank) 
          endif
      enddo
      
      RETURN
      END SUBROUTINE Stat_Crack_Connection
