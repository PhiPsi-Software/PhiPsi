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
 
      SUBROUTINE Get_Contact_State_Eles_IN(i_Contact,
     &                       c_Cracks_CalP_Aper,
     &                       tem_Elem_Conta_Sta,
     &                       Yes_Contact,
     &                       Yes_Contact_Con,num_Contact_Ele)
c     Obtain the contact status of the element (read crack opening from input variable (IN)), used for the 'simple' contact detection algorithm
c     Elem_Conta_Sta(i_E, c_C) = 1, for hydraulic fracturing analysis: crack faces without proppant are in contact and closed with each other
c                                   Non-hydraulic fracturing analysis: Crack surfaces in mutual closure contact
c                              = 2, only for hydraulic fracturing analysis, the fracture surfaces are supported by proppant
C     Yes_Contact_Con indicates whether it converges
C     Yes_Contact indicates whether there is contact

c     The following procedures together make up the penalty function contact detection algorithm:
c     Cal_Contact_Jacobin.f
c     Cal_Contact_Contact_State_Gauss.f
c     Cal_Contact_PN_and_PT.f
c     Cal_Contact_Resid.f

      !**********************
      ! Read public variable
      !**********************
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Common
      use Global_Contact
      use Global_HF
      use Global_Elem_Area_Vol
      
      !**********************
      ! Variable Declaration
      !**********************
      implicit none
      integer,intent(in)::i_Contact
      integer,intent(inout)::tem_Elem_Conta_Sta(Num_Elem,Max_Num_Cr)
      real(kind=FT),intent(in)::c_Cracks_CalP_Aper(Max_Num_Cr,
     &                                             Max_Num_Cr_CalP)
      logical,intent(out)::Yes_Contact_Con,Yes_Contact
      integer,intent(out)::num_Contact_Ele
      integer i_E,num_CalP_el,i_CalP,c_CalP_G,c_CalP_L,c_C,i_C
      real(kind=FT) c_Aperture,Aper_Tol
      logical Yes_Equal
      integer c_NN(4)
      integer Last_Elem_Conta_Sta(Num_Elem,Max_Num_Cr),
     &        New_Elem_Conta_Sta(Num_Elem,Max_Num_Cr)  
      integer c_Adj_Ele
       
      !**********************************************************************************************
      ! The crack width tolerance used to detect the contact condition of the crack surfaces; if the
      ! crack width is less than Aper_Tol (note that 0 cannot be used as the tolerance)
      ! If this value is met, the corresponding point on the crack surface is considered to be in
      ! contact!
      ! Note: This parameter has a significant impact on whether crack surface contact iterations
      ! occur!!!!!!!!!!!
      !**********************************************************************************************
      ! Aper_Tol = 1.0D-5*Ave_Elem_L_Enrich
      ! Aper_Tol = 1.0D-5*Ave_Elem_L_Enrich !This parameter has a very large impact on the convergence of
      ! contact iterations (default: 1.0D-6*Ave_Elem_L_Enrich)
      Aper_Tol = Aper_Tol_Factor*Ave_Elem_L_Enrich
      
      !******************************
      ! Contact State Initialization
      !******************************
      Last_Elem_Conta_Sta(1:Num_Elem,1:Max_Num_Cr) = 0
      New_Elem_Conta_Sta(1:Num_Elem,1:Max_Num_Cr)  = 0
      Last_Elem_Conta_Sta = tem_Elem_Conta_Sta 
      
      Yes_Contact_Con = .False.
      Yes_Contact = .False.
      
      !********************************************************************************************
      ! Unit loop (do not check the crack tip element, because the crack tip opening is very small
      ! and directly checking it is prone to errors)
      !********************************************************************************************
      do i_E = 1,Num_Elem
           c_NN    = G_NN(:,i_E)
          !c_X_NODES = G_X_NODES(i_E,:)
          !c_Y_NODES = G_Y_NODES(i_E,:)
          if(i_Contact>=1) then
              !+++++++++++++++++++++++++++++++++++++++++++++++++++++
              ! Number of calculation points in the current element
              !+++++++++++++++++++++++++++++++++++++++++++++++++++++
              num_CalP_el = Ele_NumCalP(i_E)
              !tttttttttttttttttttttttttt
              if (num_CalP_el/=0) then
              endif
              !+++++++++++++++++++++++++++++++++++++++++++++++++
              ! current element calculation loop between points
              !+++++++++++++++++++++++++++++++++++++++++++++++++
              do i_CalP = 1,num_CalP_el
                  ! Get the global number of the current calculation point
                  c_CalP_G = Ele_CalPNum(i_E,i_CalP)
                  ! Obtain the crack number corresponding to the current calculation point
                  c_C = Cracks_LocalNumCalP(c_CalP_G,1)
                  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  ! If it is the element where the crack tip is located, no judgment will be made
                  ! for now.
                  ! The determination of the crack tip element is based on the condition of the
                  ! elements adjacent to the element where the crack tip is located.
                  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  !if (any(Elem_Type(i_E,c_C) == 1))then   
                  if (Elem_Type(i_E,c_C) == 1)then 
                      exit
                  end if
                  ! Get the local computation point number corresponding to the current calculation point
                  c_CalP_L = Cracks_LocalNumCalP(c_CalP_G,2)
                  ! Extract the crack opening at the current calculation point
                  c_Aperture = c_Cracks_CalP_Aper(c_C,c_CalP_L)
                  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  ! If the crack width <= 0 or less than or equal to the proppant width, then the
                  ! current element's
                  ! Contact status marked as contacted
                  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  if ((Key_Analysis_Type == 3) .or.
     &                (Key_Analysis_Type == 4) .or.
     &                (Key_Analysis_Type == 5))  then
                      ! If the crack opening is less than or equal to a small positive value (the original program used 0)
                      if (c_Aperture<=Aper_Tol) then     
                          New_Elem_Conta_Sta(i_E,c_C) = 1
                          exit
                      end if
                  ! If not hydraulic fracturing analysis and pure water pressure analysis
                  else   
                      if (c_Aperture<=Aper_Tol) then 
                          New_Elem_Conta_Sta(i_E,c_C) = 1
                          exit
                      end if                     
                  end if
              end do
          endif
      end do
      
      !*********************************************************************************************
      ! Determination of the crack tip element (judged based on the status of the elements adjacent
      ! to
      ! the element where the crack tip is located)
      !*********************************************************************************************
      do i_E = 1,Num_Elem
          do i_C=1,num_Crack
              c_Adj_Ele = TipEle_Adjacent_Ele(i_E,i_C) 
              if (c_Adj_Ele/= 0) then
                  ! If the adjacent element of the current crack tip is closed, then the crack tip element will also
                  ! close.
                  if(New_Elem_Conta_Sta(c_Adj_Ele,i_C) == 1 )then 
                      New_Elem_Conta_Sta(i_E,i_C) = 1 
                  endif
              end if
          enddo
      end do
      
      !********************
      ! Update unit status
      !********************
      tem_Elem_Conta_Sta =  New_Elem_Conta_Sta
      ! Number of contact elements
      num_Contact_Ele = 0
      do i_E = 1,Num_Elem
          if(any(tem_Elem_Conta_Sta(i_E,1:num_Crack).eq.1).eqv..True.)
     &                       then
              num_Contact_Ele = num_Contact_Ele+1
          endif
      enddo
      
      !**********************
      ! Convergence Judgment
      !**********************
      ! If the contact status of all elements does not change between two consecutive contact iteration
      ! steps, it is considered converged.
      call Matrixes_Equal_Is_Int(New_Elem_Conta_Sta,
     &                           Last_Elem_Conta_Sta,Num_Elem,2,
     &                           Yes_Equal)
      ! Check for convergence
      if(Yes_Equal .eqv. .True.) then
          Yes_Contact_Con = .True.
      end if     
      ! Check for cracks making contact
      if(sum(Elem_Conta_Sta) >= 1) then
          Yes_Contact = .True.
      end if     
      
      RETURN
      END SUBROUTINE Get_Contact_State_Eles_IN
