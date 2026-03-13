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
 
      subroutine Cal_Support_Domain_of_Node(iNode,
     &                                  DOMAIN_Outline,
     &                                  m_DOMAIN_Outline,
     &                                  Domain_El,
     &                                  n_Domain_El)
C     This subroutine calculates the support domain of node.
      use Global_Float_Type
      use Global_Model
      
      implicit none
      integer,intent(in)::iNode
      integer,intent(out)::DOMAIN_Outline(20,2)
      integer,intent(out)::Domain_El(10)
      integer,intent(out)::m_DOMAIN_Outline,n_Domain_El
      
      integer i_E,N1,N2,N3,N4,NN(4)
      integer i,j,tem1(40,2),Uniqued_tem1(40,2),Uniqued_m
      integer Uni_Mat_Count(40)
      integer Tem_DOMAIN_Oul(40,2),cou
      
      Domain_El(1:10) =0
      DOMAIN_Outline(1:20,1:2)=0
      tem1(1:40,1:2)=0
      Tem_DOMAIN_Oul(1:40,1:2)=0
      
      n_Domain_El = 0
      do i_E = 1,Num_Elem
          N1  = Elem_Node(i_E,1)                                         
          N2  = Elem_Node(i_E,2)                                             
          N3  = Elem_Node(i_E,3)                                             
          N4  = Elem_Node(i_E,4)  
          NN  = [N1,N2,N3,N4]
          if ( ANY( NN .eq. iNode) ) then
        
              n_Domain_El = n_Domain_El+1
              Domain_El(n_Domain_El) = i_E
          end if
      end do

      do i = 1,n_Domain_El
          tem1((i-1)*4+1,:) = 
     &             [Elem_Node(Domain_El(i),1),Elem_Node(Domain_El(i),2)]
          tem1((i-1)*4+2,:) = 
     &             [Elem_Node(Domain_El(i),2),Elem_Node(Domain_El(i),3)]
          tem1((i-1)*4+3,:) = 
     &             [Elem_Node(Domain_El(i),3),Elem_Node(Domain_El(i),4)]
          tem1((i-1)*4+4,:) = 
     &             [Elem_Node(Domain_El(i),4),Elem_Node(Domain_El(i),1)]
      end do
      
      call Matrix_Sort_Row_Int(40,2,1,n_Domain_El*4,tem1)  
            
      call Matrix_Unique_Row_Int(40,2,n_Domain_El*4,tem1,
     &                             Uniqued_tem1,Uniqued_m,
     &                             Uni_Mat_Count)    
      
      j=0
      do i=1,Uniqued_m
          if(Uni_Mat_Count(i).eq.1)then
              j=j+1
              Tem_DOMAIN_Oul(j,:) = Uniqued_tem1(i,:)
          end if
      end do
      
      m_DOMAIN_Outline = j
      call Tool_Sort_by_End_to_End(20,m_DOMAIN_Outline,
     &                             Tem_DOMAIN_Oul(1:20,:),
     &                              DOMAIN_Outline,cou)
      
      return 
      end SUBROUTINE Cal_Support_Domain_of_Node                          
