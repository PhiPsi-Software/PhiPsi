 
      subroutine Cal_Support_Domain_of_Node(iNode,
     &                                  DOMAIN_Outline,
     &                                  m_DOMAIN_Outline,
     &                                  Domain_El,
     &                                  n_Domain_El)
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
