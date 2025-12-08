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
 
      SUBROUTINE Get_Material_Matrix_3D
c     Calculate the material matrix D (3D problem)

c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Dynamic
      use Global_Material
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      integer i_mat
      real(kind=FT) D1,D2,D3,c_E,c_v   
      real(kind=FT) E1,E2,E3,v12,v23,v13,G12,G23,G13,temp
      integer i_E
      
      
      
      print *,' '
      print *,'    Calculating material matrix D...'
      
      !2023-03-19. IMPROV2023031904.
      if (Key_XA  == 2) then 
          goto 1000
      endif
      
      ALLOCATE(D(num_of_Material,6,6))
      ALLOCATE(D_Comp(num_of_Material,6,6)) 
      ALLOCATE(S(num_of_Material,6,6))
      ALLOCATE(St(num_of_Material,2))
      ALLOCATE(Sc(num_of_Material,2))
      ALLOCATE(T_Alpha(num_of_Material))
      ALLOCATE(KIc(num_of_Material,2))
      ALLOCATE(E(num_of_Material,3))
      ALLOCATE(v(num_of_Material,3))
      ALLOCATE(density(num_of_Material))
      ALLOCATE(D_for_cylindrical(num_of_Material,6,6))
      
      do i_mat =1,num_of_Material
          !*******************
          !      ISO material
          !*******************
          if (Material_Type(i_mat). eq. 1  )then      
              c_E  = Material_Para(i_mat,1)
            c_v  = Material_Para(i_mat,2)
            density(i_mat)   = Material_Para(i_mat,3)          
            St(i_mat,1)      = Material_Para(i_mat,5)
            Sc(i_mat,1)      = Material_Para(i_mat,7)
              T_Alpha(i_mat)   = Material_Para(i_mat,8)
            KIc(i_mat,1)     = Material_Para(i_mat,6)  
              

              D1 = c_E*(ONE-c_v)/(ONE-TWO*c_v)/(ONE+c_v)
              D2 = c_E*c_v/(ONE-TWO*c_v)/(ONE+c_v)
              D3 = c_E/TWO/(ONE+c_v)
              
              D(i_mat,1,:)=[D1,    D2,    D2,   ZR,   ZR,  ZR]
              D(i_mat,2,:)=[D2,    D1,    D2,   ZR,   ZR,  ZR]
              D(i_mat,3,:)=[D2,    D2,    D1,   ZR,   ZR,  ZR]
              D(i_mat,4,:)=[ZR,    ZR,    ZR,   D3,   ZR,  ZR]
              D(i_mat,5,:)=[ZR,    ZR,    ZR,   ZR,   D3,  ZR]
              D(i_mat,6,:)=[ZR,    ZR,    ZR,   ZR,   ZR,  D3]
                                
              ! S matrix
              call Matrix_Inverse(D(i_mat,:,:),S(i_mat,:,:),6)                     
            
            E(i_mat,1) = c_E
            v(i_mat,1) = c_v
              
              !-----------------------------------------------------------------------------------
              ! D_for_cylindrical, used for strain calculation in cylindrical coordinate system
              ! Essentially, it is the inverse matrix of D
              ! Ref: heory_documents\024.2 Stress Tensor and Vector Transformation from Cartesian
              ! Coordinate System to Cylindrical Coordinate System 2.6-2021-09-10.html
              !-----------------------------------------------------------------------------------
              ! -- OPTION 1 --
C             D_for_cylindrical(i_mat,1,:)=
C    &                         [ONE,    -c_v,    -c_v,   ZR,   ZR,  ZR]
C             D_for_cylindrical(i_mat,2,:)=
C    &                         [-c_v,    ONE,    -c_v,   ZR,   ZR,  ZR]
C             D_for_cylindrical(i_mat,3,:)=
C    &                         [-c_v,    -c_v,    ONE,   ZR,   ZR,  ZR]
C             D_for_cylindrical(i_mat,4,:)=
C    &                         [ZR,    ZR,    ZR,TWO*(ONE+c_v), ZR,  ZR]
C             D_for_cylindrical(i_mat,5,:)=
C    &                         [ZR,    ZR,    ZR,ZR, TWO*(ONE+c_v),  ZR]
C             D_for_cylindrical(i_mat,6,:)=
C    &                         [ZR,    ZR,    ZR,ZR, ZR,  TWO*(ONE+c_v)]             
C             D_for_cylindrical(i_mat,:,:) = 
C    &                   D_for_cylindrical(i_mat,:,:)/c_E
              ! -- OPTION 2 (Calculate the inverse matrix of D, equivalent to OPTION 1)--
              call Matrix_Inverse(D(i_mat,1:6,1:6),
     &                            D_for_cylindrical(i_mat,1:6,1:6),6)  
          
     
          !***********************
          !    Composite material
          !***********************
          elseif(Material_Type(i_mat). eq. 5)then      
              c_E  = Material_Para(i_mat,1)
            c_v  = Material_Para(i_mat,2)
            density(i_mat)   = Material_Para(i_mat,3)          
            St(i_mat,1)      = Material_Para(i_mat,5)
            Sc(i_mat,1)      = Material_Para(i_mat,7)
              T_Alpha(i_mat)   = Material_Para(i_mat,8)
            KIc(i_mat,1)     = Material_Para(i_mat,6)  
              

              D1 = c_E*(ONE-c_v)/(ONE-TWO*c_v)/(ONE+c_v)
              D2 = c_E*c_v/(ONE-TWO*c_v)/(ONE+c_v)
              D3 = c_E/TWO/(ONE+c_v)
              
              D(i_mat,1,1:6)=[D1,    D2,    D2,   ZR,   ZR,  ZR]
              D(i_mat,2,1:6)=[D2,    D1,    D2,   ZR,   ZR,  ZR]
              D(i_mat,3,1:6)=[D2,    D2,    D1,   ZR,   ZR,  ZR]
              D(i_mat,4,1:6)=[ZR,    ZR,    ZR,   D3,   ZR,  ZR]
              D(i_mat,5,1:6)=[ZR,    ZR,    ZR,   ZR,   D3,  ZR]
              D(i_mat,6,1:6)=[ZR,    ZR,    ZR,   ZR,   ZR,  D3]
                                
              ! S matrix
              call Matrix_Inverse(D(i_mat,:,:),S(i_mat,:,:),6)                     
            
            E(i_mat,1) = c_E
            v(i_mat,1) = c_v
              !======================================
              ! Composite material (fiber material).
              !======================================
              D_Comp(i_mat,1:6,1:6) = ZR
              E1  = Material_Para_Added(i_mat,1)      
              E2  = Material_Para_Added(i_mat,2)    
              E3  = Material_Para_Added(i_mat,3)    
              v12 = Material_Para_Added(i_mat,4)    
              v23 = Material_Para_Added(i_mat,5)    
              v13 = Material_Para_Added(i_mat,6)    
              G12 = Material_Para_Added(i_mat,7)    
              G23 = Material_Para_Added(i_mat,8)    
              G13 = Material_Para_Added(i_mat,9)    
              temp = ONE-v12**2 -v23**2-v13**2 -TWO*v12*v23*v13
              D_Comp(i_mat,1,1) = E1*(ONE-v23**2)/temp
              D_Comp(i_mat,1,2) = E1*(v12+v13*v23)/temp
              D_Comp(i_mat,1,3) = E1*(v13+v12*v13)/temp
              D_Comp(i_mat,2,1) = E1*(v12+v13*v23)/temp
              D_Comp(i_mat,2,2) = E2*(ONE-v13**2)/temp
              D_Comp(i_mat,2,3) = E2*(v23+v12*v13)/temp
              D_Comp(i_mat,3,1) = E1*(v13+v12*v13)/temp
              D_Comp(i_mat,3,2) = E2*(v23+v12*v13)/temp     
              D_Comp(i_mat,3,3) = E3*(ONE-v12**2)/temp
              D_Comp(i_mat,4,4) = G12
              D_Comp(i_mat,5,5) = G23
              D_Comp(i_mat,6,6) = G13
              
              
              ! D_for_cylindrical, used for strain calculation in cylindrical coordinates, 2021-09-12
              ! Essentially the inverse matrix of D_Comp
              call Matrix_Inverse(D(i_mat,1:6,1:6),
     &                            D_for_cylindrical(i_mat,1:6,1:6),6)  
              
          !****************
          !         Others
          !****************
          else
              c_E  = Material_Para(i_mat,1)
            c_v  = Material_Para(i_mat,2)
            density(i_mat)   = Material_Para(i_mat,3)          
            St(i_mat,1)      = Material_Para(i_mat,5)
            Sc(i_mat,1)      = Material_Para(i_mat,7)
              T_Alpha(i_mat)   = Material_Para(i_mat,8)
            KIc(i_mat,1)     = Material_Para(i_mat,6)  
              

              D1 = c_E*(ONE-c_v)/(ONE-TWO*c_v)/(ONE+c_v)
              D2 = c_E*c_v/(ONE-TWO*c_v)/(ONE+c_v)
              D3 = c_E/TWO/(ONE+c_v)
              
              D(i_mat,1,:)=[D1,    D2,    D2,   ZR,   ZR,  ZR]
              D(i_mat,2,:)=[D2,    D1,    D2,   ZR,   ZR,  ZR]
              D(i_mat,3,:)=[D2,    D2,    D1,   ZR,   ZR,  ZR]
              D(i_mat,4,:)=[ZR, ZR, ZR,   D3,   ZR,  ZR]
              D(i_mat,5,:)=[ZR, ZR, ZR,ZR,   D3,     ZR]
              D(i_mat,6,:)=[ZR, ZR, ZR,ZR,   ZR,  D3   ]
                                
              ! S matrix
              call Matrix_Inverse(D(i_mat,:,:),S(i_mat,:,:),6)                     
            
            E(i_mat,1) = c_E
            v(i_mat,1) = c_v
              
              ! D_for_cylindrical, used for strain calculation in cylindrical coordinate system
              ! Essentially, it is the inverse matrix of D
              ! Ref: heory_documents\024.2 Stress Tensor and Vector Transformation from Cartesian Coordinate
              ! System to Cylindrical Coordinate System 2.6-2021-09-10.html
              call Matrix_Inverse(D(i_mat,1:6,1:6),
     &                            D_for_cylindrical(i_mat,1:6,1:6),6)  
              
          end if      
      end do

 1000 continue
 
      !IMPROV2023031904.
      ! Key_XA = 2. Dedicated to Xin'ao, calculate the D matrix for each element: Elem_D_XA(Num_Elem,6,6)
      if (Key_XA  == 2) then 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_E,c_v,
!$OMP&                                    D1,D2,D3) 
          do i_E =1,Num_Elem  
              c_E  = Elem_E_XA(i_E)
              c_v  = Elem_Mu_XA(i_E)
              D1 = c_E*(ONE-c_v)/(ONE-TWO*c_v)/(ONE+c_v)
              D2 = c_E*c_v/(ONE-TWO*c_v)/(ONE+c_v)
              D3 = c_E/TWO/(ONE+c_v)
              Elem_D_XA(i_E,1,1:6)=[D1,    D2,    D2,   ZR,   ZR,  ZR]
              Elem_D_XA(i_E,2,1:6)=[D2,    D1,    D2,   ZR,   ZR,  ZR]
              Elem_D_XA(i_E,3,1:6)=[D2,    D2,    D1,   ZR,   ZR,  ZR]
              Elem_D_XA(i_E,4,1:6)=[ZR,    ZR,    ZR,   D3,   ZR,  ZR]
              Elem_D_XA(i_E,5,1:6)=[ZR,    ZR,    ZR,   ZR,   D3,  ZR]
              Elem_D_XA(i_E,6,1:6)=[ZR,    ZR,    ZR,   ZR,   ZR,  D3]
          enddo
!$omp end parallel do  
      endif
      
      
      RETURN
      END SUBROUTINE Get_Material_Matrix_3D
