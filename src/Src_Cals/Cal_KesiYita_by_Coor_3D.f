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
 
      subroutine Cal_KesiYita_by_Coor_3D(Point,c_elem,Kesi,Yita,Zeta)
C     Calculate the local coordinate system (3D) based on coordinates; the derivation can be found in my notes V3-177.
      use Global_Float_Type      
      use Global_Model
      
      implicit none
      
      real(kind=FT),intent(in)::Point(3)
      integer,intent(in)::c_elem
      real(kind=FT),intent(out)::Kesi,Yita,Zeta
      
      real(kind=FT) x1,x2,x3,x4,x5,x6,x7,x8
      real(kind=FT) y1,y2,y3,y4,y5,y6,y7,y8
      real(kind=FT) z1,z2,z3,z4,z5,z6,z7,z8
      real(kind=FT) a1,a2,a3,a4,a5,a6,a7,a8
      real(kind=FT) b1,b2,b3,b4,b5,b6,b7,b8
      real(kind=FT) c1,c2,c3,c4,c5,c6,c7,c8
      
      real(kind=FT) X(3),EPS,Y(3),ONE_EIG
      integer L
      EXTERNAL FS_3D
      
      
      x1    = G_X_NODES(1,c_elem)
      x2    = G_X_NODES(2,c_elem)
      x3    = G_X_NODES(3,c_elem)
      x4    = G_X_NODES(4,c_elem)
      x5    = G_X_NODES(5,c_elem)
      x6    = G_X_NODES(6,c_elem)
      x7    = G_X_NODES(7,c_elem)
      x8    = G_X_NODES(8,c_elem)
      
      y1    = G_Y_NODES(1,c_elem)
      y2    = G_Y_NODES(2,c_elem)
      y3    = G_Y_NODES(3,c_elem)
      y4    = G_Y_NODES(4,c_elem)
      y5    = G_Y_NODES(5,c_elem)
      y6    = G_Y_NODES(6,c_elem)
      y7    = G_Y_NODES(7,c_elem)
      y8    = G_Y_NODES(8,c_elem)
      
      z1    = G_Z_NODES(1,c_elem)
      z2    = G_Z_NODES(2,c_elem)
      z3    = G_Z_NODES(3,c_elem)
      z4    = G_Z_NODES(4,c_elem)
      z5    = G_Z_NODES(5,c_elem)
      z6    = G_Z_NODES(6,c_elem)
      z7    = G_Z_NODES(7,c_elem)
      z8    = G_Z_NODES(8,c_elem)
      
      ONE_EIG = ONE/EIG
      
      a1    = ONE_EIG * (-x1+x2-x3+x4+x5-x6+x7-x8)
      a2    = ONE_EIG * ( x1-x2+x3-x4+x5-x6+x7-x8)
      a3    = ONE_EIG * ( x1-x2-x3+x4-x5+x6+x7-x8)
      a4    = ONE_EIG * ( x1+x2-x3-x4-x5-x6+x7+x8)
      a5    = ONE_EIG * (-x1+x2+x3-x4-x5+x6+x7-x8)   
      a6    = ONE_EIG * (-x1-x2+x3+x4-x5-x6+x7+x8)
      a7    = ONE_EIG * (-x1-x2-x3-x4+x5+x6+x7+x8)
      a8    = ONE_EIG * ( x1+x2+x3+x4+x5+x6+x7+x8) 
      
      b1    = ONE_EIG * (-y1+y2-y3+y4+y5-y6+y7-y8)
      b2    = ONE_EIG * ( y1-y2+y3-y4+y5-y6+y7-y8)
      b3    = ONE_EIG * ( y1-y2-y3+y4-y5+y6+y7-y8)
      b4    = ONE_EIG * ( y1+y2-y3-y4-y5-y6+y7+y8)
      b5    = ONE_EIG * (-y1+y2+y3-y4-y5+y6+y7-y8)   
      b6    = ONE_EIG * (-y1-y2+y3+y4-y5-y6+y7+y8)
      b7    = ONE_EIG * (-y1-y2-y3-y4+y5+y6+y7+y8)
      b8    = ONE_EIG * ( y1+y2+y3+y4+y5+y6+y7+y8) 
      
      c1    = ONE_EIG * (-z1+z2-z3+z4+z5-z6+z7-z8)
      c2    = ONE_EIG * ( z1-z2+z3-z4+z5-z6+z7-z8)
      c3    = ONE_EIG * ( z1-z2-z3+z4-z5+z6+z7-z8)
      c4    = ONE_EIG * ( z1+z2-z3-z4-z5-z6+z7+z8)
      c5    = ONE_EIG * (-z1+z2+z3-z4-z5+z6+z7-z8)   
      c6    = ONE_EIG * (-z1-z2+z3+z4-z5-z6+z7+z8)
      c7    = ONE_EIG * (-z1-z2-z3-z4+z5+z6+z7+z8)
      c8    = ONE_EIG * ( z1+z2+z3+z4+z5+z6+z7+z8)    
      
      X=[ZR,ZR,ZR]
      EPS=1.0D-15
      
      CALL Slove_NOLE_Gradient_Method_3D(3,EPS,X,Y,FS_3D,L,
     &                             a1,a2,a3,a4,a5,a6,a7,a8,
     &                             b1,b2,b3,b4,b5,b6,b7,b8,
     &                             c1,c2,c3,c4,c5,c6,c7,c8,
     &                             Point)
      
      Kesi=X(1)
      Yita=X(2)
      Zeta=X(3)
      return 
      end SUBROUTINE Cal_KesiYita_by_Coor_3D 
      
      SUBROUTINE Slove_NOLE_Gradient_Method_3D(N,EPS,X,Y,FS_3D,L,
     &                             a1,a2,a3,a4,a5,a6,a7,a8,
     &                             b1,b2,b3,b4,b5,b6,b7,b8,
     &                             c1,c2,c3,c4,c5,c6,c7,c8,
     &                             Point)
      use Global_Float_Type
      implicit none    
      DIMENSION X(N),Y(N)
      real(kind=FT) X,Y,F,D,S
      real(kind=FT) a1,a2,a3,a4,a5,a6,a7,a8
      real(kind=FT) b1,b2,b3,b4,b5,b6,b7,b8
      real(kind=FT) c1,c2,c3,c4,c5,c6,c7,c8
      real(kind=FT) Point(3)    
      real(kind=FT) EPS
      integer L,I,J,N
      EXTERNAL FS_3D
      
      L=1000
    5 call fs_3d(x,n,f,y,a1,a2,a3,a4,a5,a6,a7,a8,
     &           b1,b2,b3,b4,b5,b6,b7,b8,
     &           c1,c2,c3,c4,c5,c6,c7,c8,point)
      if (F>=EPS) then
        L=L-1
        if (L==0) return
        D=ZR
        do J=1,N
            D=D+Y(J)*Y(J)
        enddo
      
        S=F/D
        do I=1,N
            X(I)=X(I)-S*Y(I)
        enddo
        goto 5
      end if      
      END SUBROUTINE Slove_NOLE_Gradient_Method_3D
      
      SUBROUTINE FS_3D(X,N,F,Y,
     &              a1,a2,a3,a4,a5,a6,a7,a8,
     &              b1,b2,b3,b4,b5,b6,b7,b8,
     &              c1,c2,c3,c4,c5,c6,c7,c8,
     &              Point)
     
      use Global_Float_Type
      implicit none
      integer N
      
      real(kind=FT) Point(3)
      real(kind=FT) X(N),Y(N),F,F1,F2,F3,DF1,DF2,DF3
      real(kind=FT) a1,a2,a3,a4,a5,a6,a7,a8
      real(kind=FT) b1,b2,b3,b4,b5,b6,b7,b8
      real(kind=FT) c1,c2,c3,c4,c5,c6,c7,c8
      
      F1= a1*X(1)*X(2)*X(3)+a2*X(1)*X(2)+a3*X(1)*X(3)+a4*X(2)*X(3)
     &   +a5*X(1) + a6*X(2) + a7*X(3) + a8 -Point(1)
      F2= b1*X(1)*X(2)*X(3)+b2*X(1)*X(2)+b3*X(1)*X(3)+b4*X(2)*X(3)
     &   +b5*X(1) + b6*X(2) + b7*X(3) + b8 -Point(2)
      F3= c1*X(1)*X(2)*X(3)+c2*X(1)*X(2)+c3*X(1)*X(3)+c4*X(2)*X(3)
     &   +c5*X(1) + c6*X(2) + c7*X(3) + c8 -Point(3)
      
      F = F1*F1+F2*F2+F3*F3
      
      DF1 = a1*X(2)*X(3)+a2*X(2)+a3*X(3)+a5
      DF2 = b1*X(2)*X(3)+b2*X(2)+b3*X(3)+b5
      DF3 = c1*X(2)*X(3)+c2*X(2)+c3*X(3)+c5
      Y(1)= TWO*(F1*DF1+F2*DF2+F3*DF3)
      
      DF1 = a1*X(1)*X(3)+a2*X(1)+a4*X(3)+a6
      DF2 = b1*X(1)*X(3)+b2*X(1)+b4*X(3)+b6
      DF3 = c1*X(1)*X(3)+c2*X(1)+c4*X(3)+c6
      Y(2)= TWO*(F1*DF1+F2*DF2+F3*DF3)
      
      DF1 = a1*X(1)*X(2)+a3*X(1)+a4*X(2)+a7
      DF2 = b1*X(1)*X(2)+b3*X(1)+b4*X(2)+b7
      DF3 = c1*X(1)*X(2)+c3*X(1)+c4*X(2)+c7
      Y(3)= TWO*(F1*DF1+F2*DF2+F3*DF3)      
      RETURN
      
      end subroutine FS_3D