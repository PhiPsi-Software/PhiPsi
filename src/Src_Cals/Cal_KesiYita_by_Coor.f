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
 
      subroutine Cal_KesiYita_by_Coor(Point,c_elem,Kesi,Yita)
C     Calculate local coordinate system based on coordinates
      use Global_Float_Type      
      use Global_Model
      
      implicit none
      
      real(kind=FT),intent(in)::Point(2)
      integer,intent(in)::c_elem
      real(kind=FT),intent(out)::Kesi,Yita
      
      real(kind=FT) x1,x2,x3,x4,y1,y2,y3,y4
      real(kind=FT) a1,a2,a3,a4,b1,b2,b3,b4
      
      real(kind=FT) X(2),EPS,Y(2)
      integer L
      EXTERNAL FS
      
      
      x1    = G_X_NODES(1,c_elem)
      x2    = G_X_NODES(2,c_elem)
      x3    = G_X_NODES(3,c_elem)
      x4    = G_X_NODES(4,c_elem)

      y1    = G_Y_NODES(1,c_elem)
      y2    = G_Y_NODES(2,c_elem)
      y3    = G_Y_NODES(3,c_elem)
      y4    = G_Y_NODES(4,c_elem)

      a1    = 0.25D0 * (-x1+x2+x3-x4)
      a2    = 0.25D0 * ( x1-x2+x3-x4)
      a3    = 0.25D0 * (-x1-x2+x3+x4)
      a4    = 0.25D0 * ( x1+x2+x3+x4)

      b1    = 0.25D0 * (-y1+y2+y3-y4)
      b2    = 0.25D0 * ( y1-y2+y3-y4)
      b3    = 0.25D0 * (-y1-y2+y3+y4)
      b4    = 0.25D0 * ( y1+y2+y3+y4)    
      
      X=[ZR,ZR]
      EPS=1.0D-15
      
      CALL Slove_NOLE_Gradient_Method(2,EPS,X,Y,FS,L,
     &                             a1,a2,a3,a4,b1,b2,b3,b4,
     &                             Point)
      
      Kesi=X(1)
      Yita=X(2)
      
      return 
      end SUBROUTINE Cal_KesiYita_by_Coor 
      
      SUBROUTINE Slove_NOLE_Gradient_Method(N,EPS,X,Y,FS,L,
     &                             a1,a2,a3,a4,b1,b2,b3,b4, 
     &                              Point)
      use Global_Float_Type
      implicit none    
      DIMENSION X(N),Y(N)
      real(kind=FT) X,Y,F,D,S
      real(kind=FT) a1,a2,a3,a4,b1,b2,b3,b4,Point(2)    
      real(kind=FT) EPS
      integer L,I,J,N
      EXTERNAL FS
      

      L = 1000
    5 call FS(X,N,F,Y,a1,a2,a3,a4,b1,b2,b3,b4,Point)
    
      
      if (F .GE. EPS) then
        L=L-1
        if (L==0) return
        
        D=ZR
        do J=1,n
            D=D+Y(J)*Y(J)
        enddo
      
      
        S=F/D
        do I=1,N
            X(I)=X(I)-S*Y(I)
        enddo
        goto 5
      end if      
      
      
      
      RETURN
      END SUBROUTINE Slove_NOLE_Gradient_Method
      
      SUBROUTINE FS(X,N,F,Y,a1,a2,a3,a4,b1,b2,b3,b4,Point)
      use Global_Float_Type
      implicit none
      integer N
      real(kind=FT) Point(2)
      real(kind=FT) X(N),Y(N),F,F1,F2,DF1,DF2
      real(kind=FT) a1,a2,a3,a4,b1,b2,b3,b4
      
      F1= a4+a3*X(2)+a1*X(1)+a2*X(1)*X(2)-Point(1)
      F2= b4+b3*X(2)+b1*X(1)+b2*X(1)*X(2)-Point(2)
      F = F1*F1+F2*F2
      
      DF1 = a1+a2*X(2)
      DF2 = b1+b2*X(2)
      Y(1)= TWO*(F1*DF1+F2*DF2)
      
      DF1 = a3+a2*X(1)
      DF2 = b3+b2*X(1)
      Y(2)= TWO*(F1*DF1+F2*DF2)
      RETURN
      end subroutine FS