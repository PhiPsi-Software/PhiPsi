 
      subroutine Tool_Area_Polygon(x,y,nb,area)
       
      use Global_Float_Type
      IMPLICIT NONE
      integer,intent(in)::nb
      real(kind=FT) ,intent(in)::x(nb),y(nb)
      real(kind=FT) ,intent(out)::area
      INTEGER i, n, nm1
      real(kind=FT) a


      n = nb
      
      
      IF((abs(x(1)-x(n)) <=Tol_15) .AND. (abs(y(1)-y(n))<=Tol_15))then 
          n = n - 1
      endif
      
      SELECT CASE (n)
      CASE (:2)
          area = ZR
      CASE (3)
          area=HLF*((x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1)))
      CASE DEFAULT
          nm1 = n - 1
          a = x(1)*(y(2) - y(n)) + x(n)*(y(1) - y(nm1))
          DO  i = 2, nm1
              a = a + x(i)*(y(i+1) - y(i-1))
          END DO
          area = HLF*a
      END SELECT

      RETURN
      END subroutine Tool_Area_Polygon
