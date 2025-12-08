 
      subroutine Cal_Length_of_ELe_by_Coors_Ave(x,y,Ave_Length)

      use Global_Float_Type
      use Global_Model
      
      implicit none
      
      real(kind=FT) x,y
      real(kind=FT),intent(out)::Ave_Length
      integer OUT_Elem
      integer NN(4)
      real(kind=FT) area

      call Cal_Ele_Num_by_Coors(x,y,OUT_Elem) 
                                         
      NN  = Elem_Node(OUT_Elem,1:4)
      
      call Tool_Area_Polygon(Coor(NN,1),Coor(NN,2),4,area)
      
      Ave_Length = sqrt(area)
      
      return 
      end SUBROUTINE Cal_Length_of_ELe_by_Coors_Ave                       
