MODULE ModIntegrand
implicit none


 CONTAINS

 


FUNCTION VegasIntegrand1(yRnd,VgsWgt,res)
implicit none
integer :: VegasIntegrand1
real(8) ::  yRnd(*),res(*),VgsWgt


    res(1) = yRnd(1)*yRnd(2)


VegasIntegrand1=0
RETURN
END FUNCTION






END MODULE ModIntegrand


