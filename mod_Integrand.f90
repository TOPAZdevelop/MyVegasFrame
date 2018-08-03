MODULE ModIntegrand
implicit none


 CONTAINS

 


FUNCTION VegasIntegrand1(yRnd,VgsWgt,res)
implicit none
integer :: VegasIntegrand1
real(8) ::  yRnd(*),res(*),VgsWgt,x1,x2

    x1 = yRnd(1)
    x2 = yRnd(2)
    res(1) = 0d0

    if( dabs(x1).lt.10d0 .and. dabs(x2).lt.10d0  ) then
      res(1) = dexp(-x1**2)+dexp(-x2**2)
    endif
    res(1) = res(1)/70.89815404d0
    
   
VegasIntegrand1=0
RETURN
END FUNCTION





FUNCTION VegasIntegrand2(yRnd,VgsWgt,res)
implicit none
integer :: VegasIntegrand2
real(8) ::  yRnd(*),res(*),VgsWgt,x1,x2


    x1 = yRnd(1)
    x2 = yRnd(2)
    res(1) = 0d0

    if( dabs(x1+x2).lt.20d0 .and. dabs(x1-x2).lt.20d0  ) then
        res(1) = dexp(-(x1+x2)**2)+dexp(-(x1-x2)**2)
    endif
    res(1) = res(1)/70.89815404d0 

  
VegasIntegrand2=0
RETURN
END FUNCTION





SUBROUTINE VegasIntegrand1_DCUHRE(Ndim,yRnd,NumFun,Res)
integer :: Ndim,NumFun
real(8) :: yRnd(1:NDim),Res(1:NumFun),VgsWgt,Dummy(1:1)
VgsWgt=1d0

   Res(1) = VegasIntegrand1(yRnd,VgsWgt,Dummy)
   Res(1) = Dummy(1)

 
RETURN  
END SUBROUTINE  





SUBROUTINE VegasIntegrand2_DCUHRE(Ndim,yRnd,NumFun,Res)
integer :: Ndim,NumFun
real(8) :: yRnd(1:NDim),Res(1:NumFun),VgsWgt,Dummy(1:1)
VgsWgt=1d0

   Res(1) = VegasIntegrand2(yRnd,VgsWgt,Dummy)
   Res(1) = Dummy(1)

 
RETURN  
END SUBROUTINE  





FUNCTION VegasIntegrand1_FOAM(yRnd)
real(8) :: yRnd(*),xRnd(1:2),VegasIntegrand1_FOAM,VgsWgt,Dummy(1:1)
VgsWgt=1d0

   xRnd(1) = yRnd(1)*40d0-20d0
   xRnd(2) = yRnd(2)*40d0-20d0
   VegasIntegrand1_FOAM = VegasIntegrand1(yRnd,VgsWgt,Dummy)
   
   VegasIntegrand1_FOAM = Dummy(1)*40d0*40d0

RETURN  
END FUNCTION  





FUNCTION VegasIntegrand2_FOAM(yRnd)
real(8) :: yRnd(*),xRnd(1:2),VegasIntegrand2_FOAM,VgsWgt,Dummy(1:1)
VgsWgt=1d0

   xRnd(1) = yRnd(1)*40d0-20d0
   xRnd(2) = yRnd(2)*40d0-20d0
   VegasIntegrand2_FOAM = VegasIntegrand2(xRnd,VgsWgt,Dummy)
   VegasIntegrand2_FOAM = Dummy(1)*40d0*40d0

RETURN  
END FUNCTION  










FUNCTION VegasIntegrand3_FOAM(yRnd)
real(8) :: yRnd(*),VegasIntegrand3_FOAM

   VegasIntegrand3_FOAM = yRnd(1) * yRnd(2)**2

RETURN  
END FUNCTION  







FUNCTION VegasIntegrand2_CUBA(ndim,yRnd,ncomp,res,userdata,nvec,core)
implicit none
integer :: VegasIntegrand2_CUBA
integer ndim, ncomp, nvec, core,userdata
double precision yRnd(ndim,nvec), res(ncomp,nvec),x1,x2,VgsWgt,Dummy(1:1)



   x1 = yRnd(1,1)*40d0-20d0
   x2 = yRnd(2,1)*40d0-20d0

    res(1,1) = VegasIntegrand2((/x1,x2/),VgsWgt,Dummy)
    res(1,1) = dummy(1)*40*40
    
  
VegasIntegrand2_CUBA=0
RETURN
END FUNCTION






END MODULE ModIntegrand




