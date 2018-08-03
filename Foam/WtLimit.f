
      SUBROUTINE WtLimitFind(mmout,eeps,Effic)
*//////////////////////////////////////////////////////////////////////////////////
*// Calculates wtmax for which overflow integral is below epsilon
*// The precision of the result is limited by beam size and statistics
*// Initialization
*//         CALL WtLimitStart(id1,id2,nchx,xu)
*//                  id1,id2 = histogram idents
*//                  nchx,xu = histogram bin number and upper limit
*// For each eevent
*//         CALL WtLimitFill(Wt)
*//                  Wt  =   weight under scrutinity
*// Fianalization
*//         CALL WtLimitFind(mout,eps,Effic)
*//                  mout =  output unit number (disk)
*//                  eps  =  epsilon in the definition of WtMax
*//                  Effic=  Efficiency=<wt>/WtMax
*// Extra getters:
*//         CALL WtLimitGetAveWt(AveWt)
*//         CALL WtLimitGetWtMax(WtMax)
*//                  Remember to call WtLimitFind before!!!
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER           id1,id2,iid1,iid2,nchx,nnchx
      DOUBLE PRECISION  Wt1,eeps,eps,Effic
      INTEGER           mmout,mout
      DOUBLE PRECISION  xl,xu,xxu
      DOUBLE PRECISION  Bin, Bin1, LowEdge
      DOUBLE PRECISION  AAveWt,AveWt,WtMax
      INTEGER           ib,ibX
      DOUBLE PRECISION  GLK_hi
      DOUBLE PRECISION  sum,sumWt,WtLimit,AveWt1
      LOGICAL GLK_Exist
      INTEGER           init
      SAVE
      DATA              init/0/
*     ------------------------------------------
      mout=mmout
      IF(init.EQ.0) THEN
         WRITE(   *,*) 'WtLimitFind: WARNING!!! lack of initialization' !
         WRITE(mout,*) 'WtLimitFind: WARNING!!! lack of initialization' !
         RETURN
      ENDIF
***   CALL GLK_Print(id1 )
***   CALL GLK_Print(id2 )
      eps=eeps
      sum   = 0d0
      sumWt = 0d0
      DO ib=0,nchx+1
         Bin  = GLK_hi(id1,ib)
         Bin1 = GLK_hi(id2,ib)
         sum   = sum   +Bin
         sumWt = sumWt +Bin1
      ENDDO
      IF(sum.EQ.0d0 .OR. sumWt.EQ.0d0) THEN
         WRITE(   *,*) 'WtLimitFind: zero content of histogram !!!',sum,sumWt  !
         WRITE(mout,*) 'WtLimitFind: zero content of histogram !!!',sum,sumWt  !
         RETURN
      ENDIF
      AveWt = sumWt/sum
      DO ibX=nchx+1,1,-1
         LowEdge =xl+(ibX-1d0)*(xu-xl)/nchx
         sum   = 0d0
         sumWt = 0d0
         DO ib=0,nchx+1
            Bin  = GLK_hi(id1,ib)
            Bin1 = GLK_hi(id2,ib)
            IF(ib.GE.ibX) Bin1=LowEdge*Bin
            sum   = sum   +Bin
            sumWt = sumWt +Bin1
         ENDDO
         AveWt1 = sumWt/sum
         IF( ABS(1d0-AveWt1/AveWt) .GT. eps ) GOTO 100
      ENDDO
 100  CONTINUE
      IF(ibX.EQ.nchx+1) THEN
         WtLimit=1d200
         Effic=0d0
         WRITE(   *,*) '+++++ WtLimit undefined. Higher uper limit in histogram' !
         WRITE(mout,*) '+++++ WtLimit undefined. Higher uper limit in histogram' !
      ELSEIF(ibX.EQ.1) THEN
         WtLimit=0d0
         Effic=-1d0
         WRITE(   *,*) '+++++ WtLimit undefined. Lower uper limit or more bins' !
         WRITE(mout,*) '+++++ WtLimit undefined. Lower uper limit or more bins' !
      ELSE
         WtLimit=xl+(ibX)*(xu-xl)/nchx ! We over-estimate WtLimit, under-estimate Effic
         Effic=AveWt/WtLimit
      ENDIF
      WRITE(   *,*) '00000000000000000000000000000000000000000000000000000000000000000000000' !
      WRITE(   *,*) '00-->WtLimit: No_evt., <Wt>,  WtLimit= ',sum,AveWt,WtLimit !
      WRITE(   *,*) '00-->WtLimit: EFFICIENCY <Wt>/WtLimit= ',Effic !
      WRITE(   *,*) '00000000000000000000000000000000000000000000000000000000000000000000000' !
      WRITE(mout,*) '-----------------------------------------------------------------------' !
      WRITE(mout,*) '-->FindWtLimit: content, <Wt>,  WtLimit= ',sum,AveWt,WtLimit !
      WRITE(mout,*) '-->FindWtLimit: EFFICIENCY <Wt>/WtLimit= ',Effic !
*----------------------------------------------------
      RETURN
*----------------------------------------------------
      ENTRY WtLimitFill(Wt1)
      IF(init.EQ.0) THEN
         WRITE(*,*) 'WtLimitFill: STOP, lack of initialization' !
         STOP
      ENDIF
      CALL GLK_Fil1(id1, Wt1, 1d0)
      CALL GLK_Fil1(id2, Wt1, Wt1)
      RETURN      
*----------------------------------------------------
      ENTRY WtLimitStart(iid1,iid2,nnchx,xxu)
      init=init+1
      id1=iid1
      id2=iid2
      nchx=nnchx
      xl=0d0
      xu=xxu
      IF( GLK_Exist(id1) ) THEN
         WRITE(*,*) '------ WtLimitStart: WARNING!!!! deleting id1=', id1!
         CALL GLK_Delet(id1)
      ENDIF
      IF( GLK_Exist(id2) ) THEN
         WRITE(*,*) '------ WtLimitStart: WARNING!!!! deleting id2=', id2!
         CALL GLK_Delet(id2)
      ENDIF
      CALL GLK_Book1(iid1,'MC weight  $',nchx, xl,xu)
      CALL GLK_Book1(iid2,'MC weight  $',nchx, xl,xu)
      RETURN
*----------------------------------------------------
      ENTRY WtLimitGetAveWt(AAveWt)
      AAveWt=AveWt
      RETURN
*----------------------------------------------------
      ENTRY WtLimitGetWtMax(WtMax)
      WtMax=WtLimit
      RETURN
      END                       !!!WtLimitFind
