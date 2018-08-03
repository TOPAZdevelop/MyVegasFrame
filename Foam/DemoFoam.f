*////////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                        //
*//                 Test of Foam   (S. Jadach, November 2000)                              //
*//                                                                                        //
*//   To execute this test:       make DemoFoam
*//                               make -f Makefile DemoFoam
*//	                          make -f Makefile DemoFoamMap
*//
*////////////////////////////////////////////////////////////////////////////////////////////
      PROGRAM Main
      IMPLICIT NONE
      CALL Main2 ! 2-dim test
      END


      SUBROUTINE Main2
*////////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                        //
*//   Test of Foam with 2-dimensional testing function                                     //
*//   gmake TestFoam
*//   gmake DemoFoamMap
*//                                                                                        //
*////////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION   Density2
      EXTERNAL           Density2
      INTEGER            nBuf,nSampl,OptDrive,k
      DOUBLE PRECISION   MCresult,MCerror,MCwt,MCvector(10)
      INTEGER            nDim ! dimension of simlical subspace
      INTEGER            kDim ! dimension of hyp-cubical subspace
      INTEGER            loop,nevtot,OptEdge,eps,Effic
      DOUBLE PRECISION   Xsec,Error
*
      OPEN(16,FILE='output-DemoFoam')
*=========================================================
      nDim      = 2
      kDim      = 0
      nBuf      = 5000          ! Buffer length <5000, default = 1000
      nSampl    = 200           ! Number of MC events in exploration of a cell, default = 200
      nBuf      = 150           ! SPECIAL for ploting
      nSampl    = 1000          ! SPECIAL for ploting
      OptDrive  = 2             ! Type of Drive =0,1,2 for TrueVol,Sigma,WtMax, default = 2
      OptEdge   = 1             ! Include vertices in sampling or not, =0,1, default =1
*=========================================================
      CALL GLK_Book1(8100,'WT distribution  $', 80, 0.0d0,2d0)
*=======================================================================
      CALL FoamA_SetnDim(          nDim)
      CALL FoamA_SetkDim(          kDim)
      CALL FoamA_SetnBuf(          nBuf)
      CALL FoamA_SetnSampl(      nSampl)
      CALL FoamA_SetOptDrive(  OptDrive)
      CALL FoamA_SetOptEdge(    OptEdge)
      CALL FoamA_SetChat(1)
*------------------------------------
      CALL FoamA_Initialize(Density2)
*------------------------------------
      WRITE(*,*) ' ####################################################################'
      WRITE(*,*) ' ############## back in main program ################################'
      WRITE(*,*) ' ####################################################################'
*------------------------------------
      CALL WtLimitStart(1004,1005,1000,10d0)  
*------------------------------------
      nevtot = 200000
      nevtot = 1000000
      DO loop = 1, nevtot
*------------------------------------
         CALL FoamA_MakeEvent(Density2)    ! generate MC event
         CALL FoamA_GetMCvector(MCvector)  ! get MC event, vector
         CALL FoamA_GetMCwt(MCwt)          ! get MC weight
         IF(loop.LE.20) THEN
            WRITE(*,'(a,10f10.5)') '##### MCwt,MCvector=',MCwt,(MCvector(k),k=1,nDim) !
         ENDIF
*------------------------------------
         CALL WtLimitFill(MCwt)
         CALL GLK_Fil1(8100,MCwt,1d0)
      ENDDO
      WRITE(*,*) ' ####################################################################' !
      WRITE(*,*) ' ############## generation finished  ################################' !
      WRITE(*,*) ' ####################################################################' !
*------------------------------------
      CALL FoamA_Finalize(MCresult,MCerror)
*------------------------------------
      WRITE(*,'(a,g19.9,a,g19.9,a,f7.5)') 
     $     'MCresult= ',MCresult,' +- ',MCerror,'         RelErr= ',MCerror/MCresult   !
*------------------------------------
      eps = 5d-4
      CALL WtLimitFind(16,eps,Effic)
*------------------------------------
* Remember that nbuf=<500, because om TeX memory
      OPEN(11, FILE='demo-cell-map.txp')
      CALL FoamA_PltBegin(-11)
      IF( nBuf.LE.2500)  CALL FoamA_PltVert
      IF( nBuf.LE. 500)  CALL FoamA_PltCell
      CALL FoamA_PltEnd
*------------------------------------
      CALL FoamA_BufActPrint(16)
      CALL GLK_Print(8100)
      WRITE(*,*) ' ####### Main2 finished ######## '
      END                       ! Main2


      DOUBLE PRECISION FUNCTION Density2(Xarg)
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  Xarg(*)
      DOUBLE PRECISION  TFun2Sphere
      DOUBLE PRECISION  TFun2Diagon
      DOUBLE PRECISION  TFun2Void
*
cc      Density2  = TFun2Diagon(Xarg) ! (f_a)
      Density2  = TFun2Sphere(Xarg) ! (f_b)
cc      Density2  = TFun2Void(Xarg) ! (f_c)
      END                       !!!Density2


      DOUBLE PRECISION FUNCTION Density3(Xarg)
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  Xarg(*)
      DOUBLE PRECISION  TFun3Sphere,RFun3Sphere,TFun3Diagon,TFun3Void,RFun3Void
*----------------------------------------------
      Density3 = TFun3Diagon(Xarg) ! Cylinder along diagonal
cc      Density3 = TFun3Sphere(Xarg) ! Thin sphere
cc      Density3  =TFun3Void(Xarg)   ! Big cubical void, efficiency 0.34
      END                       !!!Density3


*///////////////////////////////////////////////////////////////////////////////////////////////////
*//               Collection of testing functions                                                 //
*//                                                                                               //
*///////////////////////////////////////////////////////////////////////////////////////////////////

      DOUBLE PRECISION FUNCTION TFun3Void(x)
*//////////////////////////////////////////////////////////////////////////////////
*//  Big empty void inside cube, only near surface thin uniform layer of density //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  x(*)
      DOUBLE PRECISION  Fun
      DOUBLE PRECISION  Thickness
      DATA   Thickness / 0.050 /             ! Thickness of the cube
*---------------------------------------------------------------
      Fun=0d0
      IF( x(1).LT.Thickness) Fun=1/Thickness
      IF( x(2).LT.Thickness) Fun=1/Thickness
      IF( x(3).LT.Thickness) Fun=1/Thickness
      IF( (1d0-x(1)).LT.Thickness) Fun=1/Thickness
      IF( (1d0-x(2)).LT.Thickness) Fun=1/Thickness
      IF( (1d0-x(3)).LT.Thickness) Fun=1/Thickness
      TFun3Void = Fun
      END

      DOUBLE PRECISION FUNCTION TFun2Void(x)
*//////////////////////////////////////////////////////////////////////////////////
*//  Big empty void inside cube, only near surface thin uniform layer of density //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  x(*)
      DOUBLE PRECISION  Fun
      DOUBLE PRECISION  Thickness
      DATA   Thickness / 0.050 /             ! Thickness of the cube walls
*---------------------------------------------------------------
      Fun=0d0
      IF( x(1).LT.Thickness) Fun=1/Thickness
      IF( x(2).LT.Thickness) Fun=1/Thickness
      IF( (1d0-x(1)).LT.Thickness) Fun=1/Thickness
      IF( (1d0-x(2)).LT.Thickness) Fun=1/Thickness
      TFun2Void = Fun
      END


      DOUBLE PRECISION FUNCTION TFun3Diagon(x)
*//////////////////////////////////////////////////////////////////////////////////
*//  Density only within cylinder along diagonal
*//  Slight parallel shift off diagonal, see vector A
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  x(*)
      DOUBLE PRECISION  Pi
      PARAMETER( Pi=3.1415926535897932d0)
      DOUBLE PRECISION  Thickness,R,Fun
      DATA   Thickness / 0.020 /             ! radius of cylinder along diagonal
      DOUBLE PRECISION  A1,A2,A3
      DATA   A1,A2,A3  / 0.15, 0.10, -0.15 / ! slight parallel shift from diagonal
*---------------------------------------------------------------
      R   = (x(1)-x(2) -A1+A2)**2 +(x(2)-x(3) -A2+A3)**2 +(x(3)-x(1) -A3+A1)**2  ! distance from diagonal
      Fun = Thickness/( R + Thickness**2)/Pi ! Breit-Wigner profile, normalized to delta
      TFun3Diagon = Fun*(5d0+x(3))           ! 20% tilt along z-axis
      END


      DOUBLE PRECISION FUNCTION TFun2Diagon(x)
*//////////////////////////////////////////////////////////////////////////////////
*//  Density only within cylinder along diagonal
*//  Slight parallel shift off diagonal, see vector A
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  x(*)
      DOUBLE PRECISION  Pi
      PARAMETER( Pi=3.1415926535897932d0)
      DOUBLE PRECISION  Thickness,R,Fun
      DATA   Thickness / 0.020 /             ! radius of cylinder along diagonal
      DOUBLE PRECISION  A1,A2
      DATA   A1,A2 / 0.15, 0.10 /            ! slight parallel shift from diagonal
*---------------------------------------------------------------
      R   = (x(1)+x(2) -1)**2            ! distance from diagonal
      Fun = Thickness/( R + Thickness**2)/Pi ! Breit-Wigner profile, normalized to delta
      TFun2Diagon = Fun
      END

      DOUBLE PRECISION FUNCTION TFun3Sphere(x)
*//////////////////////////////////////////////////////////////////////////////////
*//   3-dimensional testing function                                             //
*//   Thin sphere centered at (A1,A2,A3) with Radius and Thickness defined below //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  x(*)
      DOUBLE PRECISION  Pi
      PARAMETER( Pi=3.1415926535897932d0)
      DOUBLE PRECISION  Radius,Thickness,A1,A2,A3,R,Fun
      DATA   A1,A2,A3  / 0.25, 0.40, 0.50 / ! centre of sphere
      DATA   Radius    / 0.35  /            ! radius  of sphere
      DATA   Thickness / 0.020 /            ! thickness of sphere
*---------------------------------------------------------------
      R   = SQRT( (x(1)-A1)**2 +(x(2)-A2)**2 +(x(3)-A3)**2 )
      Fun = Thickness/( (R-Radius)**2 + Thickness**2)/Pi ! Breit-Wigner profile, normalized to delta
      Fun = Fun/(4d0*Pi*Radius**2)                       ! Normalize to one
      TFun3Sphere    = Fun
      END

      DOUBLE PRECISION FUNCTION TFun2Sphere(x)
*//////////////////////////////////////////////////////////////////////////////////
*//   3-dimensional testing function                                             //
*//   Thin sphere centered at (A1,A2,A3) with Radius and Thickness defined below //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  x(*)
      DOUBLE PRECISION  Pi
      PARAMETER( Pi=3.1415926535897932d0)
      DOUBLE PRECISION  Radius,Thickness,A1,A2,R,Fun
      DATA   A1,A2  / 0.25, 0.40 / ! centre of sphere
      DATA   Radius    / 0.35  /            ! radius  of sphere
      DATA   Thickness / 0.020 /            ! thickness of sphere
*---------------------------------------------------------------
      R   = SQRT( (x(1)-A1)**2 +(x(2)-A2)**2 )
      Fun = Thickness/( (R-Radius)**2 + Thickness**2)/Pi ! Breit-Wigner profile
      Fun = Fun/(4d0*Pi*Radius**2)                       ! Normalize to one
      TFun2Sphere    = Fun
      END


      DOUBLE PRECISION FUNCTION TFun1d1(x)
*//////////////////////////////////////////////////////////////////////////////////
*//                                                                              //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  x
      DOUBLE PRECISION  x0,gamma,pi
*---------------------------------------------------------------
      pi = 4d0*DATAN(1d0)
* breit-wigner
      x0=0.6d0
      gamma=0.001d0
      TFun1d1    = 1d0/pi*gamma/( (x-x0)**2 + gamma**2 )
      END

      DOUBLE PRECISION FUNCTION TFun1d(x)
*//////////////////////////////////////////////////////////////////////////////////
*//   sharp peak on the background of a mild function                            //
*//                                                                              //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  x
      DOUBLE PRECISION  x0,pi,gamma,lambda
*---------------------------------------------------------------
      pi = 4d0*DATAN(1d0)
* breit-wigner
      x0=0.6d0
      gamma=0.00001d0
      lambda = 3d0
      TFun1d    = 1d0/pi*gamma/( (x-x0)**2 + gamma**2 )
     $            +EXP(-lambda*x)*lambda
     $            +EXP(-lambda)
      END

      DOUBLE PRECISION FUNCTION Tfun1d2(x)
*//////////////////////////////////////////////////////////////////////////////////
*//                                                                              //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  x,eps
*---------------------------------------------------------------
* 1/x singularity
      eps =0.05d0
      Tfun1d2    = 1d0/(x+eps)
      END

      DOUBLE PRECISION FUNCTION Tfun1d3(x)
*//////////////////////////////////////////////////////////////////////////////////
*//                                                                              //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION  x,eps
*---------------------------------------------------------------
* primitive
      Tfun1d3    = 1d0-x
      END
