MODULE ModKinematics
implicit none
save


type :: Histogram
    integer :: NBins
    real(8) :: BinSize
    real(8) :: LowVal
    real(8) :: SetScale
    real(8),allocatable :: Value(:)
    real(8),allocatable :: Value2(:)
    integer,allocatable :: Hits(:)
    character :: Info*(50)
!     
    logical :: BinSmearing=.false.
    real(8) :: SmearSigma=0.1d0   
end type





integer,public :: it_sav

type(Histogram),allocatable   :: Histo(:)





!DEC$ IF(_UseMPIVegas.EQ.1)
integer,public,parameter :: NUMHISTO=40       ! this has to match the constants in pvegas_mpi.c
integer,public,parameter :: MXHISTOBINS=500
type, BIND(C) :: ReducedHistogram
    real(8) :: Value(1:MXHISTOBINS)
    real(8) :: Value2(1:MXHISTOBINS)
    integer :: Hits(1:MXHISTOBINS)
end type
type(ReducedHistogram)  :: RedHisto(1:NUMHISTO)
public :: RedHisto,getRedHisto,transferHisto,clearRedHisto
!DEC$ ENDIF



contains





!DEC$ IF(_UseMPIVegas.EQ.1)
INTEGER FUNCTION getRedHisto(TheHisto,NHisto)
implicit none
type(ReducedHistogram) :: TheHisto
integer NHisto,NBin


  do NBin=1,MXHISTOBINS
    TheHisto%Value(NBin)  = RedHisto(NHisto)%Value(NBin)
    TheHisto%Value2(NBin) = RedHisto(NHisto)%Value2(NBin)
    TheHisto%Hits(NBin)   = RedHisto(NHisto)%Hits(NBin)
  enddo

getRedHisto=0
RETURN
END FUNCTION



INTEGER FUNCTION transferHisto(TheHisto,NHisto)
use ModParameters
implicit none
type(ReducedHistogram) :: TheHisto
integer NHisto,NBin
transferHisto=0

  if( NHisto.gt.NumHistograms ) return! this is required because in pvegas we loop until max.number of histograms (NUMHISTO)
  do NBin=1,Histo(NHisto)%NBins
    Histo(NHisto)%Value(NBin)  = TheHisto%Value(NBin)
    Histo(NHisto)%Value2(NBin) = TheHisto%Value2(NBin)
    Histo(NHisto)%Hits(NBin)   = TheHisto%Hits(NBin)
  enddo

RETURN
END FUNCTION



SUBROUTINE clearRedHisto()
implicit none
integer NHisto,NBin

  do NHisto=1,NUMHISTO
  do NBin=1,MXHISTOBINS
    RedHisto(NHisto)%Value(NBin)  = 0d0
    RedHisto(NHisto)%Value2(NBin) = 0d0
    RedHisto(NHisto)%Hits(NBin)   = 0
  enddo
  enddo

RETURN
END SUBROUTINE
!DEC$ ENDIF






SUBROUTINE InitHisto()
use ModMisc
use ModParameters
implicit none
integer :: AllocStatus,NHisto,i,j

it_sav = 1



RETURN
END SUBROUTINE








SUBROUTINE EvalPhasespace_2to2(EHat,xRndPS,Mom,PSWgt)
use ModMisc
use ModParameters
implicit none
real(8) :: EHat
real(8) :: PSWgt,PSWgt2,PSWgt3,PSWgt4,PSWgt5
real(8) :: Mom(1:4,1:4),MomW(1:4),xRndPS(1:2)
real(8),parameter :: N2=2, PiWgt2 = (2d0*Pi)**(4-N2*3) * (4d0*Pi)**(N2-1)
real(8),parameter :: N3=3, PiWgt3 = (2d0*Pi)**(4-N3*3) * (4d0*Pi)**(N3-1)


!  generate PS: massless + massless --> massive(anti-top) + massive(top)
   call genps(2,Ehat,xRndPS(1:2),(/m_Top,m_Top/),Mom(1:4,3:4),PSWgt)
   PSWgt = PSWgt*PiWgt2

!  particles on the beam axis:
   Mom(1,1) =  EHat*0.5d0
   Mom(2,1) =  0d0
   Mom(3,1) =  0d0
   Mom(4,1) = +EHat*0.5d0

   Mom(1,2) =  EHat*0.5d0
   Mom(2,2) =  0d0
   Mom(3,2) =  0d0
   Mom(4,2) = -EHat*0.5d0


return
END SUBROUTINE









END MODULE
