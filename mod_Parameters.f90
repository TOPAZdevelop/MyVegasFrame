MODULE ModParameters
implicit none


integer, public :: Collider, Process, PDFSet, NLOParam
integer, public :: MasterProcess, Correction, ObsSet
integer, public :: VegasIt0,VegasIt1,VegasNc0,VegasNc1,VegasMxDim,VegasSeed
integer, public :: VegasIt0_default,VegasIt1_default,VegasNc0_default,VegasNc1_default
integer, public :: NumHistograms=0
integer(8), public, save :: EvalCounter=0
integer(8), public, save :: SkipCounter=0

logical, public :: Seed_random,WarmUp
integer(8) :: ItMx,NCall,NPrn,It
integer(8), parameter :: MXDIM=25! has to match pvegas_mpi.h
integer, public :: TheSeeds(0:2) = (/2,700470849,476/)! only used if seed_random=.false., the first entry is the total number of seeds

! PVegas (MPI) part
integer, public :: MPI_Rank
integer, parameter :: NUMFUNCTIONS=1
integer, parameter :: WORKERS=7
integer, parameter :: PDIM=0

real(8), public :: MuRen, MuFac, AvgFactor
character, public :: HistoFile*(200),FileTag*(50),DataDir*(200)
character, public :: GridFile*(200),LHEFile*(200)
integer, public :: GridIO
real(8), public :: time_start,time_end

real(8), public, parameter :: Pi = 3.1415926535897932384626433832795028842d0
real(8), public, parameter :: GeV=0.01d0


real(8), public, parameter :: GF = (1.16639d-5)/GeV**2
real(8), public            :: m_Top
real(8), public, parameter :: m_Z     = 91.1876*GeV
real(8), public, parameter :: m_W     = 80.399d0*GeV
real(8), public, parameter :: m_H     = 125.0d0*GeV
real(8), public, parameter :: g2_weak = 4d0*dsqrt(2d0)*m_W**2*GF
real(8), public, parameter :: g_weak = dsqrt(g2_weak)
!orig real(8), public, parameter :: sw = dsqrt(4.d0*Pi*alpha/g2_weak)
!orig real(8), public, parameter :: sw2 = sw**2
real(8), public, parameter :: sw2 = 1d0 - m_W**2/m_Z**2
real(8), public, parameter :: sw = dsqrt(sw2)
real(8), public, parameter :: alpha4Pi = g2_weak*sw2
real(8), public, parameter :: alpha = alpha4Pi/4d0/Pi
real(8), public, parameter :: cw = dsqrt(1d0-sw2)
real(8), public, parameter :: EL = dsqrt(4.d0*Pi*alpha)
real(8), public            :: Ga_Top(0:1)
real(8), public            :: Ga_TopExp = 1.99d0*GeV
real(8), public            :: Ga_WExp   = 2.14d0*GeV
real(8), public            :: Ga_ZExp   = 2.4952d0*GeV
real(8), public, parameter :: Vev  = 246.218458102d0*GeV!  =1.0d0/sqrt(Gf*sqrt(2.0d0))

real(8), public, parameter :: Q_up    = 2d0/3d0
real(8), public, parameter :: Q_dn    =-1d0/3d0
real(8), public, parameter :: Q_el    =-1d0
real(8), public, parameter :: Q_nu    = 0d0
real(8), public            :: Q_top
real(8), public, parameter :: Q_Wp    =+1d0
real(8), public, parameter :: Q_Wm    =-1d0
real(8), public            :: Q_in

real(8), public, parameter :: T3_up    =+1d0/2d0
real(8), public, parameter :: T3_dn    =-1d0/2d0
real(8), public, parameter :: T3_nu    =+1d0/2d0
real(8), public, parameter :: T3_el    =-1d0/2d0

real(8), public, parameter :: couplZUU_left  = -sw/cw*Q_up + 1d0/sw/cw * T3_up
real(8), public, parameter :: couplZUU_right = -sw/cw*Q_up  
real(8), public, parameter :: couplZDD_left  = -sw/cw*Q_dn + 1d0/sw/cw * T3_dn
real(8), public, parameter :: couplZDD_right = -sw/cw*Q_dn

real(8), public, parameter :: couplZEE_left  = -sw/cw*Q_el + 1d0/sw/cw * T3_el
real(8), public, parameter :: couplZEE_right = -sw/cw*Q_el

real(8), public, parameter :: couplZNN_left  = -sw/cw*Q_nu + 1d0/sw/cw * T3_nu
real(8), public, parameter :: couplZNN_right = -sw/cw*Q_nu

real(8), public :: WidthExpansion
real(8), public, parameter :: fbGeV2=0.389379d12*GeV**2

real(8), public, parameter :: Nf_light=5d0
real(8), public, parameter :: Nf_heavy=1d0
real(8), public, parameter :: Nf=Nf_light+Nf_heavy

! Jet Algorithm
integer, public            :: AlgoType           ! 1 = kT,   -1 = anti kT
integer, public, parameter :: RecombPrescr = 0   ! 0 = 4-vector addition,   1 = Ellis-Soper prescription

! PDF Set
character, public :: PDFSetString*(100)
integer, public :: LHAPDFMember

integer, parameter, public :: LHE_Up_=2
integer, parameter, public :: LHE_Dn_=1
integer, parameter, public :: LHE_Chm_=4
integer, parameter, public :: LHE_Str_=3
integer, parameter, public :: LHE_Top_=6
integer, parameter, public :: LHE_Bot_=5
integer, parameter, public :: LHE_Glu_=21
integer, parameter, public :: LHE_ElM_=11
integer, parameter, public :: LHE_MuM_=13
integer, parameter, public :: LHE_TaM_=15
integer, parameter, public :: LHE_NuE_=12
integer, parameter, public :: LHE_NuM_=14
integer, parameter, public :: LHE_NuT_=16
integer, parameter, public :: LHE_Wp_=24
integer, parameter, public :: LHE_Pho_=22
integer, parameter, public :: LHE_Z_=23
integer, parameter, public :: LHE_Hig_=25


! particle 0 = not defined
integer, public, target :: Up_  = 1
integer, public, target :: Dn_  = 2
integer, public, target :: Chm_ = 3
integer, public, target :: Str_ = 4
integer, public, target :: Top_ = 5
integer, public, target :: Bot_ = 6
integer, public, target :: Glu_ = 10
integer, public, target :: Pho_ = 11
integer, public, target :: Z0_  = 12
integer, public, target :: Wp_  = 13
integer, public, target :: HTop_= 14
integer, public, target :: Stop_= 15
integer, public, target :: Sbot_= 16
integer, public, target :: ElM_ = 20
integer, public, target :: Hig_ = 26

integer, public, target :: AUp_  = -1
integer, public, target :: ADn_  = -2
integer, public, target :: AChm_ = -3
integer, public, target :: AStr_ = -4
integer, public, target :: ATop_ = -5
integer, public, target :: ABot_ = -6
integer, public, target :: Wm_   = -13
integer, public, target :: AHTop_= -14
integer, public, target :: AStop_= -15
integer, public, target :: ASBot_= -16
integer, public, target :: ElP_  = -20


real(8), public :: Collider_Energy
real(8), public :: alpha_s, alpha_s4Pi, alpha_sOver2Pi
 

 CONTAINS



SUBROUTINE InitParameters
implicit none

END SUBROUTINE





END MODULE



