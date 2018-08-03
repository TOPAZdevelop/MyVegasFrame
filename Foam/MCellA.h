*///////////////////////////////////////////////////////////////////////////////////////
*//                                                                                   //
*//          Pseudoclass MCellA                                                       //
*//          MCell  Mega Cell MC saampler                                             //
*//                                                                                   //
*//   To get Mega cells (re)define  m_nBufMax  below                                  //
*//                                                                                   //
*///////////////////////////////////////////////////////////////////////////////////////
*
      INTEGER            m_nBinMax
      PARAMETER          (m_nBinMax=256)
*
      INTEGER        m_nBufMax
      PARAMETER  (   m_nBufMax = 100000 )         ! maximum buffer length for all cells (100k=8MB)
      INTEGER        m_KdiMax
      PARAMETER  (   m_KdiMax  = 20 )             ! maximum dimension for hypercubics
      INTEGER        m_DimMax
      PARAMETER  (   m_DimMax  = m_KdiMax )       !  maximum total dimension
      INTEGER        m_cMax
      PARAMETER  (   m_cMax    = m_nBufMax/2+1 )  ! maximum number of (active) Cells
*
      INTEGER            m_CeStat,   m_CePare,   m_CeDau1,  m_CeDau2,  m_CeSamp, m_CeVert !
      DOUBLE PRECISION   m_CeIntg,   m_CeDriv,  m_CePrim,   m_CePrCu !
      INTEGER            m_CeDivi
      INTEGER            m_LastCe,   m_LastAc,   m_nBuf,    m_nBin,    m_ActC !
      INTEGER            m_Kdim,     m_Dimen
      INTEGER            m_Chat,     m_Out,      m_nSampl,  m_Ncalls !
      INTEGER            m_OptPeek,  m_OptDrive, m_OptEdge, m_EvPerBin, m_OptRanIni, m_OptRanLux!
      DOUBLE PRECISION   m_VolTot   !
      DOUBLE PRECISION   m_Drive,    m_SumWt,    m_SumWt2,  m_NevGen,    m_WtMax,   m_WtMin !
      DOUBLE PRECISION   m_MCresult, m_MCerror,  m_MCwt,    m_MCvector !
      INTEGER            m_Ltx,      m_MagicInit
*
      COMMON /c_MCellA/   
     $ m_CeStat(m_nBufMax),           ! Cell member: (4B) status=0 inactive, =1 active
     $ m_CePare(m_nBufMax),           ! Cell member: (4B) parent cell pointer
     $ m_CeDau1(m_nBufMax),           ! Cell member: (4B) daughter1 cell pointer
     $ m_CeDau2(m_nBufMax),           ! Cell member: (4B) daughter2 cell pointer
     $ m_CeSamp(m_nBufMax),           ! Cell member: (4B) No of MC events in exploration
     $ m_CeDivi(m_nBufMax),           ! Cell member: (4B) Division ratio and edge
     $ m_CeIntg(m_nBufMax),           ! Cell member: (8B) integral estimator
     $ m_CeDriv(m_nBufMax),           ! Cell member: (8B) Drive integral estimate, for build-up
     $ m_CePrim(m_nBufMax),           ! Cell member: (8B) Primary integral estimate, MC generation
     $ m_CePrCu(0:m_cMax),            ! Active Cell member: (8B) Cumulative Primary 
     $ m_ActC(    m_cMax),            ! Active Cell member: (8B) List of all pointers to ACTIVE cells
* Actually 56B/Cell, 100k Cells = 5.6MB, factor 2 compression still possible.
     $ m_VolTot,            ! Estimate of Volume total, without error
     $ m_Drive,             ! M.C. generation Drive value of integral
     $ m_SumWt,             ! M.C. generation sum of Wt
     $ m_SumWt2,            ! M.C. generation sum of Wt**2
     $ m_NevGen,            ! M.C. generation sum of 1d0
     $ m_WtMax,             ! M.C. generation maximum wt
     $ m_WtMin,             ! M.C. generation minimum wt
     $ m_MCresult,          ! M.C. generation Final value of INTEGRAL
     $ m_MCerror,           ! M.C. generation Final walue of ERROR
     $ m_MCwt,              ! M.C. generation current event weight
     $ m_MCvector(m_DimMax),! M.C. generated vector,
     $ m_Kdim,              ! dimension of the hypercubics
     $ m_Dimen,             ! total dimension of the problem =m_Kdim
     $ m_nBuf,              ! Actual dynamic lenth of the buffer m_nBuf<m_nBufMax
     $ m_nBin,              ! Actual dynamic lenth of the buffer m_nBuf<m_nBufMax
     $ m_LastAc,            ! Last active cell
     $ m_LastCe,            ! Last cell in buffer 
     $ m_nSampl,            ! No. of sampling when dividing cell
     $ m_Ncalls,            ! No. of function calls, total
     $ m_OptPeek,           ! Flag for  random ceel choice: Peek =0,1 for maximum,
     $ m_OptDrive,          ! Flag for type of Drive =0,1,2 for TrueVol,Sigma,WtMax
     $ m_OptEdge,           ! Flag which decides whether vertices are included in the sampling
     $ m_EvPerBin,          ! Max eff. MC events/bin, (saving CPU) inactive if =0
     $ m_OptRanIni,         ! Flag default, =1 r.n. generator not initialized in MCellA
     $ m_OptRanLux,         ! Flag =-1,0,1,2,3,4 r.n. generator level
     $ m_Chat,              ! Flag for chat level in output, Chat=1 normal level
     $ m_Ltx,               ! Latex Output unit, for debug
     $ m_Out,               ! Output unit
     $ m_MagicInit          ! Magic cookie of initialization (global variable)
*
      SAVE /c_MCellA/
*//////////////////////////////////////////////////////////////////////////////////////
*//                                                                                  //
*//////////////////////////////////////////////////////////////////////////////////////
