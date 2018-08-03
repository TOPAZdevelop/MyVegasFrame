! ./NumInt VegasIt0=3 VegasNc0=100000

PROGRAM NumTwoLoop
use ModParameters
use ModKinematics
! use ifport
implicit none
real(8) :: MCresult,MCerror,chi2
!DEC$ IF(_UseMPIVegas.EQ.1)
! include 'mpif.h'
! integer ::ierror
!    call MPI_INIT(ierror)
!    call MPI_COMM_RANK(MPI_COMM_WORLD,MPI_Rank,ierror)
!DEC$ ELSE
   MPI_Rank=0
!DEC$ ENDIF


   call GetCommandlineArgs()
   call SetFileNames()
   call InitParameters()
   call InitHisto()
   call InitVegas()
   call InitRandomSeed(TheSeeds)   
   call OpenFiles()
   if( MPI_Rank.eq.0 ) then   
      call WriteParameters(6)   ! 6=stdout
   endif

   if( MPI_Rank.eq.0 ) print *, "Running"
   if( MPI_Rank.eq.0 ) call cpu_time(time_start)
   call StartPVegas()
   if( MPI_Rank.eq.0 ) then
      call cpu_time(time_end)
      write(*,"(A,1F8.1,A)") "Done in ",(time_end-time_start)," seconds."
   endif
   
   call CloseFiles()
!DEC$ IF(_UseMPIVegas .EQ.1)  
!    call MPI_FINALIZE(ierror)
!DEC$ ENDIF


END PROGRAM







SUBROUTINE GetCommandlineArgs()
use ModParameters
use ModKinematics
use ModMisc
implicit none
character :: arg*(100)
character :: env*(31)
integer :: NumArgs,NArg


   VegasNc0=-1
   VegasNc1=-1
   VegasIt0=-1
   VegasIt1=-1

#if _compiler==1
   NumArgs = NArgs()-1
#elif _compiler==2
   NumArgs = COMMAND_ARGUMENT_COUNT()
#endif   
      
   do NArg=1,NumArgs
    call GetArg(NArg,arg)
    if( arg(1:9).eq."VegasIt0=" ) then
        read(arg(10:11),*) VegasIt0
    elseif( arg(1:9).eq."VegasIt1=" ) then
        read(arg(10:11),*) VegasIt1
    elseif( arg(1:9).eq."VegasNc0=" ) then
        read(arg(10:20),*) VegasNc0
    elseif( arg(1:9).eq."VegasNc1=" ) then
        read(arg(10:20),*) VegasNc1
    endif
   enddo

   if( VegasNc0.le.0 ) call Error("VegasNc0 not set")
!    if( VegasNc1.le.0 ) call Error("VegasNc1 not set")
   if( VegasIt0.le.0 ) call Error("VegasIt0 not set")
!    if( VegasIt1.le.0 ) call Error("VegasIt1 not set")
    
return
END SUBROUTINE





SUBROUTINE SetFileNames()
use ModParameters
! use ifport
implicit none




RETURN
END SUBROUTINE




SUBROUTINE WriteParameters(TheUnit)
use ModParameters
use ModKinematics
implicit none
integer TheUnit



END SUBROUTINE







SUBROUTINE StartPVegas()
use ModIntegrand
use ModKinematics
use ModParameters
implicit none
real(8) :: MCresult,MCerror,MC_Chi2
!DEC$ IF(_UseMPIVegas .EQ.1)  
! include 'mpif.h'
!DEC$ ENDIF
integer i,init,NDim,NEVAL,IFAIL,kDim,nBuf,nSampl,OptDrive,OptEdge,loop,nregions
double precision yrange(1:2*MXDIM)
INTEGER, PARAMETER :: NW=10000000
REAL(8) :: WORK(1:NW),MCwt,MCvector(10)
integer ::  userdata, nvec, flags, mineval, maxeval, key,  spin, nbatch, gridno,nnew, nmin
double precision :: epsrel, epsabs, flatness



  ndim=2
  ncall = VegasNc0
  itmx  = VegasIt0
  nprn = 3

  do i=1,ndim
    yrange(i)=0d0
    yrange(i+ndim)=1d0
  enddo

  init=0
!DEC$ IF(_UseMPIVegas .EQ.1)  
  call ClearRedHisto()
!DEC$ ENDIF  

  yrange(1) = -20d0
  yrange(2) = -20d0
  yrange(3) = +20d0
  yrange(4) = +20d0




!DEC$ IF(_UseMPIVegas .EQ.1)  
!   call vegas_mpi(yrange(1:2*ndim),ndim,VegasIntegrand2,init,ncall,itmx,nprn,NUMFUNCTIONS,PDIM,WORKERS,MCresult,MCerror,MC_Chi2)
!DEC$ ELSE  
  call vegas_ser(yrange(1:2*ndim),ndim,VegasIntegrand2,init,ncall,itmx,nprn,NUMFUNCTIONS,PDIM,WORKERS,MCresult,MCerror,MC_Chi2)
!DEC$ ENDIF



print *, ""
call DCUHRE(ndim,NUMFUNCTIONS,yrange(1:2),yrange(3:4),ncall/2,ncall,VegasIntegrand2_DCUHRE,0d0,1d-6,0,NW,init,MCresult,MCerror,NEVAL,IFAIL,WORK)
print *, "cuhre ",MCresult,MCerror,NEVAL,IFAIL







userdata = 0 
nvec =1 
epsrel = 1d-6
epsabs = 0d0
flags = 0 
key =  9 
spin = -1

print *, ""

  call cuhre(ndim, NUMFUNCTIONS, VegasIntegrand2_CUBA, userdata, nvec, epsrel, epsabs, flags, ncall/2,ncall ,key, "", spin,nregions, NEVAL, IFAIL, MCresult,MCerror, MC_Chi2)

print *, "cuba cuhre",MCresult,MCerror,NEVAL,IFAIL





print *, ""
nbatch = 1000
gridno=1
  call vegas(ndim, NUMFUNCTIONS, VegasIntegrand2_CUBA, userdata, nvec,epsrel, epsabs, flags, 1313, ncall/2,ncall, ncall/10, ncall/10, nbatch, gridno, "", spin, NEVAL, IFAIL, MCresult,MCerror, MC_Chi2)

print *, "cuba vegas",MCresult,MCerror,NEVAL,IFAIL





print *, ""
nbatch = 1000
nnew = 1000
nmin = 2
flatness = 25d0
  call suave(ndim, NUMFUNCTIONS, VegasIntegrand2_CUBA, userdata, nvec,epsrel, epsabs, flags, 1313, ncall/2,ncall, nnew, nmin, flatness, "", spin, nregions,NEVAL, IFAIL, MCresult,MCerror, MC_Chi2)

print *, "cuba suave",MCresult,MCerror,NEVAL,IFAIL


print *, ""
! *=========================================================
      nDim      = 2
      kDim      = 0
      nBuf      = 1000          ! Buffer length <5000, default = 1000
      nSampl    = 10000           ! Number of MC events in exploration of a cell, default = 200
      OptDrive  = 2             ! Type of Drive =0,1,2 for TrueVol,Sigma,WtMax, default = 2
      OptEdge   = 1             ! Include vertices in sampling or not, =0,1, default =1
! *=========================================================
      CALL GLK_Book1(8100,'WT distribution  $', 80, 0.0d0,2d0)
! *========================================================= 
      CALL FoamA_SetnDim(          nDim)
      CALL FoamA_SetkDim(          kDim)
      CALL FoamA_SetnBuf(          nBuf)
      CALL FoamA_SetnSampl(      nSampl)
      CALL FoamA_SetOptDrive(  OptDrive)
      CALL FoamA_SetOptEdge(    OptEdge)
      CALL FoamA_SetChat(0)
! *------------------------------------
      CALL FoamA_Initialize(VegasIntegrand2_FOAM)
! *------------------------------------
      CALL WtLimitStart(1004,1005,1000,10d0)  
           
      DO loop = 1, ncall
         CALL FoamA_MakeEvent(VegasIntegrand2_FOAM)    ! generate MC event
         CALL FoamA_GetMCvector(MCvector)  ! get MC event, vector
         CALL FoamA_GetMCwt(MCwt)          ! get MC weight
         CALL WtLimitFill(MCwt)
         CALL GLK_Fil1(8100,MCwt,1d0)
      ENDDO
      CALL FoamA_Finalize(MCresult,MCerror)
!       WRITE(*,'(a,g19.9,a,g19.9,a,f7.5)') 'MCresult= ',MCresult,' +- ',MCerror,'         RelErr= ',MCerror/MCresult   !


return
END SUBROUTINE











SUBROUTINE InitVegas()
use ModKinematics
use ModParameters
implicit none




return
END SUBROUTINE








SUBROUTINE OpenFiles()
use ModParameters
implicit none
character :: filename*(200)

!    filename = trim(HistoFile)//'.dat'
!    open(unit=14,file=trim(filename),form='formatted',access= 'sequential',status='replace')            ! Histogram file

!    filename = trim(HistoFile)//'.status'
!    open(unit=15,file=trim(filename),form='formatted',access= 'sequential',status='replace')         ! Vegas status file

!    filename = trim(HistoFile)//'.tmp_histo'
!    open(unit=16,file=trim(filename),form='formatted',access= 'sequential',status='replace')         ! Histo status file


return
END SUBROUTINE



SUBROUTINE CloseFiles()
use modParameters
implicit none

!    close(14)
!    close(15)
!    close(16)


return
END SUBROUTINE




SUBROUTINE WriteHisto(TheUnit,curit,VG_CurrResult,VG_CurrError,MCresult,MCerror,Chi2,RunTime)
use ModKinematics
use ModParameters
implicit none
integer :: NBin,Hits,NHisto,SumHits,TheUnit,curit,NBin2,NBin3,NHisto2
real(8) :: BinSize,LowVal,BinVal,Value,Error,Integral
real(8) :: BinSize2,BinSize3,LowVal2,LowVal3,BinVal2,BinVal3
real(8),parameter :: ToGeV=1d2, ToPb=1d-3
real(8) :: MCresult,MCerror,RunTime,VG_CurrResult,VG_CurrError,Chi2
character :: filename*(200),arg*(500)
integer, save :: Prev_Process=-1313999


  if(TheUnit.ne.6) then 
    filename = trim(HistoFile)//'.dat'
    open(unit=TheUnit,  file=trim(filename),form='formatted',access= 'sequential',status='replace')   ! Histogram file

    
!   writing status file
    filename = trim(HistoFile)//'.status'
    if( Prev_Process.ne.Process ) then
        open(unit=TheUnit+1,file=trim(filename),form='formatted',access= 'sequential',status='replace')   ! status file
        call Get_Command(arg)
        write(TheUnit+1,'(A1,1X,A)') "#","-------------------------------------------------------------------------------------------------"
        write(TheUnit+1,'(A1,1X,A)') "#",trim(arg)
        write(TheUnit+1,'(A1,1X,A)') "#","-------------------------------------------------------------------------------------------------"
        write(TheUnit+1,'(A1,1X,A)') "#","  It.    Result              Error               Accum. result       Accum. error         Chi^2  "
        write(TheUnit+1,'(A1,1X,A)') "#","-------------------------------------------------------------------------------------------------"
        Prev_Process = Process
      else
        open(unit=TheUnit+1,file=trim(filename),form='formatted',access= 'sequential',status='old',position='append')   ! status file
        if( curit.eq.1 ) write(TheUnit+1,'(A1,1X,A)') "#","-------------------------------------------------------------------------------------------------"
    endif
    if( curit.gt.0 ) write(TheUnit+1,'(A1,1X,I3,4E20.8,F12.2)') "#",curit,VG_CurrResult,VG_CurrError,MCresult,MCerror,Chi2
  endif

  call WriteParameters(TheUnit)
  write(TheUnit,"(A,2X,1F9.2,A)") "# run time =",RunTime/60d0,"min"
  write(TheUnit,"(A,2X,1F20.10)") "# EvalCounter  =",dble(EvalCounter)
  write(TheUnit,"(A,2X,1F20.10)") "# SkipCounter  =",dble(SkipCounter)
  if( EvalCounter.gt.0 .and. dble(SkipCounter)/dble(EvalCounter) .gt. 0.02d0 ) write(TheUnit,"(A,2X)") "# **** WARNING  ****: SkipCounter is larger than 2%"
  write(TheUnit,"(A,2X,1PE20.10,2X,1PE20.5)") "#TotCS[fb]=",MCresult,MCerror
  do NHisto=1,NumHistograms
      write(TheUnit,"(A,I2,A,A)") "# Histogram",NHisto,": ",Histo(NHisto)%Info
      Integral = 0d0
      SumHits = 0
      BinSize = Histo(NHisto)%BinSize * Histo(NHisto)%SetScale
      LowVal  = Histo(NHisto)%LowVal  * Histo(NHisto)%SetScale
      do NBin=1, Histo(NHisto)%NBins
          BinVal = (LowVal+(NBin-1)*BinSize)
          Hits   = Histo(NHisto)%Hits(NBin)
          SumHits = SumHits + Hits
          
              Value  = Histo(NHisto)%Value(NBin)/BinSize/curit
              Integral = Integral + Histo(NHisto)%Value(NBin)/curit
              Error  = 1d0/(BinSize)/curit * dsqrt(dabs( Histo(NHisto)%Value2(NBin) - 1d0/curit/ncall*Histo(NHisto)%Value(NBin)**2) )
          if(Hits.ge.999999999) Hits=999999999
          write(TheUnit,"(I2,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,2X,I9,A)") NHisto,"|",BinVal,"|",Value,"|",Error,"|",Hits,"|"

       enddo
      Integral = Integral + (Histo(NHisto)%Value(0)+Histo(NHisto)%Value(Histo(NHisto)%NBins+1))/curit
      write(TheUnit,"(A,2X,1PE23.16)") "# integrated result:",Integral
      write(TheUnit,"(A,2X,1I23)") "# total number of hits:",SumHits
  enddo


  if(TheUnit.ne.6) then
    close(TheUnit)
    close(TheUnit+1)
  endif

RETURN
END SUBROUTINE







! if seed_random=true:  "Seeds" returns the seeds chosen according to system clock
! is seed_random=false: "Seeds" is used as input for initializing the random number generator
SUBROUTINE InitRandomSeed(Seeds)
use modParameters
use modMisc
implicit none
integer :: Seeds(0:12)
integer, dimension(:), allocatable :: gen_seed
integer :: n,i,sclock,SeedSize


    if( seed_random ) then 
        call random_seed()
        call random_seed(size=SeedSize)
        if(SeedSize.ne.2) call Error("ifort SeedSize has changed from 2 to ",SeedSize)
        Seeds(0) = SeedSize
        call random_seed(get=Seeds(1:SeedSize))
    else        
        call random_seed(size=n)
        if( n.ne.Seeds(0) ) call Error("Number of input seeds does not match random_seed(size=n)",n)
        call random_seed(put = Seeds(1:n))
    endif

return
END SUBROUTINE

