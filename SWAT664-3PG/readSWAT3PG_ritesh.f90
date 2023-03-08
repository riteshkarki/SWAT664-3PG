        subroutine readSWAT3PGinp

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine reads co2 concentrations
!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    
!!   
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
    !! qichun co2-
      use parm 
  
      implicit none
      integer :: icc, eof, icc_scup
      real :: SLA0C,SLA1C,tSLAC,pRxC,pRnC,CoeffCondC,rAgeC,nAgeC,pFS2C
      real :: pFS20C,aWsC,nWsC,gammaF0C,gammaF1C,tgammaFC,gammaRC,gammaN0C
      real :: gammaN1C,tgammaNC,ngammaNC,mRC,mSC,mFC,alphaCxC,molPAR_MJC,gDM_molC
      real :: thin_powerC,wSx1000C,mort_thinning_keyC,mort_age_keyC,m0C,forest_typeC
      real :: aHC, nHBC
      character (len=5) :: cnamec
      !initializing locals in the loop
      
      eof = 0
      
      do
          icc = 0
          cnamec = ""
          SLA0C = 0.
          SLA1C = 0.
          tSLAC = 0.
          pRxC = 0.
          pRnC = 0.
          CoeffCondC = 0.
          rAgeC = 0.
          nAgeC = 0.
          pFS2C = 0.
          pFS20C = 0.
          aWsC = 0.
          nWsC = 0.
          gammaF0C = 0.
          gammaF1C = 0.
          tgammaFC = 0.
          gammaRC = 0.
          gammaN0C = 0.
          gammaN1C = 0.
          tgammaNC = 0.
          ngammaNC = 0.
          mRC = 0.
          mSC = 0.
          mFC = 0.
          alphaCxC = 0.
          molPAR_MJC = 0.
          gDM_molC = 0.
          thin_powerC = 0.
          wSx1000C = 0.
          aHC = 0.
          nHBC = 0.
          mort_thinning_keyC = 0.
          mort_age_keyC = 0.
          m0C = 0.
          forest_typeC = 0.
          
          
        read (2088,*,iostat=eof) icc, cnamec
        if (eof < 0) exit
        read (2088,*,iostat=eof) SLA0C,SLA1C,tSLAC,pRxC,pRnC,CoeffCondC,rAgeC,nAgeC,pFS2C,pFS20C    

        if (eof < 0) exit
        read (2088,*,iostat=eof) aWsC,nWsC,gammaF0C,gammaF1C,tgammaFC,gammaRC,gammaN0C,gammaN1C,tgammaNC,ngammaNC

        if (eof < 0) exit
        read (2088,*,iostat=eof) mRC,mSC,mFC,alphaCxC,molPAR_MJC,gDM_molC,thin_powerC,wSx1000C,aHC,nHBC     

        if (eof < 0) exit
        read (2088,*,iostat=eof) mort_thinning_keyC,mort_age_keyC,m0C,forest_typeC


        if (eof < 0) exit

        if (icc <= 0) exit
      
      
          !Ritesh - SWAT3PG Readin Variables
          SLA0(icc) = SLA0C
          SLA1(icc) = SLA1C
          tSLA(icc) = tSLAC
          pRx(icc) = pRxC
          pRn(icc) = pRnC
          CoeffCond(icc) = CoeffCondC
          rAge(icc) = rAgeC
          nAge(icc) = nAgeC
          pFS2(icc) = pFS2C
          pFS20(icc) = pFS20C
          aWs(icc) = aWsC
          nWs(icc) = nWsC
          gammaF0(icc) = gammaF0C
          gammaF1(icc) = gammaF1C
          tgammaF(icc) = tgammaFC
          gammaR(icc) = gammaRC
          gammaN0(icc) = gammaN0C
          gammaN1(icc) = gammaN1C
          tgammaN(icc) = tgammaNC
          ngammaN(icc) = ngammaNC
          mR(icc) = mRC
          mS(icc) = mSC
          mF(icc) = mFC
          alphaCx(icc) = alphaCxC
          molPAR_MJ(icc) = molPAR_MJC
          gDM_mol(icc) = gDM_molC
          thin_power(icc) = thin_powerC
          wSx1000(icc) = wSx1000C
          aH(icc) = aHC
          nHB(icc) = nHBC
          mort_thinning_key(icc) = mort_thinning_keyC
          mort_age_key(icc) = mort_age_keyC
          m0(icc) = m0C
          forest_type(icc) = forest_typeC
      
      

      end do 
      
      !temporary override of frsd parameters using bacteria
      !frsd (icc) = 7
          !icc_scup = 7
          
         ! sla0(icc_scup) = SLA0_scup !sla0c !wdpq
         ! sla1(icc_scup) = SLA1_scup !sla1c !wgpq 
         ! prx(icc_scup) = pRx_scup !prxc   !wdlpq
         ! prn(icc_scup) = pRn_scup !prnc   !wglpq
         ! rage(icc_scup) = rAge_scup !ragec !wdps
         ! nage(icc_scup) = nAge_scup !nagec !wgps
         ! pfs2(icc_scup) = pFS2_scup !pfs2c !wdlps
         ! pfs20(icc_scup) = pFS20_scup !pfs20c !wglps
         ! gammaf0(icc_scup) = gammaF0_scup !gammaf0c !bactkdq
         ! gammaf1(icc_scup) = gammaF1_scup  !gammaf1c !thbact
         ! tgammaf(icc_scup) = tgammaF_scup  !tgammafc !wof_p
         ! gammar(icc_scup) = gammaR_scup  !gammarc   !wof_lp
         ! gamman0(icc_scup) = gammaN0_scup !gamman0c !wdpf
        !  gamman1(icc_scup) = gammaN1_scup !gamman1c !wgpf
         ! tgamman(icc_scup) = tgammaN_scup !tgammanc !wdlpf                                            
        !  alphacx(icc_scup) = alphaCx_scup !alphacxc !wglpf
          


      close (2088)
      return
      end