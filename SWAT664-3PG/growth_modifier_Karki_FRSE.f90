        subroutine growth_modifier_3PG_FRSE
        !!Ritesh - Module to simulate Evergreen forest using 3PG
        !!============================================
        !!Input variables
        !!    sol_sw(:)     |mm H20        |amount of water stored in the soil profile on
        !!                                 |current day
        !!    sol_fc(:,:)   |mm H2O        |amount of water available to plants in soil 
        !!                                 |layer at field capacity (fc - wp),Index:(layer,HRU)
        !!    sol_nly(:)    |              |number of layers in the soil profile
        !!    sol_clay(:,:) |none          |fraction clay content in soil
        !!                                 |point)
        !!            idplt |none          |land cover code from crop.dat
        !!    bioday        |kg              |biomass accumulated in the day after adjusting for stress
        !!==============================================
        !!==============================================
        !!read in parameters - defined locally at the moment
       !pRx     : Maximum fraction of NPP to roots
       !pRn     : Minimum fraction of NPP to roots
       !pFS2    : Ratio of foliage:stem partitioning at B = 2 cm
       !pFS20   : Ratio of foliage:stem partitioning at B = 20 cm
       !aWs     : Constant in the stem mass vs diameter relationship
       !nWs     : Power in the stem mass vs diameter relationshiop
       !gammaF0 : litterfall rate at t = 0 month^-1
       !gammaF1 : Maximum litterfall rate month^-1
       !tgammaF : Age at which litterfall rate has mediam value 
       !gammaF  : monthly litterfall rate
       !kg       : parameter calcualtion for calculating litterfall rate
       !gammaR  : monthly root turnover
       !gammaN0 : seedling mortality rate
       !gammaN1 : mortality rate for large t
       !tgammaN : Age at which mortality rate has median value
       !ngammaN : shape of mortality response
       !gammaN  : stem turnover rate
       !wSx1000 : max stem mass per tree @ 1000 tree per hectare
       !thinPower : Power in self-thinning rule
       !mF      : fraction mean single-tree foliage biomass lost per dead tree
       !mR      : fraction mean single-tree root biomass lost per dead tree
       !mS      : fraction mean single-tree stem biomass lost per dead tree\
       !mort_stress : Number of trees mortality due to stress
       !SLA0     : Specific Leaf Area for the calcuation of LAI at age 0 - Unit - m2/kg
       !SLA1    : Specific Leaf Area for the calculation of LAI for mature leaves - m2/kg
       !tSLA    : Age at which specific leaf area = (SLA0 + SLA1)/2 - years
       !aH      : Constant in stem-height relationship
       !nHB     : Power of DBH in stem height relationship
       !nHC     : Power of competition in the stem heigt relationship

        !!==============================================
        !! local variables
       !SWconst :     Soil water constant - relative water deficiti for 50% reduction (3PG)
       !SWpower :     Power determining shape of soil wate response
       !pFsConst:     For calculating the allometric relationship for above ground biomass partitioning
       !pFsPower:     For calculating the allometric relationship for above ground biomass partitioning
       !Sol_AWC_Max:  Maximum available soil water
       !soil_class:   soil class (1,2,3,4 for sandy, sandy loam, clayey loam, clay)
       !clay_percent: average clay percent in all soil layers
       !sol_modifier: soil water modifier
       !vpd_modifier: vapor pressure deficit modifier
       !age_modifier: tree age modifier
       !n_frac_root : root fraction for the day based on modifier
       !n_frac_foliage: 
       !f_phys      : combined modifier for calculating root biomass fraction for the day
       !FR          : fertility rating (value between 0 and 1)
       !m0          : m0 is a calibration paramter
       !m           : paramter for calculating n_frac_root
       !bio_loss    : cumulative biomass lost on the day from foliage, stem, and roots.
    
    
        use parm
        
        implicit none
        
        integer :: j, k, nlay
        real::sf=0.
       
        real :: SWconst, SWpower, Sol_AWC_Max, soil_class
        real :: clay_percent!, sol_modifier, vpd_modifier !!Making sol_modifier and vpd_modifier to global variable for analysis
        real :: age_modifier, idp !f_phys - converted to a global variable
        real :: FR
        real :: pFsPower, pFsConst
        real :: pFS
        real :: gammaF, kg
        real :: gammaN, mort_stress
        !real :: SLA
        
        !Calculation of biomass
        real :: alphaC, epsilon
        
        !Calculation of basal area
        real :: basal_area_prop !basal area, basal area proportion
        real :: stems_n_ha    !no of stems per hectare, power for self thinning rule
        real :: biom_tree_max, PI    !max stem mass per tree @ 1000 tree per ha, maximum biomass for each tree at current stocking
        real :: mort_n_thinning
        real, parameter :: accuracy = 1.d0/1000.d0
        integer :: t
        real :: n, x1, x2, fN, dfN, dN
        
        !modifying NPP:GPP ratio with age
        real :: npp_ratio
        
        !!Declaring variables for Carbon-Zhang
        real :: BLG1, BLG2, BLG3, CLG, XX
        real :: sol_min_n,  resnew_n, resnew_ne, resnew
        real :: LMF, LSF, LSLF, LSNF,LMNF
        real :: RLN, RLR 
        
       
        !initializing local variables 
        j=0.
        j=ihru
        nlay = 0.
        nlay = sol_nly(j)
        idp = 0.
        idp = idplt(j)
        bioday3PG_opt = 0.
        
        k = 0.
        SWconst = 0.
        SWpower = 0.
        Sol_AWC_Max = 0.
        soil_class = 0.
        clay_percent = 0.
        sol_modifier = 0.
        vpd_modifier = 0.
        n_frac_root = 0.
        n_frac_foliage = 0.
        n_frac_stem = 0.
        f_phys = 0. !!Made into a global variable
        SLA = 0.


        FR = 0.0

        
        !Initializing for above ground biomass partitioning 

        pFsConst = 0.
        pFsPower = 0.
        pFS = 0. 

        
        !Initialization for calculating foliage mortality (litterfall)
        gammaF = 0.
        kg = 0.
        

        
        !Initialization for calculating stem mortality based on age
        gammaN = 0.

        
        !Assigning value
        mort_stress = 0.
        mort_n_thinning = 0.
        
        !!Added by Ritesh - calculation of biomass assimilation using equations from the 3PG module
        alphaC = 0.
        epsilon = 0.

        
        !!Assigning values for calculation of mortality based on self-thinning
        basal_area_prop = 0.
        stems_n_ha = 0.
        biom_tree_max = 0.
        PI = 0.

        
        
        !Calculation of fphys -------------------------------------------------------------------------------------------
        
        !Calculation of soil modifier
        
        !maximum available soil water
        do k=1,nlay
            Sol_AWC_Max = Sol_AWC_Max + sol_fc(k,j)
        End do
        
        !Average soil  clay percent
        do k=1,nlay
            clay_percent = clay_percent + sol_clay(k,j)
        End do
        clay_percent = clay_percent/nlay/100
        
        !Assigning soil class based on clay percent 
        if (clay_percent < 0.11) then
            soil_class = 1
        else if (clay_percent > 0.10 .and. clay_percent < 0.31) then
            soil_class = 2
        else if (clay_percent > 0.30 .and. clay_percent < 0.51) then
            soil_class = 3
        else
            soil_class = 4
        end if
        
        SWconst = 0.8d0 - 0.10d0 * soil_class
        SWpower = 11.d0 - 2.d0 * soil_class
        !SWconst = 0.7
        !SWpower = 9.0
        
        !soil modifier
        sol_modifier = 1.d0 / (1.d0 + ((1.d0 -  sol_sw(j) / Sol_AWC_Max) / SWconst) ** SWpower)
        
        !Calculation of VPD modifier
        !!vpd_modifier = EXP(-2.5*VPD*10)
        vpd_modifier = EXP(-CoeffCond(idp)*vpd)
         
        !Calculation of tree age for age modifier
        !Tree age
        tree_age(j) = float(curyr_mat(j)) + float(i_mo-1)/12 !!tree age in years and month
        
        !Tree age related modifier
        if (tree_age(j) < 1) then
            age_modifier = 1
        else
            age_modifier = 1/(1+(tree_age(j)/float(mat_yrs(idp))/rAge(idp))**nAge(idp))
        end if
        
        !Calculation of total modifier for root fraction calculation
        f_phys = Min(vpd_modifier, sol_modifier) * age_modifier        !f_phys = 1
        !Calculation of FR
        !FR = 0.1
        FR = min(strsw(j), strsn(j))
        !Calculation of biomass assimilation using equations from the 3PG module -------------------------------------------------------------------------
        !alphaC = alphaCx*strstmp(j)*min(strsw(j), strsn(j))*f_phys
        !epsilon = gDM_mol*molPAR_MJ*alphaC
        !bioday3PG(j) = 0.47*epsilon*par*10 !!kg/ha 

        alphaC = alphaCx(idp)
        epsilon = gDM_mol(idp)*molPAR_MJ(idp)*alphaC 
        
        
        !modifying bioday3PG_opt for changing NPP:GPP ratio with age
        !bioday3PG_opt = 0.47*epsilon*par*10 !!kg/ha
        npp_ratio = -(9*10.**(-7)*(tree_age(j)**2)) - 0.0006*tree_age(j) + 0.4804 
        
        bioday3PG_opt = npp_ratio*epsilon*par*10 !!kg/ha
        
        !print *, tree_age(j), npp_ratio
        !Calculation of biomass fraction from previous day for use in nupF, npupF as well as CLG
        
        
       !Ritesh 
        bm_frac = 0.
        bm_frac(j) = bio_ms_3PG(j)/1000/600
       !Ritesh
        
        call nupF
        call npupF
        
        bioday3PG(j) = bioday3PG_opt*strstmp(j)*min(strsw(j), strsn(j), strsp(j))*f_phys
        
        !print *, plantn(j), npp_ratio, bioday3PG(j)
        if (bioday3PG(j) <= 0. ) then
            bioday3PG(j) = 0.
        end if
        
        !!add by zhang
        !!============
        if (cswat == 2) then
        !! NPPC_d(j) = NPPC_d(j) + bioday3PG(j) * 0.42
         !modified by Ritesh to print bioday for each day for forests 
         NPPC_d(j) = bioday3PG(j) * 0.42
        end if
        !!add by zhang
          !!============          
        !Caculation of root, stem, and leaves fraction from the daily assimilated biomass ------------------------------------------------------------------
        
        !Calculation of root fraction (n_frac_root)
        m_root = m0(idp) + (1-m0(idp))*min(strsw(j),strsn(j))
        n_frac_root(j) = (pRx(idp) * pRn(idp)) / (pRn(idp) + (pRx(idp) - pRn(idp)) * f_phys * m_root)
        
        
        !Calculation of stem fraction (n_frac_stem)
        !Partitioning coefficient pFsPower and pFsConst
        pFsPower = Log( pFS20(idp) / pFS2(idp) ) / Log( 20.d0 / 2.d0 )
        pFsConst = pFS2(idp) / 2.d0 ** pfsPower
        
        !Biomass of each tree
        bio_tree(j) = bio_stem(j)/bio_stock(j)
        
        !Average breast height diameter (dbh) for each tree
        dbh_tree(j) = (bio_tree(j)/aWs(idp))**(1/nWs(idp))
        
        !Calculation of pFS
        pFS = pFsConst*dbh_tree(j)**pFsPower
        
        !n_frac_stem
        n_frac_stem(j) = (1-n_frac_root(j))/(1+pFS)
        
        !Calculation of foliage fraction (n_frac_foliage)
        n_frac_foliage(j) = (1 - n_frac_root(j) - n_frac_stem(j))
        
        
        !Updating the biomass increase and total biomass for each tree component
        !root_incr(j) = n_frac_root(j) * bioday * reg
        root_incr(j) = n_frac_root(j) * bioday3PG(j)
        bio_root(j) = bio_root(j) + root_incr(j)
        
        !stem_incr(j) = n_frac_stem(j) * bioday * reg
        stem_incr(j) = n_frac_stem(j) * bioday3PG(j)
        bio_stem(j) = bio_stem(j) + stem_incr(j)

        !foliage_incr(j) = n_frac_foliage(j) * bioday * reg
        foliage_incr(j) = n_frac_foliage(j) * bioday3PG(j)
        bio_foliage(j) = bio_foliage(j) + foliage_incr(j)
        

        !Calculation of litterfall and root turnover -------------------------------------------------------------
        !Calculation of gammaF
         if( tgammaF(idp) * gammaF1(idp) == 0.d0 ) then
            gammaF = gammaF1(idp)
        else
            kg = 12.d0 * Log(1.d0 + gammaF1(idp) / gammaF0(idp)) / tgammaF(idp)
            gammaF = gammaF1(idp) * gammaF0(idp) / (gammaF0(idp) + (gammaF1(idp) - gammaF0(idp)) * Exp(-kg * tree_age(j)))
        end if
        
        litterfall(j) = gammaF*bio_foliage(j)
        root_turnover(j) = gammaR(idp)*bio_root(j) 
        
        !!!!printing and testing
        !if (j == 7) then

        !end if
        
       !Calculation of stem mortality based on tree age -----------------------------------------------------------
      ! mort_age_key = 1.0 !Activation key
       !Calculation of number of tree mortality due to tree age!!!
        !gammaN = gammaN1
        if (mort_age_key(idp) == 1.0) then
            if ( tgammaN(idp) /= 0.d0 ) then
                gammaN = gammaN1(idp) + (gammaN0(idp) - gammaN1(idp)) * Exp(-0.69314718 * ( tree_age(j) / tgammaN(idp)) ** ngammaN(idp)) !!LN2 = 0.69314718
            end if
            if (gammaN > 0.) then
                mort_stress = gammaN * bio_stock(j)/12.d0/100.d0
                mort_stress = Min(mort_stress, bio_stock(j)) !mort stress cannot be more than availabe no of stems
        
                foliage_loss_age(j) = mF(idp) * mort_stress * (bio_foliage(j) / bio_stock(j)) 
                stem_loss_age(j) = mS(idp) * mort_stress * (bio_stem(j) / bio_stock(j))  
                root_loss_age(j) = mR(idp) * mort_stress * (bio_root(j) / bio_stock(j))
            end if
        else 
            foliage_loss_age(j) = 0.d0
            stem_loss_age(j) = 0.d0 
            root_loss_age(j) = 0.d0
            mort_stress = 0.d0
        end if
            
        !Calculation of stem mortality based on self-thinning rule --------------------------------------------------------------------------------------    
        !self-thinning mortality activation key
        !mort_thinning_key = 1.0
        
        if (mort_thinning_key(idp) == 1.0) then
            !Calculation of stem mortality based on the self-thinning rule
            PI = 4.D0*DATAN(1.D0)
            basal_area = dbh_tree(j)**2.d0/4.d0*PI*bio_stock(j) / 10000.d0
            basal_area_prop =  basal_area/basal_area
            stems_n_ha = bio_stock(j)/basal_area_prop
        
            biom_tree_max = wSx1000(idp)*(1000.d0/stems_n_ha)**thin_power(idp)

           ! function f_get_mortality(stems_n, WS, mS, wSx1000, thinPower) result(mort_n)
            !f_get_mortality( stems_n_ha(i), biom_stem(i) / basal_area_prop(i) ,mS(i), wSx1000(i), thinPower(i) ) 
            

            if (biom_tree_max < bio_tree(j)) then    
                !actual calculation
                n = stems_n_ha/1000.d0
                x1 = 1000.d0 * mS(idp) * (bio_stem(j)/1000/basal_area_prop) / stems_n_ha
                t = 0

                do
                    t = t + 1
                    if (n <= 0.d0) exit !added in 3PG+

                    x2 = wSx1000(idp) * n ** (1.d0 - thin_power(idp))
                    fN = x2 - x1 * n - (1.d0 - mS(idp)) * (bio_stem(j)/1000/basal_area_prop)
                    dfN = (1.d0 - thin_power(idp)) * x2 / n - x1
                    dN = -fN / dfN
                    n = n + dN

                    if (abs(dN) <= accuracy .Or. t >= 5) exit

                end do

                mort_n_thinning = stems_n_ha - 1000.d0 * n
            
                if (mort_n_thinning < bio_stock(j)) then
                
                    foliage_loss_thinning(j) = mF(idp) * mort_n_thinning * (bio_foliage(j) / bio_stock(j))
                    root_loss_thinning(j) = mR(idp) * mort_n_thinning * (bio_root(j) / bio_stock(j))
                    stem_loss_thinning(j) = mS(idp) * mort_n_thinning * (bio_stem(j) / bio_stock(j))

                else
                    foliage_loss_thinning(j) = 0.d0
                    root_loss_thinning(j) = 0.d0
                    stem_loss_thinning(j) = 0.d0
                end if
                !print*, biom_tree_max, bio_tree(j), mort_n_thinning, stems_n_ha, bio_stock(j)
            end if

        else
            !foliage_loss_thinning(j) = 0.d0
            !root_loss_thinning(j) = 0.d0
            !stem_loss_thinning(j) = 0.d0
            !mort_n_thinning = 0.d0
        end if
            
        
        !Calculating the total biomass lost each day due to turnover, age-related mortality, and self-thinning----------------------------------------------------  
        !total biomass loss for foliage
        bio_loss_foliage(j) = litterfall(j) + foliage_loss_age(j) +  foliage_loss_thinning(j)
        
        !total biomass loss for root
        bio_loss_root(j) = root_turnover(j) + root_loss_age(j) + root_loss_thinning(j)
        
       ! print*,i,j, tgammaF*gammaF1,gammaF,litterfall(j), bio_loss_foliage(j)
        
        !total biomass loss for stem
        bio_loss_stem(j) = stem_loss_age(j) + stem_loss_thinning(j)
        
        !!combining all biomass losses for the day
        bio_loss(j) = bio_loss_foliage(j) + bio_loss_root(j) + bio_loss_stem(j)
        
        !!Updating root, foliage, and stem biomass and tree count at the end of each day-------------------------------------------------------------------------
        bio_foliage(j) = bio_foliage(j) - bio_loss_foliage(j)
        bio_root(j) = bio_root(j) - bio_loss_root(j)
        bio_stem(j) = bio_stem(j) - bio_loss_stem(j) 
        bio_stock(j) = bio_stock(j) - mort_stress - mort_n_thinning
        !!!Calculation of biomass loss on foliage, stem, and root fraction due to daily mortality and tree mortality
       
         !Total 3PG Biomass
        bio_ms_3PG(j) = bio_root(j) + bio_stem(j) + bio_foliage(j)      
        AGB(j) = bio_stem(j) + bio_foliage(j)
        
        !Calculation of LAI based on 3PG---------------------------------------------------------------------------------------------------------
        SLA(j) = SLA1(idp) + (SLA0(idp) - SLA1(idp)) * Exp(-0.69314718*(tree_age(j)/tSLA(idp))**2.0)
        LAI_3PG(j) = SLA(j)*bio_foliage(j)/1000*0.1
        
        !!Providing the calculated 3PG biomass and LAI to bio_ms and laiday values to connect to all other sub routines
        
        bio_ms(j) = bio_ms_3PG(j)
        laiday(j) = LAI_3PG(j)
        
        !! calculate fraction of total biomass that is in the roots
        rwt(j) = bio_root(j)/bio_ms(j)  
        
        !! calculate new canopy height

        cht(j) = aH(idp)*dbh_tree(j)**nHB(idp)

        
        select case (1)  !! Junyu Qi 02/26/2022
        case(1)    !! Junyu Qi 02/26/2022
        !!!Residue calculation - copied from dormant-yang.
        !orgc_f = 0.
        BLG1 = 0.
        BLG2 = 0.
        BLG3 = 0.
        CLG = 0.
        sol_min_n = 0.
        resnew_n = 0.
        resnew_ne = 0.
        LMF = 0.
        LSF = 0.
        LSLF = 0.
        LSNF = 0.
        LMNF = 0.
        XX = 0.
        RLN = 0.
        RLR = 0.
        
        !idorm(j) = 1
        resnew = 0.
        resnew = bio_loss(j)
          
            
        !!add by zhang
        !!===================
        if (cswat == 2) then 
            rsdc_d(j) = rsdc_d(j) + resnew*0.42
        end if
        !!add by zhang
        !!===================

            !!insert new biomss by zhang
            !!=============================
        if (cswat == 2 .and. resnew > 0.) then
	        !!all the lignin from STD is assigned to LSL, 
	        !!add STDL calculation
	        !!
	        !sol_LSL(k,ihru) = sol_STDL(k,ihru)
	        !CLG=BLG(3,JJK)*HUI(JJK)/(HUI(JJK)+EXP(BLG(1,JJK)-BLG(2,JJK)*&HUI(JJK))
	        ! 52  BLG1 = LIGNIN FRACTION IN PLANT AT .5 MATURITY
            ! 53  BLG2 = LIGNIN FRACTION IN PLANT AT MATURITY
            !CROPCOM.dat BLG1 = 0.01 BLG2 = 0.10
            !SUBROUTINE ASCRV(X1,X2,X3,X4)
            !EPIC0810
            !THIS SUBPROGRAM COMPUTES S CURVE PARMS GIVEN 2 (X,Y) POINTS.
            !USE PARM
            !XX=LOG(X3/X1-X3)
            !X2=(XX-LOG(X4/X2-X4))/(X4-X3)
            !X1=XX+X3*X2
            !RETURN
            !END 
            !HUI(JJK)=HU(JJK)/XPHU               
                
            BLG1 = 0.01/0.10
            BLG2 = 0.99
                !!   !! SWAT-3PG-----------
            if(ifor==2) then
                               
                !if(idplt(j) == 6 .or. idplt(j) == 7 .or. idplt(j) == 8) then
                if (idc(idplt(j)) == 7) then
                
                    BLG3 = 0.2 
                else
                    BLG3 = 0.1 
                endif
                
            else
                BLG3 = 0.1 
            end if 
                !!   !! SWAT-3PG-----------
            XX = log(0.5/BLG1-0.5)
            BLG2 = (XX -log(1./BLG2-1.))/(1.-0.5)
            BLG1 = XX + 0.5*BLG2
            !CLG=BLG3*phuacc(j)/(phuacc(j) +  EXP(BLG1-BLG2*phuacc(j))) !! Question - Need to understand what this is for and how the use of phuacc impacts this eqn
            CLG=BLG3*bm_frac(j)/(bm_frac(j) +  EXP(BLG1-BLG2*bm_frac(j)))!Ritesh - changed phuacc to bm_frac
	        !if (k == 1) then
		    !sf = 0.05
	        !else
		    !sf = 0.1
	        !end if	

            !kg/ha  
	        sol_min_n = 0.	
	        sol_min_n = (sol_no3(1,j)+sol_nh3(1,j))
	          	          
	          
                
            !!resnew = bio_ms(j) * bio_leaf(idplt(j)) !!Original 
                
            !resnew = bio_loss !!Modified Ritesh - litterfall is calculated using 
          !  resnew_n = resnew * plantn(j) !!Does this line need to be added     	  
           
            !! added by Junyu Qi, 02/26/2022
             !! plantn(j) is not a fraction, use pltfr_n(j) instead.   
            !!! Better to use different N fraction for different component of biomass, leaf, stem ,root
             resnew_n = resnew * pltfr_n(j)   !bio_leaf(idplt(j))    !! added by Junyu Qi, 02/26/2022
	       
	       
	        if (resnew > 10.) then   
        	    resnew_ne = resnew_n + sf * sol_min_n        
        	else
        	    resnew_ne = resnew_n
        	end if
        	!Not sure 1000 should be here or not!
        	!RLN = 1000*(resnew * CLG/(resnew_n+1.E-5))
        	RLN = (resnew * CLG/(resnew_n+1.E-5))
        	RLR = MIN(.8, resnew * CLG/1000/(resnew/1000+1.E-5))
        	    
        	LMF = 0.85 - 0.018 * RLN
        	if (LMF <0.01) then
        	    LMF = 0.01
        	else
        	    if (LMF >0.7) then
        	        LMF = 0.7
        	    end if
        	end if      	  
	          !if ((resnew * CLG/(resnew_n+1.E-5)) < 47.22) then
		        !    LMF = 0.85 - 0.018 * (resnew * CLG/(resnew_n+1.E-5))
	          !else
		        !    LMF = 0.
	          !end if 	

	        LSF =  1 - LMF  
        	  
	        sol_LM(1,j) = sol_LM(1,j) + LMF * resnew
	        sol_LS(1,j) = sol_LS(1,j) + LSF * resnew
        	  
	        !here a simplified assumption of 0.5 LSL
	        LSLF = 0.0
	        LSLF = CLG          
	          
	        sol_LSL(1,j) = sol_LSL(1,j) + RLR* LSF * resnew	          
	        sol_LSC(1,j) = sol_LSC(1,j) + 0.42*LSF * resnew  
	          
	        sol_LSLC(1,j) = sol_LSLC(1,j) + RLR*0.42*LSF * resnew
	        sol_LSLNC(1,j) = sol_LSC(1,j) - sol_LSLC(1,j)              
                
            !X3 = MIN(X6,0.42*LSF * resnew/150) 
                
	        if (resnew_ne >= (0.42 * LSF * resnew /150)) then
		        sol_LSN(1,j) = sol_LSN(1,j) + 0.42 * LSF * resnew / 150
		        sol_LMN(1,j) = sol_LMN(1,j) + resnew_ne - (0.42 * LSF * resnew / 150) + 1.E-25
	        else
		        sol_LSN(1,j) = sol_LSN(1,j) + resnew_ne
		        sol_LMN(1,j) = sol_LMN(1,j) + 1.E-25
	        end if	
        	
	        !LSNF = sol_LSN(1,j)/(sol_LS(1,j)+1.E-5)	
        	  
	        sol_LMC(1,j) = sol_LMC(1,j) + 0.42 * LMF * resnew	
	        !LMNF = sol_LMN(1,j)/(sol_LM(1,j) + 1.E-5)           
            if (resnew > 10.) then               !!Junyu Qi
            !update no3 and nh3 in soil 
                sol_no3(1,j) = sol_no3(1,j) * (1-sf)
                sol_nh3(1,j) = sol_nh3(1,j) * (1-sf)
            end if

        end if   !!Addedy by Ritesh- End
            !!insert new biomss by zhang
            !!===========================


            
        sol_rsd(1,j) = sol_rsd(1,j) + resnew
        sol_rsd(1,j) = Max(sol_rsd(1,j),0.)
        sol_fon(1,j) = resnew * pltfr_n(j) + sol_fon(1,j)
        sol_fop(1,j) = resnew * pltfr_p(j) + sol_fop(1,j)
        !bio_hv(icr(j),j) = bio_ms(j) + bio_hv(icr(j),j)   !!Harvested biomass
        !bio_yrms(j) = bio_yrms(j) + bio_ms(j) / 1000.     !!
        !bio_ms(j) = bio_ms(j) * (1. - bio_leaf(idplt(j)))
        plantn(j) = plantn(j) - resnew * pltfr_n(j)
        plantp(j) = plantp(j) - resnew * pltfr_p(j)
        strsw(j) = 1.
        !!laiday(j) = alai_min(idplt(j))- Original
        !laiday(j) = 0. !!Modified by Ritesh - LAI for deciduous is set to 0 with 3PG.
        !phuacc(j) = 0.
        !laimxfr(j) = 0.        !Sue White - dormancy
        ncrops(icr(j),j) = ncrops(icr(j),j) + 1

        !! this biomixing process will reduce residue and orgN in the first soil layer  !! Junyu Qi 02/26/2022
       !! if (biomix(j) > .001) call newtillmix (j,biomix(j))            !! !! Junyu Qi 02/26/2022

       ! print*,i, bio_loss ,resnew_n ,sol_rsd(1,j) ,sol_fon(1,j),sol_LSN(1,j),   sol_LMN(1,j)  !! Junyu Qi 02/26/2022
		      
          

        end select    !! Junyu Qi 02/26/2022    
                
                
          !! beginning of perennial (pasture/alfalfa) dormant period
    end subroutine
    
    


