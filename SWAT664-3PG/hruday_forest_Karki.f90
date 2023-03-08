      subroutine hruday_F

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine writes daily HRU output to the output.hru file

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    aird(:)       |mm H2O        |amount of water applied to HRU on current
!!                                 |day

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ii          |none          |counter
!!    j           |none          |HRU number
!!    pdvas(:)    |varies        |array to hold HRU output values
!!    pdvs(:)     |varies        |array to hold selected HRU output values
!!                               |when user doesn't want to print all
!!    sb          |none          |subbasin number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm
       
      implicit none
      integer :: j, sb, ii, iflag, nlay, k
      real, dimension (mhruo) :: pdvas, pdvs
      character (len=4) :: cropname
      integer:: ly, idplant, icl
      real :: soil_nitrate
    
      j = 0
      j = ihru
      nlay = 0.
      nlay = sol_nly(j)
      k = 0


      sb = hru_sub(j)
      iflag = 0
      do ii = 1, itoth
        if (ipdhru(ii) == j) iflag = 1
      end do
      if (iflag == 0) return
             
      !Average soil  clay percent
      soil_nitrate = 0. 
        do k=1,nlay
            soil_nitrate = soil_nitrate + sol_no3(k,j)
        End do  
     
       pdvas = 0.
       pdvs = 0.
       
       pdvas(1) = cht(j)!bio_stock(j) !tree/ha
       pdvas(2) = dbh_tree(j)  !Avg DBH of each tree         
       pdvas(3) = bio_tree(j)/1000 !average biomass of each tree tons
       pdvas(4) = etday!basal_area   !basal_area
       pdvas(5) = bio_loss(j)/1000!bio_lossbioday3PG(j) !tree mortality due to thinning                            
       pdvas(6) = bio_loss_root(j)/1000      ! daily loss in root          
       pdvas(7) = bio_loss_stem(j)/1000       ! daily loss in stem                                
       pdvas(8) = (bio_loss_foliage(j) + frsd_dormancy_lf(j))/1000    ! daily loss in foliage   
       pdvas(9) = bio_root(j)/1000            ! root biomass     
       pdvas(10) = bio_stem(j)/1000           ! stem biomass
       pdvas(11) = bio_foliage(j)/1000        ! foliage biomass  
                   !!- (bio_loss_foliage(j) + bio_loss_stem(j) + bio_loss_root(j))!total biomass balance
       pdvas(12) = LAI_3PG(j)!SurQ_DON_0(j)/1000!LAI_3PG(j)!(bio_root(j) + bio_stem(j) + bio_foliage(j)) &
                    !/ 1000 !(bioday*reg) - ((bio_root(j) + bio_stem(j) + bio_foliage(j))-(bio_root_ini(j) + bio_stem_ini(j) + bio_foliage_ini(j)) &
                  !! + (bio_loss_foliage(j) + bio_loss_stem(j) + bio_loss_root(j) + frsd_dormancy_lf(j)))!total biomass balance
       pdvas(13) = NPPC_d(j)!bioday3PG(j)*0.42!LatQT_DON_0(j)/1000!laiday(j)!tree_age(j)!foliage_incr(j) - ((bio_foliage(j) - bio_foliage_ini(j)) + bio_loss_foliage(j) + frsd_dormancy_lf(j))                    !foliage biomass balance
       pdvas(14) = bio_ms_3PG(j)/1000!nh4_min(j)/1000!bio_ms_3PG(j)/1000 !3PG total biomass
       pdvas(15) = AGB(j)/1000
       pdvas(16) = rspc_d(j)
       pdvas(17) = NPPC_d(j) - rspc_d(j)




      call xmon 
          
      !if (j == 7) then
      !    print*, bio_root(j), bio_root_ini(j), root_incr(j), bio_loss_root(j)
      !    print*, bio_stem(j), bio_stem_ini(j), bio_foliage_ini(j), (bioday*reg)
      !            !! - (bio_loss_foliage(j) + bio_loss_stem(j) + bio_loss_root(j))
      !end if 
      
      
      idplant = idplt(j)
      if (idplant > 0) then
        cropname = cpnm(idplant)
      else
        cropname = "NOCR"
      endif
      
      
      write (399,1199) cropname, j, subnum(j),        &
          hruno(j), sb, nmgt(j), i_mo, icl(iida), iyr, hru_km(j),     &
         (pdvas(ii), ii = 1, 17)
1199  format (a4,i5,1x,a5,a4,i5,1x,i4,1x,i2,1x,i2,1x,i4,1x,e10.3,       &
        17f10.3)
      

      return

      end
