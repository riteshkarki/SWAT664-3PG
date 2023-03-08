      subroutine header_forest

!!    ~ ~ ~ PURPOSE ~ ~ ~                                               
!!    This subroutine defines header titles for the different output files
     !! added by Ritesh Karki
!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~                                    
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    hedb(:)     |NA            |column titles in subbasin output file
!!    hedr(:)     |NA            |column titles in reach output file
!!    hedrsv(:)   |NA            |column titles in reservoir output file
!!    heds(:)     |NA            |column titles in HRU output file
!!    hedwtr(:)   |NA            |column titles in HRU impoundment output 
!!                               |file
!!    icolb(:)    |none          |space number for beginning of column in
!!                               |subbasin output file
!!    icolr(:)    |none          |space number for beginning of column in
!!                               |reach output file
!!    icolrsv(:)  |none          |space number for beginning of column in
!!                               |reservoir output file
!!    icols(:)    |none          |space number for beginning of column in
!!                               |HRU output file
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~                        

      use parm

     !------biomass_massbalance_ritesh---------
            
       hedbm = (/"  TreeHt_m","       DBH"," AvjBmEch_T","      ETmm",  &       !4
                 " biolosday"," root_loss","  stem_loss","  fol_loss",  &       !8
                 "  bio_root","  bio_stem"," bio_folage","   LAI_3PG",  &
                 "  NPPC_day","  3PGBmass","       AGB" ," soil_resp",  &
                 "       NEP"  /)
       
       !hedbm = (/" bioday3PG"," OrgNRes_d","  OrgPRes_d"," OrgCRes_d",  &       !4
       !          "  bio_loss","  nplntday","  plantfr_n","   n_defic",  &       !8
       !          "  bio_root","  bio_stem"," bio_folage","   LAI_3PG",  &
       !          "  tree_age","  3PGBmass"," SWATBmass"  /)

       !hedbm = (/" bioday3PG"," OrgNRes_d","  orgn_fert"," org_grazf",  &       !4
       !          "  no3_immo","    ab_no3"," absorb_nh3"," solc_orgn",  &       !8
       !          " solcOrgIn","  sedOrgnO","   PerQ_DON"," Surq_DonO",  &
       !          " LatQTDonO","   nh4_min","  SWATBMass"  /)
      return
      end                                                            