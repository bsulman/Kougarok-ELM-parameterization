#!/bin/bash

# Exit on errors
set -e

OLMT_dir=/home/b0u/models/OLMT

pftvars="'LEAFC','FROOTC','WOODC','LIVECROOTC','TLAI','CPOOL','STORVEGC','STORVEGN','LEAFN','FROOTN','LIVECROOTN','LIVESTEMC','LIVESTEMN','DEADCROOTC','DEADCROOTN','DEADSTEMC','DEADSTEMN','TOTVEGC','TOTVEGN','GPP','NPP','HTOP','AGNPP','BGNPP','LEAF_MR','FROOT_MR','LIVESTEM_MR','LIVECROOT_MR','MR','VCMAX25TOP','GR','AVAILC','PLANT_CALLOC','EXCESS_CFLUX','XSMRPOOL_RECOVER','XSMRPOOL','CPOOL','LEAFC_XFER_TO_LEAFC','FROOTC_XFER_TO_FROOTC','DOWNREG','INIT_GPP','SMINN_TO_NPOOL','PLANT_PDEMAND','PLANT_NDEMAND','SMINP_TO_PPOOL','FROOT_PROF','CROOT_PROF'"

cmd="python2 ${OLMT_dir}/site_fullrun.py --site AK-KM64 --sitegroup NGEEArctic --dailyvars --var_list_pft ${pftvars} --nopftdyn --surfdata_grid --nyears_ad_spinup 200 --nyears_final_spinup 600 --tstep 1 --cpl_bypass --spinup_vars --machine cades --compiler gnu --mpilib openmpi --batch_build --model_root /home/b0u/models/E3SM --caseroot /home/b0u/cases"
#cmd="${cmd} --gswp3"
cmd="${cmd} --metdir /home/b0u/driver_data/atm_forcing.datm7.GSWP3_daymet.1x1pt_kougarok-NGEE/cpl_bypass_full" 
#cmd="${cmd} --site_forcing kougarok-NGEE"

# Run 1: Downsampled surface properties and PFTs from original E3SM
${cmd} --caseidprefix Kougarok_daymet_E3SMpfts 

# Run 2: E3SM PFTs but with varying soil depths
# So, need local domain and surface properties specification but original NGEE-Arctic PFTs, same distribution in each cell
# I set this up with a netcdf file manually set to the right PFTs, but might be better to autogenerate it from the Run 1 file?
${cmd} --caseidprefix Kougarok_daymet_E3SMpfts_soilthickness \
--var_soilthickness --np=6 \
--domainfile=/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/share/domains/domain.clm/domain.lnd.51x63pt_kougarok-NGEE_TransA_navy.nc \
--surffile=/home/b0u/Kougarok_param_edits/param_files/surfdata_Kougarok_downscaled_PFTs_soildepths.nc

# Run 3: E3SM PFTs but with different PFT distributions across the gradient
# Local domain file, local surface file with E3SM PFT definitions but varying distributions based on original/Arctic PFT crosstalk
# Check PFT defs in surffile compared with downscaled E3SM PFTs
#   Downscaled E3SM uses evergreen_needleleaf_boreal_tree, not broadleaf_evergreen_shrub
${cmd} --caseidprefix Kougarok_daymet_E3SMpfts_communities \
--var_soilthickness --np=6 \
--domainfile=/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/share/domains/domain.clm/domain.lnd.51x63pt_kougarok-NGEE_TransA_navy.nc \
--surffile=/home/b0u/Kougarok_param_edits/param_files/surfdata_Kougarok_defaultPFTs_all-shrubs-decid-boreal.nc
# /home/b0u/Kougarok_param_edits/param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c190604default.nc

# Run 4: Kougarok transect surface data and Arctic PFTs
# Local domain file, local surface file with new PFTs and new PFT (obs) distributions and parameters
${cmd} --caseidprefix Kougarok_daymet_Arcticpfts \
--var_soilthickness --np=6 \
--maxpatch_pft=12 \
--nopftdyn \
--domainfile=/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/share/domains/domain.clm/domain.lnd.51x63pt_kougarok-NGEE_TransA_navy.nc \
--surffile=/home/b0u/Kougarok_param_edits/param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c190604-sub12_updated_2019-06-17.nc \
--mod_parm_file=/home/b0u/Kougarok_param_edits/clm_params_updated.nc

