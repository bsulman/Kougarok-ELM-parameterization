#!/bin/bash

# Exit on errors
set -e
set -u

run_1=false
run_2=false
run_3=false
run_4=false
name="Kougarok"
met="Daymet"
dry_run=false
while getopts '1234an:m:d' OPTION; do
    case "$OPTION" in
        1)
            run_1=true
            ;;
        2)
            run_2=true
            ;;
        3)
            run_3=true
            ;;
        4)
            run_4=true
            ;;
        a)
            run_1=true
            run_2=true
            run_3=true
            run_4=true
            ;;
        d)
            dry_run=true
            ;;
        n)
            name=$OPTARG
            ;;
        m)
            met=$OPTARG
            ;;
        ?)
            echo "Options are:"
            echo " -1: E3SM grid cell"
            echo " -2: E3SM PFTs with Kougarok soil depths"
            echo " -3: E3SM PFTs with Kougarok community distributions"
            echo " -4: Arctic PFTs"
            echo " -a: Do all four runs"
            echo " -n <name>: Use run prefix <name>. Default is 'Kougarok'"
            echo " -m <meteorology>: SNAP or Daymet (Default Daymet)"
            echo " -d: Dry run (just print commands)"
            exit 1
            ;;
    esac
done
shift $((OPTIND -1))

echo " Do run 1: $run_1"
echo " Do run 2: $run_2"
echo " Do run 3: $run_3"
echo " Do run 4: $run_4"
echo " Dry run?  $dry_run"
echo " Run prefix will be: $name"

OLMT_dir=/home/b0u/models/OLMT

pftvars="'LEAFC','FROOTC','WOODC','LIVECROOTC','TLAI','CPOOL','STORVEGC','STORVEGN','LEAFN','FROOTN','LIVECROOTN','LIVESTEMC','LIVESTEMN','DEADCROOTC','DEADCROOTN','DEADSTEMC','DEADSTEMN','TOTVEGC','TOTVEGN','GPP','NPP','HTOP','AGNPP','BGNPP','LEAF_MR','FROOT_MR','LIVESTEM_MR','LIVECROOT_MR','MR','VCMAX25TOP','GR','AVAILC','PLANT_CALLOC','EXCESS_CFLUX','XSMRPOOL_RECOVER','XSMRPOOL','CPOOL','LEAFC_XFER_TO_LEAFC','FROOTC_XFER_TO_FROOTC','DOWNREG','INIT_GPP','SMINN_TO_NPOOL','PLANT_PDEMAND','PLANT_NDEMAND','SMINP_TO_PPOOL'" #,'FROOT_PROF','CROOT_PROF'"

cmd="python2 ${OLMT_dir}/site_fullrun.py --site AK-KM64 --sitegroup NGEEArctic --dailyvars --var_list_pft ${pftvars} --nopftdyn --surfdata_grid --nyears_ad_spinup 200 --nyears_final_spinup 600 --tstep 1 --cpl_bypass --spinup_vars --machine cades --compiler gnu --mpilib openmpi --batch_build --model_root /home/b0u/models/E3SM --caseroot /home/b0u/cases"

if [[ ${met} = "Daymet" ]]; then
    cmd="${cmd} --metdir /home/b0u/driver_data/atm_forcing.datm7.GSWP3_daymet.1x1pt_kougarok-NGEE/cpl_bypass_full" 
elif [[ ${met} = "SNAP" ]]; then
    cmd="${cmd} --metdir /home/b0u/driver_data/SNAP_Kougarok/NCAR-CCSM4/precip_corrected/"
else
    echo "Error: met must be either 'Daymet' or 'SNAP'"
    exit 1
fi

#cmd="${cmd} --gswp3"
#cmd="${cmd} --metdir /home/dmricciuto/models/SNAP_met/1x1pt_US-Kou/"
#cmd="${cmd} --site_forcing kougarok-NGEE"

# Change forcing to RCP8.5
cmd="${cmd} --co2_file fco2_datm_rcp8.5_1765-2500_c110919.nc --aero_rcp85 --ndep_rcp85"

if [[ $dry_run = "true" ]]; then
    cmd="echo ${cmd}"
fi

# Run 1: Downsampled surface properties and PFTs from original E3SM
if [[ $run_1 = "true" ]]; then 
    ${cmd} --caseidprefix ${name}_E3SMpfts 
fi

# Run 2: E3SM PFTs but with varying soil depths
# So, need local domain and surface properties specification but original NGEE-Arctic PFTs, same distribution in each cell
# I set this up with a netcdf file manually set to the right PFTs, but might be better to autogenerate it from the Run 1 file?
if [[ $run_2 = "true" ]]; then
${cmd} --caseidprefix ${name}_E3SMpfts_soilthickness \
--var_soilthickness --np=6 \
--domainfile=/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/share/domains/domain.clm/domain.lnd.51x63pt_kougarok-NGEE_TransA_navy.nc \
--surffile=/home/b0u/Kougarok_param_edits/param_files/surfdata_Kougarok_downscaled_PFTs_soildepths.nc
fi

# Run 3: E3SM PFTs but with different PFT distributions across the gradient
# Local domain file, local surface file with E3SM PFT definitions but varying distributions based on original/Arctic PFT crosstalk
# Check PFT defs in surffile compared with downscaled E3SM PFTs
#   Downscaled E3SM uses evergreen_needleleaf_boreal_tree, not broadleaf_evergreen_shrub
if [[ $run_3 = "true" ]] ; then
${cmd} --caseidprefix ${name}_E3SMpfts_communities \
--var_soilthickness --np=6 \
--domainfile=/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/share/domains/domain.clm/domain.lnd.51x63pt_kougarok-NGEE_TransA_navy.nc \
--surffile=/home/b0u/Kougarok_param_edits/param_files/surfdata_Kougarok_defaultPFTs_all-shrubs-decid-boreal.nc
# /home/b0u/Kougarok_param_edits/param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c190604default.nc
fi

# Run 4: Kougarok transect surface data and Arctic PFTs
# Local domain file, local surface file with new PFTs and new PFT (obs) distributions and parameters
if [[ $run_4 = "true" ]]; then
${cmd} --caseidprefix ${name}_Arcticpfts \
--var_soilthickness --np=6 \
--maxpatch_pft=12 \
--nopftdyn \
--domainfile=/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/share/domains/domain.clm/domain.lnd.51x63pt_kougarok-NGEE_TransA_navy.nc \
--surffile=/home/b0u/Kougarok_param_edits/param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c190604-sub12_updated_2020-10-08.nc \
--mod_parm_file=/home/b0u/Kougarok_param_edits/clm_params_updated.nc
fi


