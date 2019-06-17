#!/bin/bash

# Exit on error
set -e
# verbose (print commands)
set -x

modeltype="ELMuserpft_defaultparams_noP_adspinuptest"
sitename="Kougarok"
# For spinup:
compset="ICB1850CNPRDCTCBC"
# For historical:
# compset="ICB20TRCNPRDCTCBC"

# Create case in new directory
cd /home/b0u/models/E3SM/cime/scripts
./create_newcase --case /home/b0u/cases/${modeltype}_${sitename}_${compset} --res CLM_USRDAT --mach cades --compiler gnu --compset $compset --project ccsi --walltime 24:00:00

cd /home/b0u/cases/${modeltype}_${sitename}_${compset}
./xmlchange SAVE_TIMING=FALSE,CLM_USRDAT_NAME=1x1pt_kougarok-NGEE,LND_DOMAIN_PATH=/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/share/domains/domain.clm,LND_DOMAIN_FILE=domain.lnd.51x63pt_kougarok-NGEE_TransA_navy.nc

./xmlchange NTASKS=6
./xmlchange NTASKS_ESP=1

# If doing cold start
./xmlchange CLM_ACCELERATED_SPINUP=on,CLM_FORCE_COLDSTART=on

# comment out PFLOTRAN_INC and PFLOTRAN_LIB


# Run case setup. It will fail because some things are missing from the nml file
./case.setup || echo "case.setup failed. Continuing, as this was expected"

# Add things to nml file
cat >> user_nl_clm << EOF
fsurdat = '/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/lnd/clm2/surfdata_map/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c181115default.nc'
paramfile = '/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/lnd/clm2/paramdata/clm_params_c180524.nc'
!fsurdat = '/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/lnd/clm2/surfdata_map/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c181115-sub12.nc'
!paramfile = '/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/lnd/clm2/paramdata/clm_params_c180524-sub12.nc'
nyears_ad_carbon_only = 25
spinup_mortality_factor = 10
metdata_type = 'gswp3v1_daymet'
metdata_bypass = '/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/atm/datm7/atm_forcing.datm7.GSWP3_daymet.1x1pt_kougarok-NGEE/cpl_bypass_full'
aero_file = '/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/ACME_inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_1850_mean_1.9x2.5_c090803.nc'
CO2_file = '/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/ACME_inputdata/atm/datm7/CO2/fco2_datm_1765-2007_c100614.nc'

! History file controls
hist_empty_htapes = .false. ! If true, turns off default history output (of everything at subgrid level. If false, first file is defaults (need to account for that in subsequent controls
hist_fincl2 = 'LEAFC','FROOTC','WOODC','LIVECROOTC','TLAI','CPOOL','STORVEGC','STORVEGN','LEAFN','FROOTN','LIVECROOTN','LIVESTEMC','LIVESTEMN',
              'DEADCROOTC','DEADCROOTN','DEADSTEMC','DEADSTEMN','TOTVEGC','TOTVEGN','GPP','NPP','HTOP','AGNPP','BGNPP'
              'LEAF_MR','FROOT_MR','LIVESTEM_MR','LIVECROOT_MR','MR','VCMAX25TOP','GR','AVAILC','PLANT_CALLOC','EXCESS_CFLUX','XSMRPOOL_RECOVER','XSMRPOOL','CPOOL'
hist_dov2xy = .true., .false. ! True means subgrid-level output, false means patch (PFT) level output
hist_nhtfrq = 0, -24 ! Output frequency. 0 is monthly, -24 is daily, 1 is timestep
hist_mfilt  = 365 365 ! History file writing frequency: number of points from hist_nhtfrq
hist_avgflag_pertape = 'A', 'A', 'A' ! A for averaging, I for instantaneous
EOF


# ./xmlchange CLM_BLDNML_OPTS="-maxpft 12" --append 
./xmlchange CLM_BLDNML_OPTS="-bgc bgc -nutrient cn -nutrient_comp_pathway rd  -soil_decomp ctc -methane -nitrif_denitrif  -maxpft 12"


./xmlchange STOP_OPTION=nyear,REST_OPTION=nyear,REST_N=20
./xmlchange STOP_N=150


echo "Two more things to do by hand:"
echo " 1. Comment out PFLOTRAN_INC and PFLOTRAN_LIB in env_mach_specific.xml"
echo " 2. Add "-DCPL_BYPASS" to CPPDEFS in Macros.make"
echo "Then run case.build and case.submit to start run"
