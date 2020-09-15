#!/bin/bash

#prefix="Kougarok_daymet_rhizometurnovercode_fixed"
prefix=$1
today=`date +'%Y%m%d'`

#bash archive_sim.bash ${prefix}_E3SMpfts_AK-KM64_ICB20TRCNPRDCTCBC ${prefix}_E3SMpfts_$today
#bash archive_sim.bash ${prefix}_E3SMpfts_soilthickness_AK-KM64_ICB20TRCNPRDCTCBC ${prefix}_E3SMpfts_soilthickness_$today
#bash archive_sim.bash ${prefix}_E3SMpfts_communities_AK-KM64_ICB20TRCNPRDCTCBC ${prefix}_E3SMpfts_communities_$today
bash archive_sim.bash ${prefix}_Arcticpfts_AK-KM64_ICB20TRCNPRDCTCBC ${prefix}_Arcticpfts_$today

