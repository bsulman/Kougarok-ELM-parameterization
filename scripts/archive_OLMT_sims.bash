#!/bin/bash

today=`date +'%Y%m%d'`

bash archive_sim.bash Kougarok_daymet_E3SMpfts_AK-KM64_ICB20TRCNPRDCTCBC E3SMpfts_$today
bash archive_sim.bash Kougarok_daymet_E3SMpfts_soilthickness_AK-KM64_ICB20TRCNPRDCTCBC E3SMpfts_soilthickness_$today
bash archive_sim.bash Kougarok_daymet_E3SMpfts_communities_AK-KM64_ICB20TRCNPRDCTCBC E3SMpfts_communities_$today
bash archive_sim.bash Kougarok_daymet_Arcticpfts_gramrhizomes_AK-KM64_ICB20TRCNPRDCTCBC Arcticpfts_$today

