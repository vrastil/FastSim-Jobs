#!/bin/bash

if [[ $HOSTNAME == koios* ]]; then
  echo "Use this script from localhost"
else
  SRC="koios:/home/users/vrastil/GIT/FastSim/jobs/output/report/"
  DEST="/home/michal/Documents/GIT/FastSim/report/clanek/"
  OTHER="-azmh --stats --info=progress2"
  rsync $OTHER $SRC $DEST
fi
