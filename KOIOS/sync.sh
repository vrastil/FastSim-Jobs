#!/bin/bash

SRC="koios:/home/users/vrastil/GIT/FastSim-Jobs/output/"
DEST="/home/vrastil/Documents/GIT/FastSim/jobs/output/"

USAGE="
Usage:
\tSRC = $SRC
\tDEST = $DEST
cmd opt:
\tall:\tsync ALL files from SRC to DEST
\tdat:\tsync all *.dat files from SRC to DEST
\tpng:\tsync all *.png files from SRC to DEST
\tgif:\tsync all *.gif files from SRC to DEST
\tres:\tsync results files (results/*) from SRC to DEST
\tsupp:\tsync pwr_diff files (pwr_diff/*) from SRC to DEST
\tstack:\tsync stacked files (STACK_*/*) from SRC to DEST
\t+test:\tperform a trial run with no changes made (--dry-run)
"

if [[ $# != 1 ]] && [[ $# != 2 ]]; then
  printf "$USAGE"
  exit
fi

OTHER="-azmh --stats --info=progress2"
INCLUDE="--include=*.log --include=*/"
EXCLUDE="--exclude=*"

if [ $# == 2 ]; then
  if [ $2 != "test" ]; then
    printf "$USAGE"
    exit
  else
    echo Test run:
    OTHER+=" -n"
  fi
fi

if [ $1 == "all" ]; then
  echo "Syncing all files..."
  EXCLUDE=""
elif [ $1 == "dat" ]; then
  echo "Syncing *.dat files..."
  INCLUDE+=" --include=*.dat"
elif [ $1 == "png" ]; then
  echo "Syncing *.png files..."
  INCLUDE+=" --include=*.png"
elif [ $1 == "gif" ]; then
  echo "Syncing *.gif files..."
  INCLUDE+=" --include=*.gif"
elif [ $1 == "res" ]; then
  echo "Syncing results/* files..."
  INCLUDE+=" --include=results/*"
elif [ $1 == "supp" ]; then
  echo "Syncing pwr_diff/* files..."
  INCLUDE+=" --include=pwr_diff/*"
elif [ $1 == "stack" ]; then
  echo "Syncing STACK_*/*** files..."
  INCLUDE+=" --include=STACK_*/***"
else
  printf "$USAGE"
  exit
fi

rsync $OTHER $INCLUDE $EXCLUDE $SRC $DEST
