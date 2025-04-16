#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

# Tamanini cases
$QFDS $DEBUG $QUEUE -d $INDIR Tamanini_D01_Q30_DS_05.fds
$QFDS $DEBUG $QUEUE -p 8 -d $INDIR Tamanini_D01_Q30_DS_10.fds
$QFDS $DEBUG $QUEUE -p 16 -d $INDIR Tamanini_D01_Q30_DS_20.fds
$QFDS $DEBUG $QUEUE -d $INDIR Tamanini_D38_Q30_DS_05.fds
$QFDS $DEBUG $QUEUE -p 8 -d $INDIR Tamanini_D38_Q30_DS_10.fds
$QFDS $DEBUG $QUEUE -p 16 -d $INDIR Tamanini_D38_Q30_DS_20.fds
$QFDS $DEBUG $QUEUE -d $INDIR Tamanini_D38_Q62_DS_05.fds
$QFDS $DEBUG $QUEUE -p 8 -d $INDIR Tamanini_D38_Q62_DS_10.fds
$QFDS $DEBUG $QUEUE -p 16 -d $INDIR Tamanini_D38_Q62_DS_20.fds

echo FDS cases submitted
