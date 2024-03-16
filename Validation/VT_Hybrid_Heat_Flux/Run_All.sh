#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR wasson50_1p0.fds
$QFDS $DEBUG $QUEUE -d $INDIR wasson50_1p5.fds
$QFDS $DEBUG $QUEUE -d $INDIR wasson50_2p0.fds
$QFDS $DEBUG $QUEUE -d $INDIR wasson90_1p0.fds
$QFDS $DEBUG $QUEUE -d $INDIR wasson90_1p5.fds
$QFDS $DEBUG $QUEUE -d $INDIR wasson90_2p0.fds

echo FDS cases submitted
