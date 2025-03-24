#!/bin/csh

if ( $#argv == 0 ) then
    echo "Usage: auto-check.sh config-[arch].sh"
    exit(1)
endif

echo "================================================="
date
echo " _______________________"
hostname

source $argv[1]

cd $M3DC1_CODE_DIR
rm unstructured/regtest/check.log

cd unstructured/regtest
./check $BATCH_SUFFIX > check.log
cat check.log

if ("$MAILTO" != "") then
    echo "Sending mail to $MAILTO"
    sendmail $MAILTO < check.log
endif
