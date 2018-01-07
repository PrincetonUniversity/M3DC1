#!/bin/bash

echo "calls to int*(: " `grep "int[12345](" *.f90 | grep -v "!" | grep -v "function" | wc -l`
echo "calls to intx*(: " `grep "intx[12345](" *.f90 | grep -v "!" | grep -v "function" | wc -l`
echo "calls to int*( in metricterms_new.f90: " `grep "int[12345](" metricterms_new.f90 | grep -v "!" | grep -v "function" | wc -l`
echo "calls to intx*( in metricterms_new.f90: " `grep "intx[12345](" metricterms_new.f90 | grep -v "!" | grep -v "function" | wc -l`

grep -n "(:,:," *.f90 | grep "int[12345]" | grep -v "intx[12345](" | grep -v "function"
grep -n "intx[12345](" *.f90 | grep -v "(:,:," | grep -v "function" | grep -v "temp[abcdef]"
grep -n "intx[12345](" *.f90 | grep -v "(e(" | grep -v "(mu79(" | grep -v "(trialx(" | grep -v "function" | grep -v "(trial(" | grep -v "(temp[abcdef],"
grep -n "temp[abcdef]" *.f90 | grep "int[12345]"
grep -n "temp79a" *.f90 | grep "tempa"
grep -n "temp79b" *.f90 | grep "tempb"
grep -n "temp79c" *.f90 | grep "tempc"
grep -n "temp79d" *.f90 | grep "tempd"
grep -n "temp79e" *.f90 | grep "tempe"
grep -n "temp79f" *.f90 | grep "tempf"
