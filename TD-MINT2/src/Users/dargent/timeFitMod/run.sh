#!/bin/bash
#
#for i in {1..50}; do for j in {1..50};  do ./timeFit $i $j < signal_new.txt ; done ; done
for i in {1..50}; do ./timeFit $i $i < signal_new.txt ; done
