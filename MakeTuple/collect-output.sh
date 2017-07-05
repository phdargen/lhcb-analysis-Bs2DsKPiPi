#!/bin/bash
#
#./collect_output FILE

#set the enviroment

#LbLogin -c x86_64-slc6-gcc49-opt
#lb-run LHCbDIRAC lhcb-proxy-init

FILE=$1

echo "Collecting files from $FILE..."
counter=1

while read line
do
   name='b2dhhh'
   outfile=$name\_$counter.root
   echo $outfile
   if [ -e /auto/data/dargent/BsDsKpipi/data16/$outfile ]
   then
       echo "   -> File existing, moving to next line....."
   else
       lb-run LHCbDIRAC dirac-dms-get-file $line
       mv b2dhhh.root /auto/data/dargent/BsDsKpipi/data16/$outfile
   fi
   echo "   Done."
   let counter=$counter+1
done < $FILE

echo "Done."


