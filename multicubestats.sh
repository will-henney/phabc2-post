##
## Run globstats for all (or some) timesteps of a given model
##
## 14 Jul 2008 WJH - modified from mhdfitsemiss.sh
## 16 Jan 2008 WJH - modified for MHD models
## 01 Dec 2005 WJH - initial version
##
## Example usage:
##
##        sh multiglobstats.sh glob11-128x
##

runid=$1

# optionally restrict lowest and highest time
if [ $# -ge 2 ]; then
    i1=$2
else
    i1=0
fi
if [ $# -ge 3 ]; then
    i2=$3
else
    i2=9999
fi
if [ $# -ge 4 ]; then
    iskip=$4
else
    iskip=1
fi


execdir=$(dirname $0)

for i in $(seq $i1 $iskip $i2); do
    timeid=$(printf '%4.4i' $i)
    if [ -f ${runid}"-dd"${timeid}".fits" ]; then
	if [ ! -f ${runid}"-globstats"${timeid}".dat" ]; then
	    echo "Calculating cube statistics for $runid $i"
	    printf '%s\n' $runid $i | ${execdir}/cubestats > /dev/null 
	    printf '%s\n' $runid $i | ${execdir}/cubevstats > /dev/null 
	    printf '%s\n' $runid $i | ${execdir}/cubedenstats > /dev/null 
	fi
    fi
done
echo "No more files...done."
exit 0
