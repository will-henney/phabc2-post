##
## Make emissivity datacube FITS files
##
## 16 Jan 2008 WJH - modified for MHD models
## 01 Dec 2005 WJH - initial version
##
## Example usage:
##
##        sh mhdfitsemiss.sh t50m32e01 Halpha
##

runid=$1
emtype=$2

# optionally restrict lowest and highest time
if [ $# -ge 3 ]; then
    i1=$3
else
    i1=1
fi
if [ $# -ge 4 ]; then
    i2=$4
else
    i2=9999
fi
if [ $# -ge 5 ]; then
    iskip=$5
else
    iskip=1
fi


execdir=$(dirname $0)

for i in $(seq $i1 $iskip $i2); do
    timeid=$(printf '%4.4i' $i)
    echo ${runid}"-dd"${timeid}".fits"
    if [ -f ${runid}"-dd"${timeid}".fits" ]; then
	if [ ! -f ${runid}"-te"${timeid}".fits" ]; then
	    echo "Calculating temperature cube for $runid $i"
	    printf '%s\n' $runid $i | ${execdir}/cubet > /dev/null
	fi
	if [ ! -f ${runid}"-e-"${emtype}${timeid}".fits" ]; then
	    echo "Calculating $emtype emissivity cube for $runid $i"
	    printf '%s\n' $runid $i $emtype | ${execdir}/cubeemiss > /dev/null
	fi
	# We only calculate the cubes, not the maps
	# The latter is done by makerotmap via makemovie.py
    fi
done
echo "No more files...done."
exit 0
