#! /bin/bash
dir=$(dirname $0)
t1=$1
t2=$2
for line in Halpha N26584 O35007; do
    for suff in dm0 zerob-e01; do
	sh $dir/mhdfitsemiss.sh t50m32-weak-$suff $line $t1 $t2
    done
done
