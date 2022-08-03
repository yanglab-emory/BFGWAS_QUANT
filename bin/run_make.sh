#!/usr/bin/bash

#######################################################################
### Input Arguments for Run Makefile Bash Script
#######################################################################
# --wkdir: Working directory with writing access
# --mkfile: Path to make file
# --njob: Number of jobs to run in parallele
#######################################################################

VARS=`getopt -o "" -a -l \
wkdir:,mkfile:,njob: \
-- "$@"`

if [ $? != 0 ]
then
    echo "Terminating....." >&2
    exit 1
fi

eval set -- "$VARS"

while true
do
    case "$1" in
	    --wkdir|-wkdir) wkdir=$2; shift 2;;
        --mkfile|-mkfile) mkfile=$2; shift 2;;
        --njob|-njob) njob=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done


cd $wkdir

# echo Run make clean
# make -f ${mkfile} clean

echo Run makefile $mkfile with $njob parallel jobs

make -k -C ${wkdir} -f ${mkfile} -j ${njob} > ${wkdir}/make.output 2> ${wkdir}/make.err

exit


