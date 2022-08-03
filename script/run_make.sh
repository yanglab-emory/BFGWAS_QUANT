#!/usr/bin/bash

j=$1
wkdir=$2
study=$3
i=$SGE_TASK_ID
mkfile=${wkdir}/${study}_BFGWAS.mk
cd $wkdir

echo Run make clean
make -f ${mkfile} clean


echo Run make with $mkdir and j=$j parallel jobs

make -k -C ${wkdir} -f ${mkfile} -j ${j} > ${wkdir}/make.output 2> ${wkdir}/make.err

exit
