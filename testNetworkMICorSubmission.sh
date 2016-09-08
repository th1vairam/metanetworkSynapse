#!/bin/bash

. ./config.sh

#number of cores to reserve for job
nthreads=1

#path to csv file with annotations to add to file on Synapse
annotationFile="$outputpath/annoFile.txt"

#path to csv file with provenance to add to file on synapse
provenanceFile="$outputpath/provenanceFile.txt"

# prefixes
c3="c3net"
mrnet="mrnet"
wgcna="wgcna"

qsub -v s3=$s3,dataFile=$dataFile,pathv=$pathv,numberCore=$nthreads,outputpath=$outputpath,s3b=$s3b,parentId=$parentId,annotationFile=$annotationFile,provenanceFile=$provenanceFile -pe mpi $nthreads -S /bin/bash -V -cwd -N $c3 -e "$outputpath/${c3}error.txt" -o "$outputpath/${c3}out.txt" $pathv/networkScripts/c3net.sh

qsub -v s3=$s3,dataFile=$dataFile,pathv=$pathv,numberCore=$nthreads,outputpath=$outputpath,s3b=$s3b,parentId=$parentId,annotationFile=$annotationFile,provenanceFile=$provenanceFile -pe mpi $nthreads -S /bin/bash -V -cwd -N $mrnet -e "$outputpath/${mrnet}error.txt" -o "$outputpath/${mrnet}out.txt" $pathv/networkScripts/mrnet.sh

qsub -v s3=$s3,dataFile=$dataFile,pathv=$pathv,numberCore=$nthreads,outputpath=$outputpath,s3b=$s3b,parentId=$parentId,annotationFile=$annotationFile,provenanceFile=$provenanceFile -pe mpi $nthreads -S /bin/bash -V -cwd -N $wgcna -e "$outputpath/${wgcna}error.txt" -o "$outputpath/${wgcna}out.txt" $pathv/networkScripts/wgcna.sh

#s3=$s3 dataFile=$dataFile pathv=$pathv c3net=1 mrnet=1 wgcnaTOM=1 sparrowZ=0 lassoCV1se=0 ridgeCV1se=0 genie3=0 tigress=0 numberCore=$nthreads outputpath=$outputpath s3b=$s3b parentId=$parentId annotationFile=$annotationFile provenanceFile=$provenanceFile bucket=$bucket $pathv/pushNetS3.sh
