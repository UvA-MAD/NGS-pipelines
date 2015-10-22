#!/bin/bash
if [ ! -z $1 ] ; then
    ORG=$1
else
    echo "Download hairpin sequences and annotation from Mirbase"
    echo "arguments: "
    echo " 1 : mandatory 3-letter code for organism"
    echo " 2 : optional version of mirbase default is 21" 
    echo " 3 : optional destination folder of fa and gff3 file " 
    exit; 
fi

#Mirbase version 
VERSION=${2:-21}
RESULT_DIR=${3:-.}

MIRBASE_FTPDIR="ftp://mirbase.org/pub/mirbase/${VERSION}/"

HAIRPINFA="/zfs/datastore0/group_root/MAD-RBAB/05_Reference-db/external/Mirbase/Mirbase${VERSION}_hairpin.fa"

#download hairpin sequence data if needed 
if [ ! -f ${HAIRPINFA}  ] ; then 
    echo "Downloading sequences for Mirbase $VERSION"
    curl ${MIRBASE_FTPDIR}hairpin.fa.gz | gunzip > ${HAIRPINFA}
else 
   echo "Found sequencedata for Mirbase ${VERSION} in 05_Reference-db/external/Mirbase, skipping download of hairpins"
fi

if [ ! -f ${HAIRPINFA}  ] ; then 
    echo "Could not download Mirbase $VERSION."
    exit
fi

#download gff3 annotation if needed 

ORG_ANN=${RESULT_DIR}/Mirbase${VERSION}_${ORG}.gff3 

if [ ! -f ${ORG_ANN}  ] ; then 
   curl ${MIRBASE_FTPDIR}/genomes/${ORG}.gff3  > ${ORG_ANN}  
   if [ $? -eq 78 ] ; then 
           echo "Annotation for ${ORG} not found on mirbase"
           rm ${ORG_ANN}
           exit
   fi 
else 
   echo "Found gff file for ${ORG} skipping download"
fi

if [ ! -f ${ORG_ANN}  ] ; then 
    echo "Could not download gff3 file for $ORG"
    exit
fi

#Mirbase hairpin collection contains sequences for all organisms
#Extract sequences for organism based on the start of the name. 
#Switch to state p=1 which prints sequence data if name is found. 
#b=1 is used to handle newline printing before printing id of next sequence
AWK_SCR='BEGIN {b=1; p=0} $1 ~ /^>'${ORG}'/ {if (b==0) {print "" }; b=0; p=1; printf("%s %s\n",$1,$2); next}  $1 ~ /^>/ {p=0; } p==1  { printf("%s",$0)} END {print}' 
awk  "$AWK_SCR" ${HAIRPINFA}  > ${RESULT_DIR}/Mirbase${VERSION}_${ORG}.fa
