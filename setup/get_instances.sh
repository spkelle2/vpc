#!/usr/bin/env bash

if [ -z $1 ]; then
  if [ -z "$VPC_DIR" ]
  then
    if [ ! -z "${REPOS_DIR}" ]
    then
      echo "Please define VPC_DIR (the root vpc dir, possibly ${REPOS_DIR}/vpc):"
    else
      echo "Please define VPC_DIR (the root vpc dir):"
    fi
    read VPC_DIR
    if [ -z "$VPC_DIR" ]
      then echo "Need to define VPC_DIR. Exiting."
      exit
    fi
  fi
  DATA_DIR="${VPC_DIR}/data/instances/original"

  echo "VPC_DIR is set to $VPC_DIR"
else
  DATA_DIR="${1}"
fi

echo "Instances will be downloaded to $DATA_DIR"

START_DIR=$PWD
cd $DATA_DIR

INSTANCES=miplib2017
echo "Downloading ${INSTANCES}"
wget -O ${INSTANCES}.zip 'https://miplib.zib.de/downloads/collection.zip'
mkdir -p ${DATA_DIR}/${INSTANCES}
unzip -q ${INSTANCES}.zip -d ${DATA_DIR}/${INSTANCES}
mv ${DATA_DIR}/${INSTANCES}/revised-submissions/*/instances/*.mps.gz ${DATA_DIR}/${INSTANCES}
#mv ${DATA_DIR}/${INSTANCES}/mas74.mps.gz ${DATA_DIR}/miplib2017/mas074.mps.gz
#mv ${DATA_DIR}/${INSTANCES}/mas76.mps.gz ${DATA_DIR}/miplib2017/mas076.mps.gz

# In MIPLIB 2010, unzipping will create desired directory
INSTANCES=miplib2010
echo "Downloading ${INSTANCES}"
wget -o ${INSTANCES}.zip 'https://miplib2010.zib.de/download/miplib2010-1.1.3-benchmark.zip'
unzip -q ${INSTANCES}.zip -d ${DATA_DIR}

INSTANCES=miplib2003
echo "Downloading ${INSTANCES}"
wget https://miplib2010.zib.de/miplib2003/download/miplib2003.tar
mkdir -p ${DATA_DIR}/${INSTANCES}
tar -xf ${INSTANCES}.tar --directory ${DATA_DIR}/${INSTANCES}

INSTANCES=miplib3
echo "Downloading ${INSTANCES}"
wget https://miplib2010.zib.de/miplib3/miplib3.tar.gz
tar -xf ${INSTANCES}.tar.gz --directory ${DATA_DIR}
rm ${DATA_DIR}/${INSTANCES}/miplib.cat
rm ${DATA_DIR}/${INSTANCES}/mps_format
rm ${DATA_DIR}/${INSTANCES}/references
cd ${DATA_DIR}/${INSTANCES}
for f in *; do mv "$f" "$f.mps"; done
for f in *; do gzip "$f"; done
cd ${DATA_DIR}

INSTANCES=miplib2
echo "Downloading ${INSTANCES}"
wget -O ${INSTANCES}.tar.gz 'https://miplib2010.zib.de/miplib2/miplib.tar.gz'
mkdir ${DATA_DIR}/${INSTANCES}
tar -xf ${INSTANCES}.tar.gz --directory ${DATA_DIR}/${INSTANCES}
mv ${DATA_DIR}/${INSTANCES}/miplib/* ${DATA_DIR}/${INSTANCES}
rm ${DATA_DIR}/${INSTANCES}/miplib.cat
rm ${DATA_DIR}/${INSTANCES}/mps_format
rm ${DATA_DIR}/${INSTANCES}/references
cd ${DATA_DIR}/${INSTANCES}
rm -r miplib
for f in *; do mv "$f" "$f.mps"; done
for f in *; do gzip "$f"; done
#mv *.mps.gz ..
#cd ..
cd ${DATA_DIR}

#cp ${DATA_DIR}/miplib2/stein27_nocard.mps.gz ${DATA_DIR}/miplib3/
#cp ${DATA_DIR}/miplib2/stein45_nocard.mps.gz ${DATA_DIR}/miplib3/

INSTANCES=coral
echo "Downloading ${INSTANCES}"
wget -O ${INSTANCES}.tar https://coral.ise.lehigh.edu/wp-content/uploads/mip-instances/instances/ALL_INSTANCE.tar
mkdir -p ${DATA_DIR}/${INSTANCES}
tar -xf ${INSTANCES}.tar --directory ${DATA_DIR}/${INSTANCES}
mv ${INSTANCES}/mcf2.mps.bz2 ${INSTANCES}/danoint.mps.bz2

echo "Finished downloading instances to $DATA_DIR"
cd ${START_DIR}
