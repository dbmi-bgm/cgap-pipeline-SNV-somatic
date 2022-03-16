#!/bin/bash

filename=$1

vcf-validator $filename

if [[ $? -eq 0 ]];
  then
      echo -e "quickcheck\tOK" > ${filename%%.*}_integrity_check
  else
      echo -e "quickcheck\tFAILED" > ${filename%%.*}_integrity_check
fi
