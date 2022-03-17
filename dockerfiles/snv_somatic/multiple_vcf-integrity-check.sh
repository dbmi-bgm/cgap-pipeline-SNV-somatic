#!/bin/bash

filename=$1 || exit 1

#need to create isolate filename (no path, no extension)
no_path_filename=${filename##*/} || exit 1
#no_extension_filename=${no_path_filename%%.*} || exit 1

vcf-validator $filename

if [[ $? -eq 0 ]];
  then
      echo -e "quickcheck\tOK" > integrity_check_${no_path_filename}
  else
      echo -e "quickcheck\tFAILED" > integrity_check_${no_path_filename}
fi
