#!/bin/bash

# variables from command line
input_vcf=$1
reference=$2
regionfile=$3
# data sources
vep_tar_gz=$4
clinvar_gz=$5
dbnsfp_gz=$6
fordownload_tar_gz=$7
spliceai_snv_gz=$8
spliceai_indel_gz=$9
gnomad_gz=${10}
gnomad_gz2=${11}
CADD_snv=${12}
CADD_indel=${13}
phylop100bw=${14}
phylop30bw=${15}
phastc100bw=${16}
# parameters
nthreads=${17}
version=${18} # 101
assembly=${19} # GRCh38

# self variables
directory=VCFS/

# rename with version
dbnsfp=dbNSFP4.1a.gz

# rename dbNSFP
ln -s $dbnsfp_gz $dbnsfp
ln -s ${dbnsfp_gz}.tbi ${dbnsfp}.tbi
ln -s ${dbnsfp_gz%.*}.readme.txt dbnsfp.readme.txt

# unpack data sources
tar -xzf $vep_tar_gz
tar -xzf ${vep_tar_gz%%.*}.plugins.tar.gz
tar -xzf $fordownload_tar_gz

# setting up output directory
mkdir -p $directory

# command line VEP
# plugins
plugin_entscan="--plugin MaxEntScan,fordownload"
plugin_dbnsfp="--plugin dbNSFP,${dbnsfp},phyloP100way_vertebrate_rankscore,GERP++_RS,GERP++_RS_rankscore,SiPhy_29way_logOdds,SiPhy_29way_pi,PrimateAI_score,PrimateAI_pred,PrimateAI_rankscore,CADD_raw_rankscore,Polyphen2_HVAR_pred,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_score,SIFT_pred,SIFT_converted_rankscore,SIFT_score,REVEL_rankscore,REVEL_score,Ensembl_geneid,Ensembl_proteinid,Ensembl_transcriptid"
plugin_spliceai="--plugin SpliceAI,snv=${spliceai_snv_gz},indel=${spliceai_indel_gz}"
plugin_CADD="--plugin CADD,${CADD_snv},${CADD_indel}"

plugins="--dir_plugins VEP_plugins --plugin SpliceRegion,Extended --plugin TSSDistance $plugin_entscan $plugin_dbnsfp $plugin_spliceai $plugin_CADD"

# customs
cutsom_clinvar="--custom ${clinvar_gz},ClinVar,vcf,exact,0,ALLELEID,CLNSIG,CLNREVSTAT,CLNDN,CLNDISDB,CLNDNINCL,CLNDISDBINCL,CLNHGVS,CLNSIGCONF,CLNSIGINCL,CLNVC,CLNVCSO,CLNVI,DBVARID,GENEINFO,MC,ORIGIN,RS,SSR"
custom_gnomad="--custom ${gnomad_gz},gnomADg,vcf,exact,0,AC,AC-XX,AC-XY,AC-afr,AC-ami,AC-amr,AC-asj,AC-eas,AC-fin,AC-mid,AC-nfe,AC-oth,AC-sas,AF,AF-XX,AF-XY,AF-afr,AF-ami,AF-amr,AF-asj,AF-eas,AF-fin,AF-mid,AF-nfe,AF-oth,AF-sas,AF_popmax,AN,AN-XX,AN-XY,AN-afr,AN-ami,AN-amr,AN-asj,AN-eas,AN-fin,AN-mid,AN-nfe,AN-oth,AN-sas,nhomalt,nhomalt-XX,nhomalt-XY,nhomalt-afr,nhomalt-ami,nhomalt-amr,nhomalt-asj,nhomalt-eas,nhomalt-fin,nhomalt-mid,nhomalt-nfe,nhomalt-oth,nhomalt-sas"
custom_gnomad2="--custom ${gnomad_gz2},gnomADe2,vcf,exact,0,AC,AN,AF,nhomalt,AC_oth,AN_oth,AF_oth,nhomalt_oth,AC_sas,AN_sas,AF_sas,nhomalt_sas,AC_fin,AN_fin,AF_fin,nhomalt_fin,AC_eas,AN_eas,AF_eas,nhomalt_eas,AC_amr,AN_amr,AF_amr,nhomalt_amr,AC_afr,AN_afr,AF_afr,nhomalt_afr,AC_asj,AN_asj,AF_asj,nhomalt_asj,AC_nfe,AN_nfe,AF_nfe,nhomalt_nfe,AC_female,AN_female,AF_female,nhomalt_female,AC_male,AN_male,AF_male,nhomalt_male,AF_popmax"
custom_phylop100="--custom ${phylop100bw},phylop100verts,bigwig,exact,0"
custom_phylop30="--custom ${phylop30bw},phylop30mams,bigwig,exact,0"
custom_phastcons100="--custom ${phastc100bw},phastcons100verts,bigwig,exact,0"

customs="$cutsom_clinvar $custom_gnomad $custom_gnomad2 $custom_phylop100 $custom_phylop30 $custom_phastcons100"

basic_vep="--sift b --polyphen b --ccds --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --uniprot --tsl --appris --gene_phenotype --pubmed --var_synonyms --variant_class --mane"

# options and full command line
options="--fasta $reference --assembly $assembly --use_given_ref --offline --cache_version $version --dir_cache . $basic_vep --force_overwrite --vcf"

command="tabix -h $input_vcf {} > {}.sharded.vcf || exit 1; if [[ -e {}.sharded.vcf ]] || exit 1; then if grep -q -v '^#' {}.sharded.vcf; then vep -i {}.sharded.vcf -o ${directory}{}.vep.vcf $options $plugins $customs || exit 1; fi; fi; rm {}.sharded.vcf || exit 1"


# runnning VEP in parallel
echo "Running VEP"
cat $regionfile | xargs -P $nthreads -i bash -c "$command" || exit 1

# merging the results
echo "Merging vcf files"
array=(${directory}*.vep.vcf)

IFS=$'\n' sorted=($(sort -V <<<"${array[*]}"))
unset IFS

grep "^#" ${sorted[0]} > combined.vep.vcf

for filename in ${sorted[@]};
  do
    if [[ $filename =~ "M" ]]; then
      chr_M=$filename
    else
      grep -v "^#" $filename >> combined.vep.vcf
      rm -f $filename
    fi
  done

if [[ -v  chr_M  ]]; then
  grep -v "^#" $chr_M >> combined.vep.vcf
  rm -f $chr_M
fi

# compress and index output vcf
bgzip combined.vep.vcf || exit 1
tabix -p vcf combined.vep.vcf.gz || exit 1