    vcf_dir= $proj_dir/vcf/$prefix
    if [-d ]; then  
    mkdir -p
    freebayes \
        --bam $bam_file \
        --vcf $vcf_file \
        --fasta-reference $shared_dir/genomes/mg1655ref.fa \
        --pvar 0.0001 \
        --ploidy 2 \
        --min-alternate-count 4 \
        --no-ewens-priors \
        --use-mapping-quality \
        --no-marginals