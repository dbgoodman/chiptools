java -Xmx2G -jar /opt/snpEff/snpEff.jar eff \
    -i vcf -o vcf \
    -s ${path}.eff.html \
    -ud 50 \
    -no-downstream \
    ecoli <(perl -pe "s/Chromosome/U00096/g" $path) \
