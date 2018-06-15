# install BLAST
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-*+-x64-linux.tar.gz
tar -zxvf ncbi-blast-*+-x64-linux.tar.gz
mv ncbi-blast*/bin/makeblastdb software
mv ncbi-blast*/bin/blastp software/
rm -r ncbi-blast-*

# install ORFfinder
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz
gzip -d ORFfinder.gz
mv ORFfinder software/

