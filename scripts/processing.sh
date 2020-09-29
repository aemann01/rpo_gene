#######################
#RPO GENE PRIMER DESIGN
#######################
# tblastn ncbi nt using protein reference sequences from UniProt
cd rpoA 
cat uniprot_rpoa_reviewed_yes_taxonomy_Bacteria.fasta | parallel --gnu --block 25k --recstart '>' --pipe \
	tblastn -evalue 1e-10 \
	-max_target_seqs 5 \
	-outfmt 6 \
	-db /data2/refdb/ncbi/nt \
	-query - > blast.out &
cd ../rpoB
cat uniprot-bacteria+rpob-filtered-reviewed_yes.fasta | parallel --gnu --block 50k --recstart '>' --pipe \
	tblastn -evalue 1e-10 \
	-max_target_seqs 5 \
	-outfmt 6 \
	-db /data2/refdb/ncbi/nt \
	-query - > blast.out &
cd ../rpoC
cat uniprot-rpoc_taxonomy_Bacteria_AND_reviewed_yes.fasta | parallel --gnu --block 50k --recstart '>' --pipe \
	tblastn -evalue 1e-10 \
	-max_target_seqs 5 \
	-outfmt 6 \
	-db /data2/refdb/ncbi/nt \
	-query - > blast.out &

# RPOC PROCESSING
cd rpoC/
# pull entries that had 80% identity to protein sequence, alignment of at least 1kb
awk -F"\t" '$3>=80.0 && $4>=1000' blast.out > good.hits
# get rid of rpoB sequences
grep -v "RPOB" good.hits > temp
mv temp good.hits
# format for bedtools to pull rpoC regions
awk '{if($10>$9){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' good.hits > good.hits.bed
# slice database
bedtools getfasta -fi /data2/refdb/ncbi/nt -bed good.hits.bed > good.hits.fa
# sort by length
vsearch --sortbylength good.hits.fa --output good.hits.sort.fa
# remove sequence duplicates
vsearch --derep_fulllength good.hits.sort.fa \
	--output good.hits.derep.fa \
	--uc good.hits.derep.uc \
	--sizeout
# cluster at 99% identity
vsearch --cluster_fast good.hits.derep.fa \
	--id 0.99 \
	--centroids good.hits.clust.fa \
	--profile good.hits.clust.profile.fa \
	--sizein \
	--sizeout \
	--threads 0
# add in chloroplast outgroup (zea mays rpoc1)
cat good.hits.clust.fa zea_mays.rpoc1.fa > temp
mv temp good.hits.clust.fa
# align sequences
mafft --auto good.hits.clust.fa > good.hits.align.fa
# change sequence headers so trimal doesn't cut them out and raxml doesn't freak out
sed -i 's/:/_/' good.hits.align.fa
sed -i 's/;/_/g' good.hits.align.fa
# trim
trimal -in good.hits.align.fa \
	-out good.hits.trimal.fa \
	-resoverlap 0.80 \
	-seqoverlap 80 \
	-gt 0.5 \
	-st 0.001 
# reference tree
rm *tre
raxmlHPC-PTHREADS -T 60 \
	-m GTRCAT \
	-c 25 \
	-e 0.001 \
	-p 31415 \
	-f a \
	-N 100 \
	-x 02938 \
	-n tre \
	-s good.hits.trimal.fa
# fix headers in alignment file
sed -i 's/;/ /g' good.hits.align.fa
rm *untrim.tre
# untrimmed reference tree
raxmlHPC-PTHREADS -T 30 \
	-m GTRCAT \
	-c 25 \
	-e 0.001 \
	-p 31415 \
	-f a \
	-N 100 \
	-x 02938 \
	-n untrim.tre \
	-s good.hits.align.fa
# get taxonomy from accession IDs
grep ">" good.hits.trimal.fa | sed 's/>//' | awk -F"_" '{print $1}' > batchez_query.ids
IFS=$'\n'; for line in $(cat batchez_query.ids); 
	do esearch -db nuccore -query $line | \
	elink -target taxonomy | \
	efetch -format native -mode xml | \
	grep -m 1 ScientificName | \
	sed 's/<ScientificName>//' | \
	sed 's/<\/ScientificName>//' | \
	sed 's/^[ \t]*//g' | \
	sed 's/ /_/g'; done > taxonomy.txt
# make annotations file
grep ">" good.hits.trimal.fa | sed 's/>//' | paste - taxonomy.txt > annotations.txt
sed -i 1i"accession\ttaxonomy" annotations.txt
# still have two distinct clades -- need to figure out if there are actually two copies, run panaroo
# download from NCBI Assembly database
# search query: ("Streptococcus"[Organism] OR streptococcus[All Fields]) AND ((latest[filter] OR "latest refseq"[filter]) AND "complete genome"[filter] AND (all[filter] NOT "derived from surveillance project"[filter] AND all[filter] NOT anomalous[filter]) AND "genbank has annotation"[Properties])

#######################
# PANAROO STREPTOCOCCUS
#######################
# install panaroo
# conda install -c conda-forge -c bioconda -c defaults panaroo
# install updated version mash 
# install prokka
# conda install -c conda-forge -c bioconda -c defaults prokka
# upgrade tbl2asn
# get dependency glib-2.14
# mkdir ~/src/glibc_install; cd ~/src/glibc_install 
# wget http://ftp.gnu.org/gnu/glibc/glibc-2.14.tar.gz
# tar zxvf glibc-2.14.tar.gz
# cd glibc-2.14
# mkdir build
# cd build
# ../configure --prefix=/opt/glibc-2.14
# make -j4
# cd /opt/glibc-2.14/etc
# sudo sh -c "echo '/usr/local/lib' >> ld.so.conf" 
# sudo sh -c "echo '/opt/lib' >> ld.so.conf"
# cd ~/src/glibc_install/glibc-2.14/build
# sudo make install
# export LD_LIBRARY_PATH="/opt/glibc-2.14/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
# now download the newest version of tbl2asn
# cd ~/src
# wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz
# gzip -d linux64.tbl2asn.gz
# chmod +x linux64.tbl2asn
# mv linux64.tbl2asn ~/src/miniconda3/bin/tbl2asn 
# # first need to run prokka on all strep genomes
# cd ~/rpo_gene/panaroo/ncbi-genomes-2020-09-14
# export LC_ALL=C
# # max 20 cores per job
# ls *fna | sed 's/.fna//' | while read line; do prokka --force --prefix $line $line.fna --cpus 20 2> prokka.log; done &

# fix mash install to 2.2
# wget https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar
# tar xvf mash-Linux64-v2.2.tar
# cd mash-Linux64-v2.2
# mv ./mash ~/src/miniconda3/bin/mash

# quality control checks 
# get reference database
# wget https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh
panaroo-qc -i /data2/vince/oral/ncbi-genbank/strep/genomic/str_smu/smu-id/fna/*/*gff \
	-o results/ \
	-t 15 \
	--graph_type all \
	--ref_db refseq.genomes.k21s1000.msh
# I'm not sure how to interpret the results, continue on with panaroo pipeline, come back to this
panaroo -i /data2/vince/oral/ncbi-genbank/strep/genomic/str_smu/smu-id/fna/*/*gff \
	-o results/ \
	--mode strict \
	-a pan \
	-t 40
# now want to pull the locus tags that are annotated by prokka for each rpo gene
# get locus tags for each
grep "rpoA" gene_presence_absence.csv | tr "," "\n" | grep "smu" > rpoA_locus.ids
grep "rpoB" gene_presence_absence.csv | tr "," "\n" | grep "smu" > rpoB_locus.ids
grep "rpoC" gene_presence_absence.csv | tr "," "\n" | grep "smu" > rpoC_locus.ids
# how many of each?
wc -l *ids
 # 172 rpoA_locus.ids
 # 172 rpoB_locus.ids
 # 172 rpoC_locus.ids
# looks like we have a representative for each across same number of genomes, only a single locus tag per each....
# now want to pull those locus tags and build a tree from them to see if the sequences are different
# pull fasta for each gene
cat rpoA_locus.ids | while read line; do grep -w $line gene_data.csv | awk -F"," '{print ">rpoA_"$1"|"$2"|"$4"\n"$6}' ; done > rpoA_locus.fa
cat rpoB_locus.ids | while read line; do grep -w $line gene_data.csv | awk -F"," '{print ">rpoB_"$1"|"$2"|"$4"\n"$6}' ; done > rpoB_locus.fa
cat rpoC_locus.ids | while read line; do grep -w $line gene_data.csv | awk -F"," '{print ">rpoC_"$1"|"$2"|"$4"\n"$6}' ; done > rpoC_locus.fa
# align
ls rpo*fa | sed 's/.fa//' | parallel --gnu 'mafft {}.fa > {}.align.fa'
# align as one tree for funzies
cat rpoA_locus.fa rpoB_locus.fa rpoC_locus.fa > all_rpo.fa
mafft all_rpo.fa > all_rpo.align.fa
# build trees
ls *align.fa | sed 's/.align.fa//' | parallel --gnu 'fasttree < {}.align.fa > {}.tre'

# talked to vince, want to try to cluster our rpoC reads against the full oral genome db at 50% ID
# on hillary
cd /home/allie/rpo_gene/rpoC
vsearch \
	--usearch_global /data2/vince/oral/ncbi-genbank/homd-oral2/genomic/no-plasmid-subsample-char/n4584-vsearch/n4584.fnn \
	--db good.hits.clust.fa \
	--otutabout n4584_rpoC_50PID_otutable.txt \
	--strand both \
	--id 0.50 \
	--threads 60 \
	--matched n4584_rpoC_50PID.fa \
	--uc n4584_rpoC_50PID.uc

# while this is running, do primer prospector analyze primers using short amplicons
analyze_primers.py -f ~/refdb/oral_bac_genomes/n4584.fnn -P primersF.txt
analyze_primers.py -f ~/refdb/oral_bac_genomes/n4584.fnn -P primersR.txt
get_amplicons_and_reads.py \
	-f ~/refdb/oral_bac_genomes/n4584.fnn \
	-i Univf_rpoB_deg_n4584_hits.txt:Univr_rpoB_deg_n4584_hits.txt \
	-t 1.5
# sort and dereplicate
vsearch --sortbylength Univf_Univr_amplicons.fasta \
	--output Univf_Univr_amplicons.sort.fa
vsearch --derep_fulllength Univf_Univr_amplicons.sort.fa \
	--outpout Univf_Univr_amplicons.derep.fa \
	--sizeout
# align and build into tree
mafft Univf_Univr_amplicons.derep.fa > Univf_Univr_amplicons.align.fa
# format headers for raxml
sed -i 's/;/ /g' Univf_Univr_amplicons.align.fa
# build tree
raxmlHPC-PTHREADS -T 8 \
	-m GTRCAT \
	-c 25 \
	-e 0.001 \
	-p 31415 \
	-f a \
	-N 100 \
	-x 02938 \
	-n tre \
	-s Univf_Univr_amplicons.align.fa
# while this is running, get list of hits, find out which taxa are not represented
grep ">" Univf_Univr_amplicons.fasta | sed 's/>//' | awk -F"|" '{print $1}' > primer_hits.txt
diff temp primer_hits.txt > primer_diff.txt
# there's quite a bit of drop out, try to optimze primers?
optimize_primers.py -i Univf_rpoB_deg_n4584_hits.txt
optimize_primers.py -i Univr_rpoB_deg_n4584_hits.txt
# analyze optimized primers
analyze_primers.py -f ~/refdb/oral_bac_genomes/n4584.fnn -P primersF_opt.txt &
analyze_primers.py -f ~/refdb/oral_bac_genomes/n4584.fnn -P primersR_opt.txt &
# get amplicons from optimized primers
get_amplicons_and_reads.py \
	-f ~/refdb/oral_bac_genomes/n4584.fnn \
	-i rpoBf_opt_n4584_hits.txt:rpoBr_opt_n4584_hits.txt \
	-t 1.5
# sort and dereplicate
vsearch --sortbylength rpoBf_rpoBr_amplicons.fasta \
	--output rpoBf_rpoBr_amplicons.sort.fa
vsearch --derep_fulllength rpoBf_rpoBr_amplicons.sort.fa \
	--output rpoBf_rpoBr_amplicons.derep.fa \
	--sizeout
# align and build into tree
mafft rpoBf_rpoBr_amplicons.derep.fa > rpoBf_rpoBr_amplicons.align.fa
# format headers for raxml
sed -i 's/;/ /g' rpoBf_rpoBr_amplicons.align.fa
# build tree
raxmlHPC-PTHREADS -T 8 \
	-m GTRCAT \
	-c 25 \
	-e 0.001 \
	-p 31415 \
	-f a \
	-N 100 \
	-x 02938 \
	-n opt.tre \
	-s rpoBf_rpoBr_amplicons.align.fa
# what species/strains do we gain with optimized primers?
grep ">" Univf_Univr_amplicons.fasta | sed 's/>//' | awk -F"|" '{print $1}' | sed 's/_/./' | awk -F"_" '{print $1}' | sort | uniq  > temp1
grep ">" rpoBf_rpoBr_amplicons.fasta | sed 's/>//' | awk -F"|" '{print $1}' | sed 's/_/./' | awk -F"_" '{print $1}' | sort | uniq  > temp2
diff temp1 temp2 > optimized_primer_added_species.txt
# get annotations for both trees
grep ">" rpoBf_rpoBr_amplicons.fasta | sed 's/>//' > temp1
cat temp1 | awk -F"|" '{print $1}' | while read line; do grep -w $line ~/refdb/oral_bac_genomes/ID-n4584-1.txt ; done > temp2
paste temp1 temp2 > annotations.txt 
# add headers in nano
# what taxa are we losing completely with both?
cat annotations.txt | sort > temp1
cat ~/refdb/oral_bac_genomes/ID-n4584-1.txt | sort > temp2
diff temp1 temp2 | grep ">" | sed 's/> //' > lost_taxa.txt
# try clustering really lax threshold amplicons (3.0) to get as many strains/species as possible
get_amplicons_and_reads.py \
	-f ~/refdb/oral_bac_genomes/n4584.fnn \
	-i rpoBf_opt_n4584_hits.txt:rpoBr_opt_n4584_hits.txt \
	-t 3.0 \
	-o amplicons_for_vsearch
# sort and dereplicate
cd amplicons_for_vsearch
vsearch --sortbylength rpoBf_rpoBr_amplicons.fasta --output rpoBf_rpoBr_amplicons.sort.fa
# what threshold to trim really long sequences?
Rscript plot_GC_length.R rpoBf_rpoBr_amplicons.derep.fa
# dereplicate
vsearch --derep_fulllength rpoBf_rpoBr_amplicons.sort.fa \
	--output rpoBf_rpoBr_amplicons.derep.fa \
	--sizeout \
	--minseqlength 300 \
	--maxseqlength 800
# align and build tree
mafft rpoBf_rpoBr_amplicons.derep.fa > rpoBf_rpoBr_amplicons.align.fa 
# format headers for raxml
sed -i 's/;/ /g' rpoBf_rpoBr_amplicons.align.fa
# build tree
raxmlHPC-PTHREADS -T 8 \
	-m GTRCAT \
	-c 25 \
	-e 0.001 \
	-p 31415 \
	-f a \
	-N 100 \
	-x 02938 \
	-n trim.tre \
	-s rpoBf_rpoBr_amplicons.align.fa
# are the lost taxa actually lost (i.e., no representatives of that species in the tree?)
awk '{print $2}' lost_taxa.txt | sed 's/_/./' | awk -F"_" '{print $1}' | sort | uniq | sed 's/\./_/' > lost_species_temp
awk -F"\t" '{print $3}' annotations.txt | sed 's/_/./' | awk -F"_" '{print $1}' | sort | uniq | sed 's/\./_/' > found_taxa_temp
diff lost_species_temp found_taxa_temp
diff lost_species_temp found_taxa_temp | grep "<" | sed 's/< //' > no_representatives.txt
#####################
# PLACEMENT TREE TEST
#####################

# what if we tried something like the eukref pipeline (which is unfortunately currently down)?
# can use those reads that formed a nice tree and blasted them against nt at low seq similarity
# then place them into our nice tree and see if we can't expand the references?
# working off our 99% clustered, untrimmed reads as query
cat good.hits.clust.fa | parallel --gnu --block 50k --recstart '>' --pipe \
	blastn -evalue 1e-10 \
	-max_target_seqs 5 \
	-perc_identity 0.60 \
	-outfmt 6 \
	-db /data2/refdb/ncbi_nt_blast/nt \
	-query - > good.hits.clust.blast.out 
# pull regions of interest
awk '{if($10>$9){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' \
	good.hits.clust.blast.out > for.placement.bed
# slice database
bedtools getfasta -fi /data2/refdb/ncbi_nt_blast/nt -bed for.placement.bed >  for.placement.hits.fa
# sort by length
vsearch --sortbylength for.placement.hits.fa --output for.placement.hits.sort.fa
# remove sequence duplicates
vsearch --derep_fulllength for.placement.hits.sort.fa \
	--output for.placement.hits.derep.fa \
	--uc for.placement.hits.derep.uc \
	--sizeout
# check for chimeric sequences in newly downloaded sequences
vsearch --uchime_denovo for.placement.hits.derep.fa --sizein --sizeout --nonchimeras for.placement.nochim.fa
# then add to untrimmed alignment using sina
sina -i good.hits.align.fa --prealigned -o good.hits.align.arb
sina -i for.placement.nochim.fa -r good.hits.align.arb -o for.placement.sina.fa -p 2
# make placement tree using good tree
# fix headers in ref to match tree
sed -i 's/ /_/g' good.hits.align.fa
cat for.placement.sina.fa good.hits.align.fa > queryalign_plus_refalign.fa
sed -i 's/:/_/g' queryalign_plus_refalign.fa
sed -i 's/;/_/g' queryalign_plus_refalign.fa
sed -i 's/,//g' queryalign_plus_refalign.fa
# remove sequences with same name and seq
cat queryalign_plus_refalign.fa |\
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |\
sort -t $'\t' -k1,1 -u |\
tr "\t" "\n" > fix.align.fa
# make tree with reduced file generated by raxml
rm *const.tre
~/raxmlHPC-AVX-v8/raxml -f a -N 100 -G 0.2 -m GTRCAT -n const.tre -s fix.align.fa -g ../RAxML_bipartitions.tre -T 4 -x 25734 -p 25793














