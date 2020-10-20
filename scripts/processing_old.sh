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

# talked to 












# two distinct clades? pull just streptococcus and realign to see what's up
mkdir strep_only
grep "Streptococcus" annotations.txt | awk '{print $1}' > strep_only/strep.ids
cd strep_only
seqtk subseq ../good.hits.trimal.fa strep.ids > strep.trimal.fa
# realign everything
mafft strep.trimal.fa > strep.realign.fa
# very different, same approx pid and length, all but a couple are annotated as rpoc
# split into two clades, run primer prospector on both
sed -i 's/:/ /' good.hits.clust.fa
seqtk subseq good.hits.clust.fa clade-a.ids > clade-a.fa
seqtk subseq good.hits.clust.fa clade-b.ids > clade-b.fa
# which taxa have multiple hits per record?
awk -F"\t" '{print $2, $9, $10}' good.hits | sort | uniq -c | \
	sed 's/^[ \t]*//g' | grep "^1" | awk -F" " '{print $2}' | \
	sort | uniq -c | sed 's/^[ \t]*//g' | grep -v "^1" | \
	awk '{print $2}' | while read line; 
	do grep $line good.hits ; 
	done > multiple_hits_same_record.txt
# which genera are found in either clade a or b?
cat clade-a.ids | while read line; do grep $line annotations.txt ; done > clade-a.tax
cat clade-b.ids | while read line; do grep $line annotations.txt ; done > clade-b.tax
awk '{print $2}' clade-a.tax | awk -F"_" '{print $1}' | sort | uniq > tempa
awk '{print $2}' clade-b.tax | awk -F"_" '{print $1}' | sort | uniq > tempb
diff tempa tempb > diff_a-b.txt
grep -E ">|<" diff_a-b.txt | sed 's/>//' | sed 's/<//' | sed 's/^[ \t]*//g' > temp
mv temp diff_a-b.txt
# what are the annotations for those with multiple, non-overlapping hits?








IFS=$'\n'; for line in $(cat batchez_query.ids); 
	do esearch -db nuccore -query $line | \
	elink -target taxonomy | \
	efetch -format native -mode xml | \
	grep -m 1 ScientificName | \
	sed 's/<ScientificName>//' | \
	sed 's/<\/ScientificName>//' | \
	sed 's/^[ \t]*//g' | \
	sed 's/ /_/g'; done > taxonomy.txt



# realign
mafft clade-a.fa > clade-a.align.fa
mafft clade-b.fa > clade-b.align.fa














# long branch sequences removed from alignment and tree -- final file: all.final.fa
# now run primer prospector on this alignment
generate_primers_denovo.py -i all.final.fa \
	-o rpob_primers.txt \
	-a ../rpoB_alignment.fasta \
	-p 0.98
# sort primers
sort_denovo_primers.py -i rpob_primers.txt -a 900:5000
# get scores for primers
# initialize primer queries
# grep -v "Est" amplicon_len_pairs.txt | sort | uniq | sed 's/^/analyze_primers.py -f all.fa -p /' | sed 's/\t/ -s /' | sed 's/$/ -o primers/'
analyze_primers.py -f all.fa -p 1285r -s TCGGWAGTTYTCCCYTTGGG -o primers &
analyze_primers.py -f all.fa -p 1320r -s AGTTCACTCGMYTTTCCCAA -o primers &
analyze_primers.py -f all.fa -p 1321r -s GAGTTCACTCGMYTTTCCCA -o primers &
analyze_primers.py -f all.fa -p 1322r -s AGAGTTCACTCGMYTTTCCC -o primers &
analyze_primers.py -f all.fa -p 1323r -s CAGAGTTCACTCGMYTTTCC -o primers &
analyze_primers.py -f all.fa -p 1523r -s CRTGCGCCGCCATYTTTCCC -o primers &
analyze_primers.py -f all.fa -p 1524r -s CCRTGCGCCGCCATYTTTCC -o primers &
analyze_primers.py -f all.fa -p 1570r -s CGGKCCTYKWGWAGGCATTC -o primers &
analyze_primers.py -f all.fa -p 1571r -s ACGGKCCTYKWGWAGGCATT -o primers &
analyze_primers.py -f all.fa -p 1572r -s CACGGKCCTYKWGWAGGCAT -o primers &
analyze_primers.py -f all.fa -p 1605r -s CCRAYRTTCATCGAGGACCC -o primers &


analyze_primers.py -f all.fa -p 1606r -s GCCRAYRTTCATCGAGGACC -o primers &
analyze_primers.py -f all.fa -p 1612r -s TCAACTGCCRAYRTTCATCG -o primers &
analyze_primers.py -f all.fa -p 1613r -s YTCAACTGCCRAYRTTCATC -o primers &
analyze_primers.py -f all.fa -p 1614r -s GYTCAACTGCCRAYRTTCAT -o primers &
analyze_primers.py -f all.fa -p 1690r -s RKCYATCTCTKGCCCTCAAA -o primers &
analyze_primers.py -f all.fa -p 1755r -s AACGGCRTGATYTTTCTCAC -o primers &
analyze_primers.py -f all.fa -p 1756r -s GAACGGCRTGATYTTTCTCA -o primers &
analyze_primers.py -f all.fa -p 1757r -s TGAACGGCRTGATYTTTCTC -o primers &
analyze_primers.py -f all.fa -p 1758r -s CTGAACGGCRTGATYTTTCT -o primers &
analyze_primers.py -f all.fa -p 1786r -s CYTTACCCCAGGYTGTGGTA -o primers &
analyze_primers.py -f all.fa -p 1787r -s GCYTTACCCCAGGYTGTGGT -o primers &
analyze_primers.py -f all.fa -p 1788r -s GGCYTTACCCCAGGYTGTGG -o primers &
analyze_primers.py -f all.fa -p 1794r -s ACCAATGGCYTTACCCCAGG -o primers &
analyze_primers.py -f all.fa -p 1795r -s CACCAATGGCYTTACCCCAG -o primers &
analyze_primers.py -f all.fa -p 1796r -s CCACCAATGGCYTTACCCCA -o primers &
analyze_primers.py -f all.fa -p 1822r -s TCAGCCCAACTCCATTCCCA -o primers &
analyze_primers.py -f all.fa -p 1823r -s TTCAGCCCAACTCCATTCCC -o primers &
analyze_primers.py -f all.fa -p 1824r -s CTTCAGCCCAACTCCATTCC -o primers &
analyze_primers.py -f all.fa -p 1825r -s GCTTCAGCCCAACTCCATTC -o primers &
analyze_primers.py -f all.fa -p 1826r -s AGCTTCAGCCCAACTCCATT -o primers &
analyze_primers.py -f all.fa -p 1827r -s TAGCTTCAGCCCAACTCCAT -o primers &
analyze_primers.py -f all.fa -p 1828r -s CTAGCTTCAGCCCAACTCCA -o primers &
analyze_primers.py -f all.fa -p 1829r -s CCTAGCTTCAGCCCAACTCC -o primers &
analyze_primers.py -f all.fa -p 1830r -s CCCTAGCTTCAGCCCAACTC -o primers &
analyze_primers.py -f all.fa -p 1831r -s GCCCTAGCTTCAGCCCAACT -o primers &
analyze_primers.py -f all.fa -p 1832r -s MGCCCTAGCTTCAGCCCAAC -o primers &
analyze_primers.py -f all.fa -p 1833r -s GMGCCCTAGCTTCAGCCCAA -o primers &
analyze_primers.py -f all.fa -p 1834r -s TGMGCCCTAGCTTCAGCCCA -o primers &
analyze_primers.py -f all.fa -p 1835r -s RTGMGCCCTAGCTTCAGCCC -o primers &
analyze_primers.py -f all.fa -p 269f -s TATYCCTWACCGGGTCTGGT -o primers &
analyze_primers.py -f all.fa -p 525f -s GGAAYCGTCGRTCGGTGGGA -o primers &
analyze_primers.py -f all.fa -p 525f -s TGGAAYCGTCGRTCGGTGGG -o primers &
analyze_primers.py -f all.fa -p 525f -s YTGGAAYCGTCGRTCGGTGG -o primers &
analyze_primers.py -f all.fa -p 591f -s GWWCARTTCCARTTYATGGA -o primers &
analyze_primers.py -f all.fa -p 632f -s CCACAARCGCGYTTCGCTGG -o primers &
analyze_primers.py -f all.fa -p 695f -s CAYTAYGGCGRTTGCCATGA -o primers &
analyze_primers.py -f all.fa -p 868f -s RTGTCGTGCCCKATYCCTTC -o primers &
analyze_primers.py -f all.fa -p 869f -s TGTCGTGCCCKATYCCTTCT -o primers &
analyze_primers.py -f all.fa -p 870f -s GTCGTGCCCKATYCCTTCTG -o primers &
analyze_primers.py -f all.fa -p 894f -s YGAYGAYKCAACCGYGCTAT -o primers &
analyze_primers.py -f all.fa -p 895f -s GAYGAYKCAACCGYGCTATG -o primers &




# get amplicons







# place these into the tree with tagged primer set to make sure they're actually amplifying 
