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

# RPOA PROCESSING

# RPOB PROCESSING

# RPOC PROCESSING












# now need to pull those entries that had 90% identity to the protein sequence, and an alignment of at least 1kb
awk -F"\t" '$3>=90.0 && $12>=2000' homd_blast.out > homd_blast.90pid2k.out
awk -F"\t" '$3>=90.0 && $12>=2000' n4k_blast.out > n4k_blast.90pid2k.out
# get list of fasta that will need to be rev comp none in 4k
awk '{if($9>$10){print($2)}}' homd_blast.90pid2k.out > homd_blast.revcomp.ids
# list of ids that are in correct orientation
awk '{if($10>$9){print($2)}}' homd_blast.90pid2k.out > homd_blast.normal.ids
# format for bedtools
awk '{if($10>$9){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' homd_blast.90pid2k.out > homd_blast.90pid2k.bed
awk '{if($10>$9){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' n4k_blast.90pid2k.out > n4k_blast.90pid2k.bed
# slice genomes
bedtools getfasta -fi ~/refDB/homd_genomes/ALL_genomes.fna \
	-bed homd_blast.90pid2k.bed > homd_blast.90pid2k.fa
bedtools getfasta -fi ~/refDB/vince_oral_genomes/n4584.fnn \
	-bed n4k_blast.90pid2k.bed > n4k_blast.90pid2k.fa
# fix headers for seqtk
sed -i 's/:/\t/' homd_blast.90pid2k.fa
# reverse complement
seqtk subseq homd_blast.90pid2k.fa homd_blast.revcomp.ids | seqtk seq -r > homd_blast.revcomp.fa
# normal complement
seqtk subseq homd_blast.90pid2k.fa homd_blast.normal.ids > homd_blast.normal.fa
# clean up
rm homd_blast.90pid2k.fa
cat homd_blast.normal.fa homd_blast.revcomp.fa > homd_blast.90pid2k.fa
rm homd_blast.normal.fa homd_blast.revcomp.fa
rm *ids
# cat together both datasets
cat n4k_blast.90pid2k.fa homd_blast.90pid2k.fa > all.fa
# remove sequence duplicates
vsearch --derep_fulllength all.fa --output all.derep.fa --sizeout
# sort by length
vsearch --sortbylength all.derep.fa --output all.sort.fa
# cluster at 99%
vsearch --cluster_fast all.sort.fa --id 0.99 --centroids all.clust.fa --sizein --sizeout
# add in chloroplast outgroup
cat all.clust.fa arabidopsis_rpob.fa > temp
mv temp all.clust.fa
# align sequences
mafft --auto all.clust.fa > all.align.fa
# trim 
trimal -in all.align.fa \
	-out all.trimal.fa \
	-resoverlap 0.90 \
	-seqoverlap 90 \
	-gt 0.5 \
	-st 0.001 \
	-scc \
	-sident \
	1> all.trimal.out
# fix headers for raxml
sed -i 's/;/ /' all.trimal.fa
# manually removed poorly aligned duplicates in aliview, saved as clean copy
# Sequence names of taxon 102 and 243 are identical, they are both called SEQF1300_CP000076.1
# Sequence names of taxon 106 and 244 are identical, they are both called SEQF2332_CP003041.1
# Sequence names of taxon 388 and 692 are identical, they are both called SEQF2672_CP006936.2
# Sequence names of taxon 391 and 406 are identical, they are both called SEQF1130_CP001658.1
# Sequence names of taxon 391 and 519 are identical, they are both called SEQF1130_CP001658.1
# Sequence names of taxon 406 and 519 are identical, they are both called SEQF1130_CP001658.1
# Sequence names of taxon 419 and 500 are identical, they are both called SEQF2206_AENR01000042.1
# Sequence names of taxon 422 and 504 are identical, they are both called SEQF2204_GL545268.1
# ERROR: Found 8 taxa that had equal names in the alignment, exiting...
# build a reference tree
~/raxmlHPC-AVX-v8/raxml -T 4 \
	-m GTRCAT \
	-c 25 \
	-e 0.001 \
	-p 31415 \
	-f a \
	-N 100 \
	-x 02938 \
	-n all.tre \
	-s all.clean.fa
# get taxonomy info for tree annotations
grep ">SEQ" all.clean.fa | sed 's/>//' | \
	awk -F"_" '{print $1}' | \
	while read line; do grep -w $line ~/refDB/homd_genomes/SEQID_info.txt ; \
	done > temp
grep ">" all.clean.fa | grep -v "SEQ" | \
	sed 's/>//' | awk -F"|" '{print $1}' | \
	sed 's/NC.*//' | while read line; do \
	grep -w $line ~/refDB/vince_oral_genomes/ID-n4584-1.txt ; done > temp2
# remove wonky sequences from alignment
# SEQF1302_AM181176.4
# SEQF1303_CP000304.1
# SEQF1645_FM211192.1
# SEQF1836_GG696773.1
# SEQF2196_CM001025.1
# SEQF2333_AGSL01000031.1
# SEQF2352_AJMR01000041.1
# SEQF2381_CM001513.1
# SEQF2382_CP003677.1
# SEQF2383_AFEH01000043.1
# SEQF2423_HG916826.1
# SEQF2513_AJXE01000041.1
# SEQF2521_CP003725.1
# esc_col_97|ECTT12B_5153
# myc_tub_93|Z534_00164
# myc_tub_99|X376_03205
# pse_flu_11|pse_flu_11_04670
# pse_flu_12|pse_flu_12_32765
# pse_flu_14|pse_flu_14_05555
# pse_flu_16|pse_flu_16_08460
# pse_flu_1|SRM1_05185
# pse_flu_22|pse_flu_22_28050
# pse_flu_23|pse_flu_23_06855
# pse_flu_27|pse_flu_27_02895
# pse_flu_29|NL64_10740
# pse_flu_34|PFLU4_43860
# pse_flu_37|pse_flu_37_20515
# pse_flu_38|A1D17_17590
# pse_flu_39|AO356_14615
# pse_flu_40|TK06_24250
# pse_flu_41|AO353_07030
# pse_flu_43|pse_flu_43_06290
# pse_flu_45|RY26_30130
# pse_flu_47|VD17_07105
# pse_flu_49|AO066_26680
# pse_flu_55|RU10_25275
# pse_flu_57|UG46_24200
# pse_flu_58|pse_flu_58_14350
# pse_flu_59|PF1751_v1c49570
# pse_flu_5|pse_flu_5_12300
# pse_flu_64|AWV77_16650
# pse_flu_66|AN403_604
# pse_flu_67|NX10_25670
# pse_flu_68|QS95_25800
# pse_flu_70|RL74_14730
# pse_flu_75|PSF113_5301
# pse_flu_79|B723_01820
# pse_flu_84|Pfl01_5085
# pse_flu_85|pse_flu_85_24490
# pse_flu_87|PflQ2_4978
# pse_flu_9|pse_flu_9_09615
# pse_pse_2|AU05_14060
# sta_aur_99|T933_01327
# myc_gen_2|MG_341
# myc_pne_55|C680_02920
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
