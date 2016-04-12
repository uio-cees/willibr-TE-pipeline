#!/bin/sh
#SBATCH --job-name=pipeline
#SBATCH --account=nn9244k
#SBATCH --time=96:0:0
#SBATCH --mem-per-cpu=3000M 
#SBATCH --cpus-per-task=16
#SBATCH --output=repeats_slurm_out/slurm-%j-%N.out

if [ -n "$SLURM_JOB_ID" ]; then
    # running in a slurm job
    source /cluster/bin/jobsetup
fi

module purge
 
TMP=$LOCALTMP
TMPDIR=
 
mkdir -p repeats
 
species=${1%}
origdir=$PWD
 
mkdir -p repeats/${species}
echo ${species}
 
mkdir ${TMPDIR}/${species}_repeats
cd ${TMPDIR}/${species}_repeats
 
cp $origdir/repeats_worker_script.slurm ${TMPDIR}/${species}_repeats
cp $origdir/repeats_master_pipeline.slurm ${TMPDIR}/${species}_repeats
 
cp $origdir/$1 ${TMPDIR}/${species}_repeats

#------------------------------Pointing to paths-----------------------------------------#

#	The repository contains databases and scripts used during the pipeline.

export REPO=

#---Content of the repository:

#	Databases:

#	eukaryotic-tRNAs.fa #http://lowelab.ucsc.edu/GtRNAdb/download.html
#	uniprot_sprot17nov2015.fasta #http://www.ebi.ac.uk/uniprot/database/download.html
#	RepeatPeps.lib #Comes with RepeatMasker distributions
#	TEseeker.DNAtransposons.fasta #http://www3.nd.edu/~teseeker/download.html

#	Python scripts (tested in python2/2.7.9):

#	change_headers_to_seqN.py
#	reprint.tPSI.lib.py
#	reprint.ltrharvest.lib.py
#	reprint.ltrharvest.lib.new.py
#	reprint.filtered.lib.py
#	sequence_cleaner.py
#	N50.py
#	reprint.genome.py

#	LTRharvest/LTRdigest filtering scripts (courtesy of Ning Jiang): Will hopefully be published soon.

#	Misc:
#	assemblathon_stats.pl : http://korflab.ucdavis.edu/datasets/Assemblathon/Assemblathon2/Basic_metrics/assemblathon_stats.pl
#	filter_protein_match.lua: https://github.com/satta/ltrsift/blob/master/filters/filter_protein_match.lua 

#	The only parameter given the pipeline is a FASTA file of the genome.

export GENOME=$1

#	Point to a modified version of TransposonPSI with dual core usage.

export tPSI= # path to transposonPSI.pl

#	Point to CD-HIT -- Can probably change this with using the module system

export CDHIT= #path to cd-hit-est

#-----------------------------Loading modules--------------------------------------------#

module load genometools/1.5.7 
module load repeatmodeler/1.0.8
module load blast/2.2.26
module load perlmodules/5.10_2
module load hmmer/3.1b1
module load muscle/3.8.31
module load python2/2.7.9
module load blast+/2.2.29
module load repeatmasker/4.0.6
module load bedtools/2.17.0
module load cd-hit/4.6.4

#------------------------------Genome index table----------------------------------------#
#	NB! Genome FASTA headers must not contain any special characters (only letters and/or 
#	numbers).
 
#	Renaming headers (to SEQ0, SEQ1, SEQ2 etc.)
python2 $REPO/reprint.genome.py -i $GENOME > $GENOME.seq
 
#	Printing headers of the original genome file
grep ">" $GENOME > halfoftable1
 
#	Printing headers of the renamed genome file
grep ">" $GENOME.seq > halfoftable2
 
#	Creating a file with two columns containing original header name and new header name
paste halfoftable1 halfoftable2 > $GENOME.table
# 
#	Pointing the main variable, $GENOME to the version with renamed headers.

GENOME=$GENOME.seq

#---------------------------------Harvesting repeats-------------------------------------#

#	Starting LTRharvest/LTRdigest, RepeatModeler and TransposonPSI

arrayrun 1-3 repeats_worker_script.slurm $REPO $GENOME $tPSI $CDHIT 

wait


#-------------------------------tPSI FASTA extraction------------------------------------#

# In order to extract sequences from the TransposonPSI file, files had to be altered
cut -f 1,4,5,9 $GENOME.TPSI.allHits.chains.bestPerLocus.gff3 | sed 's/ID=.*Target.//g' | \
sed 's/;.*//g' | column -t | awk '{print $4}' \
 	> tPSI.classes
 
 paste tPSI.classes $GENOME.TPSI.allHits.chains.bestPerLocus.gff3 | column -s $'\t' -t \
 	> tmp4
 
awk '{print $2,$3,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' tmp4 \
 	> tmp5

sed 's/ /\t/g' tmp5 \
 	> tmp6
 
#	The actual extraction done using bedtools (getfasta)
 
bedtools getfasta -name -fi $GENOME -bed tmp6 -fo $GENOME.tPSI.fasta
 
#	The TransposonPSI library is classified according to protein homology
python2 $REPO/reprint.tPSI.lib.py -i $GENOME.tPSI.fasta | fold -w 60 \
 	> $GENOME.tPSI.classified.fasta
 
#################################
#	Removal of excess files  
# 
rm tmp* tPSI.classes *tPSI.fasta *TPSI.allHits.chains.gff3  tPSI.log
rm *TPSI.allHits.chains.bestPerLocus *TPSI.allHits.chains *TPSI.allHits
################################
 
#---------------------------------Clustering---------------------------------------------#
# Clustering of tPSI elements are done using CD-HIT-EST with 80% cutoff.
 
$CDHIT -i $GENOME.tPSI.classified.fasta -o $GENOME.tPSI.classified.fasta.nr.80.n5 \
-c 0.8 -T 16 -n 5

#--------------------------------Merging of libraries------------------------------------#

cd RM*
sed 's/ .*//g' consensi.fa.classified > consensi.fa.classified.stripped
cd ..
cat $GENOME.retrotransposons.TRIMs.lib.classified.reclassified \
RM*/consensi.fa.classified.stripped \
$GENOME.tPSI.classified.fasta.nr.80.n5 | fold -w 60 > $GENOME.repeats.fasta

#---------------------------Removing false positives-------------------------------------#

#	Removing false transposons by checking against SwissProt-UniProt. TIME: ~1 hour 

blastx -query $GENOME.repeats.fasta \
-db $REPO/uniprot_sprot17nov2015.fasta -num_threads 10 -evalue 1e-5 -max_target_seqs 50 \
-outfmt 6 -out $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50

blastx -query  $GENOME.repeats.fasta \
-db $REPO/RepeatPeps.lib -num_threads 10 -evalue 1e-5 -max_target_seqs 50 \
-outfmt 6 -out $GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50

sort -k1,1 -k12,12nr -k11,11n \
$GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50 | sort -u -k1,1 \
--merge > $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.highest_scoring

sort -k1,1 -k12,12nr -k11,11n \
$GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50 | sort -u -k1,1 \
--merge > $GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.highest_scoring

cut -f 1 $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.highest_scoring \
> $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.first_column

cut -f 1 $GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.highest_scoring \
> $GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.first_column

echo "Num entries only hit to RepeatPeps" 
comm -13 $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.first_column \
$GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.first_column \
> $GENOME.repeats.only_repeatpeps
	
cat $GENOME.repeats.only_repeatpeps |wc -l

echo "Num entries only hit to UniProt"

comm -23 $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.first_column \
$GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.first_column \
> $GENOME.repeats.only_uniprot
	
cat $GENOME.repeats.only_uniprot |wc -l

echo "Num entries hit to both"

comm -12 $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.first_column \
$GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.first_column | wc -l

python2 $REPO/reprint.filtered.lib.py -i $GENOME.repeats.fasta -l \
$GENOME.repeats.only_uniprot | fold -w 60 \
	> $GENOME.repeats.srt.nr.no_uniprot

#--------------------------BLASTX for classification-------------------------------------#

blastx -query $GENOME.repeats.srt.nr.no_uniprot \
-db $REPO/RepeatPeps.lib -num_threads 10 -evalue 1e-5 -max_target_seqs 1 -outfmt \
'6 bitscore length pident qseqid sacc' | awk '{print $5," ",$4," ",$1," ",$2," ",$3}' \
	> $GENOME.ClassificationScheme

python2 $REPO/classifier.py -i \
$GENOME.repeats.srt.nr.no_uniprot -c $GENOME.ClassificationScheme \
	> $GENOME.classified_repeats.fasta

awk '{print $2}' $GENOME.classified_repeats.fasta | grep "#" > list_to_fastagrep

fastagrep -v -f list_to_fastagrep $GENOME.repeats.srt.nr.no_uniprot \
	> $GENOME.non_classified_repeats.fasta

cat $GENOME.non_classified_repeats.fasta $GENOME.classified_repeats.fasta \
	> $GENOME.denovolib

#-----------------------------HMMR classification of "Unknowns"--------------------------#

#	Retrieving the FASTA sequences of unknown repeats from RepeatModeler

fastagrep -p "#Unknown" $GENOME.denovolib | fastagrep -p ">rnd" \
	> $GENOME.repmod.unknowns.lib

#	Running nhmmscan to search for homology in Dfam.

gt splitfasta -numfiles 10 $GENOME.repmod.unknowns.lib

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.1 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.1 1> $GENOME.DFAM.out.1 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.2 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.2 1> $GENOME.DFAM.out.2 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.3 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.3 1> $GENOME.DFAM.out.3 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.4 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.4 1> $GENOME.DFAM.out.4 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.5 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.5 1> $GENOME.DFAM.out.5 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.6 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.6 1> $GENOME.DFAM.out.6 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.7 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.7 1> $GENOME.DFAM.out.7 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.8 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.8 1> $GENOME.DFAM.out.8 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.9 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.9 1> $GENOME.DFAM.out.9 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.10 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.10 1> $GENOME.DFAM.out.10

cat $GENOME.repmod.unknowns.dfam_results* > $GENOME.repmod.unknowns.all.dfam_results

#	Running a set of python scripts to produce annotated versions of the libraries

python2 $REPO/dfammer.py -i $GENOME.repmod.unknowns.all.dfam_results \
	> $GENOME.repmod.unknowns.dfam_tophits

python2 $REPO/dfammer2.py -i $GENOME.repmod.unknowns.dfam_tophits \
	> $GENOME.repmod.unknowns.annotated.list

awk '{print $2}' $GENOME.repmod.unknowns.annotated.list \
	> $GENOME.repmod.unknowns.annotated.list.to_fastagrep

python2 $REPO/classifier.py -i $GENOME.repmod.unknowns.lib -c \
$GENOME.repmod.unknowns.annotated.list \
	> $GENOME.repmod.unknowns.dfam_classified.lib

#	Extracting all sequences except the unknowns not classified by HMM model homology

fastagrep -v -f $GENOME.repmod.unknowns.annotated.list.to_fastagrep $GENOME.denovolib \
	> $GENOME.denovolib.no_classifiable_unknowns

#	Removing bottom empty line

sed -i '/^$/d' $GENOME.repmod.unknowns.dfam_classified.lib 

cat $GENOME.repmod.unknowns.dfam_classified.lib \
$GENOME.denovolib.no_classifiable_unknowns \
	> $GENOME.denovolib.classified

#----------------------------------RepeatMasker------------------------------------------#

#	Creating a total repeat library

cat $GENOME.denovolib.classified $REPO/RepBase_for_RepeatMasker.20150807.eukaryotes \
> $GENOME.total.repeat.library

#	Making folders for each library

mkdir $GENOME.repmask.denovo
mkdir $GENOME.repmask.total
mkdir $GENOME.repmask.repbase

#	RepeatMasking and creating summaries

RepeatMasker -lib $GENOME.denovolib.classified -gff -a -s -pa 16 -dir \
$GENOME.repmask.denovo/ $GENOME

/cluster/software/VERSIONS/repeatmasker-4.0.6/util/buildSummary.pl -species teleostei \
$GENOME.repmask.denovo/*.out \
	> $GENOME.repmask.denovo.summary

RepeatMasker -lib $GENOME.total.repeat.library -gff -a -s -pa 16 -dir \
$GENOME.repmask.total/ $GENOME 

/cluster/software/VERSIONS/repeatmasker-4.0.6/util/buildSummary.pl -species teleostei \
$GENOME.repmask.total/*.out \
	> $GENOME.repmask.total.summary

RepeatMasker -species eukaryota -gff -a -s -pa 16 -dir \
$GENOME.repmask.repbase/ $GENOME 

/cluster/software/VERSIONS/repeatmasker-4.0.6/util/buildSummary.pl -species teleostei \
$GENOME.repmask.repbase/*.out \
	> $GENOME.repmask.repbase.summary

#--------------------------------RepeatLandscape-----------------------------------------#

#	The createRepeatLandscape.pl script requires genome size estimates to create pie 
#	charts. 

#	Getting genome size

assemblathon_stats.pl -csv $GENOME | grep "Total size of scaffolds" | awk '{print $5}' \
> genomesize

sh $REPO/landscaping.sh $GENOME

#---------------------------------Result parsing-----------------------------------------#

cd $GENOME.repmask.total/
sh $REPO/get.rm.master.sh $GENOME
cd ..

cd $GENOME.repmask.denovo/
sh $REPO/get.rm.master.sh $GENOME
cd ..

cd $GENOME.repmask.repbase
sh $REPO/get.rm.master.sh $GENOME
cd ..

################################
#	Removal of excess files  

rm list_of_CACTAS
rm $GENOME.denovolib.classified.2
rm $GENOME.csv
rm $GENOME.fai
rm $GENOME.sds $GENOME.md5 $GENOME.esq $GENOME.des $GENOME.csv $GENOME.suf $GENOME.prj
rm $GENOME.llv $GENOME.lcp $GENOME.N50* formatdb.log error.log *.masked
rm $GENOME.TRIM* $GENOME.retro99*  $GENOME.retro85* *.dgt *.sort *Passed*
rm -r retro85 TRIM85 retro99 TRIM99 *fasta_files
rm $GENOME.retrotransposons.result99 $GENOME.retrotransposons.outinner99*
rm $GENOME.retrotransposons.out99 *retrotransposons.gff99
rm $GENOME.retrotransposons.result85 $GENOME.retrotransposons.outinner85*
rm $GENOME.retrotransposons.out85 $GENOME.retrotransposons.gff85 *withdomains*
rm $GENOME.retrotransposons.TRIMs.lib $GENOME.retrotransposons.TRIMs.lib.masked
rm *DFAM* $GENOME.repmod.unknowns.lib $GENOME.repmod.unknowns.all.dfam_results
rm $GENOME.repmod.unknowns.dfam_tophits $GENOME.repmod.unknowns.annotated.list
rm $GENOME.repmod.unknowns.annotated.list.to_fastagrep 
rm $GENOME.repmod.unknowns.dfam_classified.lib $GENOME.denovolib.no_classifiable_unknowns
rm $GENOME.tPSI.classified.fasta *80* $GENOME.repeats.fasta 
rm $GENOME.repeats.srt.nr.no_uniprot *unknowns* $GENOME.denovolib
rm list_of_CACTAS
rm $GENOME.denovolib.classified.2
rm $GENOME.ClassificationScheme $GENOME.classified_repeats.fasta list_to_fastagrep
rm $GENOME.non_classified_repeats.fasta
rm *blastx*
rm *repeats.only_uniprot
rm *repeats.only_repeatpeps
rm $GENOME.retrotransposons.TRIMs.lib.classified
rm unaligned.fa
rm halfoftable1					
rm halfoftable2	
################################

#------------------------------------End of pipeline-------------------------------------#

rsync -ravz * ${origdir}/repeats/${species} 
