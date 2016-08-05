make-jetset
===========

This is a pipeline to generate probe set scores and gene mappings 
used in the *jetset* R package.

The finished *jetset* R package itself is available from 
[CRAN](https://cran.r-project.org/web/packages/jetset/index.html)
 or from
[GitHub](https://github.com/aroneklund/jetset)

The algorithm used here is described in this publication:
> Jetset: selecting an optimal microarray probe set to represent a gene  
> Qiyuan Li, Nicolai J. Birkbak, Balazs Gyorffy, Zoltan Szallasi, and Aron C. Eklund  
> BMC Bioinformatics 2011 Dec 15; 12:474  
> http://www.biomedcentral.com/1471-2105/12/474/abstract  


Prerequisites
=============


You need to have:

1. Linux/unix/MacOS
2. [R](https://www.r-project.org)
3. [NCBI Standalone BLAST](http://www.ncbi.nlm.nih.gov/books/NBK52640/).
  I am using blastall 2.2.25, and I haven't checked whether other versions work as well.


How to proceed
==============

1. Obtain the probe sequences in FASTA format. I need to find my notes on where this came
  from, but fwiw my file "probe.hgu133a.fa" starts like this:

        >1007_s_at:467:181:3330:A
        CACCCAGCTGGTCCTGTGGATGGGA
        >1007_s_at:531:299:3443:A
        GCCCCACTGGACAACACTGATTCCT
        >1007_s_at:86:557:3512:A
        TGGACCCCACTGGCTGAGAATCTGG

2. Get the latest RefSeq fasta files from NCBI:

        wget "ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.*.rna.fna.gz"
        date +%F > human.rna.version
        gunzip -c human.*.rna.fna.gz >human.rna.fna

3. Create a version with simpler headers (RefSeq ID only):

        perl -pe 's/^>.*ref\|(.._[\d\.]+)\|.*$/>\1/' <human.rna.fna >human.rna2.fna

4. Create a file containing the length of every peptide in RefSeq:

        perl getFastaLengths.pl <human.rna2.fna >human.rna.len

5. Create a database for BLAST searching:

        formatdb -i human.rna2.fna -p F -n refseq.human.rna

6. Create a RefSeq ID to Entrez Gene ID lookup table (based on the R package 
   "org.Hs.eg.db"). Note that this script will obtain the latest version of org.Hs.eg.db:

        Rscript makeREFSEQ2EG.R

7. Run the BLAST searches (I use 4 cores at a time):

        blastall -p blastn -d refseq.human.rna -i probe.hgu95av2.fa    -o blastresult.hgu95av2.refseq.txt    -F F -m 8 -e 1 -W 8 -a 4
        blastall -p blastn -d refseq.human.rna -i probe.hgu133a.fa     -o blastresult.hgu133a.refseq.txt     -F F -m 8 -e 1 -W 8 -a 4
        blastall -p blastn -d refseq.human.rna -i probe.hgu133plus2.fa -o blastresult.hgu133plus2.refseq.txt -F F -m 8 -e 1 -W 8 -a 4
        blastall -p blastn -d refseq.human.rna -i probe.u133x3p.fa     -o blastresult.u133x3p.refseq.txt     -F F -m 8 -e 1 -W 8 -a 4

8. Come back in about 11 hours (on my machine, at least).

9. Perform the score calculations based on BLAST results (these are single-threaded, so I do these simultaneously):

        Rscript calculateJetset.R hgu95av2
        Rscript calculateJetset.R hgu133a
        Rscript calculateJetset.R hgu133plus2
        Rscript calculateJetset.R u133x3p

10. Come back in about an hour (on my machine, at least).

11. Optionally, reduce the size of the RData files (in R):

        library(tools)
        resaveRdaFiles('.')

12. At this point in my usual pipeline, I add the resulting .RData files into the jetset
  package.  This step is not described here.

12. Optionally, create CSV versions of the files. *Note that this makes CSV files from
   whichever version of jetset is currently in your R installation.*

        Rscript makeCSVfiles.R


Ta-dah!


