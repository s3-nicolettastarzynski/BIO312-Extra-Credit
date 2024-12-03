# Bioinformatics Final Repository

This repository provides a comprehensive, step-by-step guide to the independent analysis of dynein genes performed throughout the semester. The analysis aims to determine the evolutionary history of these genes, exploring their phylogeny, functional evolution, and genetic variation. Each repository section corresponds to a specific laboratory session or analytical approach used during the course. All necessary code, explanations of used software, and parameters implemented are included.

**1. Initial Gene Exploration and Finding Homologs (Lab 3)**

To begin exploring the dynein gene family, one member of the gene family was identified and used as a query to other members of the gene family. The StartingGene was NP_079039.4. Using the GenBank gene page, the StartingGene was used as an accession, where this particular gene was identified as dynein axonemal intermediate chain 4 from *H. sapiens.*

The first step in determining the evolution of the dynein gene family was to identify homologs using BLAST KEY, and filtering the BLAST output to only include high-scoring homologs by requiring an e-value less than 1e-30. BLAST will allow you to find homologous sequences to the query protein.

**Instructions**

1. Create a new working directory within the home directory to store files related to dynein analysis.

```
mkdir ~/lab03-$MYGIT/dyneins
```

2. Navigate to this folder using the cd command and verify that you are in the correct directory using pwd.

```
cd ~/lab03-$MYGIT/dyneins

pwd 
```

3. Download the StartingGene query protein from NCBI

```
ncbi-acc-download -F fasta -m protein "NP_079039.4"
```

4. Perform the BLAST search and view the condensed output  

```
blastp -db ../allprotein.fas -query NP_079039.4.fa -outfmt 0 -max_hsps 1 -out dyneins.blastp.typical.out

less -S dyneins.blastp.typical.out
```

5. For easier interpretation, the BLAST search can also be requested in a tabular output, if preferred.

```
blastp -db ../allprotein.fas -query NP_079039.4.fa  -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out dyneins.blastp.detail.out
```

6. Require the e-value to be less than 1e-30 to filter out low-scoring homologs

```
awk '{if ($6< 1e-30)print $1 }' dyneins.blastp.detail.out > dyneins.blastp.detail.filtered.out
```

7. Count the total number of hits in the BLAST results after implementing the e-value

```
wc -l dyneins.blastp.detail.filtered.out 
```

8. Determine the number of Paralogs found in each species

```
grep -o -E "^[A-Z]\.[a-z]+" dyneins.blastp.detail.filtered.out  | sort | uniq -c
```

**2. Sequence Alignment (Lab 4)**

Utilizing the filtered results from the BLAST search obtained earlier, a multiple sequence alignment can be performed to predict homologous positions across dynein homologs. This alignment will also predict which positions may contain insertions or deletions when compared to the other homologous sequences. The average percent identity of all sequences in the alignment is further determined using T-COFFEE and AlignBuddy.

**Instructions**

1. Create a new folder to store all files for Lab 4 and navigate to this folder

```
mkdir ~/lab04-$MYGIT/dyneins
cd ~/lab04-$MYGIT/dyneins
```

2. Use SeqKit to obtain the sequences from the filtered BLAST output file containing dynein homologs

```
seqkit grep --pattern-file ~/lab03-$MYGIT/dyneins/dyneins.blastp.detail.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/dyneins/dyneins.homologs.fas 
```

3. Perform the global multiple sequence alignment with MUSCLE

```
muscle -align ~/lab04-$MYGIT/dyneins/dyneins.homologs.fas -output ~/lab04-$MYGIT/dyneins/dyneins.homologs.al.fas
```

4. View the Alignment using Alv

```
alv -kli  ~/lab04-$MYGIT/dyneins/dyneins.homologs.al.fas | less -RS
```

5. Save the Alignment as a PDF

```
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R  ~/lab04-$MYGIT/dyneins/dyneins.homologs.al.fas
```

6. Determine the Total Length of the Alignment

```
alignbuddy  -al  ~/lab04-$MYGIT/dyneins/dyneins.homologs.al.fas
```

7. Determine the Total Length of the Alignment After Removing Columns Containing Gaps

```
alignbuddy -trm all  ~/lab04-$MYGIT/dyneins/dyneins.homologs.al.fas | alignbuddy  -al
```

8. Calculate The Total Length of the Alignment After Removing Invariant Positions

```
 alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/dyneins/dyneins.homologs.al.fas | alignbuddy  -al
 ```

9. Calculate the Average Percent Identity Using T-Coffee

```
 t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/dyneins/dyneins.homologs.al.fas -output sim
 ```

10. Calculate the Average Percent Identity Using AlignBuddy

```
 alignbuddy -pi ~/lab04-$MYGIT/dyneins/dyneins.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
     END{ print(100*sum/num) } '
```

**3. Investigating the Phylogenetic Relationships of Dynein Genes using IQ-Tree (Lab 5)**

Using the global multiple sequence alignments obtained using MUSCLE, the optimal phylogenetic tree for dynein homologs can be reconstructed using IQ-TREE. This software allows us to predict the best amino exchange and rate of heterogeneity models for the global multiple sequence alignment. For the dynein family, no distinct outgroups are identified at this moment; hence, the tree is rooted at the midpoint. A midpoint-rooted tree will allow us to identify the oldest divergence event in the dynein gene phylogeny.

**Instructions**

1. Create a directory in Lab5 for the dynein gene family and navigate to it

```
mkdir ~/lab05-$MYGIT/dyneins
cd ~/lab05-$MYGIT/dyneins
```

2. Remove any sequence containing a duplicate label and input the file into the directory created

```
sed 's/ /_/g'  ~/lab04-$MYGIT/dyneins/dyneins.homologs.al.fas | seqkit grep -v -r -p "dupelabel" >  ~/lab05-$MYGIT/dyneins/dyneins.homologsf.al.fas
```

3. Use IQ-TREE to find the maximum likehood tree estimate with ultrafast bootstrap support levels.

```
iqtree -s ~/lab05-$MYGIT/dyneins/dyneins.homologsf.al.fas -bb 1000 -nt 2
```

4. Display the Newick formatted tree file

```
nw_display ~/lab05-$MYGIT/dyneins/dyneins.homologsf.al.fas.treefile
```

5. Graphically Display the Unrooted Phylogenetic Tree with a Small Text Size (0.4) and Appropriate Label Lengths for Optimal Visualization (15)

```
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R  ~/lab05-$MYGIT/dyneins/dyneins.homologsf.al.fas.treefile ~/lab05-$MYGIT/dyneins/dyneins.homologsf.al.fas.treefile.pdf 0.4 15
```

6. Reroot the Tree at the Midpoint Using GoTree

```
gotree reroot midpoint -i ~/lab05-$MYGIT/dyneins/dyneins.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile
```

7. View the Rooted Tree in a Text File

```
nw_order -c n ~/lab05-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile  | nw_display -
```

8. Convert the Rooted Tree File into an Svg Image for Better Visualization

```
nw_order -c n ~/lab05-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s  >  ~/lab05-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile.svg -
```

9. Convert the Svg File into a PDF file

```
convert  ~/lab05-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile.svg  ~/lab05-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile.pdf
```

**4. Reconciling a Gene and Species Tree (Lab 6)**

Once the phylogenetic tree was constructed, the gene and species trees need to be reconciled using the Notung software to predict reconciliation cost, duplication, and loss events. This reconciliation further allows us to distinguish orthologs and paralogs from the identified dynein gene homologs obtained via BLAST search. Using the Thirdkind software, the gene-species reconciliation tree is created as a graphic to investigate changes in dynein gene copies across different vertebrate lineages. In this image, evolutionary events can be distinguished as various shapes. The tree's root is denoted as a triangle, whereas loss events, gene duplications, and speciation events are depicted as x’s, squares, and circles, respectively.

**Instructions**

1. Create a New Folder to Store Files for Lab 6 and Navigate to it

```
mkdir ~/lab06-$MYGIT/dyneins 
cd ~/lab06-$MYGIT/dyneins
```

2. Copy the Midpoint-Rooted Gene Tree from Lab 5 into the Lab 6 Folder

```
cp/lab05-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile ~/lab06-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile 
```

3. Perform the Gene-Species Reconciliation using Notung from the Midpoint-Rooted Tree File

```
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-$MYGIT/species.tre -g ~/lab06-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/dyneins/
```

4. Generate a RecPhyloXML Object and View the Gene-Within-Species Tree via Thirdkind

```
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile.rec.ntg --include.species
```

5. Create a Gene-Reconciliation-Within Species Tree Reconciliation Graphic using Thirdkind

```
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile.rec.ntg.xml -o  ~/lab06-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile.rec.svg 
```

6. Convert this graphic to a PDF

```
convert  -density 150 ~/lab06-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile.rec.svg ~/lab06-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile.rec.pdf
```

**5. Protein Domain Prediction (Lab 8)**

To investigate the conservation of protein domains within the dynein gene family and identify domains within homologous dynein sequences, RPS-BLAST and Pfam software. To exclusively obtain only strong matches in domains, the e-value required will be set to 0.001. From the Pfam software, the predicted domains will be plotted onto the midpoint-rooted phylogenetic tree obtained from IQ-TREE using the ggtools and drawproteins packages in R script. This will allow us to determine the evolutionary relationships of protein domains present in dynein genes and individually visualize the domains present across various dynein homologs.

1. Copy the Raw Unaligned Dynein Gene Sequence and Remove Stop Codons

```
sed 's/*//' ~/lab04-$MYGIT/dyneins/dyneins.homologs.fas > ~/lab08-$MYGIT/dyneins/dyneins.homologs.fas
```

2. Run the RPS-BLAST with a Required e-value of 0.001

```
rpsblast -query ~/lab08-$MYGIT/dyneins/dyneins.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/dyneins/dyneins.rps-blast.out  -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .001
```

3. Create a Copy of the Phylogenetic Gene Tree from Lab5

```
cp ~/lab05-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile ~/lab08-$MYGIT/dyneins
```

4. Plot the Pfam Domain Predictions from RPS-BLAST Next to Their Relative Protein on the Phylogenetic Tree

```
Rscript  --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/dyneins/dyneins.homologsf.al.mid.treefile ~/lab08-$MYGIT/dyneins/dyneins.rps-blast.out ~/lab08-$MYGIT/dyneins/dyneins.tree.rps.pdf
```

5. View the Annotations in the Output File

```
mlr --inidx --ifs "\t" --opprint  cat ~/lab08-$MYGIT/dyneins/dyneins.rps-blast.out | tail -n +2 | less -S
```

6. Determine the Pfam Domains with an Annotation of Zero

```
cut -f 1 ~/lab08-$MYGIT/dyneins/dyneins.rps-blast.out | sort | uniq -c
```

7. Determine the Most Common Pfam Domain Annotation

```
cut -f 6 ~/lab08-$MYGIT/dyneins/dyneins.rps-blast.out | sort | uniq -c
```

8. Determine the Longest Annotated Protein Domain Found

```
awk '{a=$4-$3;print $1,'\t',a;}' ~/lab08-$MYGIT/dyneins/dyneins.rps-blast.out |  sort  -k2nr
```

9. Determine the e-values of the Protein Domain Hits

```
cut -f 1,5 -d $'\t' ~/lab08-$MYGIT/dyneins/dyneins.rps-blast.out
```
