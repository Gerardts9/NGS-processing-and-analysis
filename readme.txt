1.

# We use the edgeR.analysis.R R function to perform the differential expression analysis between brain and liver. We use the same matrix and metadata used to compare brain and heart but now we select coefficient 3 to specify this tissue comparison:

edgeR.analysis.R --input_matrix ../quantifications/encode.mouse.gene.expected_count.idr_NA.tsv \
                 --metadata /tutorial/data/gene.quantifications.index.tsv \
                 --fields tissue \
                 --coefficient 3 \
                 --output brain_X_liver


# We select only the significant genes using the last column of the file we just generated. Then, using the second column and a threshold, depending on whether it takes positive or negative values, we select genes in both directions and save them into a file:

awk '$NF<0.01 && $2<-10{print $1"\tover_brain_X_liver"}' edgeR.cpm1.n2.brain_X_liver.tsv > edgeR.0.01.over_brain_X_liver.txt

awk '$NF<0.01 && $2>10 {print $1"\tover_liver_X_brain"}' edgeR.cpm1.n2.brain_X_liver.tsv > edgeR.0.01.over_liver_X_brain.txt


# Count the differentially expressed genes in each group:

wc -l edgeR.0.01.over*.txt


# First, we select, from the gencode file, the gene id and gene type attributes that belong to genes using column 3. Then, a python script joins these two columns with the column generated with the previous code. Everything is saved into a new file with 3 columns that will be used to create the heatmap:

awk '$3=="gene"{ match($0, /gene_id "([^"]+).+gene_type "([^"]+)/, var); print var[1],var[2] }' OFS="\t" /tutorial/refs/gencode.vM4.gtf \ 
| join.py --file1 stdin \ 
	  --file2 <(cat edgeR.0.01.over*.txt) \ 
| sed '1igene\tedgeR\tgene_type' > gene.edgeR.tsv


# We select the first column of the gene.edgeR.tsv file and all lines except for number 1 that belongs to the column name. Then, we use the TPMs from the encode.mouse.gene.TPM.idr_NA.tsv file and we introduce this into the R script to create a heatmap. In the command we specify to introduce several things into the plot such as the metadata, the dendrograms and the color palettes for each of the parts that need them:

cut -f1 gene.edgeR.tsv \
| tail -n+2 \
| selectMatrixRows.sh - ../quantifications/encode.mouse.gene.TPM.idr_NA.tsv \
| ggheatmap.R --width 5 \
              --height 8 \
              --col_metadata /tutorial/data/gene.quantifications.index.tsv \
              --colSide_by tissue \
              --col_labels labExpId \
              --row_metadata gene.edgeR.tsv \
              --merge_row_mdata_on gene \
              --rowSide_by edgeR,gene_type \
              --row_labels none \
              --log \
              --pseudocount 0.1 \
              --col_dendro \
              --row_dendro \
              --matrix_palette /tutorial/palettes/palDiverging.txt \
              --colSide_palette /tutorial/palettes/palTissue.txt \
              --output heatmap.brain_X_liver.pdf


2. 

# We prepare a file with a list of all ensembl genes in the annotation using column 10 and removing values after the . to use the as universe for the enrichment analysis.

awk '{split($10,a,/\"|\./); print a[2]}' /tutorial/refs/gencode.vM4.gtf | sort -u > universe.txt


# We select the ENSEMBL genes id from the two files generated with the differential expression analysis and then we use them to perform the gene ontology analysis with the GO_enrichment.R function. We set some parameters such as the gene universe and, in our case, we will use biological process (BP) as the GO category used:

awk '{split($1,a,"."); print a[1]}' edgeR.0.01.over_brain_X_liver.txt \
| GO_enrichment.R --universe universe.txt \
                  --genes stdin \
                  --categ BP \
                  --output edgeR.over_brain_X_liver \
                  --species mouse

awk '{split($1,a,"."); print a[1]}' edgeR.0.01.over_liver_X_brain.txt \
| GO_enrichment.R --universe universe.txt \
                  --genes stdin \
                  --categ BP \
                  --output edgeR.over_liver_X_brain \
                  --species mouse


# We use the GO files generated to extract the Gene Ontology terms and their p-values to use them in the REVIGO web page:

awk 'NR==1{$1="% "$1}{print $1,$2}' edgeR.over_brain_X_liver.BP.tsv
awk 'NR==1{$1="% "$1}{print $1,$2}' edgeR.over_liver_X_brain.BP.tsv


3. 

# The first file we need to create is a list of all transcripts' ids from the gencode file of protein coding genes. We use column 3 to select transcripts and we only select protein coding from the gene type column:

awk '$3=="transcript" && $0~/gene_type "protein_coding"/{ match($0, /transcript_id "([^"]+)/, id); print id[1] }' /tutorial/refs/gencode.vM4.gtf |sort -u > protein_coding_transcript_IDs.txt


# Now we create an annotation file only for exons using, again, the whole annotation file and the file recently created: 

cat /tutorial/refs/gencode.vM4.gtf |awk '$3=="exon"' |grep -Ff protein_coding_transcript_IDs.txt > exon-annot.gtf


# We use the bash script selectMatrixRows.sh to filter the TPM matrix to only contain protein-coding transcripts:

selectMatrixRows.sh protein_coding_transcript_IDs.txt /tutorial/quantifications/encode.mouse.transcript.TPM.idr_NA.tsv > pc-tx.tsv


# We create a transcript expression matrix file for each of the 3 tissues using a for loop:

for tissue in Brain Heart Liver; do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 pc-tx.tsv > expr.${tissue}.tsv done


# We extract alternative splicing events from the exon annotation files of the following types: Skipping Exon (SE), Retained Intron (RI), Mutually Exclusive Exons (MX) and Alternative First exons (AF, included in FL with alternative last exon):

suppa.py generateEvents -i exon-annot.gtf -e SE RI MX FL -o localEvents -f ioe



# We obtained the following: 12,254 of Skipping exon, 2763 of retained intron, 1054 of mutually exclusive exons and 20128 of alternative first exon.

wc -l localEvents*ioe


# We create a bar plot with the different alternative splicing types used in this exercise and including also A3 and A5 ones:

wc -l localEvents*ioe | grep -v total | awk '{split($2,a,"_"); print a[2]"\t"$1}' | sed '1iEvent\tNumber' | /tutorial/teaching-utils/ggbarplot.R -i stdin --header -o test.pdf -t Barplot_for_our_data


# We compute PSI for each of the types and for each of the tissues using, first, suppa.py and, second, a loop and the selectMatrixColumns bash script. To do this, we substitute the event variable for each type of splicing:

# SE:
event=SE; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}
event=SE; for tissue in Brain Heart Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done

#RI:
event=RI; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}
event=RI; for tissue in Brain Heart Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done

#MX:
event=MX; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}
event=MX; for tissue in Brain Heart Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done

#AF:
event=AF; suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file pc-tx.tsv -o PSI-${event}
event=AF; for tissue in Brain Heart Liver ;do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi;done


# We compute the differential splicing using suppa.py between brain and liver for the four different types of splicing we are analyzing obtaining a dpsi file that we will use to create the heatmap in the future:

event=SE; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}
event=RI; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}
event=MX; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}
event=AF; suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi --tpm expr.Brain.tsv expr.Liver.tsv -c -gc -o DS.${event}

# We create a heatmap for each of the four types of splicing we are working with. In order to do this, first, we create a txt file with the ENSEMBL IDs of the transcripts we want to incorporate into the heatmap, removing nan values and using columns 2 and 3 of the .dpsi file. Then, we create a tsv file adding the PSI values as columns at the txt file using the selectMatrixRows.sh script. Finally, we use the R function ggheatmap.R and the tsv file just created to plot the heatmap.

event=SE; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.4 || $2<-0.4) && $3<0.05{print $1}' DS.${event}.dpsi > input-SE-heatmap.txt
selectMatrixRows.sh input-SE-heatmap.txt DS.SE.psivec > input-SE-heatmap.tsv
ggheatmap.R -i input-SE-heatmap.tsv --row_dendro -o test-heatmap-SE.pdf --matrix_palette /tutorial/palettes/palSequential.txt --matrix_fill_limits "0,1" -B 8 –col_dendro

event=RI; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.2 || $2<-0.2) && $3<0.05{print $1}' DS.${event}.dpsi > input-RI-heatmap.txt
selectMatrixRows.sh input-RI-heatmap.txt DS.RI.psivec > input-RI-heatmap.tsv
ggheatmap.R -i input-RI-heatmap.tsv --row_dendro -o test-heatmap-RI.pdf --matrix_palette /tutorial/palettes/palSequential.txt --matrix_fill_limits "0,1" -B 8 --col_dendro

event=MX; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.2 || $2<-0.2) && $3<0.1{print $1}' DS.${event}.dpsi > input-MX-heatmap.txt
selectMatrixRows.sh input-MX-heatmap.txt DS.MX.psivec > input-MX-heatmap.tsv
ggheatmap.R -i input-MX-heatmap.tsv --row_dendro -o test-heatmap-MX.pdf --matrix_palette /tutorial/palettes/palSequential.txt --matrix_fill_limits "0,1" -B 8 --col_dendro

event=AF; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.5 || $2<-0.5) && $3<0.05{print $1}' DS.${event}.dpsi > input-AF-heatmap.txt
selectMatrixRows.sh input-AF-heatmap.txt DS.AF.psivec > input-AF-heatmap.tsv
ggheatmap.R -i input-AF-heatmap.tsv --row_dendro -o test-heatmap-AF.pdf --matrix_palette /tutorial/palettes/palSequential.txt --matrix_fill_limits "0,1" -B 8 --col_dendro


# 4:

# We use the function bedtools intersect to obtain a tsv file containing the common Chip-seq peaks between brain and liver. We use the narrow peaks files of these two tissues and then we chose the columns we want to include. We select only the rows where the last column that belongs to intersection is different to 0:

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoLiver.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Brain_coordinates","Brain_peak","Liver_coordinates","Liver_peak","intersection"}$NF!=0{print $1":"$2"-"$3,$4,$11":"$12"-"$13,$14,$NF}' > common-peaks-Brain-liver.tsv


# We use the same function and files but choose different columns to create the bed file:

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoLiver.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Brain_coordinates","Brain_peak","Liver_coordinates","Liver_peak","intersection"}$NF!=0{print $1,$2,$3,$4";"$14}' > common-peaks-Brain-liver.bed


# To create the files containing specific peaks for each tissue, we use the same structure as before but selecting only peaks with no intersection:

bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoLiver.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Brain_coordinates","Brain_peak"}$NF==0{print $1":"$2"-"$3,$4}' > Brain-specific-peaks.tsv
bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoLiver.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Brain_coordinates","Brain_peak"}$NF==0{print $1,$2,$3,$4}' > Brain-specific-peaks.bed
bedtools intersect -a /tutorial/results/CHIPembryoLiver.narrowPeak -b /tutorial/results/CHIPembryoBrain.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Liver_coordinates","Liver_peak"}$NF==0{print $1":"$2"-"$3,$4}' > Liver-specific-peaks.tsv
bedtools intersect -a /tutorial/results/CHIPembryoLiver.narrowPeak -b /tutorial/results/CHIPembryoBrain.narrowPeak -wao|awk 'BEGIN{FS=OFS="\t"; print "Liver_coordinates","Liver_peak"}$NF==0{print $1,$2,$3,$4}' > Liver-specific-peaks.bed


# Finally, to create the barplot we select the files that we just created and count their lines to know the number of elements we have in each group. Then, we use R function ggbarplot.R, introducing the values we obtained and some other information such as the palette to form it:

wc -l *.bed | grep -v total | awk '{split($2,a,"-"); print a[1]"\t"$1}' | sed '1iEvent\tNumber' | /tutorial/teaching-utils/ggbarplot.R -i stdin --header -o Barplot.pdf -t Barplot_Exercise_4 --palette_fill /tutorial/palettes/palTissue.txt --fill_by 1


# 5:

# We create a sorted  bed file with all the protein coding genes like the one used during the hands-on:

awk 'BEGIN{FS=OFS="\t"}$3=="gene" && $0~/gene_type "protein_coding"/ && $7=="+"{ match($0, /transcript_id "([^"]+)/, id); print $1,$5-200,$5+200,id[1],$7 }$3=="gene" && $0~/gene_type "protein_coding"/ && $7=="-"{ match($0, /transcript_id "([^"]+)/, id); print $1,$5-200,$5+200,id[1],$7 }' /tutorial/refs/gencode.vM4.gtf |sort -u > protein-coding-genes-200up_downTSS.bed


# Subset the overexpressed genes in each of the two tissues using the bed file we just created:

cut -f1 ../analysis/edgeR.0.01.over_liver_X_brain.txt |sort -u| grep -Ff - protein-coding-genes-200up_downTSS.bed > genes-over-liver-TSS.bed
cut -f1 ../analysis/edgeR.0.01.over_brain_X_liver.txt |sort -u| grep -Ff - protein-coding-genes-200up_downTSS.bed > genes-over-brain-TSS.bed 


# Select the specific methylation peaks and the common ones:

cat ../chip-analysis/Brain-specific-peaks.bed | tail -n+2 | sort -u > Brain-peaks.bed 
cat ../chip-analysis/Liver-specific-peaks.bed | tail -n+2 | sort -u > Liver-peaks.bed 
cat ../chip-analysis/Common-peaks-Brain-liver.bed | tail -n+2 | sort -u > Common-peaks.bed


# Finally, we use bedtools intersect to obtain overexpressing genes and with methylation peaks in each of the different conditions proposed:

bedtools intersect -a genes-over-brain-TSS.bed -b Brain-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > brain-peaks-over-brain.bed
bedtools intersect -a genes-over-liver-TSS.bed -b Liver-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > liver-peaks-over-liver.bed
bedtools intersect -a genes-over-brain-TSS.bed -b Liver-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > liver-peaks-over-brain.bed
bedtools intersect -a genes-over-liver-TSS.bed -b Brain-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > brain-peaks-over-liver.bed
bedtools intersect -a genes-over-liver-TSS.bed -b Common-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > common-peaks-over-liver.bed
bedtools intersect -a genes-over-brain-TSS.bed -b Common-peaks.bed -wao|awk 'BEGIN{FS=OFS="\t"}$NF!=0{print $1,$2,$3,$4,$11,$12,$13,$14,$NF}' > common-peaks-over-brain.bed


# We count then number of lines in the files:

wc -l *brain.bed
wc -l *liver.bed


