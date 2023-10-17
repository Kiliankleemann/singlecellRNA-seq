CONDA_SUBDIR=linux-64 conda create -n scTE_pipeline   # create a new environment
conda activate scTE_pipeline
conda env config vars set CONDA_SUBDIR=linux-64
conda deactivate
conda activate scTE_pipeline
conda install star
conda install scTE






##scTE pipeline
gtf2bed < /media/kilian/DATA/Human_PBMC_NEUTRO/reference/GRCh38_GENCODE_rmsk_TE.gtf > /media/kilian/DATA/Human_PBMC_NEUTRO/reference/GRCh38_GENCODE_rmsk_TE.bed

scTE_build -te /media/kilian/DATA/Human_PBMC_NEUTRO/reference/hsflnil1_8438_rm.bed -gene /media/kilian/DATA/Human_PBMC_NEUTRO/reference/GRCh38.refGene.gtf -o /media/kilian/DATA/Human_PBMC_NEUTRO/reference/scTE_hg38


#Generate STAR genome
STAR  --runMode genomeGenerate --runThreadN 16 --genomeDir STAR_index_hg38 --genomeFastaFiles reference/GRCh38.p13.genome.fa  --sjdbGTFfile reference/GRCh38.refGene.gtf


#Running STARsolo for alignment - output unfiltered BAM files for scTE pipeline
 --outSAMattributes NH HI AS nM CR CY UR UY \
 --readFilesCommand zcat \
 --outFilterMultimapNmax 100 \
 --winAnchorMultimapNmax 100 \
 --outMultimapperOrder Random \
 --runRNGseed 777 \
 --outSAMmultNmax 1 \
 --twopassMode Basic



#Make samplelist of all fastq files
echo 'Detected fastq samples'
find ./ -name "*_1.fastq.gz" -maxdepth 1 -type f -exec basename "{}" \; |  cut -d '_' -f1 | sort -u > sample_list.txt
cat sample_list.txt

#Mulimapping with EM algorithm
cat sample_list.txt | while read sample; do
	STAR --runMode alignReads \
	--soloType Droplet \
	--genomeDir STAR_index_hg38 \
	--readFilesIn fastq/${sample}_4.fastq.gz fastq/${sample}_3.fastq.gz \
	--soloCBwhitelist /media/kilian/DATA/CSF_aging_TE_analysis/reference/3M-february-2018.txt \
	--outSAMtype BAM Unsorted \
	--outSAMattributes NH HI AS nM CR CY UR UY \
	--readFilesCommand zcat \
	--outFilterMultimapNmax 100 \
	--winAnchorMultimapNmax 100 \
	--outMultimapperOrder Random \
	--runRNGseed 777 \
	--outSAMmultNmax 1 \
	--twopassMode Basic \
	--soloUMIlen 10 \
	--soloMultiMappers EM \
	--runThreadN 16 \
	--outFileNamePrefix Alignment1/ \
	--outReadsUnmapped Fastx \
	--outSAMunmapped Within \
	--soloBarcodeReadLength 0
done

#Mulimapping with UNIFORM algorithm
cat sample_list.txt | while read sample; do
	STAR --runMode alignReads \
	--soloType Droplet \
	--genomeDir STAR_index_hg38 \
	--readFilesIn fastq/${sample}_4.fastq.gz fastq/${sample}_3.fastq.gz \
	--soloCBwhitelist /media/kilian/DATA/CSF_aging_TE_analysis/reference/3M-february-2018.txt \
	--outSAMtype BAM Unsorted \
	--outSAMattributes NH HI AS nM CR CY UR UY \
	--readFilesCommand zcat \
	--outFilterMultimapNmax 100 \
	--winAnchorMultimapNmax 100 \
	--outMultimapperOrder Random \
	--runRNGseed 777 \
	--outSAMmultNmax 1 \
	--twopassMode Basic \
	--soloUMIlen 10 \
	--soloMultiMappers Uniform \
	--runThreadN 12 \
	--outFileNamePrefix Alignment_uniform/ \
	--outReadsUnmapped Fastx \
	--outSAMunmapped Within \
	--soloBarcodeReadLength 0
done


#scTE
mkdir scTE_hd5ad_files_alignment_uniform
cat sample_list.txt | while read sample; do
 scTE -i ${sample}/Aligned.out.bam -o scTE_hd5ad_files_alignment_uniform/${sample} -x mm10.exclusive.idx --hdf5 True -CB CR -UMI UR -p 8
done
