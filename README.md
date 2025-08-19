# Applied Genomics – Slide Recap

Full slide deck: [Google Drive](https://drive.google.com/file/d/1tV58Ldxbase2jpN9sgMsIyVQvUfpziG5/view?usp=sharing)

This document summarizes the key concepts from **README1–12**, integrating branches of genetics, molecular tools, sequencing technologies, genome assembly, annotation, population applications, and applied genomics.

---

## 1. Branches of Genetics

* **Classical genetics**: Mendelian inheritance, pedigree analysis, laws of segregation & independent assortment.
* **Molecular genetics**: DNA as genetic material, PCR, Sanger sequencing.
* **Population genetics**: Hardy–Weinberg equilibrium, allele/genotype frequencies, LD, inbreeding, ROH.
* **Quantitative genetics**: phenotype = genetics + environment, variance decomposition, heritability, QTL.
* **Genomics**: whole-genome study, WGS, metadata, multi-omics.

---

## 2. Molecular Proof of Genes

* DNA demonstrated as genetic material (Griffith → Avery → Hershey–Chase).
* **PCR**: exponential amplification of DNA.
* **Sanger sequencing**: ddNTP-based chain termination, gold standard.
* Foundation for NGS and genomics.

---

## 3. Population-Level Genetics

* **HWE**: expected genotype frequencies under random mating.
* **Allele frequency** calculation from genotype counts.
* **Inbreeding coefficient (F)**: probability alleles are identical by descent.
* **LD (r², D′)**: non-random association of alleles at loci.
* **ROH**: signatures of inbreeding or bottlenecks.

---

## 4. Quantitative Traits

* Phenotypes shaped by **G + E + G×E**.
* Variance partitioned into VA, VD, VI.
* **Heritability**: H² (broad), h² (narrow, predictive).
* **QTL mapping & GWAS** link traits to genomic regions.

---

## 5. Sequencing Technologies

* **Sanger**: accurate, low throughput.
* **Illumina**: short, accurate reads, cluster generation, paired-end.
* **Ion Torrent**: detects H+ release, ionograms.
* **454 Roche & ABI SOLiD**: early NGS, discontinued.
* **PacBio (SMRT, HiFi)**: long accurate reads.
* **Nanopore (MinION, PromethION)**: ultra-long reads, real-time.
* Costs dropped faster than Moore’s Law.

---

## 6. Genome Assembly

* **Shotgun sequencing**: random fragmentation.
* Algorithms: Greedy, OLC, de Bruijn graphs.
* **k-mers**: basis of assembly and genome size estimation.
* **Scaffolding** with mate-pairs or long reads.
* **N50, coverage, BUSCO**: assess assembly quality.

---

## 7. Genome Annotation

* **Repeat annotation**: RepeatMasker, Dfam, TEannot.
* **Gene models**:

  * Ab initio (AUGUSTUS).
  * Homology-based (BLAST, Exonerate).
  * Transcript evidence (RNA-seq).
  * Integrative tools (MAKER, BRAKER2).
* Functional annotation: GO, KEGG, domains.

---

## 8. Data Formats & Tools

* **FASTQ**: raw reads + quality.
* **SAM/BAM**: alignments, CIGAR strings.
* **VCF**: SNPs, indels, SVs.
* **GFF/GTF, BED**: annotations & intervals.
* QC: FastQC, trimming (Trimmomatic).
* Alignment: BWA-MEM; Variant calling: GATK.
* Visualization: IGV, Ensembl.

---

## 9. Population Applications

* **Genotyping**: SNP arrays, GBS, RAD, Re-GBS.
* **Structural variants**: SNPs, CNVs, inversions, translocations.
* **GWAS pipeline**: QC → population structure → association → Manhattan plots.
* **MAS/MAB**: genomic selection in breeding.

---

## 10. Specialized Sequencing

* **aCGH**: CNV detection via hybridization.
* **Pool-seq**: allele frequencies in populations.
* **Targeted sequencing**: amplicon/hybrid capture panels.
* **Methyl-seq / Bisulfite**: DNA methylation at base resolution.
* **RNA-seq**: transcriptome profiling, isoforms, expression.

---

## 11. Applied Genomics

* **Comparative genomics**: synteny, phylogeny, pathogen evolution.
* **Epigenomics**: DNA/histone modifications, chromatin accessibility.
* **Transcriptomics**: bulk and single-cell RNA-seq.
* **Functional genomics**: CRISPR, RNAi, gene perturbation.
* **Applications**:

  * Medicine (precision genomics).
  * Agriculture (genomic breeding).
  * Industry (synthetic biology).
  * Environment (metagenomics).

---

## ✅ Final Takeaway

Genetics has evolved from Mendelian inheritance to **whole-genome, multi-omics approaches**. Modern genomics integrates **sequencing technologies, assembly, annotation, population analyses, and applied biotechnology**. Together, these tools enable us to link **genotype → phenotype → application** across medicine, agriculture, and environmental sciences.
