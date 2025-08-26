# Applied Genomics – Slide Recap

Full slide deck: [Google Drive](https://drive.google.com/file/d/1tV58Ldxbase2jpN9sgMsIyVQvUfpziG5/view?usp=sharing)

This document summarizes the key concepts from **README1–12**, integrating branches of genetics, molecular tools, sequencing technologies, genome assembly, annotation, population applications, and applied genomics.

---

**Genetics** is the discipline that studies **biological diversity at different levels: genomes, populations, individuals and organisms***.
**Diversity** is the **driving-force** of genetics.


## 1. Branches of Genetics

* **Classical genetics (transmission/formal)**: Mendelian inheritance, pedigree analysis, laws of segregation & independent assortment.
  -->heredity and segregation of traits.
     - focuses on the **principles of heredity** and how traits are transmitted.
     - **chromosomes** are the physical carriers of heredity, containing multiple genes.
     - **crossing over** during meiosis reshuffles paired chromosomes; recombination frequency depends on distance betweeen loci.
     - **first genetic maps** were developed in Drosophila later replaced by genome sequencing.
     - loci were studied by phenotype changes-->Only **loci with different alleles** could be detected in classical genetics.
     - Mendel formulated the **laws of heredity** from Pisum sativum experiments-->
              -**Genotype** (allele combination) vs **Phenotype** (observable traits: genotype + environment)
              -**Dominant** allele masks recessive; **Codominance**-->intermediate phenotypes possible
              -**1) Law of segregation**: each gamete carries one allele form each locus; classic 3:1 phenotype ratio in F2 generation.
              -**2) Law of Independent Assortment**: alleles for different traits assort independently if loci are on different chromosomes; 9:3:3:1 in dihybrid crosses.
     - **Pedigrees**: graphical representation of inheritance across generations.
         -Symbols: o female and square male rombous unknown sex; empty=healthy filled=affected and half-filled=carrier
         -**PLINK** is the software for genotype-phenotype data management;
            input: text file (one raw per sample and multiple fields); columns=family, individual ID, parents, sex, phenptype.
            facilitates high-throughput genotyping analysis and genome-wide association studies (GWAS).
     - **cytogenetics**: visulization of chromosomes during metaphase (condensed); involves staining techniques to identify structure and number.
     - **chromosomal map**: represent gene positions and distances based on recombination frequency; closely linked genes-->low recombination.
     - 1 cM = 1% recombination
     - **Sex chromosomes** and detrmination: in mammals XX is the female and XY the male, PARR (pseudoautosomal regions) allow delimited recombination between X and Y.
         in other system X0 (insect): XX female and X male; ZW (birds): ZZ male and ZW female; haplo-dyploid (bees): fertilized diploid female, unfertilized haploid male;           temperature-dependent sex in reptiles.
       
* **Molecular genetics**: DNA as genetic material, PCR, Sanger sequencing.
  -->DNA and gene-level mechanisms.
  
* **Population genetics**: Hardy–Weinberg equilibrium, allele/genotype frequencies, LD, inbreeding, ROH.
  -->variability within and across populations.
  -LD (linkage disequilibrium): non-random association of alleles at different loci in a population. Stronger LD when loci are physically closer tend to be inherited          together.
  
* **Quantitative genetics**: phenotype = genetics + environment, variance decomposition, heritability, QTL.
  -->polygenic traits
  
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
