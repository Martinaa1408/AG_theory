# Applied Genomics – Slide Recap

Full slide deck: [Google Drive](https://drive.google.com/file/d/1tV58Ldxbase2jpN9sgMsIyVQvUfpziG5/view?usp=sharing)

This document summarizes the key concepts from **README1–12**, integrating branches of genetics, molecular tools, sequencing technologies, genome assembly, annotation, population applications, and applied genomics.

---

**Genetics** is the discipline that studies **biological diversity at different levels: genomes, populations, individuals and organisms**.
**Diversity** is the **driving-force** of genetics.

---

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
     - **Sex chromosomes** and determination: in mammals XX is the female and XY the male, PARR (pseudoautosomal regions) allow delimited recombination between X and Y.
         in other system X0 (insect): XX female and X male; ZW (birds): ZZ male and ZW female; haplo-dyploid (bees): fertilized diploid female, unfertilized haploid male;           temperature-dependent sex in reptiles.
       
---

* **Molecular genetics**: DNA as genetic material, PCR, Sanger sequencing.
  
  -->focuses on the chemical nature of genes, their structure, function, and how genetic information is encoded, replicated, expressed and regulated.
  
     - **DNA replication, transcription, and translation**.
     - **gene regulation**.
     - **post-transcriptional and post-translational processes**.
     - the technical foundations of molecular genetics include: **recombinant DNA technology**, **DNA sequencing**, **library construction**, **PCR amplification**,               ***hybridization techniques**, **gel electrophoresis**.
     - **SANGER SEQUENCING: first generation sequencing, an early but fundamental sequencing method**.
           Principle:
            - based on controlled DNA synthesis with chain-terminating nucleotides.
            - requires: DNA template (to be sequenced), a primer (defines start point); dNTPs normal nucleotides, ddNTPs (modified nucleotides, terminate synthesis) and                            DNA polymerase enzyme
           Process:
           1- run four reactions (one per nucleotide);
           2- incorporate normal or modified nucleotides randomly during synthesis;
           3- modified nucletoides block extension, generating DNA fragments of different lengths.
           4- fragments are separated by electrophoresis which sorts by size;
           5- detection is based on light or color emitted by the final ddNTP;
           6- by analyzing fragment lengths and terminal nucleotides the DNA sequence is reconstructed.
       -->key features: uses linear amplification millions of fragments of varying lengths, highly accurate but low throughput compared to modern methods, provided the               foundation for today's sequencing technologies.
       
     - **PCR (polymerase chain amplification)**: method for amplification. Its goal is to take a specific tiny segment of DNA and make billions copies of it. This process           is exponential amplification, each cycle double the number of copies.
           Process:
             - denaturation --> annealing --> extension
           Principle:
             - requires: template DNA, primers, Taq DNA polymerase, dNTPs and buffer
       
     --> link: PCR is first used to isolate and massively amplify the specific DNA region of interest, then Sanger is used to determine the precise genetic code of that             amplified PCR product.

---

* **Population genetics**: Hardy–Weinberg equilibrium, allele/genotype frequencies, LD, inbreeding, ROH.
  
  -->variability within and across populations, infact studies the genetic composition of populations (group of individual of the same species) and howthis composition          changes over space and time

     - focus on **allele frequencies** in populations.
     - allele frequencies can be used to infer **genotype frequencies** under certain assumptions.
     - **HWE (hardy weinberg equilibrium)**: a central model (null model) in population genetics to predict expected genotype frequencies from allele frequencies.
       Assumptions: the population must satisfy the **following conditions**:
         - diploid organisms
         - only sexual reproduction occurs
         - non overlapping generations
         - random mating
         - infinite population size
         - equal allele frequencies in sexes
         - no migration, selection, mutation or gene flow
      If all assumptions hold, allele and genotype frequencies remain constant across generations.
       **Allele frequencies = Number of copies of a specific allele/total number of alleles at that locus in the population**
       **Genotype frequencies = Number of individuals with a given genotype/Total number of individuals in the population**
      Ex. biallelic locus: let p= frequency of the dominant allele and let q= frequency of the recessive allele, since only two alleles exist: p + q = 1.
      The expected genotype frequencies follow: **(p+q)^2= p^2 + 2pq + q^2 = 1**.
        where p^2 = homozygous dominant (AA), 2pq = Aa heterozygous (Aa), q^2 = homozygous recessive (aa).
      Interpretation: if the observed population frequencies differ from HWE expectations,this indicates that evolutionary forces may be acting; the model is not a               description of real populations but a reference baseline to study deviations.

     - **LD (linkage disequilibrium)**: non-random association of alleles at different loci in a population. Stronger LD when loci are physically closer tend to be                inherited together.
       Influenced by:
        - recombination rate (low recombination=high LD).
        - population size (small effective size= high LD).
        - inbreeding and bottlenceks (increase LD).
     - **Inbreeding coefficients**: Fped (pedigree-based probability that two alleles are identical by descent); FROH (fraction of genome covered by ROH); Fis (inbreeding         coefficient measuring excess homozygosity relative to HWE.
     - LD (r², D′): non-random association of alleles at loci.
     - **ROH (runs of homozigosity)**: longs continuous stretches of homozygous genotypes in the genome; indicate inbreeding both chromosomes inherited from a common              ancestor. Long ROH = recent inbreeding and short ROH = ancient inbreeding.
     - **centiMorgan (cM)**: genetic distance unit based on recombination frequency. 1 cM = 1% recombination (chance of crossover per generation). Loci close together show        high LD.
     - **F_ST**= measures genetic differentiation between populations
       Var(p)/p(1-p) ranges: 0 no differentiation and 1 populations are fixed for different alleles.

LINK-->HWE is the theoretical baseline for LD, ROH, and inbreeding measures; LD is shaped by recombination measured in cM and mantained by inbreeding or drift; cM quantifies recombination which breaks down LD over time; ROH provide a genomic measure of inbreeding; all capture inbreeding but from different sources (pedigree, genomic data, frequency-based); F_ST extends the concept of inbreeding and HWE deviations to the population level.

[HWE] baseline equilibrium-->(deviations)-->[inbreeding]--[ROH] genomic signal-->(quantified by Fped,FROH,FIS)-->(reduces effective recombination)-->[LD high]-->(linked to)-->[cM]-->(across populations)-->[F_st] differentiation between populations

Population Applications

* **Genotyping**: SNP arrays, GBS, RAD, Re-GBS.
* **Structural variants**: SNPs, CNVs, inversions, translocations.
* **GWAS pipeline**: QC → population structure → association → Manhattan plots.
* **MAS/MAB**: genomic selection in breeding.

---
       
* **Quantitative genetics**: phenotype = genetics + environment, variance decomposition, heritability, QTL.
  
  -->polygenic traits, studies traits influenced by many genes and environmental factors (ex. stature, eye color, running ability, survival to adulthood)
     Traits can be continuous (height), discrete (litter size), or binary (survival yes/no).
    -key concept: **phenotype (P) = Genotype effect (G) + Environment effect (E)**
      component: **genotype effect** = **addictive effect (A)** sum of contributions of alleles across loci; **dominance effect (D)** interaction between alleles at the                     same locus; **interaction effect or epistatic (I)** gene-gene interactions.
                 **environment effect** = permanent environmental effect and temporary environmental effect.
     ---> thus: Var (P) = Var (G) + Var (E)
    - **Heritability** is the proportion of phenotypic variance explained by genetic variance. (h^2=Var(G)/Var(P))
      low heritability: < 0.1
      medium heritability: 0.1 - 0.4
      high geritability: > 0.4
      To estimate heritability, genetic realtionships between individuals must be considered.
    - Application: quantitative genetics is essential to understand variation and covariation among relatives in natural and managed populations, study the dynamics of           evolutionary change, improve animals and plants and investigate complex disease in humans.
      
---
     
**Genomics**: whole-genome study, WGS, metadata, multi-omics.
-->what is the **genome**? the entire genetic content of an organism (DNA sequence).
-->what is **genomics**? the study of all genes in an organism and their interactions. Is a multidisciplinary field combining: genetics, molecular biology, robotics and computing.

**1990-2001 Human Genome Project (HGP)**: sequenced 3 billion base pairs, high quality (<1 error per 10k bases)

**Omics expansion**: functional genomics (gene expression and function), transcriptomics (RNA profiles), proteomics (protein complement), metabolomics (metabolites), phenomics (phenotype data)

**diversity of genomes**: 3 domains of life: **bacteria** (simple but share some genes with humans), **eukaryota** (complex genomes-->nuclear + organellar), **archea** (extreme environments, thermophiles, helophiles and methogens).

**types of genomes**: viruses (DNA or RNA genomes), prokaryotes (circular genomes + plasmids no reference genome concept), eukaryotes (nuclear genome the main DNA, mitochondrial 60 kb, chloroplast).

**data and databases**: genebank repository of dna sequences with metadata (species, source, annotations); WGS project randomly sequence dna fragments and assemble later; comparative genomics compare annotated and unannotated genomes to infer gene functions.

**big data in genomics**: NGS has made genomic data astronomical: Terabyte (TB) = 10^12 bytes, Petabytes (PB) = 10^15 bytes, Exabytes (EB) = 10^18 bytes, Zettabytes (ZB) = 10^21 bytes (1 million TB)
Challenges: Storage, Transfer, Analyses.

**metadata**=contextual information accompanying genomic data.

---

## Branches of Genetics – Conceptual Links

| Branch                  | Core Idea / Remember Link                                                                                                      |
|--------------------------|-------------------------------------------------------------------------------------------------------------------------------|
| **Classical genetics**   | **Mendel → Laws (Segregation, Independent Assortment)** → inheritance of traits → chromosomes as carriers → crossing-over → genetic maps → pedigrees → sex determination systems |
| **Molecular genetics**   | **DNA → Information flow (Replication → Transcription → Translation)** → gene regulation → PCR amplifies → Sanger reads sequence → foundation of modern sequencing & recombinant DNA |
| **Population genetics**  | **HWE (p² + 2pq + q² = 1) → baseline equilibrium** → deviations = evolutionary forces → LD (non-random alleles) → recombination (cM) → inbreeding (Fped, FROH, Fis) → ROH = genomic inbreeding signal → F_ST = differentiation among populations |
| **Quantitative genetics**| **Phenotype (P) = Genotype (G) + Environment (E)** → variance decomposition (A+D+I+E) → heritability (h² = VarG/VarP) → QTL detect loci → polygenic traits (continuous, discrete, binary) → applied to breeding & complex disease |
| **Genomics**             | **Genome = all DNA** → Human Genome Project (HGP, 1990–2001) → sequencing technologies (NGS) → big data challenges (TB → ZB) → multi-omics (functional, transcriptomics, proteomics, metabolomics, phenomics) → comparative genomics & databases (GenBank, WGS) |

---

## 2. Sequencing Technologies

* **NGS next generation sequencing** is the main technology in modern genomics, enabling the production and analysis of massive amounts of sequencing data; requires consideration of error rates throughput costs and timelines and provides unbiased sequencing of millions of DNA fragments simultaneously.
  -many DNA targets sequenced in parallel
  -millions of fragments analyzed in a single run
  -uses clonal PCR amplification to boost signal detection
  -sequences are unknown in advance identified only after computational analysis.
  
* **Sanger (1st generation)**: accurate, low throughput; one target DNA-->one PCR-->one sequencing reaction.
  -each capillary= one lane max 58 sequences/run
  -targeted: sequence is know in advance (via specific primers)
  -expensive, accurate, low throughput.
  
LINK: from 'we know what we sequence (SANGER)' to 'we discover what we sequenced (NGS)'.

**data production**: modern sequencers up to 48 samples/run; run duration is 1.5 hours and read length 1.5 kb per sample.
Every sequencing experiments involves 3 stages: Planning, data production and analyses.
In sanger era the data production was the most expensive and limiting factor and few sequences produced long timelines and limited insights
In present NGS data production is cheap and fast planning is crucial to optimize study deisgn and analysis is the bottleneck.

Ex. HGP 100,000,000 per genome today--> <1000 per genome, cost per 1 Mb DNA <0,01 
-->NGS has dramatically outspaced MOOre's law in reducing costs.

---
* **Ion Torrent sequencing system**: detects H+ release, ionograms.
  Ion torrent (now Thermo Fisher Ion S5 system) is a short-read NGS technology based on semiconductor sequencing.
  
---
* **Illumina**: short, accurate reads, cluster generation, paired-end.
---
* **454 Roche & ABI SOLiD**: early NGS, discontinued.
---
* **PacBio (SMRT, HiFi)**: long accurate reads.
---
* **Nanopore (MinION, PromethION)**: ultra-long reads, real-time.
---


## 3. Genome Assembly

* **Shotgun sequencing**: random fragmentation.
* Algorithms: Greedy, OLC, de Bruijn graphs.
* **k-mers**: basis of assembly and genome size estimation.
* **Scaffolding** with mate-pairs or long reads.
* **N50, coverage, BUSCO**: assess assembly quality.

---

## 4. Genome Annotation

* **Repeat annotation**: RepeatMasker, Dfam, TEannot.
* **Gene models**:

  * Ab initio (AUGUSTUS).
  * Homology-based (BLAST, Exonerate).
  * Transcript evidence (RNA-seq).
  * Integrative tools (MAKER, BRAKER2).
* Functional annotation: GO, KEGG, domains.

---

## 5. Data Formats & Tools

* **FASTQ**: raw reads + quality.
* **SAM/BAM**: alignments, CIGAR strings.
* **VCF**: SNPs, indels, SVs.
* **GFF/GTF, BED**: annotations & intervals.
* QC: FastQC, trimming (Trimmomatic).
* Alignment: BWA-MEM; Variant calling: GATK.
* Visualization: IGV, Ensembl.

---

## 6. Specialized Sequencing

* **aCGH**: CNV detection via hybridization.
* **Pool-seq**: allele frequencies in populations.
* **Targeted sequencing**: amplicon/hybrid capture panels.
* **Methyl-seq / Bisulfite**: DNA methylation at base resolution.
* **RNA-seq**: transcriptome profiling, isoforms, expression.

---

## 7. Applied Genomics

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
