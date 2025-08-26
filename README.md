# Applied Genomics – Slide Recap

Full slide deck: [Google Drive](https://drive.google.com/file/d/1tV58Ldxbase2jpN9sgMsIyVQvUfpziG5/view?usp=sharing)

This document summarizes the key concepts from **README1–12**, integrating branches of genetics, molecular tools, sequencing technologies, genome assembly, annotation, population applications, and applied genomics.

---

**Genetics** is the discipline that studies **biological diversity at different levels: genomes, populations, individuals and organisms**.
**Diversity** is the **driving-force** of genetics.

---

## 1. Branches of Genetics

## **Classical genetics (transmission/formal)**: Mendelian inheritance, pedigree analysis, laws of segregation & independent assortment.

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

## **Molecular genetics**: DNA as genetic material, PCR, Sanger sequencing.
  
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

## **Population genetics**: Hardy–Weinberg equilibrium, allele/genotype frequencies, LD, inbreeding, ROH.
  
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

Population applications:  
- **Genotyping (SNP arrays, GBS, RAD)** → captures genome-wide variants; essential for diversity studies, structure, and as input for GWAS. Arrays = known variants (cheap, reproducible), GBS/RAD = discover SNPs in non-model species.  
- **Structural variants (CNVs, inversions, translocations)** → large-scale changes often affect gene dosage or regulation; critical in cancer genomics, evolution, and adaptation. Detected with arrays, short reads (depth, discordant pairs), or more effectively with long reads.  
- **GWAS** → links genotypes to phenotypes without prior hypotheses. Pipeline: QC → correction for structure (PCA, kinship) → association → multiple testing correction → visualization (Manhattan). Strength: unbiased trait discovery; Limits: rare variants hard to detect, risk of stratification bias.  
- **MAS / MAB** → bridges population data to breeding. MAS selects individuals with markers linked to traits; MAB/genomic selection uses genome-wide SNPs + predictive models to estimate breeding values (GEBVs), speeding up improvement of polygenic traits.  
- **Conceptual link to population genomics** → LD underlies SNP associations; ROH and inbreeding shape background signals; F_ST reflects differentiation and thus resolution of associations across populations.  

---
       
## **Quantitative genetics**: phenotype = genetics + environment, variance decomposition, heritability, QTL.
  
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

## **NGS next generation sequencing** is the main technology in modern genomics, enabling the production and analysis of massive amounts of sequencing data; requires consideration of error rates throughput costs and timelines and provides unbiased sequencing of millions of DNA fragments simultaneously.
  -many DNA targets sequenced in parallel
  -millions of fragments analyzed in a single run
  -uses clonal PCR amplification to boost signal detection
  -sequences are unknown in advance identified only after computational analysis.
  
## **Sanger (1st generation)**: accurate, low throughput; one target DNA-->one PCR-->one sequencing reaction.
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

| Aspect                  | **Sanger Sequencing (1st Gen)**                              | **NGS (Next Generation Sequencing)**                            |
|--------------------------|-------------------------------------------------------------|-----------------------------------------------------------------|
| **Concept**             | *We know what we sequence* (target defined by primers)       | *We discover what we sequenced* (massively parallel, unbiased)  |
| **Throughput**          | Low: one DNA → one PCR → one sequencing reaction; max ~96 capillaries/run (~58 seqs/run) | High: millions of fragments sequenced in parallel in a single run |
| **Input / Target**      | Specific known sequence (targeted, primer-based)             | Random, genome-wide, no prior knowledge needed                  |
| **Signal amplification**| PCR product sequenced directly                               | Clonal amplification (emulsion PCR, bridge PCR) → boosts signal |
| **Data production**     | Few sequences, expensive, time-consuming                     | Dozens of samples/run, fast (~1.5h), cost-efficient, high yield |
| **Read length**         | ~800–1000 bp per read                                        | Short reads (100–400 bp) or long reads (PacBio/Nanopore)        |
| **Error rate**          | Very low (high accuracy)                                     | Higher error rate (depends on platform), computational correction needed |
| **Main limitation**     | Costly, low throughput, limited discovery power              | Analysis bottleneck: huge data → requires planning & computation |
| **Era bottleneck**      | Sanger: **data production** (slow, costly, few insights)     | NGS: **analysis** (planning crucial, data cheap & fast)         |

---

## **Ion Torrent sequencing system**: detects H+ release, ionograms.
  Ion torrent (now Thermo Fisher Ion S5 system) is a **short-read NGS technology based on semiconductor sequencing**.
  Instead of using fluorescence or optics, the system detects pH changes caused by nucleotide incorporation.
  Sequencing occurs in tiny wells on a silicon ship, each acting as an independent reaction chamber.
  Each well= one DNA fragment per thousands of clonal copies (to amplify signal)

  -->**chips**: the chip is the core of the ion torrent system; each well = one sequencing reaction with a pH sensor to detect ions; dna fragments are loaded into wells         and sequenced in parallel; different chip models= different throughput scalability (510 2,5 Gb output 2.3 million reads).
  -->**sequencing principle**: dna fragments are engineered with universal primer sites, one nucleotide type is flowed across the chip at a time. If complementary the           polymerase adds it to the DNA strand releasing H+ ions detected as pH change and pyrophosphate, each pH change=incorporation event converted into digital signal.           Signal intensity is proportional to the number of identical nucleotides added (homopolymers)
  -->**workflow**: DNA fragmentation (sonication or enzymatic)-->library prep (adapters + barcodes)-->emulsion PCR(clonal amplification)-->chip loading (dna loaded beads        into wells)-->sequencing(ion detection via pH sensors)-->data output (torrent suite software for qc and read statistics).
  -->**library and template prep**: dna fragmentation require sonication (physical) or enzymatic digestion, then library preparation require fragments engineered with           adapters + universal primer sites, barcoding allows multiplexing samples from multiple individuals and clonal amplification performed by emulsion PCR (one fragment         DNA per micro-droplet), creates thousand of copies of the same DNA fragment and the goal is 1 dna fragment per sphere avoids polyclonals (mixed signals).
  -->**flow and sequencing**: sequencing occus by programmed flows of nucleotides; homopolymer regions (AAA, TTT) signal proportionality problem higher error rate; output       is displayed in **ionograms** similar to elctropherograms.
  -->**strengths**: fast and relatively cheap (no fluorescently modified nucleotides and block), simple chemistry, highly scalable with different chips and support             barcoding (cost-efficient multiplexing)
  -->**limitation**: error-prone in homopolymers regions (cannot distinguish 3 vs 4 identical bases), require careful DNA quantification to avoid polyclonals (multiple         template per sphere and empty wells, chips are single-use, read-length limited to 200-400 bp.
  
| Aspect                  | Key Points / Keywords                                                                 |
|--------------------------|---------------------------------------------------------------------------------------|
| **Technology**           | Short-read NGS, **Sequencing by Synthesis (SBS)** → detects **H⁺ ions (pH change)** instead of fluorescence/optics |
| **Chip**                 | Silicon chip with wells = independent reaction chambers; each well = 1 DNA fragment (thousands of clonal copies); pH sensor detects incorporation; different chip models = scalable throughput (e.g. 510 chip ≈ 2.5 Gb, ~2.3M reads) |
| **Sequencing principle** | DNA fragments with universal primers; nucleotides flowed sequentially; if incorporated → polymerase adds base → H⁺ released → pH change detected; **signal intensity ∝ # of incorporated bases (homopolymers issue)**; output = **ionograms** |
| **Workflow**             | DNA fragmentation (sonication/enzymatic) → library prep (adapters + barcodes) → **emulsion PCR (clonal amplification)** → chip loading (DNA beads into wells) → sequencing (pH detection, SBS method) → data analysis (Torrent Suite software) |
| **Library/Template prep**| DNA fragments + adapters + universal primer sites; barcoding for multiplexing; clonal amplification by emulsion PCR (1 fragment per droplet → 1 bead with thousands of identical copies); avoid polyclonals/empty wells |
| **Flow & sequencing**    | Programmed flows of single nucleotide types → incorporation → H⁺ release → signal; **homopolymer regions (AAA/TTT)** cause signal proportionality error; output visualized as ionograms (similar to electropherograms) |
| **Strengths**            | Fast, relatively cheap, simple chemistry (no fluorescent dyes/blockers), scalable with different chips, barcoding enables multiplexing |
| **Limitations**          | Homopolymer error (3 vs 4 bases indistinguishable), requires accurate DNA quantification, risk of polyclonals, single-use chips, short read length (200–400 bp) |

-->**Electropherogram**= fluorescence colors → high accuracy → low throughput
-->**Ionogram** = pH intensity → high speed/parallelism → homopolymer error risk

---

## **Illumina SBS**: short, accurate reads, cluster generation, paired-end.
  Illumina is the **leading NGS platform** in genomics today. It is based on **Sequencing by Synthesis (SBS)** where nucleotides are incorporated one at a time detected by   fluorescence and then reset for the next cycle.
  -produces short-reads (100-300 bp)
  -extremely high accuracy (Q30-99.9%)
  -scalable with different machines and flow cells
  -->**workflow**: 1-library prep: DNA fragmentation by physical sonification or enzymatic digestion; adapters are added on DNA fragments enable binding to flow cell            oligos, contain barcodes for multiplexing and provide primer binding sites.
     2-cluster generation (bridge amplification): performed on a flow cell (glass slide coated with oligonucleotides) each fragment binds the surface and forms a bridge,        PCR amplification occurs creating thousand of clonal copiesd of each fragments and each cluster = one DNA fragment amplified, ordered clusters in modern flow cells --      >higher accuracy and resolution.
     3-sequencing by synthesis: reversible terminator nucleotides each base has a fluorescent label + chemical blocking group, only one nucleotide is incorporated per           cycle. Workflow per cycle: -added 4 fluorescently labeled nucleotides; -DNA polymerase incorporates the complementary base; imaging system captures the fluorescent         color of each cluster; -blocking group is removed and synthesis can continue; -cycle repeats. Output: fluorescence signals translated into read sequence.
  -->**chemistry and improvements**: 4 channel SBS classic method-->4 fluorochromes one per nucleotide, 2 channel SBS uses only 2 dyes for 4 bases reduced cost/time, 1          channel SBS uses one dye with chemical tricks (on/off states) to distinguish nucleotides. The advantages of SBS chemistry is to prevent homopolymers errors, high           precision because reaction is stopped and read at each cyle. Trade-off: chemistry is more expansive due to modified nucleotides, sequencing is slower because of stop-      read-reset-cycles.
  -->**accuracy and read length**: error rate extremely low, gold standard for WSG, exome sequencing, RNA-Seq and small variant detection (SNPs, indel).
  -->**advantages**: very high accuracy, wide range of applications, scalable throughput (different flow cells and instruments), multiplexing via barcodes.
  -->**limitationss**: more expansive chemistry (modified nucleotides), short-reads compared to long-read technologies, slower per cycle (due to imaging and chemical reset      steps).
  -->fluorescent label attach to the base, and blocking group in the 3'.

| Aspect                  | Key Points / Keywords                                                                 |
|--------------------------|---------------------------------------------------------------------------------------|
| **Technology**           | Leading NGS platform, **Sequencing by Synthesis (SBS)** → nucleotides incorporated one at a time, detected by **fluorescence**, reset each cycle |
| **Read type**            | Short reads (100–300 bp), **paired-end** option, Q30 ≥ 99.9% (extremely accurate)     |
| **Workflow**             | 1. **Library prep** → DNA fragmentation (sonication/enzymatic) + adapters (barcodes, primer sites) <br> 2. **Cluster generation** → bridge amplification on flow cell, clonal copies per cluster <br> 3. **SBS sequencing** → reversible terminator nucleotides, fluorescence detection per cycle, blocking group removed, repeat |
| **Chemistry**            | Classic: **4-channel SBS** (4 dyes, one/base) <br> 2-channel SBS: 2 dyes → 4 bases <br> 1-channel SBS: 1 dye with on/off chemical states <br> **Advantage**: prevents homopolymer errors, precise base-by-base reading <br> **Trade-off**: chemistry expensive, slower cycles (stop → read → reset) |
| **Accuracy & length**    | Very low error rate, gold standard for WGS, exome, RNA-Seq, SNP/indel detection       |
| **Advantages**           | Ultra-high accuracy, wide applications, scalable throughput (machines & flow cells), supports multiplexing via barcodes |
| **Limitations**          | Expensive chemistry (modified nucleotides), short reads (assembly harder than long-read tech), slower cycles due to imaging & reset |

-->**Illumina Paired End Sequencing**: Illumina sequencing is usually short-read sequencing. To overcome the limitation of sequencing only small fragments paired-end sequencing allow reading both ends of a DNA fragment.
Ex. 1000 bp fragmented DNA illumina paired end sequence first 200 bp and last 200 bp the middle 600 bp remain unsequenced directly but are physically linked.
2 short reads separated by an insert size (unsequenced gap).
Applications: genome assembly (de novo, resequencing), structural variation (reveal translocations, inversion, gene fusions), transcriptomics (splicing events),improved reconstruction of missing sequences.
-->**Mate-Pair Sequencing (Illumina special library prep)** → Designed for longer insert sizes (2–20 kb). DNA is circularized, fragmented, and sequenced at the junctions. As a result, the two reads come from far apart in the original genome (kb distance), not just hundreds of bp like paired-end.
Applications: de novo genome assembly scaffolding, resolving repetitive regions, detecting large structural rearrangements, mapping across long gaps.
-->Link between them:
Paired-end = short insert (~200–600 bp) → precise, high-resolution, local context.
Mate-pair = long insert (2–20 kb) → long-range linking information, helps scaffold assemblies and detect large structural variation.

---

## **454 Roche & ABI SOLiD**: early NGS, discontinued.
  **ABI SOLiD**: sequencing by ligation system, It identifies the sequence by matching and ligating pre-made probes with ligase.
    Its key defining feature was its use of DNA ligase instead od DNA polymerase to determine the sequence, and its unique two-base encoding system for color space. The        throughput is impressive but have complex data analysis and short-reads (50-75 bp).
    -->Step: emPCR; bind beads to a glass slide; hybridize and ligate fluorescent probes; detect color from ligation.
  **454 Roche**: SBS method, It build the DNA strand by adding nucleotides with polymerase.
    Ion Torrent is to faster and significantly cheaper as it eliminated the need for expansive enzymatic reagents and complex optical imaging systems.
    -->Step: emPCR; load beads into Picotiterplate wells; flow nucleotides; detect light from incorporation.
    -->the detection method is pyrosequencing PPi (a biochemical reaction that produces light)


| Platform        | Principle / Chemistry                                                                 | Workflow Steps                                                                                          | Features & Read Length           | Pros                                                                 | Cons                                                                 |
|-----------------|---------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|----------------------------------|----------------------------------------------------------------------|----------------------------------------------------------------------|
| **ABI SOLiD**   | **Sequencing by Ligation**: uses **DNA ligase** + fluorescent probes; unique **two-base encoding (color-space)** | 1. emPCR (clonal amplification on beads) <br> 2. Beads immobilized on slide <br> 3. Fluorescent probes hybridized & ligated <br> 4. Color signals detected | Short reads **50–75 bp**         | High throughput, ligation-based accuracy (error-correction via color-space) | Complex data analysis (color-space), short reads, platform discontinued |
| **454 Roche**   | **Sequencing by Synthesis (SBS)** via **pyrosequencing** (detects light from PPi release) | 1. emPCR (DNA on beads) <br> 2. Load beads into **PicoTiterPlate wells** <br> 3. Flow nucleotides sequentially <br> 4. Incorporation → PPi release → luciferase reaction → light detected | Reads up to **400–700 bp**       | Longer reads than Illumina at the time; fast compared to Sanger; scalable with picotiterplate | Expensive reagents (luciferase/enzymes), homopolymer errors, high cost per base, discontinued |

---

**Long-read sequencing**: while Illumina dominates short-read sequencing two main technologies provide long reads (up to tens of kb or more): ONT and PAcBio.
Long read sequencing helps resolve: repetitive regions, structural variations, complex genome assemblies, phasing of haplotypes.

## **PacBio (SMRT, HiFi)**: long accurate reads.
  -->Principle: a single DNA molecule is immobilized in a well (Zero-Mode Waveguide, ZMW). DNA polymerase incorporates fluorescently labeled nucleotides. Each                   incorporation event emits light captured by a camera. Interpulse duration, time between signals, helps determine the base.
  -->**CCS Circular Consensus Sequencing**: DNA fragments are circularized SMRTbell templates; sequenced multiple times as polymerase loops around the circle; errors are        random corrected by overlapping multiple passes.
  -->Features: read length typically 10-25 kb, can reach 50 kb. Error correction: high accuracy with CCS 8HiFi reads). Bias reduction: no issues with GC-rich regions.           Phasing: can separate maternal and paternal haplotypes.
  -->Pros: accurate long reads (HIFi reads with Q30+); excellent for strctural variants and complex genomes; NO GC-bias; phasing of alleles possible.
  -->Cons: high cost (3-4x Illumina per genome); lower throughput than short-read systems; library prep is complex and requires high-quality non degraded DNA; larger lab-       based instruments.
  -->fluorescent label on the phosphate chain (not on the base)
  
| Aspect                  | Key Points / Keywords                                                                 |
|--------------------------|---------------------------------------------------------------------------------------|
| **Principle**            | **Single-Molecule Real-Time (SMRT)**: DNA immobilized in **Zero-Mode Waveguide (ZMW)** wells; polymerase incorporates fluorescent nucleotides → light pulses detected; interpulse duration aids base calling |
| **CCS / HiFi**           | **Circular Consensus Sequencing (CCS)**: DNA circularized into SMRTbell templates → polymerase makes multiple passes → random errors corrected → **HiFi reads (Q30+, high accuracy)** |
| **Read length**          | Typical: **10–25 kb**; up to **50 kb** possible                                      |
| **Unique features**      | Random error correction with CCS; no GC-bias; **haplotype phasing** (maternal vs paternal separation) |
| **Pros**                 | Accurate **long reads** (HiFi Q30+); excellent for structural variants & complex genomes; unbiased (GC-rich ok); haplotype phasing |
| **Cons**                 | High cost (3–4× Illumina); lower throughput; requires high-quality DNA; complex library prep; large lab instruments |



## **Nanopore (MinION, PromethION)**: ultra-long reads, real-time.
  -->Principle: sequencing native DNA/RNA strands (no synthesis step). A single stranded DNA passes through a nanopore embedded in a membrane. Electric current runs across      the pore each nucleotide or combination of bases causes a specific voltage disruption. Volatge shifts are translated into base calls.
  -->Features: read length up to 10-20 kb in practice, theoretical maximum >1 Mb. Portability devices like MinION are samll, USB-powered and usable in the field. Direct         sequencing of RNA (no cDNA needed).
  -->Pros: very cheap instrument, portable usable on-site, ultra long reads possible, sequencing of RNA directly
  -->cons: high error rate (5% historically 20%), throughput less stable, library prep still required and not suitable for degraded samples.

| Aspect                  | Key Points / Keywords                                                                 |
|--------------------------|---------------------------------------------------------------------------------------|
| **Principle**            | Native **DNA/RNA strands** pass through nanopore in a membrane; electric current across pore → each nucleotide causes **voltage disruption** → signal translated into bases |
| **Read length**          | Practical: **10–20 kb**; Theoretical: **>1 Mb (ultra-long reads)**                     |
| **Devices**              | **MinION** (USB, portable, field use) <br> **PromethION** (high-throughput, lab-scale) |
| **Unique features**      | Direct RNA sequencing (no cDNA needed); real-time data streaming                       |
| **Pros**                 | Very low-cost device, portable (field sequencing), ultra-long reads, RNA direct sequencing |
| **Cons**                 | Higher error rates (~5–20%), throughput less stable, library prep still required, not good for degraded samples |


---

## 3. Raw data processing and file formats

**NGS** produces raw signal signal data that must be converted into reads with per-base quality scores.
Step include:
*-1 Basecalling*: conversion of raw signals into nucleotides
*-2 Read filtering*: remove poor-quality sequences
    Not all reads are equally useful. Filtering is so crucial: remove polyclonal reads, discarded reads<25 bases too short for reliable alignment, trim low-quality regions     (sliding windows/moving windows) and adjust parameters for low-complexity libraries.
*-3 Alignment*: mapping reads to a reference genome: once sequencing data (FASTQ) are aligned against a reference genome teh results are stored in SAM and BAM files.
*-4 Variant Calling*: identifying SNPs, indels or structural variants: After sequencing, alignment, QC, and filtering, the next step is variant calling. The goal is to         identify genetic differences (SNPs, indels, SVs) compared to a reference genome. Main software: GATK (Genome Analysis Toolkit) 

---

* **FASTA** format: plain text with nucleotide sequence (ACGT..) + identifier. No quality information.
  
* **FASTQ**: raw reads + quality-->FASTA + Phred quality scores per base
    Four lines per entry: Sequence ID (header); Nucleotide sequence; Separator (+); Quality scores (ASCII symbols)
  **quality scores (Phred scale)**: Q=-10log10(e) where Q stay for quality score and e for probability of incorrect base call.
  ex. Q20-->1 error in 100 bases (99% accuracy); Q30-->1 error in 1000 bases (99.9% accuracy)
  Ion torrent error rates Q20 due to homopolymers improvinf towards Q30.
  FASTQ-->unaligned reads.

  **Quality control (QC)** is a fondamental step in NGS workflows: ensures that the raw data (FASTQ...) are reliable; identifies biases introduced by library preparation,    sequencing platform or sample quality; prevent False positives in downstream variant discovery.
  Why matters: poor-quality reads, GC biases misinterpretation of gene regions, duplicates..
  **Tool for QC**: FastQC java-based, modular, HTML report, and Prinseq efficient for trimming and filtering.

  **FASTQC** -->provides modular analyses: import data from BAM, SAM, FASTQ, quick overview of potential problems, summary graphs, export results and tables and work         offline.
  -Basic statistic (sequence count, length distribution)
  -Per-Base sequence quality (signal): boxplots per position along the read, aggregated phred scores show accuracy, Q>=20 and green zone good orange warning; quality drops    towards the end of reads
  -Per-sequence quality score:in the X Phred in the Y the number of reads, shows distribution of read-level quality, reads with mean Q<20 should be discarded.
  -Per-base sequence content: shows frequency of ATCG at each position; stable percentages across read length; fluctuantions at the start/end: sequencing errors or library    prerp bias.
  -GC content: distribution of GC content per read; compared against theoretical distribution of the target genome; deviations indicate biases in library prep and             contaminants; if multiple peaks-->contaminants or mixed DNA sources.
  -Sequence duplication levels: most reads are unique (expected); high duplication bias in library prep; amplicon sequencing duplication expected; duplicated reads can        bias varaint calls, must be removed; after deduplication most reads should be unique.

  **Trimming**-->removing low-quality bases increases dataset reliability. 2 main approaches:
  -Threshold-based trimming: define Q threshold (Q20), removes bases below threshold until high quality base reached
  -Window-based trimming (preferred): define window size (5 nt), calculate average quality in the window, if average<threshold you identify trim region, retains more nt by    smoothing local fluctuations.

  *PIPELINE OF QC*-->run FASTQC; apply trimming/filtering; remove low-quality bases and duplicated reads; re-run FastQC; only then proceed with alignment and variant         discovery
  
---

* **SAM/BAM**: alignments, CIGAR strings.
  **BAM** stay for binary alignment map. compressed binary format with both reads + alignment information; requires index file; not human readable. It is the binary          version of SAM.
  **SAM** stay for sequence alignment map: a text based format, human-readable.
  A SAM file has 2 sections:
  **a) Header**
  Lines start with @.
  Contains metadata about alignment:
  Reference sequence IDs & lengths.
  Alignment program & command line (@PG).

  **b) Alignment Section**
  Each read is represented by a row with 11 mandatory fields:
  QNAME → read identifier.
  
  FLAG → integer encoding alignment properties.
         The FLAG encodes properties of each read as integers.
         Examples: 4 → read is unmapped. 256 → secondary alignment. 1024 → PCR duplicate.
         Multiple values can be combined.
  
  RNAME → reference sequence name (chromosome).
  POS → starting position on reference.
  MAPQ → mapping quality score.
  
  CIGAR → compact description of alignment.
          CIGAR = Compact Idiosyncratic Gapped Alignment Report. Describes how each read aligns to the reference.
          Common symbols: M → match/mismatch (aligned). I → insertion (present in read, absent in reference). D → deletion (absent in read, present in reference).
                          S → soft clipping (part of the read not aligned).
  
  RNEXT → mate/next read info.
  PNEXT → position of mate/next read.
  TLEN → observed template length.
  SEQ → read sequence.
  QUAL → base quality (Phred).

  
  BAM and SAM are standard formats in genomics and widely used in downstream analysis.

  *PIPELINE*:start form FASTQ file-->align read to reference genome with BWA-MEM-->the output is the BAM file (compressed alignment)-->BAM can be converted back to SAM.
  
  **Quality control** continues after alignment:
  -Mapping Quality (MAPQ)-->Probability that read is correctly mapped.
   MAPQ = 0 → read could map to multiple locations (e.g., repetitive elements, assembly errors).
  Reads with low MAPQ often discarded.
  -Duplicate Removal-->Duplicates come from PCR amplification.
  They don’t add information, but inflate coverage.
  Can create false positives in variant calling.
  Tools: Picard for duplicate marking/removal.
  PCR-free library preparation avoids this issue (requires more input DNA).
  
---

**Calling a variant** requires careful evaluation of multiple factors: Base call quality of each supporting base. Mapping quality (MQ) of aligned reads. Sequencing depth → minimum number of reads supporting the variant. Proximity to indels or homopolymer runs (error-prone regions).
Single-sample vs multi-sample calling. More samples = higher confidence.

**variant calling strategies**: 
a) Per-sample calling
   Each sample analyzed separately.
   Produces one VCF file per sample.
b) Joint calling--> All samples analyzed together. Produces one merged VCF file.
   Advantages:
   Increases statistical power.
   Better detection of heterozygous variants.
   Reduces false negatives from low coverage.
   Recommended: joint calling for population studies.
  
* **VCF**: SNPs, indels, SVs.
  VCF stay for variant call format: text format for variants (SNPs, indels, structural differences)
  A VCF file has two sections:
  a) Header
  Lines start with ##.
  Define metadata, filters, INFO fields, FORMAT descriptors.
  Last header line starts with #CHROM → column names.
  
  b) Variant Records
  Each row represents one variant.
  Columns include:
  CHROM → Chromosome ID.
  POS → Position of variant.
  ID → Variant ID (e.g., dbSNP rs number, or . if novel).
  REF → Reference allele.
  ALT → Alternative allele(s).
  QUAL → Phred-scaled variant quality score.
  FILTER → PASS/FAIL (e.g., “PASS”, “q10”).
  INFO → Additional annotations (DP = depth, AF = allele frequency, etc.).
  FORMAT → Defines per-sample genotype fields (e.g., GT, DP, GQ).
  10+. Samples → Genotypes and metrics for each sample.
  Genotype Codes
  0|0 → homozygous reference.
  0|1 or 1|0 → heterozygous.
  1|1 → homozygous alternate.
  If multiple alternate alleles exist: 2, 3, etc.

  Visualization: **IGV (Integrative Genomics Viewer)** → gold standard for inspecting alignments.
  Inputs: reference genome (FASTA), alignment (BAM), and variants (VCF).
  Features:
  Visualize reads aligned at a given locus.
  Detect true variants vs sequencing errors.
  Depth bar shows number of supporting reads.
  Color-coding distinguishes reference vs alternate alleles.

  Other Tools for Variant Annotation:
  **Ensembl Genome Browser**-->Allows visualization of genomes, genes, transcripts, proteins, and variants.
  Example: searching the KMO gene gives: Chromosomal location. Gene description. Transcript IDs and number. Encoded protein + UniProt link. Variant table for known           polymorphisms (dbSNP, Ensembl).

  **Variant Effect Predictor (VEP)**-->Provided by Ensembl.
  Takes VCF as input.
  Outputs: Genomic location of variants. Variant type (missense, synonymous, frameshift…). Functional consequences on gene/protein.
  **dbSNP**-->Repository of known SNPs and small indels. Variants deposited with IDs (rsXXXX). Can be cross-referenced to validate if a variant is novel.


  **Types of variants**-->
  **Small Variants**: SNPs (Single Nucleotide Polymorphisms), Substitution of one base with another.Must occur in ≥1% of the population → polymorphism.
  If <1%, it is considered a mutation.
  Indels: Small insertions or deletions.
  Substitutions: Multiple bases replaced by different ones.

  **Structural Variants**--> CNVs (Copy Number Variants) → duplicated or deleted regions.
  Inversions → DNA region flipped in orientation.
  Translocations → block of DNA moved to another chromosome or region.

---
  
* **GFF/GTF, BED**: annotations & intervals.
  BED stay for browser extensible data: tab-delimited text format defining genomic features.
  **GFF3/GTF** = rich annotation format (genes, transcripts, exons, CDS, regulatory elements).
  -->Purpose: Stores annotations (genes, exons, CDS, regulatory features, etc.) on a reference genome.
     Format: 9 mandatory tab-delimited columns
  seqid → chromosome or scaffold name
  source → annotation source (e.g. Ensembl, maker, augustus)
  type → feature type (gene, exon, CDS, mRNA, repeat…)
  start → start coordinate (1-based)
  end → end coordinate
  score → numerical value (or “.” if not used)
  strand → + or – (strand of feature)
  phase → for CDS: 0, 1, or 2 (frame of translation start)
  attributes → semicolon-separated key=value pairs (e.g., ID=gene1;Name=BRCA1)
  
  **BED** = simple intervals (coordinates, lightweight for browsers like UCSC/IGV).
  -->Purpose: Defines genomic intervals (regions of interest).
  Format: Very light, 3 mandatory columns, max 12 optional.
  Optional (up to 12 fields): name, score, strand, thickStart, thickEnd, itemRgb, blockCount, etc.

Raw sequencing output can be huge (terabytes); filtering and compression reduce size dramatically (2.5 TB raw data-->30 GB FASTQ); effcient storage and file format choice are critical for downstream analysis.

| Format      | Purpose / Content                                     | Structure / Key Fields                                                                 | Example (simplified)                                   |
|-------------|-------------------------------------------------------|----------------------------------------------------------------------------------------|--------------------------------------------------------|
| **FASTA**   | Stores raw nucleotide or protein sequences (no quality) | `>identifier` + sequence lines                                                          | `>seq1` <br> `ATGCCGTA...`                             |
| **FASTQ**   | Raw reads + **per-base quality scores**               | 4 lines/entry: <br> 1. `@ID` <br> 2. sequence <br> 3. `+` <br> 4. quality (ASCII, Phred) | `@read1` <br> `ATGCC` <br> `+` <br> `IIIII`            |
| **SAM**     | Text format: aligned reads to reference genome        | Header (`@SQ`, `@PG` …) + alignment section (11 mandatory fields: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL) | `read1  0  chr1  100  60  50M  *  0  0  ATGCC...  IIIII` |
| **BAM**     | Binary compressed version of SAM                     | Same fields as SAM but **binary (indexed)** for speed & storage                          | (binary, not human-readable)                          |
| **VCF**     | Variant calls: SNPs, indels, SVs                     | Header (`##` metadata, `#CHROM` cols) + records (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, samples) | `chr1  123456  rs123  A  G  99  PASS  DP=30  GT:DP  0/1:15` |
| **GFF3**    | Genome annotations: genes, exons, CDS, regulatory features | 9 tab-delimited cols: seqid, source, type, start, end, score, strand, phase, attributes | `chr1 Ensembl gene 1000 5000 . + . ID=gene1;Name=BRCA1` |
| **BED**     | Genomic intervals (lightweight)                      | 3 required cols: chrom, start, end (up to 12 optional: name, score, strand, color, blocks) | `chr1  1000  5000  BRCA1`                              |

---

**Before sequencing: DNA Quality Assessment & Sequencing Coverage**-->

1. Quality Assessment of DNA
Before starting NGS library preparation, DNA quality must be carefully evaluated.
DNA purity and integrity affect downstream steps such as library prep, amplification, and sequencing efficiency.

**Purity Ratios (Spectrophotometry)**
A260/280 ratio
Expected: ~1.8
High values → RNA contamination
Low values → protein contamination, residual phenol, or very low DNA concentration (<10 ng/µl)
Influenced by pH (acidic solution ↓ ratio; basic solution ↑ ratio)
A260/230 ratio
Expected: 2.0 – 2.2
Low values → carbohydrate carryover (plants), residual phenol, guanidine, glycogen
High values → improper blank solution

**DNA Integrity**
Checked by agarose gel electrophoresis:
High-quality DNA → sharp band, high molecular weight (upper gel region)
Degraded DNA → smeared band or lower fragments
Low integrity does not prevent sequencing, but long-read sequencing (Nanopore, PacBio) requires intact long fragments.

2. Sample Requirements for NGS
Companies providing sequencing services usually require:
Condition: genomic DNA
Quantity:
≥ 1 µg (small fragment library)
1.5 µg (PCR-free library)
Concentration: ≥ 12.5 ng/µl
Purity: A260/280 between 1.8 – 2.0
If requirements fail, sequencing may still be possible, but with increased risk of poor results.

3. Sequencing Coverage
**Coverage (breadth)** → % of target bases sequenced at least “X” times.
Depth of coverage → average number of times each base is sequenced.
Depth of Coverage} = LN/G

Where:
L = read length
N = number of reads
G = haploid genome length
Expressed as “X-fold” (e.g., 10x, 20x, 40x).
Examples
5x coverage → each base sequenced ~5 times (low reliability, missing variants).
20x coverage → sufficient for variant discovery (SNPs/indels).
40x coverage → high-quality genome resequencing (each base seen ≥40 times).
Genome assembly → requires much higher depth and often long-read technologies.

4. Why Coverage Matters
Variant discovery → ≥20x needed for reliable SNP/indel detection.
Clinical genomics → high coverage required to avoid missing pathogenic variants.
De novo assembly → high depth + long reads needed to resolve repeats and structural variants.

---

## 4. Genome Assembly

Genome assembly is the process of reconstructing the complete genomic sequence of an organism from millions of sequencing reads. Modern sequencing platforms produce vast numbers of short or long fragments, which must be computationally pieced together to rebuild the genome.

### Shotgun Sequencing
The most common strategy is **shotgun sequencing**, where DNA is randomly fragmented into smaller pieces that are sequenced in parallel. This provides unbiased coverage of the genome but creates the challenge of correctly reassembling overlapping fragments.

### Assembly Algorithms
Different computational strategies are used depending on read length and sequencing technology:

| Feature                | **Overlap–Layout–Consensus (OLC)**                                | **de Bruijn Graph**                                       |
|-------------------------|-------------------------------------------------------------------|-----------------------------------------------------------|
| Input                   | Long reads (Sanger, PacBio, Nanopore)                            | Short reads (Illumina)                                    |
| Principle               | Reads compared pairwise for overlaps → graph of overlaps → consensus | Reads decomposed into **k-mers**; overlaps represented as graph edges |
| Pros                    | Accurate with long reads; good for smaller datasets              | Efficient with huge short-read datasets; scales well       |
| Cons                    | Computationally expensive (all-vs-all overlaps); not efficient for short reads | Struggles with sequencing errors, repeat resolution depends on *k* size |

### k-mers and Genome Size Estimation
**k-mers** (substrings of length *k*) are the foundation of short-read assembly. Counting the frequency of all k-mers across sequencing reads produces a distribution that can be used to estimate:
- **Genome size**:  
 Genome size = Total k-mers/Peak depth
- **Heterozygosity**: heterozygous genomes show two peaks (diploid k-mer distribution).  
- **Repeat content**: repetitive sequences distort the k-mer curve and increase multiplicity.  

### C-Value (Absolute Genome Size)
The **C-value** is the amount of DNA contained in a haploid nucleus (1C, e.g., in gametes). It represents the absolute genome size, measured in base pairs or picograms of DNA.  
- Example: the human haploid genome has a **C-value ≈ 0.978 × 10⁹ bp (~3.2 pg DNA)**.  
- C-value is species-specific and important in planning sequencing projects (coverage and depth requirements).  

### Scaffolding
Once contigs (continuous assembled sequences) are generated, **scaffolding** uses additional information to connect and order them:
- **Paired-end reads**: link short fragments separated by a few hundred bp  
- **Mate-pair libraries**: provide longer inserts (2–20 kb) for spanning repeats  
- **Long reads (PacBio, Nanopore)**: bridge large gaps and improve contiguity  

### Pipeline of Correct Genome Assembly
1. **DNA extraction & QC** → high-quality, high-molecular-weight DNA is critical  
2. **Sequencing** → short reads (Illumina) for accuracy; long reads (PacBio/Nanopore) for continuity  
3. **Preprocessing** → adapter trimming, quality filtering, error correction  
4. **Contig assembly** → via OLC (long reads) or de Bruijn graphs (short reads)  
5. **Scaffolding** → use mate-pairs, long reads, Hi-C for chromosome-level assembly  
6. **Polishing** → correct errors using high-accuracy reads (Illumina short reads on top of long reads)  
7. **Quality assessment** → N50, coverage, BUSCO gene completeness  
8. **Annotation** → identify genes, repeats, and functional elements  

### Assembly Quality Metrics

| Metric       | Definition & Details                                                                                                      |
|--------------|---------------------------------------------------------------------------------------------------------------------------|
| **N50**      | A contiguity statistic. All contigs/scaffolds are ordered from longest to shortest, and lengths are cumulatively summed until 50% of the total assembly size is reached. The length of the contig at this point is the **N50**. <br> Example: if genome = 3 Gb, and adding contigs from longest to shortest passes 1.5 Gb at contig length = 5 Mb, then **N50 = 5 Mb**. |
| **Coverage** | Average sequencing depth per base (e.g. 20×, 40×). Calculated as: <br> \[ \text{Coverage} = \frac{L \times N}{G} \] where *L* = read length, *N* = number of reads, *G* = genome length. Higher coverage improves accuracy and completeness. |
| **BUSCO**    | **Benchmarking Universal Single-Copy Orthologs**: searches for a set of evolutionarily conserved orthologous genes expected in a lineage. Results are reported as percentages of: <br> - **Complete (C)** → full-length ortholog found <br> - **Single-copy (S)** → present once <br> - **Duplicated (D)** → multiple copies detected <br> - **Fragmented (F)** → partial genes <br> - **Missing (M)** → not found <br> Example: BUSCO = 95% Complete (90% S, 5% D), 3% Fragmented, 2% Missing → indicates a very complete assembly. |

**In summary**, genome assembly integrates sequencing technologies, graph-based algorithms, k-mer analysis for genome size and C-value estimation, scaffolding strategies, and rigorous quality metrics (N50, coverage, BUSCO) to transform fragmented reads into a biologically meaningful genome sequence.

---

## 5. Genome Annotation

Once a genome has been assembled, the next step is **annotation**, the process of identifying and describing the functional elements within the sequence. Annotation transforms a raw collection of contigs and scaffolds into a biologically meaningful map of genes, repeats, regulatory regions, and other features.  

---

### Repeat Annotation
A large proportion of most eukaryotic genomes is composed of **repetitive elements**, including transposable elements (TEs), tandem repeats, and low-complexity regions. Detecting and masking repeats is crucial because they can cause false gene predictions and complicate downstream analyses.  

- **RepeatMasker** → the most widely used tool, screens DNA sequences for interspersed repeats and low complexity sequences using curated repeat libraries (e.g., RepBase, Dfam).  
- **Dfam** → curated database of transposable element families, based on profile HMMs.  
- **TEannot (REPET pipeline)** → specialized in de novo transposable element discovery and annotation.  

Repeat annotation ensures that repetitive sequences are catalogued, masked when necessary, and separated from protein-coding genes. This prevents false gene calls and improves downstream prediction.  

**Table – Repeat Annotation Tools**

| Tool/DB        | Principle                              | Output                           | Applications                       |
|----------------|----------------------------------------|----------------------------------|------------------------------------|
| **RepeatMasker** | Uses libraries (RepBase, Dfam) to find interspersed repeats | Annotated repeat-masked genome   | Genome masking, TE annotation      |
| **Dfam**       | Profile HMMs for TE families           | TE classification and consensus  | Reference for repeat detection     |
| **TEannot**    | De novo discovery + homology           | Custom TE library + annotation   | Non-model species, new TE families |

---

### Gene Prediction and Gene Models
The central task of genome annotation is to define **gene models**, which describe the structure of genes, including exons, introns, untranslated regions (UTRs), and coding sequences (CDS).  

Several strategies exist, often combined for best accuracy:  

1. **Ab initio prediction** → uses intrinsic sequence signals (start/stop codons, splice sites, codon bias).  
   - **AUGUSTUS** is the gold-standard tool; trained on reference genomes, it predicts complete gene structures from sequence alone.  
2. **Homology-based prediction** → aligns known proteins/transcripts from related species to infer gene models. Tools: **BLAST, Exonerate**.  
3. **Transcript evidence** → integrates RNA-seq data; mapped reads reveal expressed exons, introns, and splice junctions.  
4. **Integrative approaches** → pipelines like **MAKER, BRAKER2** combine ab initio predictions with transcript and homology evidence, creating consensus models.  

This **multi-layered approach** produces more reliable gene models, capturing both coding and non-coding RNAs.  

**Table – Ab initio vs Extrinsic/Integrative Approaches**

| Approach / Tool   | Input Data                        | Strengths                                            | Limitations                          |
|-------------------|-----------------------------------|------------------------------------------------------|--------------------------------------|
| **Ab initio (AUGUSTUS)** | DNA sequence only                  | Detects novel genes, independent of external data     | False positives in repetitive regions, needs training |
| **Homology-based** | Known proteins/transcripts        | High accuracy for conserved genes, useful in non-models | Misses species-specific/novel genes   |
| **Transcript evidence** | RNA-seq reads, ESTs               | Captures real expression, splicing isoforms           | Limited to expressed genes, condition-dependent |
| **Integrative (MAKER, BRAKER2)** | Combines ab initio + homology + RNA-seq | Highest reliability, consensus models, widely used    | Computationally intensive, needs multiple data types |

---

### Functional Annotation
After structural annotation, the next step is to assign **biological meaning** to genes and transcripts. This is achieved through cross-referencing with biological databases:  

- **Gene Ontology (GO)** → assigns terms for *Molecular Function, Biological Process, Cellular Component*.  
- **KEGG (Kyoto Encyclopedia of Genes and Genomes)** → maps genes into pathways (metabolism, signaling).  
- **Protein domains (Pfam, InterPro)** → detect conserved motifs, functional domains, and evolutionary relationships.  

Functional annotation enables downstream interpretation, linking raw sequence to biology, pathways, and phenotypes.  

---

### Genome Annotation – Pipeline Overview

1. **Repeat Annotation** → identify and mask repeats (RepeatMasker, Dfam, TEannot).  
2. **Structural Annotation** → predict gene models (ab initio, homology, transcript evidence, integrative pipelines).  
3. **Functional Annotation** → assign GO terms, KEGG pathways, Pfam/InterPro domains.  
4. **Curation** → manual review of key genes, correction of misannotations.  

---

### In Summary
Genome annotation transforms raw assemblies into functional blueprints of organisms. It requires **repeat masking**, **robust gene prediction**, and **functional annotation**. High-quality annotation is critical for:  
- **Comparative genomics** (ortholog/paralog detection, synteny).  
- **Transcriptomics** (RNA-seq alignment and quantification).  
- **Functional studies** (gene discovery, pathway analysis).  
- **Applied biotechnology** (breeding, engineering, synthetic biology).  

### Comparative Table – Gene Annotation Approaches

| Approach / Tool              | Input Data                                    | Output                          | Best Use Case                                         | Limitations |
|-------------------------------|-----------------------------------------------|---------------------------------|------------------------------------------------------|-------------|
| **Ab initio (e.g., AUGUSTUS)** | DNA sequence only                             | Predicted gene models (exons, CDS, introns, UTRs) | Detects **novel genes** without prior data; useful in non-model species | High false positives; requires species-specific training |
| **Homology-based (BLAST, Exonerate)** | Assembled genome + protein/cDNA sequences from related species | Gene structures aligned to known orthologs | Reliable for **conserved genes**; good for cross-species annotation | Misses lineage-specific genes; limited by quality of reference database |
| **Transcript evidence (RNA-seq, ESTs)** | RNA-seq reads or EST libraries mapped to genome | Expressed gene models, splice isoforms | Defines **real transcription evidence**, detects isoforms | Limited to expressed genes under sampled conditions |
| **Integrative pipelines (MAKER, BRAKER2)** | Genome + ab initio + homology + transcript data | Consensus, high-confidence gene models | **Gold standard**: combines multiple evidence sources; best for reference genomes | Computationally intensive; requires multiple datasets |

---

## 6. Specialized Sequencing

In addition to whole-genome sequencing, a range of **specialized approaches** targets specific biological questions. These methods focus on defined regions or molecular layers and complement standard WGS/RNA-seq.

---

### aCGH (Array Comparative Genomic Hybridization)

**Principle**: Detects **Copy Number Variations (CNVs)** using hybridization on a microarray. Test and reference DNA are labeled with different fluorescent dyes, co-hybridized to probes on a chip, and the fluorescence ratio indicates gains or losses.  

**Pipeline**:
1. **DNA extraction** (test + reference).  
2. **Labeling** → test (red) and control (green).  
3. **Co-hybridization** on **microarray chip** with thousands of probes.  
4. **Washing & scanning** → measure fluorescence intensity.  
5. **Signal analysis** → log₂ ratio (test/control).  

**CNV signal interpretation**:
- **Ratio ~0** → normal copy number.  
- **Positive log₂ ratio** → duplication/gain.  
- **Negative log₂ ratio** → deletion/loss.  

**Microarray “chip pipeline”**:  
DNA/RNA → labeling → hybridization on chip probes → scanner reads fluorescent signals → computational normalization → intensity plots → interpretation.  

**Table – aCGH Overview**

| Feature        | Details                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Principle**  | Test vs reference DNA → labeled → hybridized on microarray → intensity ratio = CNV |
| **Resolution** | Typically 20–100 kb (depends on probe density)                          |
| **Applications** | CNV detection in cancer genomics, constitutional disorders, prenatal screening |
| **Strengths**  | Genome-wide CNV detection, relatively cheap, robust technology          |
| **Limitations**| Cannot detect balanced rearrangements (inversions/translocations); limited resolution compared to NGS |

---

### Pool-seq

**Principle**: DNA from many individuals is pooled and sequenced together to estimate **allele frequencies**.  

- **Equimolar DNA pool**: each individual contributes the same DNA amount to avoid bias.  
- Provides allele frequency spectra, selective sweeps, and F_ST between populations, but no individual genotypes.  

**Table – Pool-seq Overview**

| Feature        | Details                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Principle**  | Pool DNA from multiple individuals → sequence together                  |
| **Output**     | Allele frequency spectra, population-level diversity                    |
| **Applications** | Population genomics, selective sweep detection, allele frequency estimation |
| **Strengths**  | Cost-effective, fast, scalable                                          |
| **Limitations**| No individual genotypes, sensitive to unequal DNA contributions         |

---

### Targeted Sequencing

**Principle**: Focuses sequencing on specific regions instead of the whole genome.  

- **Amplicon sequencing**: PCR amplifies specific loci.  
- **Hybrid capture panels**: probes enrich predefined genomic regions (e.g., cancer or exome panels).  

**Table – Targeted Sequencing Overview**

| Feature        | Details                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Principle**  | Focused sequencing of selected loci (amplicons or capture panels)        |
| **Applications** | Cancer gene panels, rare disease diagnosis, pharmacogenomics          |
| **Strengths**  | High depth of coverage, cost-efficient, tailored to specific questions  |
| **Limitations**| Misses variants outside target, capture bias, design required in advance |

---

### Exome (Whole-Exome Sequencing — WES)

**Principle**: Captures and sequences protein-coding regions (~1–2% of genome).  

**Pipeline**:
1. DNA fragmentation.  
2. Library prep with adapters.  
3. **Hybrid capture** with biotinylated probes against exons.  
4. Wash & pull-down with streptavidin beads.  
5. PCR amplification.  
6. Sequencing → variant calling.  

**Applications**: Mendelian diagnostics, high-depth coding variant discovery, trio analysis.  

**Table – WES Overview**

| Feature        | Details                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Principle**  | Hybrid capture of coding regions (~20,000 genes, ~1–2% of genome)       |
| **Applications** | Clinical diagnostics, Mendelian disorders, coding SNP/indel discovery |
| **Strengths**  | Cheaper than WGS, deeper coverage, focuses on known disease-causing regions |
| **Limitations**| Misses regulatory/non-coding variants, capture bias, uneven coverage    |

---

### Methyl-seq / Bisulfite Sequencing

**Principle**: Bisulfite converts **unmethylated C → U** (read as T), while methylated C remains unchanged. Sequencing reveals methylation at **base resolution**.  

**Pipeline**:
1. High-quality DNA extraction.  
2. Bisulfite treatment (C → U if unmethylated).  
3. Library preparation (adapters, PCR).  
4. Sequencing.  
5. Alignment to reference genome.  
6. **Methylation calling**: compare C vs T reads → % methylation per CpG.  

**Table – Methyl-seq Overview**

| Feature        | Details                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Principle**  | Sodium bisulfite converts unmethylated C to U (→ T in sequencing); methylated C unchanged |
| **Resolution** | Single base (CpG methylation profiles)                                  |
| **Applications** | Epigenetics, imprinting, cancer methylome, developmental biology     |
| **Strengths**  | High resolution, genome-wide, quantitative                             |
| **Limitations**| DNA damage from bisulfite, incomplete conversion artifacts, requires high input |

---

### RNA-seq

**Principle**: Profiles the transcriptome by sequencing cDNA derived from RNA. Captures isoforms, splicing, and expression levels.  

**Pipeline**:
1. RNA extraction & quality check (RIN score).  
2. **Library prep**: poly(A) selection or rRNA depletion; RNA fragmented and reverse-transcribed into cDNA; adapters added.  
3. Sequencing (usually Illumina).  
4. **QC & trimming** (FastQC).  
5. **Alignment** to reference genome (STAR, HISAT2) or transcriptome (Salmon, Kallisto).  
6. **Quantification** of expression (counts/TPM/FPKM).  
7. Downstream: differential expression (DESeq2, edgeR), isoform analysis, pathway enrichment.  

**Table – RNA-seq Overview**

| Feature        | Details                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Principle**  | Sequence cDNA derived from RNA (polyA selected or rRNA-depleted)        |
| **Output**     | Expression profiles (counts, TPM, FPKM); isoforms; fusion transcripts   |
| **Applications** | Differential expression, splicing analysis, gene fusion discovery     |
| **Strengths**  | Genome-wide transcriptome profiling, detects novel isoforms             |
| **Limitations**| Sensitive to RNA quality, batch effects, alignment complexity           |

---

### Genome-Wide Association Studies (GWAS)

**Principle**: Detects statistical associations between SNPs and phenotypes across large populations.  

**Pipeline**:
1. **Genotyping & QC** (call rate, Hardy–Weinberg equilibrium, minor allele frequency).  
2. **Population structure** correction (PCA/MDS, kinship matrices).  
3. **Association testing** (linear/logistic regression, mixed models).  
4. **Multiple testing correction** (Bonferroni, FDR).  
5. **Visualization** → Manhattan plot, QQ plot.  
6. **Interpretation** → fine-mapping, candidate gene identification, pathway analysis.  

**Manhattan plots**: each SNP plotted by genomic position (x-axis) vs –log₁₀(p) (y-axis); significant peaks highlight candidate loci.  

**Bonferroni correction**: very conservative (α/N); often supplemented by FDR to retain power.  

**Challenges**:  
- Biased toward common variants with moderate effect.  
- Rare variants require sequencing and gene-based burden tests.  
- Population stratification can inflate false positives.  
- GWAS often limited to European cohorts (ancestry bias).  

**Links to Population Genomics**:  
- **LD**: fundamental to detect associations (tag SNPs).  
- **ROH / inbreeding**: impact homozygosity and recombination patterns.  
- **PCA/MDS**: visualize and correct for population structure.  
- **Effective population size**: influences LD decay and GWAS resolution.  

**Table – GWAS Overview**

| Feature        | Details                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Principle**  | Association study testing SNPs across genome vs phenotype               |
| **Pipeline**   | QC → structure correction → association → multiple testing → visualization → interpretation |
| **Output**     | Manhattan plots, QQ plots, fine-mapped loci                             |
| **Strengths**  | Genome-wide, unbiased, powerful for common variants                     |
| **Limitations**| Rare variants hard to detect, ancestry bias, population stratification issues |

---

**In summary**, specialized methods (aCGH, Pool-seq, targeted panels, **WES**, methyl-seq, RNA-seq) extend genomic analysis beyond WGS. **GWAS** builds on these data to link genetic variation with traits, but requires careful handling of QC, population structure, and multiple testing to produce biologically valid results.

