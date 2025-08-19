# README10 – Population Genomics Applications

## 1. Genotyping Technologies

### SNP Arrays

* High-density chips with thousands to millions of SNP markers.
* Platforms: Illumina Infinium, Affymetrix Axiom.
* Pros: reproducible, cost-effective for known SNPs.
* Cons: limited to pre-selected variants.

### GBS (Genotyping-by-Sequencing)

* Reduced-representation sequencing method using restriction enzymes.
* Captures random genome-wide markers.
* Pros: flexible, scalable.
* Cons: missing data, biased by enzyme choice.

### RAD-seq (Restriction site Associated DNA sequencing)

* Similar to GBS but uses stricter fragment selection.
* Useful for non-model organisms.

### Re-GBS (Refined GBS)

* Optimization of GBS with improved reproducibility and reduced bias.

---

## 2. Structural Variant Detection

* **SNPs**: single-base substitutions.
* **CNVs**: copy number variations.
* **Inversions & translocations**: large structural changes.
* Methods:

  * Read depth (coverage differences).
  * Read pair (discordant insert size/orientation).
  * Split read (partial alignments).
  * Assembly-based approaches.

---

## 3. GWAS (Genome-Wide Association Studies)

* Goal: identify SNPs associated with traits/diseases.
* **Pipeline**:

  1. Input genotype data (.ped, .map, or VCF).
  2. Quality control (filtering missingness, MAF, HWE deviations).
  3. Correct for population structure (PCA, MDS, mixed models).
  4. Statistical testing (logistic/linear regression).
  5. Multiple testing correction (Bonferroni, FDR).
  6. Visualization (Manhattan and QQ plots).

---

## 4. GWAS Output

* **Manhattan plot**:

  * x-axis = genomic position.
  * y-axis = –log10(p-value).
  * Peaks highlight genomic regions associated with traits.
* **QQ plot**: expected vs observed p-values, used to check inflation.

---

## 5. Inbreeding, LD, and ROH in GWAS

* **Inbreeding coefficient (F)**: used to correct for relatedness.
* **LD (linkage disequilibrium)**: helps in tagging variants and fine-mapping.
* **ROH (runs of homozygosity)**: detect autozygosity, consanguinity.

---

## 6. QTL Mapping

* QTL (Quantitative Trait Loci) mapping: identifies genomic regions associated with trait variation.
* Requires controlled crosses or pedigrees.
* Uses linkage information rather than population-wide association.

---

## 7. Marker-Assisted Selection (MAS) & Marker-Assisted Backcrossing (MAB)

* **MAS**: selecting individuals based on genetic markers linked to desirable traits.
* **MAB**: using markers to accelerate introgression of traits from donor to recipient lines.
* Applications: plant/animal breeding, disease resistance, yield improvement.

---

## ✅ Summary

* Genotyping technologies: SNP arrays (fixed), GBS/RAD (flexible), Re-GBS (optimized).
* GWAS pipeline: genotype → QC → structure correction → association testing → visualization.
* Structural variation detection relies on read depth, split reads, and assembly.
* MAS/MAB integrate genomic data into breeding programs for applied genetics.
