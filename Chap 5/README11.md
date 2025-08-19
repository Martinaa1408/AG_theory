# README11 – Specialized Sequencing Approaches

## 1. aCGH (Array Comparative Genomic Hybridization)

* Detects **copy number variations (CNVs)**.
* Principle:

  * Label test and reference DNA with different fluorophores.
  * Hybridize to microarray with genomic probes.
  * Measure fluorescence ratio → CNV detection.
* Applications: cancer genomics, genetic diagnostics.

---

## 2. Pool-Seq (Pooled Sequencing)

* Sequencing DNA from pooled individuals.
* Provides allele frequency estimates at population level.
* Advantages: cost-effective, scalable.
* Limitations: no individual genotypes, sensitive to pooling errors.

---

## 3. Targeted DNA Sequencing

* Focused sequencing of specific genomic regions.
* Methods:

  * **Amplicon-based**: PCR amplify selected loci.
  * **Hybrid capture**: probes capture targeted regions.
* Applications: clinical panels, candidate gene studies.

---

## 4. DNA Methylation Sequencing

* **Methyl-seq**: genome-wide profiling of DNA methylation.
* **Bisulfite sequencing**:

  * Sodium bisulfite converts unmethylated cytosines → uracil.
  * Methylated cytosines remain unchanged.
  * Sequencing reveals methylation status at single-base resolution.
* Applications: epigenetics, cancer, imprinting studies.

---

## 5. RNA Sequencing (RNA-seq)

* High-throughput sequencing of transcriptomes.
* **Workflow**:

  1. RNA extraction.
  2. Conversion to cDNA.
  3. Library preparation (fragmentation + adapters).
  4. Sequencing (Illumina, ONT, PacBio Iso-Seq).
  5. Data analysis: read alignment, quantification, differential expression.
* Applications: gene expression profiling, isoform detection, fusion transcripts, eQTL studies.

---

## ✅ Summary

* **aCGH**: CNV detection with hybridization arrays.
* **Pool-seq**: population-level allele frequency estimation.
* **Targeted sequencing**: focused approach for diagnostics and specific traits.
* **Methyl-seq/bisulfite**: base-resolution DNA methylation.
* **RNA-seq**: transcriptome-wide expression and isoform analysis.
