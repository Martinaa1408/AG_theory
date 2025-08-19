# README11 – Specialized Sequencing Approaches

## 1. aCGH (Array Comparative Genomic Hybridization)
<img width="540" height="618" alt="Screenshot 2025-08-19 172753" src="https://github.com/user-attachments/assets/a2b3faec-fb8d-4a81-9b75-6edcdcbff638" />

* Detects **copy number variations (CNVs)**.
* Principle:

  * Label test and reference DNA with different fluorophores.
  * Hybridize to microarray with genomic probes.
  * Measure fluorescence ratio → CNV detection.
* Applications: cancer genomics, genetic diagnostics.

---

## 2. Pool-Seq (Pooled Sequencing)
<img width="493" height="622" alt="Screenshot 2025-08-19 172424" src="https://github.com/user-attachments/assets/e25317a2-5539-47e4-a3c1-6566d53f591f" />

* Sequencing DNA from pooled individuals.
* Provides allele frequency estimates at population level.
* Advantages: cost-effective, scalable.
* Limitations: no individual genotypes, sensitive to pooling errors.

---

## 3. Targeted DNA Sequencing
<img width="591" height="637" alt="Screenshot 2025-08-19 172439" src="https://github.com/user-attachments/assets/cb6f4d47-316d-4121-91cb-81fed7925f1e" />


* Focused sequencing of specific genomic regions.
* Methods:

  * **Amplicon-based**: PCR amplify selected loci.
  * **Hybrid capture**: probes capture targeted regions.
* Applications: clinical panels, candidate gene studies.

---

## 4. DNA Methylation Sequencing
<img width="541" height="660" alt="Screenshot 2025-08-19 172501" src="https://github.com/user-attachments/assets/fe0f0ac0-d7ac-4d36-801c-53a4c1f7b338" />


* **Methyl-seq**: genome-wide profiling of DNA methylation.
* **Bisulfite sequencing**:

  * Sodium bisulfite converts unmethylated cytosines → uracil.
  * Methylated cytosines remain unchanged.
  * Sequencing reveals methylation status at single-base resolution.
* Applications: epigenetics, cancer, imprinting studies.

---

## 5. RNA Sequencing (RNA-seq)
<img width="459" height="392" alt="Screenshot 2025-08-19 172517" src="https://github.com/user-attachments/assets/ae8a0043-94b4-475c-916a-2beb581fad54" />


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
