# README9 – Data Formats & Bioinformatics Tools

## 1. Raw Sequencing Data

### FASTQ Format
<img width="583" height="619" alt="Screenshot 2025-08-19 171644" src="https://github.com/user-attachments/assets/420eaac8-34c1-4ff9-890d-2f784f046afa" />

<img width="536" height="674" alt="Screenshot 2025-08-19 171707" src="https://github.com/user-attachments/assets/d051eec3-2761-497a-b128-133344c45776" />


* Stores raw reads + quality scores.
* Each entry = 4 lines:

  1. `@` identifier
  2. Sequence (A, C, G, T, N)
  3. `+` separator
  4. Quality string (ASCII-encoded, Phred scores)
* **Phred quality score (Q)**:
<img width="664" height="455" alt="Screenshot 2025-08-19 171741" src="https://github.com/user-attachments/assets/3ae6d091-f8f3-4d29-8770-3326dd46d14c" />


  * Q = -10 log₁₀(p error)
  * Q20 = 99% accuracy, Q30 = 99.9% accuracy.

---

## 2. Alignment Data


### SAM/BAM Format
<img width="599" height="628" alt="Screenshot 2025-08-19 171838" src="https://github.com/user-attachments/assets/fbf0fbc6-87c0-40a4-b763-271bd79132a6" />

<img width="560" height="718" alt="Screenshot 2025-08-19 171851" src="https://github.com/user-attachments/assets/e176eaa6-9264-4f19-b8f9-e568ba25855d" />

* **SAM (Sequence Alignment/Map)**: text-based alignment format.
* **BAM**: binary compressed version of SAM.
* Key fields:

  * QNAME: read ID
  * FLAG: read properties (paired, mapped, etc.)
  * RNAME: reference sequence
  * POS: alignment position
  * MAPQ: mapping quality
  * CIGAR: alignment string (e.g., 76M = 76 matches)
* **CIGAR strings**:

  * M = match/mismatch
  * I = insertion
  * D = deletion
  * S = soft-clipped bases

---

## 3. Variant Data

<img width="543" height="578" alt="Screenshot 2025-08-19 171909" src="https://github.com/user-attachments/assets/a890862c-452b-46b7-be00-d91bd09ee154" />

### VCF (Variant Call Format)

<img width="576" height="590" alt="Screenshot 2025-08-19 171926" src="https://github.com/user-attachments/assets/61c43226-4213-4498-a1fd-a9677ccd6bac" />

* Stores genetic variants (SNPs, indels, structural variants).
* Structure:

  * **Header**: metadata, filters, references.
  * **Body**: one line per variant.

    * CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, samples.
* Example:

  ```
  1   1234567   rs123   A   G   50   PASS   DP=100;AF=0.5
  ```
* Widely used in GWAS, population genomics, clinical genetics.

---

## 4. Annotation & Genomic Coordinates

### GFF/GTF

* **GFF (General Feature Format)**, **GTF (Gene Transfer Format)**.
* Store gene models and functional elements.
* Fields: seqname, source, feature, start, end, score, strand, phase, attributes.

### BED

* Simple format for genomic regions.
* Columns: chrom, start, end, (optional) name, score, strand.
* Useful for interval queries.

### Flat Files

* Legacy text-based representations of genomic features.
* Less used, replaced by standardized formats.

---

## 5. Quality Control & Processing

* **FastQC**: assess read quality, GC content, adapter contamination.
* **Trimming tools**: Cutadapt, Trimmomatic.
* **Filtering**: remove low-quality or short reads.

---

## 6. Key Tools in Data Processing

* **Alignment**: BWA-MEM, Bowtie2, Minimap2.
* **Variant calling**: GATK, bcftools, FreeBayes.
* **Visualization**:

  * **IGV (Integrative Genomics Viewer)**: browse alignments, variants, annotations.
    <img width="573" height="760" alt="Screenshot 2025-08-19 172402" src="https://github.com/user-attachments/assets/078e04cd-f867-4058-a955-a8b3002facbf" />

  * **Ensembl genome browser**: gene annotation, comparative genomics.

---

## ✅ Summary

* **FASTQ** stores raw reads, **SAM/BAM** store alignments, **VCF** stores variants.
* **GFF/GTF and BED** describe genomic features and coordinates.
* Quality control (FastQC, trimming) ensures reliable downstream analyses.
* Core tools: BWA-MEM for alignment, GATK for variant calling, IGV/Ensembl for visualization.
