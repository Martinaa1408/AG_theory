# README9 – Data Formats & Bioinformatics Tools

## 1. Raw Sequencing Data

### FASTQ Format

* Stores raw reads + quality scores.
* Each entry = 4 lines:

  1. `@` identifier
  2. Sequence (A, C, G, T, N)
  3. `+` separator
  4. Quality string (ASCII-encoded, Phred scores)
* **Phred quality score (Q)**:

  * Q = -10 log₁₀(p error)
  * Q20 = 99% accuracy, Q30 = 99.9% accuracy.

---

## 2. Alignment Data

### SAM/BAM Format

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

### VCF (Variant Call Format)

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
  * **Ensembl genome browser**: gene annotation, comparative genomics.

---

## ✅ Summary

* **FASTQ** stores raw reads, **SAM/BAM** store alignments, **VCF** stores variants.
* **GFF/GTF and BED** describe genomic features and coordinates.
* Quality control (FastQC, trimming) ensures reliable downstream analyses.
* Core tools: BWA-MEM for alignment, GATK for variant calling, IGV/Ensembl for visualization.
