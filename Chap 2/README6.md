# README6 – Sequencing Technologies

## 1. Sanger Sequencing vs NGS

### Sanger Sequencing

* Developed in 1977 by Frederick Sanger.
* Uses **dideoxynucleotides (ddNTPs)** for chain termination.
* Produces highly accurate reads (\~1000 bp).
* Low throughput, expensive for large genomes.

### Next-Generation Sequencing (NGS)

* Emerged in mid-2000s.
* Massively parallel sequencing of millions of DNA fragments.
* Short reads (50–300 bp) but very high throughput.
* Revolutionized genomics (faster, cheaper, scalable).

### Moore’s Law and Sequencing Costs

* Sequencing costs dropped **faster than Moore’s Law**.
* From \$3B for the Human Genome Project → <\$1000 per human genome today.

---

## 2. Illumina Sequencing

* **Principle**: Sequencing by synthesis with reversible terminators.
* **Flow cell**: glass slide with lanes coated with adapters.
* **Cluster generation**: bridge amplification forms dense clusters of identical DNA fragments.
* **Sequencing**: fluorescently labeled nucleotides incorporated, imaged cycle by cycle.
* **Paired-end sequencing**: both ends of DNA fragments sequenced, improves mapping.
* **Output**: FASTQ files with base calls and quality scores.

---

## 3. Ion Torrent Sequencing

* **Principle**: semiconductor sequencing, detects H+ ions released when a nucleotide is incorporated.
* **Signal output**: ionograms (peaks of pH changes).
* **Advantages**: fast, no fluorescence/optics.
* **Limitations**: homopolymer errors (AAAA runs).
* **Chips**: 510, 520, 530 vary in throughput and read capacity.

---

## 4. Roche 454 Pyrosequencing (Historical)

* One of the first NGS platforms (2005).
* Principle: pyrophosphate release detected via luciferase (light signal).
* Read length: \~700 bp (longer than Illumina at the time).
* Limitations: homopolymer errors, expensive reagents.
* Discontinued in 2016.

---

## 5. ABI SOLiD Sequencing (Historical)

* Sequencing by ligation with fluorescent probes.
* Two-base encoding ensured high accuracy.
* Limitations: short reads, complex data processing.
* Discontinued but important for early NGS development.

---

## 6. Long-Read Sequencing

### Pacific Biosciences (PacBio)

* **SMRT (Single Molecule Real-Time) sequencing**.
* Uses zero-mode waveguides (ZMWs).
* Produces **long reads (10–50 kb)**.
* **HiFi/CCS (Circular Consensus Sequencing)** improves accuracy (>99%).
* Useful for structural variants, genome assembly, isoform sequencing.

### Oxford Nanopore Technologies (ONT)

* DNA passes through nanopore, changes in electrical current measured.
* Portable sequencers (MinION, GridION, PromethION).
* Read length: can exceed 1 Mb.
* Cost-effective, real-time sequencing.
* Error rate higher than Illumina, but improving.

---

## 7. Cost per Megabase (Approximate, 2024)

* **Illumina short-read**: \~\$0.01–0.05 / Mb.
* **ONT**: \~\$0.05–0.10 / Mb.
* **PacBio HiFi**: \~\$0.15–0.20 / Mb.
* **Sanger**: \~\$500 / Mb (too costly for WGS).

---

## ✅ Summary

* **Sanger**: gold standard, low throughput.
* **Illumina**: dominant NGS platform, short accurate reads.
* **Ion Torrent**: fast, low-cost, but homopolymer issues.
* **Roche 454 & ABI SOLiD**: historical, now discontinued.
* **PacBio & ONT**: long-read technologies, crucial for structural genomics.
* Costs have plummeted, enabling population-scale genomics.
