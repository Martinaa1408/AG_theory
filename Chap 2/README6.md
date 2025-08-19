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

<img width="543" height="658" alt="Screenshot 2025-08-19 170716" src="https://github.com/user-attachments/assets/08f24524-35ee-4987-95b7-d7ecae5b43f3" />


### Moore’s Law and Sequencing Costs

* Sequencing costs dropped **faster than Moore’s Law**.
* From \$3B for the Human Genome Project → <\$1000 per human genome today.

---

## 2. Illumina Sequencing

<img width="664" height="456" alt="Screenshot 2025-08-19 171204" src="https://github.com/user-attachments/assets/733faf0a-668c-4627-b18f-3ca99f6676cb" />


* **Principle**: Sequencing by synthesis with reversible terminators.
* **Flow cell**: glass slide with lanes coated with adapters.
* **Cluster generation**: bridge amplification forms dense clusters of identical DNA fragments.
* **Sequencing**: fluorescently labeled nucleotides incorporated, imaged cycle by cycle.
* **Paired-end sequencing**: both ends of DNA fragments sequenced, improves mapping.
* **Output**: FASTQ files with base calls and quality scores.

<img width="583" height="663" alt="Screenshot 2025-08-19 171222" src="https://github.com/user-attachments/assets/5bf04b93-02fe-45d5-828e-9d94de6caeff" />
<img width="536" height="686" alt="Screenshot 2025-08-19 171233" src="https://github.com/user-attachments/assets/2971f8de-b2ef-47bb-8366-037df519a72c" />
<img width="563" height="626" alt="Screenshot 2025-08-19 171257" src="https://github.com/user-attachments/assets/1a4d8664-8349-46ff-ad8a-56ff8b7776b9" />

---

## 3. Ion Torrent Sequencing

<img width="622" height="591" alt="Screenshot 2025-08-19 170905" src="https://github.com/user-attachments/assets/ffce8050-6282-4b51-a5e5-b38fa27b1136" />

<img width="589" height="599" alt="Screenshot 2025-08-19 170942" src="https://github.com/user-attachments/assets/42f576a1-7347-4942-87c8-64c2280d3565" />

* **Principle**: semiconductor sequencing, detects H+ ions released when a nucleotide is incorporated.
* **Signal output**: ionograms (peaks of pH changes).
* **Advantages**: fast, no fluorescence/optics.
* **Limitations**: homopolymer errors (AAAA runs).
* **Chips**: 510, 520, 530 vary in throughput and read capacity.
  
<img width="547" height="549" alt="Screenshot 2025-08-19 170809" src="https://github.com/user-attachments/assets/6e38b2e6-9ff8-4b1b-b812-b328d643a5f9" />

<img width="681" height="558" alt="Screenshot 2025-08-19 170842" src="https://github.com/user-attachments/assets/51c4b6af-a724-435d-bf5b-2fb0a65b9686" />


---

## 4. Roche 454 Pyrosequencing (Historical)

* One of the first NGS platforms (2005).
* Principle: pyrophosphate release detected via luciferase (light signal).
* Read length: \~700 bp (longer than Illumina at the time).
* Limitations: homopolymer errors, expensive reagents.
* Discontinued in 2016.

<img width="574" height="632" alt="Screenshot 2025-08-19 171126" src="https://github.com/user-attachments/assets/f30c3927-aed8-4c6e-b0f0-67c0e02e1bcc" />

---

## 5. ABI SOLiD Sequencing (Historical)

* Sequencing by ligation with fluorescent probes.
* Two-base encoding ensured high accuracy.
* Limitations: short reads, complex data processing.
* Discontinued but important for early NGS development.

<img width="604" height="580" alt="Screenshot 2025-08-19 171132" src="https://github.com/user-attachments/assets/41341886-9909-4e8e-a92f-80f2d1a7cf04" />


---

## 6. Long-Read Sequencing

<img width="443" height="580" alt="Screenshot 2025-08-19 171534" src="https://github.com/user-attachments/assets/1f4bb9fe-2c75-4f31-9949-548d87716dc2" />


### Pacific Biosciences (PacBio)

<img width="563" height="643" alt="Screenshot 2025-08-19 171442" src="https://github.com/user-attachments/assets/8963e223-b295-4a3d-b831-9b73978b2786" />
<img width="529" height="687" alt="Screenshot 2025-08-19 171518" src="https://github.com/user-attachments/assets/726315d5-f439-4a63-a4a9-892eaaaa2a91" />
<img width="509" height="670" alt="Screenshot 2025-08-19 171526" src="https://github.com/user-attachments/assets/bade749b-84e5-4df9-8c51-8cb3d8f2a95a" />

* **SMRT (Single Molecule Real-Time) sequencing**.
* Uses zero-mode waveguides (ZMWs).
* Produces **long reads (10–50 kb)**.
* **HiFi/CCS (Circular Consensus Sequencing)** improves accuracy (>99%).
* Useful for structural variants, genome assembly, isoform sequencing.

### Oxford Nanopore Technologies (ONT)

<img width="594" height="760" alt="Screenshot 2025-08-19 171414" src="https://github.com/user-attachments/assets/0e709db4-333d-4e5a-bdba-2a570d12fdfe" />
<img width="530" height="691" alt="Screenshot 2025-08-19 171431" src="https://github.com/user-attachments/assets/77ff4351-e2b6-4fd2-86cd-5aeb05159858" />

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

<img width="597" height="545" alt="Screenshot 2025-08-19 171539" src="https://github.com/user-attachments/assets/83867911-19d9-4bd8-9cd2-a52c89f45b10" />

---

## ✅ Summary

* **Sanger**: gold standard, low throughput.
* **Illumina**: dominant NGS platform, short accurate reads.
* **Ion Torrent**: fast, low-cost, but homopolymer issues.
* **Roche 454 & ABI SOLiD**: historical, now discontinued.
* **PacBio & ONT**: long-read technologies, crucial for structural genomics.
* Costs have plummeted, enabling population-scale genomics.
