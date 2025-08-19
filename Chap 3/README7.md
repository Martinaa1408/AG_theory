# README7 – Genome Assembly

## 1. Shotgun Sequencing and Assembly Pipeline

* **Shotgun sequencing**: DNA fragmented randomly, sequenced, then computationally assembled.
* **Pipeline**:

  1. DNA extraction.
  2. Library preparation (fragmentation + adapters).
  3. Sequencing (short-read or long-read).
  4. Read mapping or de novo assembly.
  5. Assembly polishing and validation.

---

## 2. Genome Size and k-mers

* **Genome size estimation** can be performed using k-mer frequency analysis.
* **k-mers**: subsequences of length k.

  * Provide information on genome size, heterozygosity, and repeat content.
* Histogram of k-mers helps identify sequencing depth and errors.

---

## 3. Assembly Algorithms

### Greedy Algorithm

* Joins reads iteratively based on best overlap.
* Simple but error-prone in repetitive regions.

### OLC (Overlap–Layout–Consensus)

* Construct overlap graph of all reads.
* Layout defines path through graph.
* Consensus sequence built from aligned reads.
* Used in early Sanger-based assemblies.

### de Bruijn Graph (DBG)

* Reads broken into k-mers.
* Nodes = k-mers, edges = overlaps.
* Efficient for short reads (Illumina).
* Handles large datasets but sensitive to errors and repeats.


---

## 4. Mate-Pair Sequencing and Scaffolding

* **Mate-pair sequencing**: generates pairs of reads separated by large inserts (2–10 kb).
* Helps resolve repeats and connect contigs.
* **Scaffolding**: ordering and orienting contigs using mate-pairs or long reads.
* Produces larger contiguous sequences (scaffolds).

  
<img width="577" height="562" alt="Screenshot 2025-08-19 172151" src="https://github.com/user-attachments/assets/27876f6a-1dab-42c4-92f9-a952bfa4eead" />

---

## 5. Gap Filling

* Gaps remain after scaffolding.
* Filled using:

  * Additional sequencing (long reads).
  * Local assembly.
  * PCR-based finishing.

---

## 6. Assembly Quality Metrics

* **N50**: length such that 50% of the genome is in contigs/scaffolds of this size or longer.
* **Coverage (depth)**: average number of times each base is sequenced.

  * Coverage = (total bases sequenced) / (genome size).
* **C-value**: DNA content per haploid genome, used as reference for genome size.


---

## 7. Assembly Completeness

* **BUSCO (Benchmarking Universal Single-Copy Orthologs)**:

  * Searches for expected conserved genes.
  * Categories: Complete, Duplicated, Fragmented, Missing.
  * Provides an estimate of assembly and annotation quality.

<img width="560" height="511" alt="Screenshot 2025-08-19 172211" src="https://github.com/user-attachments/assets/1491530c-27b3-4b57-acc1-19c52c868af7" />

---

## ✅ Summary

* Genome assembly transforms raw reads into contiguous sequences.
* Algorithms: Greedy, OLC, and de Bruijn graphs.
* Mate-pairs and scaffolding improve contiguity.
* Metrics like N50, coverage, and BUSCO assess assembly quality.
* Hybrid approaches (short + long reads) often yield the best assemblies.
