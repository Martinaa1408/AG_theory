# README8 – Genome Annotation

## 1. Introduction to Genome Annotation

* **Genome annotation**: identifying functional elements in assembled sequences.
* Two main types:

  1. **Structural annotation**: locating genes, exons, introns, repeats, regulatory regions.
  2. **Functional annotation**: assigning biological meaning (e.g., protein function, pathways).

---

## 2. Repeat Annotation

* Large portions of eukaryotic genomes are repetitive.
* **Tools**:

  * **RepeatMasker**: screens DNA for interspersed repeats and low-complexity sequences using databases like Dfam.
  * **Dfam**: database of repetitive DNA families.
  * **TEdenovo/TEannot**: de novo identification and annotation of transposable elements.
* Importance: repeats complicate assembly and annotation, but are biologically significant.

---

## 3. Structural Annotation of Genes

### Gene Components

* **ORFs (Open Reading Frames)**: potential protein-coding regions.
* **Exons**: coding segments.
* **Introns**: non-coding segments spliced out of RNA.
* **UTRs**: untranslated regions at mRNA ends.
* **Pseudogenes**: nonfunctional gene copies due to mutations.

### Approaches

1. **Ab initio prediction**

   * Uses algorithms to predict genes from DNA sequence alone.
   * Tools: **AUGUSTUS**, GeneMark.
   * Pros: works without external data.
   * Cons: prone to false positives.

2. **Homology-based prediction**

   * Aligns known genes/proteins from related species.
   * Tools: BLAST, Exonerate.
   * Pros: accurate for conserved genes.
   * Cons: limited for novel or lineage-specific genes.

3. **Evidence-based (transcriptomic)**

   * Uses RNA-seq or ESTs to confirm expressed genes.
   * Provides direct evidence of exons and splice sites.

4. **Integrative approaches**

   * Combine ab initio, homology, and transcript evidence.
   * Tools: **MAKER**, **BRAKER2**, **EVM (Evidence Modeler)**.

---

## 4. Functional Annotation

* Assigns biological meaning to predicted genes.
* Methods:

  * Protein domain detection (Pfam, InterPro).
  * Homology searches (BLASTp vs UniProt/NR).
  * Gene Ontology (GO) terms.
  * Pathway mapping (KEGG, Reactome).

---

## 5. Quality Metrics

* **Annotation edit distance (AED)**: measures agreement between prediction and evidence.
* **BUSCO**: evaluates completeness of annotated gene set.

---

## ✅ Summary

* Annotation = structural (gene models, repeats) + functional (biological roles).
* Tools: RepeatMasker/Dfam for repeats, AUGUSTUS for ab initio, BLAST for homology, MAKER/BRAKER2 for integration.
* Functional annotation assigns roles via homology, domains, and ontology.
* Accurate annotation is critical for downstream genomics, transcriptomics, and comparative studies.
