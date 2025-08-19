# README3 – Population Genetics Basics

## 1. Hardy–Weinberg Equilibrium (HWE)

* **Principle**: In a large, randomly mating population with no selection, migration, or mutation, allele and genotype frequencies remain constant.
* **Equation**:

  * Allele frequencies: p + q = 1
  * Genotype frequencies: p² (AA) + 2pq (Aa) + q² (aa) = 1
* **Uses**:

  * Test for evolutionary forces acting on populations.
  * Compare observed vs expected genotype frequencies.

**Example**: If allele A frequency p = 0.7, then q = 0.3.

* Expected genotypes: AA = 0.49, Aa = 0.42, aa = 0.09.

<img width="500" height="670" alt="Screenshot 2025-08-19 170518" src="https://github.com/user-attachments/assets/fa920efe-d3b8-4939-a1e7-96b3f3a7b019" />

<img width="500" height="670" alt="Screenshot 2025-08-19 170518" src="https://github.com/user-attachments/assets/214d5083-8195-45bf-9062-06eec067fedf" />

---

## 2. Allele and Genotype Frequencies

* **Allele frequency**: proportion of a specific allele among all alleles at a locus.
* **Genotype frequency**: proportion of individuals with a specific genotype.
* **Calculation**:

  * Allele A frequency = (2 × #AA + #Aa) / (2 × total individuals).

---

## 3. Inbreeding and the Inbreeding Coefficient (F)

* **Inbreeding**: mating between related individuals.
* **Effect**: increases homozygosity, decreases heterozygosity.
* **Inbreeding coefficient (F)**: probability that two alleles are identical by descent (IBD).

  * Example: offspring of first cousins → F = 1/16 = 0.0625.
* Consequences: inbreeding depression (loss of fitness).

---

## 4. Runs of Homozygosity (ROH)

* Long stretches of homozygous genotypes in the genome.
* Indicators of recent inbreeding or population bottlenecks.
* Used to estimate genomic inbreeding levels.
  
<img width="554" height="596" alt="Screenshot 2025-08-19 173113" src="https://github.com/user-attachments/assets/a7f74236-7eef-4514-a08a-2dd117cfa1e0" />

---

## 5. Linkage Disequilibrium (LD)

* **Definition**: non-random association of alleles at two or more loci.
* **Causes**: linkage, small population size, selection, demographic events.
* **Measures**:

  * **D′**: historical recombination (D′ = 1 = complete LD).
  * **r²**: correlation between loci (r² = 1 = perfect predictability).
* Important in GWAS for tagging SNPs and mapping traits.

---

## 6. Population Structure

* Populations often show stratification (subdivision by ancestry).
* Detected by **PCA or MDS** (multidimensional scaling).
* Essential correction in GWAS to avoid false positives.

---

## ✅ Summary

* HWE provides baseline expectations for allele and genotype frequencies.
* Inbreeding increases homozygosity; quantified by coefficient F and ROH.
* LD measures how alleles at different loci associate; key for mapping.
* Population structure correction (PCA/MDS) is critical for genomic studies.
