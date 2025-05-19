# 🧬Benchmarking Nanopore-Based RNA Modification Detection

This repository contains code and data processing workflows for the study:

**"Benchmarking Nanopore-Based RNA Modification Detection: Comparative Analysis of m6A Methylation and Pseudouridine Profiling with GLORI and BID-seq."**

---

## 📂 Datasets

This project compares two major types of RNA modifications — **m6A (N6-methyladenosine)** and **ψ (Pseudouridine)** — using multiple detection methods.

### m6A Datasets

- `m6A_cov10.bed`: Nanopore-based m6A prediction (Dorado v0.8.2)
- `GLORI.chrALL.tsv`: GLORI-based absolute methylation measurements ([Liu et al., 2024](https://www.nature.com/articles/s41587-022-01487-9))

### ψ Datasets

- `psi_cov10.bed`: Nanopore-based pseudouridine prediction (Dorado v0.8.2)(Dai et al., 2023)
- `BID_seq_HEK293T.bed `: BID-based pseudouridine measurements ([Dai et al., 2023](https://www.nature.com/articles/s41587-022-01505-w))

### Genome Reference

Sequence context (5-mer) extraction was based on the GRCh38 genome build:
/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa


