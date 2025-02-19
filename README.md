# **GenBankRefs**

## **First-Time Setup**
When running for the first time, when prompted **"run blast(y/n)"**, enter **"y"** to build the BLAST database.

---

## **📖 Terminology**
### **GenBank**
- **GenBank Record** → A single GenBank sequence entry.
- **GenBank Reference** → A reference within a GenBank record (each record may contain multiple references).
- **Record Features** → Metadata of a GenBank record (excluding gene information).
- **Record Genes** → Gene ranges and translated protein sequences.

### **AI Integration**
- **AI Extraction** → AI-assisted full-text review (currently using **GPT-4o**) to classify literature relevance and extract relevant sequence metadata and GenBank accession number.

---

## **📂 Repository Structure and File Descriptions**

### **Main Directories**
| Directory | Description |
|-----------|------------|
| **`ReferenceData/`** | Stores BLAST results and GenBank reference sequences (nucleotide & protein) for each gene in each virus. |
| **`phylogenetic/`** | R scripts for generating **phylogenetic trees** for each virus. |
| **`OutputData/`** | Contains processed outputs for each virus, including GenBank data, AI classification results, and metadata. |
| **`viruses/`** | Class methods and functions for handling **virus-specific** data, including sequence selection for **phylogenetic analysis**. |
| **`PubMed/`** | Stores **AI-extracted** PubMed references, manually curated references, and hard-linked matches. |

---

### **Processed Output Files (`OutputData/{virus}/`)**
| File | Description |
|------|------------|
| **`{virus}_GenBankGenes.xlsx`** | Maps GenBank accession numbers to gene names and CDS regions. |
| **`{virus}.db`** | The generated **SQLite database** containing all processed **GenBank & PubMed** information. |
| **`{virus}_fixed_PubID.csv`** | Fixes publication IDs, mapping **internal IDs to PMID/DOI**. |
| **`{virus}_fixed_refID.csv`** | Fixes reference IDs, linking **GenBank submission sets** to their respective metadata. |
| **`{virus}_datalog_pubmed_workflow.txt`** | Reviewer & GPT title-abstract classification workflow: In short: if **Reviewer 1** marked an article as "likely," **GPT** further classified it as "yes" or "no" for full-text sequence inclusion. If disagreement remained, **Reviewer 2 (R2)** reviewed it to resolve. |
| **`{virus}_datalog_compare_matched.txt`** | Statistics on **matched GenBank & PubMed sequences**. |
| **`{virus}_genbank.txt`** | Statistics on **all GenBank submissions**. |
| **`{virus}_pubmed.txt`** | Statistics on **all AI-extracted Literature Information**. |

---

### **PubMed Extraction Files (`PubMed/`)**
| File | Description |
|------|------------|
| **`ReferenceSummary_{Date}.xlsx`** | AI extraction results from **full-text PubMed papers** after initial review. |
| **`ReferenceSummary_Genbank_{Date}.xlsx`** | **GenBank submissions** with PMIDs **not included** in initial PubMed review. |
| **`ReferenceSummary_PubMed_Missing_{Date}.xlsx`** | **GenBank submissions missing PMIDs**, additional papers found via manual search. |
| **`Reference_Hardlink_{Date}.xlsx`** | Manually linked papers without extracted accessions (matched manually). |

---

### **Core Python Scripts**
| Script | Description |
|--------|------------|
| **`Main.py`** | The **main workflow** that orchestrates data processing. |
| **`genbank_records.py` & `GenBankFunctions.py`** | **Helper functions** for **processing, cleaning, and extracting GenBank data**. |
| **`DataFrameLogic.py`** | Dataframe **helper functions**, including **combining references from the same submission set** and **merging features with references**. |
| **`summarize_genbank.py`** | Generates `{virus}_datalog_genbank.txt`, summarizing GenBank data by **submission sets, references, sequences, and isolates**. |
| **`summarize_pubmed.py`** | Generates `{virus}_datalog_pubmed.txt` and `{virus}_datalog_pubmed_workflow.txt`. Compares **Reviewer 1 and GPT classifications**, calculates **R1 scores**, and tallies all **PubMed feature data extracted by AI**. |
| **`AI_match_paper.py`** | AI attempts to **find corresponding papers** for GenBank submissions, using **authors, submission title, journal, and metadata**. Outputs results to `GenBankRefs/OutputData/{virus}/excels/AI_cache.csv`. |
| **`database.py`** | **Creates an SQLite database** using: **Three** dataframes from GenBank (**Reference, Features, Gene**) and **Two** dataframes from PubMed (complete dataset & matched subset). |
| **`match_pubmed_GB.py`** | Matches **GenBank references** to **PubMed records** via **PMID, accession, or title** (AI-assisted). |
| **`Utilities.py` & `bioinfo.py`** | General **helper functions** for the project. |
| **`chord_diagram.py`** | Generates **chord diagrams** (currently not included in workflow). |

---

## **📌 Notes**
- **AI extraction** improves review efficiency but is **not 100% accurate**—manual review is recommended. AI matching did not give meaningful results.
- The **chord diagram** generator is currently **disabled** but can be reactivated if needed.
- Ensure that **BLAST** is run at least once to populate sequence databases.

---
