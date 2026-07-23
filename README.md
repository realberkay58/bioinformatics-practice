# 🏥 Clinical Decision Support System (CDSS) for Genomic & Transcriptomic Analytics

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![Streamlit](https://img.shields.io/badge/Streamlit-App-red.svg)](https://streamlit.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A state-of-the-art, web-based **Clinical Decision Support System (CDSS)** designed for translational bioinformatics and precision oncology. This platform bridges the gap between raw transcriptomic RNA-Seq data and actionable clinical insights, offering real-time multi-source data integration, machine learning-driven diagnostics, and interactive genomic visualizations.

---

## ✨ Key Features & Architecture

* **🧬 Genomic Master Dataset Explorer:** High-performance tabular viewer with interactive gene selection and instant cross-referencing capabilities.
* **🩺 Clinical Variant Integration (MyVariant / ClinVar):** Real-time API querying to retrieve pathogenicity statuses (Pathogenic, Benign, VUS) and official NCBI ClinVar records for selected biomarkers.
* **🤖 AI Diagnostics & Predictive Insights:** Machine Learning classification (Random Forest) to evaluate patient risk levels, diagnostic confidence scores, and feature importance (Biomarker Impact Weights).
* **📊 Advanced Expression Analytics:** 
  * **Volcano Plot:** Dynamic scatter visualization highlighting statistical significance vs. fold-change, with custom target markers (`🎯`).
  * **Cohort Expression Heatmap:** Clear separation between Control and Patient cohorts using customized expression matrices.
* **🌐 Global Cross-Validation (NCBI GEO):** Live integration with NCBI's Gene Expression Omnibus (GEO) to fetch external cohort metadata for independent result validation.
* **🔬 Molecular Data & Structural Biology:** 
  * **Evolutionary Tree (Phylogeny):** Automated ortholog alignment and phylogenetic tree generation using MUSCLE and BioPython.
  * **3D Protein Structure Viewer:** Interactive 3D visualization powered by 3Dmol.js, mapping amino acid positions and structural contexts directly from UniProt/PDB.
* **🎯 Ultimate Sticky UI/UX:** Custom CSS-injected modern glassmorphism sticky navigation tabs and card-based responsive panels preventing text overlap during deep scrolling.

---

## 🛠️ Tech Stack

* **Frontend & UI:** Streamlit, Plotly, Custom HTML/CSS (Glassmorphism & Sticky Headers)
* **Bioinformatics Core:** Biopython, GEOparse, SciPy, Pandas, NumPy
* **APIs & Databases:** NCBI Entrez, Ensembl REST API, MyVariant.info, UniProt / PDB

---

## 🚀 Installation & Local Setup

Follow these steps to run the application locally on your machine:

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/YOUR_USERNAME/YOUR_REPOSITORY_NAME.git](https://github.com/YOUR_USERNAME/YOUR_REPOSITORY_NAME.git)
   cd YOUR_REPOSITORY_NAME

   pip install -r requirements.txt

   streamlit run app.py

├── app.py                      # Main Streamlit application & UI layout
├── analiz.py                   # RNA-Seq differential expression processing pipeline
├── results/                    # Generated datasets (CSV) and dynamic assets (trees)
├── data_raw/                   # Raw GEO data cache
└── requirements.txt            # Project dependencies

How to Use:
Select a Gene: Choose a target gene from the sidebar or the Master Dataset tab.

Review the Header: Track the real-time AI Confidence Score, Primary Biomarker impact, and Patient Risk status fixed at the top.

Analyze Results: Scroll through Clinical Variants, AI Reports, Expression Analytics, Global GEO Searches, and Molecular Structures based on the selected gene.