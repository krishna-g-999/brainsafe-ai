# BrainSafe AI - Brain Health Compound Explorer

## Overview
A scientifically accurate Streamlit web app for exploring 108 nutraceuticals, supplements, phytochemicals, flavonoids, vitamins, amino acids, and drugs and their links to brain health. Part of the SAI-Net Translational Module at Sri Sathya Sai Institute of Higher Learning.

## Features
- SAI-Net logo embedded in the header alongside the title (base64 encoding)
- Fuzzy / case-insensitive compound search (difflib)
- Brain Health Radar chart (Plotly) — 7 dimensions with reference ring and value annotations
- Interactive Pathway & Metabolite Network graph with directed arrows and color-coded edge types
- BBB penetration ratings, Neuroprotective scores (0-100)
- Bioactivity profile table with progress bars
- Molecular pathway and metabolite/biomarker chips
- Neurodegenerative disease relevance table (ALS, Alzheimer's, Parkinson's, Huntington's)
- Browse all compounds grouped by type (tab)
- About SAI-Net page with research credits
- Professional navy/gold design — Inter font, no decorative elements

## Compound Database (108 entries)
### Flavonoids
- Flavonols: Quercetin, Kaempferol, Fisetin, Myricetin, Galangin, Isorhamnetin, Rutin, Taxifolin
- Flavones: Apigenin, Luteolin, Baicalein, Wogonin, Nobiletin
- Flavanones: Hesperetin, Naringenin, Eriodictyol
- Flavan-3-ols: EGCG, Catechin, Epicatechin
- Isoflavones: Daidzein, Genistein
- Anthocyanins: Cyanidin, Delphinidin
- Stilbenoids: Resveratrol, Pterostilbene

### Vitamins & Cofactors
- B1 (Thiamine), B2 (Riboflavin), B3 (Niacin), B5 (Pantothenate), B6 (Pyridoxine), B7 (Biotin), B9 (Folate), B12
- Vitamins C, D, E

### Supplements & Nutraceuticals
- CoQ10, Ubiquinol (reduced CoQ10), Idebenone (synthetic analog)
- NR (Nicotinamide Riboside), PQQ, Creatine, Spermidine, Spermine
- Magnesium L-Threonate, Citicoline, Alpha-GPC, Huperzine A
- Palmitoylethanolamide (PEA), Saffron, Gotu Kola, Honokiol, Sesamin, Oleocanthal

### Amino Acids
- Taurine, Glycine, N-Acetylcysteine, L-Carnitine, Acetyl-L-Carnitine
- Cystine, L-Tyrosine, L-Tryptophan, 5-HTP, L-Glutamine, L-Arginine, L-Methionine
- Ornithine, GABA, Carnosine, SAMe
- Glutathione (GSH)

### Herbs & Botanicals
- Ashwagandha, Ginkgo Biloba, Bacopa Monnieri, Lion's Mane, Rhodiola
- Ursolic Acid, Berberine, Sulforaphane, Allicin, Piperine

### Fatty Acids & Lipids
- Omega-3 (DHA/EPA), Oleic Acid, Alpha-Lipoic Acid
- BHB (Beta-Hydroxybutyrate), Phospholipids

### Inositol & Small Molecules
- Inositol, Caffeine, Melatonin

### Approved Drugs (from DrugBank / ChEMBL)
- Metformin, Memantine, Galantamine, Vinpocetine, Lithium Orotate, Donepezil, Rivastigmine

## Network Graph
- Concentric ring layout: compound (center) → pathways (ring r=1.85) → metabolites (ring r=3.6) → brain regions + diseases (outer ring r=5.3)
- Directed edge arrows using plotly `add_annotation` with arrowhead=2
- Edge color-coding by relationship type:
  - Blue (#3B82F6): Compound → Pathway
  - Green (#10B981): Pathway → Metabolite
  - Red (#EF4444): Metabolite → Disease
  - Amber (#F59E0B): Pathway → Brain Region
- Dashed concentric ring guides with ring labels

## Tech Stack
- Python 3.11
- Streamlit
- Plotly (radar chart, network graph with annotations)
- NetworkX (graph utilities)
- difflib (fuzzy matching)
- Local JSON data — no API, no DB, no ML
- Inter font (Google Fonts)

## Radar Chart Dimensions
1. Antioxidant Activity
2. Anti-Inflammatory
3. Mitochondrial Support
4. Aggregation Modulation
5. Cognitive Enhancement
6. Neurogenesis Support
7. Synaptic Plasticity

## Scoring Formula
raw = antioxidant*3 + anti_inflammatory*3 + mitochondrial_support*2 + aggregation_modulation*2
score = min(100, (raw / 40) * 100)
green >= 70, orange 40-69, red < 40

## Files
- app.py — main Streamlit app
- compounds.json — full compound database (108 entries)
- scorer.py — neuroprotective scoring formula
- sai_net_logo.png — SAI-Net logo (base64 embedded in header)
- .streamlit/config.toml — Streamlit server config

## Research Team
- PI: Prof. Venketesh Sivaramakrishnan, SSSIHL
- Developer: Krishnasalini Gunanathan, SSSIHL
