# BrainSafe AI: An Interactive Translational Informatics Platform for Neuroprotective Compound Profiling in Neurodegenerative Diseases

---

## Authors

Krishnasalini Gunanathan¹, Venketesh Sivaramakrishnan¹*

¹ Department of [Your Department], Sri Sathya Sai Institute of Higher Learning, Puttaparthi, Andhra Pradesh, India

*Corresponding author: [prof.email@sssihl.edu.in]

---

## Abstract

Neurodegenerative diseases — including Alzheimer's disease (AD), Parkinson's disease (PD), amyotrophic lateral sclerosis (ALS), and Huntington's disease (HD) — collectively affect hundreds of millions of individuals worldwide, yet effective disease-modifying therapies remain limited. A growing body of peer-reviewed evidence implicates dietary flavonoids, vitamins, nutraceuticals, and approved pharmacological agents in neuroprotection through mechanisms including antioxidant activity, anti-inflammatory signalling, mitochondrial support, and protein aggregation modulation. However, no publicly accessible interactive tool currently enables researchers or educators to explore these compounds systematically across multiple neuroprotective dimensions.

We present **BrainSafe AI**, a Streamlit-based web application developed as the translational educational module of SAI-Net (Structure-Activity Intelligence Network) at Sri Sathya Sai Institute of Higher Learning. The platform hosts a curated database of 128 compounds — spanning flavonoids, polyphenols, vitamins, amino acids, nutraceuticals, and approved drugs — each annotated with seven neuroprotective scores, blood-brain barrier (BBB) penetration status, mechanistic pathway data, downstream metabolite profiles, brain region targeting, and relevance scores for four neurodegenerative diseases. The platform provides (i) seven-dimensional radar charts for compound neuroprotection visualisation, (ii) interactive directed-edge network graphs mapping compound-to-pathway-to-metabolite-to-disease relationships, and (iii) a class-based inference engine that generates estimated neuroprotective profiles for compounds not present in the database. BrainSafe AI is freely accessible at [INSERT DEPLOYED URL] and the source code is available at https://github.com/krishna-g-999/brainsafe-ai.

**Keywords:** neuroprotection, flavonoids, nutraceuticals, neurodegenerative diseases, translational informatics, Alzheimer's disease, Parkinson's disease, ALS, Huntington's disease, network pharmacology

---

## 1. Motivation and Significance

Neurodegenerative diseases impose an escalating global burden, with Alzheimer's disease alone projected to affect over 130 million individuals by 2050 (WHO, 2023). Despite significant investment in drug discovery, the clinical pipeline for these conditions remains sparse, and approved disease-modifying therapies are largely absent for most conditions in this spectrum. Complementary and integrative approaches — particularly the neuroprotective potential of naturally occurring compounds — have attracted substantial research interest, with thousands of peer-reviewed studies (PubMed/PMC, 2020–2026) reporting mechanistic evidence for compounds ranging from flavonoids such as quercetin and fisetin to approved agents such as donepezil and memantine.

A critical translational gap exists between this distributed literature and the ability of researchers, clinicians, and students to access and compare compound profiles interactively. Existing databases such as ChEMBL and DrugBank provide rich chemical and pharmacological annotations, but lack (i) neuroprotection-specific scoring frameworks, (ii) visual multi-dimensional profiling, and (iii) disease-relevant network representations tailored to neurodegeneration.

BrainSafe AI addresses this gap by providing a unified, interactive, and publicly accessible platform for systematic exploration of neuroprotective compounds. The platform was developed within the SAI-Net translational module framework, which applies structure-activity intelligence and network pharmacology principles to neurodegenerative disease research at Sri Sathya Sai Institute of Higher Learning.

---

## 2. Software Description

### 2.1 Architecture and Technology Stack

BrainSafe AI is implemented in Python 3 using the Streamlit framework (v1.x), enabling rapid deployment of an interactive web interface without requiring client-side JavaScript. The application is hosted on Replit's autoscale infrastructure and is accessible via a public URL. The core components are:

- **`app.py`** — Main application module containing the user interface, search logic, radar chart generation, and network graph construction
- **`compounds.json`** — Curated compound database (128 entries) with structured annotations
- **`scorer.py`** — Composite neuroprotective scoring module

The application requires no installation by end users and runs entirely in a standard web browser.

### 2.2 Compound Database

The compound database comprises 128 entries curated from three primary sources:

| Source | Scope | Compounds Derived |
|--------|-------|------------------|
| PubMed/PMC (2020–2026) | Peer-reviewed mechanistic literature | Scoring, pathways, disease relevance |
| ChEMBL | Bioactivity and pharmacological database | Approved drugs, bioactivity annotations |
| DrugBank | Approved drug database | Clinical drugs, BBB data |

The 128 compounds span the following classes: flavonoids (quercetin, luteolin, fisetin, apigenin, kaempferol, myricetin, rutin, icariin, baicalin, chrysin, naringin, hesperidin, and others), polyphenols and curcuminoids (curcumin, resveratrol, pterostilbene, ellagic acid, rosmarinic acid, ferulic acid, caffeic acid, carnosic acid, thymoquinone, and others), vitamins and cofactors (vitamins B1, B2, B3, B6, B9, B12, C, D, E, K2, nicotinamide, pantothenate, and others), amino acids and derivatives (glycine, taurine, L-carnitine, NAC, ergothioneine, glutathione, creatine, and others), nutraceuticals (coenzyme Q10, ubiquinol, idebenone, PEA, sesamin, oleic acid, DHA, EPA, and others), and approved pharmacological agents (donepezil, rivastigmine, galantamine, memantine, rasagiline, selegiline, valproic acid, and others).

Each compound entry contains the following structured fields:

- **BBB penetration** — categorical (High / Moderate / Low / Very Low), derived from literature and ChEMBL annotations
- **Seven neuroprotective dimensions** — scored 0–10: antioxidant activity, anti-inflammatory activity, mitochondrial support, protein aggregation modulation, cognitive enhancement, neurogenesis support, synaptic plasticity modulation
- **Composite neuroprotective score** — weighted composite (0–100) emphasising antioxidant, anti-inflammatory, mitochondrial, and aggregation dimensions
- **Molecular pathways** — list of activated or inhibited signalling pathways (e.g., Nrf2/GSH, NF-κB, AMPK, SIRT1, BDNF/TrkB, PI3K/Akt, mTOR, CREB)
- **Downstream metabolites and biomarkers** — e.g., GSH, ATP, BDNF, IL-6, TNF-α, acetylcholine
- **Brain regions targeted** — e.g., hippocampus, prefrontal cortex, substantia nigra, cerebellum
- **Disease relevance** — four-level classification (High / Med / Low / None) for ALS, Alzheimer's disease, Parkinson's disease, and Huntington's disease

### 2.3 Composite Neuroprotective Scoring

The composite neuroprotective score (NPS, 0–100) is computed as a weighted sum emphasising the four mechanistic dimensions with the strongest evidence base in the neurodegeneration literature:

```
NPS = min(100, (antioxidant × 3 + anti_inflammatory × 3 + mitochondrial × 2 + aggregation_modulation × 2) / 40 × 100)
```

This weighting reflects the predominance of oxidative stress and neuroinflammation as shared pathological mechanisms across all four neurodegenerative disease spectra considered.

### 2.4 Visualisation Modules

**Radar Charts.** A seven-dimensional radar (spider) chart is rendered using Plotly for each compound, displaying scores across: Antioxidant, Anti-Inflammatory, Mitochondrial, Aggregation Mod., Cognitive, Neurogenesis, and Synaptic Plasticity axes. A reference ring at score 5 is displayed to contextualise relative strengths and weaknesses.

**Directed-Edge Network Graphs.** An interactive network graph is generated for each compound using a concentric ring layout:
- Centre node: the query compound
- Inner ring (radius 1.85 units): molecular pathways activated or inhibited
- Middle ring (radius 3.6 units): downstream metabolites and biomarkers
- Outer ring (radius 5.3 units): targeted brain regions and disease nodes

Edges are colour-coded by relationship type: blue (compound → pathway), green (pathway → metabolite), red (metabolite → disease), amber (pathway → brain region). Arrowheads indicate directionality of influence.

### 2.5 Search and Class-Based Inference Engine

The platform supports two search modes:

1. **Dropdown selection** — direct lookup of any of the 128 database compounds
2. **Free-text search** — fuzzy string matching (difflib SequenceMatcher, similarity threshold 0.65) to handle spelling variants and partial names

For queries that do not match any database entry at the 0.65 similarity threshold, the platform activates a **class-based inference engine**. The query name is classified into one of nine chemical classes (Flavonoid, Polyphenol, Amino Acid, Alkaloid, Terpenoid, Fatty Acid, Vitamin, Approved Drug, General Nutraceutical) using heuristic name-pattern matching. A class-level template profile is then returned with typical score ranges, representative pathways, metabolites, and brain regions for that class. Estimated reports are clearly flagged with a visible "ESTIMATED" badge and an explanatory disclaimer, distinguishing them from database-derived reports.

---

## 3. Illustrative Examples

### 3.1 Quercetin — High-confidence Database Entry

Searching for "quercetin" returns the full database report: NPS = 97/100, BBB = Moderate, highest scores in antioxidant (9.5) and anti-inflammatory (9.0) dimensions. The network graph maps connections through Nrf2/GSH and NF-κB pathways to downstream metabolites (GSH, IL-6, TNF-α) and disease nodes (Alzheimer's disease: High relevance, Parkinson's disease: High relevance).

### 3.2 Donepezil — Approved Drug

Donepezil (ChEMBL ID: DB00843) returns: NPS = 72/100, BBB = High, primary mechanism through AChE inhibition pathway, high relevance for Alzheimer's disease, low relevance for ALS and Huntington's disease. The radar chart reveals a distinct cognitive enhancement peak (8.5) with moderate scores in other dimensions, reflecting its selective cholinergic mechanism.

### 3.3 An Unlisted Compound — Class Inference

Querying a compound not present in the database (e.g., a novel flavanone) triggers the class-based inference engine. The system identifies the compound as a Flavonoid class based on the suffix pattern and returns a class-template profile with typical flavonoid ranges, common Nrf2/NF-κB pathway assignments, and an "ESTIMATED" flag, enabling preliminary orientation while clearly communicating the non-compound-specific nature of the estimate.

---

## 4. Impact

BrainSafe AI serves three primary audiences:

**Research community.** The structured 128-compound database with consistent seven-dimensional scoring enables systematic cross-compound comparisons that are not feasible through manual literature review. Researchers can rapidly identify compounds with specific mechanistic profiles — for example, compounds with high mitochondrial support and aggregation modulation scores for ALS-focused studies.

**Educators and students.** The interactive visualisations (radar charts and network graphs) provide intuitive representations of complex neuropharmacological relationships suitable for graduate-level neuroscience and pharmacology education.

**Translational and clinical context.** By integrating approved drugs alongside nutraceuticals and flavonoids within the same scoring framework, the platform facilitates direct comparison of natural and pharmacological agents across shared mechanistic dimensions — supporting evidence-based discussions around adjunctive or complementary approaches.

The platform is particularly relevant in the context of growing interest in dietary and lifestyle interventions for brain health, where a translational tool bridging basic science literature and accessible public health education has been absent.

---

## 5. Conclusions

BrainSafe AI provides an open-access, interactive, and scientifically grounded platform for neuroprotective compound profiling. The curated 128-compound database, seven-dimensional scoring framework, directed-edge network visualisations, and class-based inference engine collectively address a clear gap in translational neuroscience informatics. The platform is freely accessible and designed to support research, education, and evidence-based public health communication in the neurodegenerative disease domain.

Future development directions include: (i) expansion of the compound database beyond 128 entries; (ii) integration of SAI-Net multiomics analytical APIs for real-time pathway predictions; (iii) incorporation of protein docking and structure-activity relationship modules; and (iv) user-contributed compound annotations with editorial review.

---

## Conflict of Interest

The authors declare no conflict of interest.

## Funding

[Add any funding sources or state: This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.]

---

## References

1. World Health Organization. (2023). *Dementia*. WHO Global Health Estimates.
2. GBD 2019 Dementia Forecasting Collaborators. (2022). Estimation of the global prevalence of dementia in 2019 and forecasted prevalence in 2050. *The Lancet Public Health*, 7(2), e105–e125.
3. Bhullar, K. S., et al. (2022). Polyphenols and neurological disorders. *Frontiers in Nutrition*, 9, 829500.
4. Fišar, Z. (2023). Pitfalls in the quest for brain biomarkers of neurodegenerative diseases. *International Journal of Molecular Sciences*, 24(3), 2393.
5. Guo, C., et al. (2022). Oxidative stress, mitochondrial dysfunction, and the mitochondria theory of aging. *Cells*, 11(7), 1115.
6. Jenwitheesuk, A., et al. (2020). Melatonin regulates aging and neurodegeneration through energy metabolism, epigenetics, autophagy and circadian rhythm pathways. *International Journal of Molecular Sciences*, 21(16), 5700.
7. Knekt, P., et al. (2002). Flavonoid intake and risk of chronic diseases. *American Journal of Clinical Nutrition*, 76(3), 560–568.
8. Gao, X., et al. (2012). Habitual intake of dietary flavonoids and risk of Parkinson disease. *Neurology*, 78(15), 1138–1145.
9. Abubakar, M. B., et al. (2023). Flavonoids as potential therapeutic agents against neurodegenerative diseases. *Molecules*, 28(6), 2624.
10. Mendez-David, I., et al. (2021). Towards a biology of resilience: Brain and behavioral mechanisms contributing to healthy aging. *Neuroscience and Biobehavioral Reviews*, 124, 1–8.
11. Wishart, D. S., et al. (2022). DrugBank 5.0: a major update to the DrugBank database for 2018. *Nucleic Acids Research*, 46(D1), D1074–D1082.
12. Mendez, D., et al. (2019). ChEMBL: towards direct deposition of bioassay data. *Nucleic Acids Research*, 47(D1), D930–D940.
13. Pushpakom, S., et al. (2019). Drug repurposing: progress, challenges and recommendations. *Nature Reviews Drug Discovery*, 18(1), 41–58.
14. Nakagawa, S., & Deli, M. A. (2020). Blood–brain barrier: An overview. *Neuroscience & Biobehavioral Reviews*, 118, 81–98.

---

## Software Availability

| Item | Description |
|------|-------------|
| Name | BrainSafe AI |
| Permanent DOI | [To be assigned via Zenodo after GitHub release] |
| Version | 1.0 |
| Licence | MIT |
| Developer | Krishnasalini Gunanathan |
| Programming language | Python 3 |
| Software requirements | Web browser (no installation required) |
| Live application | [INSERT DEPLOYED URL] |
| Source code repository | https://github.com/krishna-g-999/brainsafe-ai |
| Support email | [your email] |
