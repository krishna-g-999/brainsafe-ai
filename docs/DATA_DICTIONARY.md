# Data dictionary: compound_library.csv

One row per unique compound (keyed by InChIKey of the standardised parent structure). Empty cells mean *not measured*, never imputed.

| Column | Type | Definition |
|---|---|---|
| inchikey | text | Standard InChIKey of the canonical parent (salt-stripped, largest fragment). Primary key. |
| canonical_smiles | text | RDKit canonical SMILES of the parent structure. |
| name | text | Compound name where known (from the clinical reference set). |
| BBB_label | 0/1 | Measured active(1)/inactive(0) call for blood-brain-barrier permeability; pChEMBL>=6 active, <5 inactive. |
| BBB_pchembl | float | Median measured pChEMBL for blood-brain-barrier permeability. |
| AChE_label | 0/1 | Measured active(1)/inactive(0) call for acetylcholinesterase (Alzheimer); pChEMBL>=6 active, <5 inactive. |
| AChE_pchembl | float | Median measured pChEMBL for acetylcholinesterase (Alzheimer). |
| BChE_label | 0/1 | Measured active(1)/inactive(0) call for butyrylcholinesterase (Alzheimer); pChEMBL>=6 active, <5 inactive. |
| BChE_pchembl | float | Median measured pChEMBL for butyrylcholinesterase (Alzheimer). |
| BACE1_label | 0/1 | Measured active(1)/inactive(0) call for beta-secretase 1 (Alzheimer); pChEMBL>=6 active, <5 inactive. |
| BACE1_pchembl | float | Median measured pChEMBL for beta-secretase 1 (Alzheimer). |
| GSK3B_label | 0/1 | Measured active(1)/inactive(0) call for glycogen synthase kinase-3 beta (neuroprotection); pChEMBL>=6 active, <5 inactive. |
| GSK3B_pchembl | float | Median measured pChEMBL for glycogen synthase kinase-3 beta (neuroprotection). |
| MAO_A_label | 0/1 | Measured active(1)/inactive(0) call for monoamine oxidase A (depression); pChEMBL>=6 active, <5 inactive. |
| MAO_A_pchembl | float | Median measured pChEMBL for monoamine oxidase A (depression). |
| MAO_B_label | 0/1 | Measured active(1)/inactive(0) call for monoamine oxidase B (Parkinson); pChEMBL>=6 active, <5 inactive. |
| MAO_B_pchembl | float | Median measured pChEMBL for monoamine oxidase B (Parkinson). |
| D2_label | 0/1 | Measured active(1)/inactive(0) call for dopamine D2 receptor; pChEMBL>=6 active, <5 inactive. |
| D2_pchembl | float | Median measured pChEMBL for dopamine D2 receptor. |
| A2A_label | 0/1 | Measured active(1)/inactive(0) call for adenosine A2A receptor; pChEMBL>=6 active, <5 inactive. |
| A2A_pchembl | float | Median measured pChEMBL for adenosine A2A receptor. |
| HT2A_label | 0/1 | Measured active(1)/inactive(0) call for serotonin 5-HT2A receptor; pChEMBL>=6 active, <5 inactive. |
| HT2A_pchembl | float | Median measured pChEMBL for serotonin 5-HT2A receptor. |
| SERT_label | 0/1 | Measured active(1)/inactive(0) call for serotonin transporter; pChEMBL>=6 active, <5 inactive. |
| SERT_pchembl | float | Median measured pChEMBL for serotonin transporter. |
| hERG_label | 0/1 | Measured active(1)/inactive(0) call for hERG channel (cardiac safety); pChEMBL>=6 active, <5 inactive. |
| hERG_pchembl | float | Median measured pChEMBL for hERG channel (cardiac safety). |
| antioxidant_pIC50 | float | Median measured DPPH radical-scavenging pIC50. |
| n_endpoints_measured | int | Number of endpoints with a measured value for this compound. |
| role | text | training_label (has a measured endpoint) or reference (clinical structure only). |
| sources | text | Semicolon-separated measured-data sources contributing to this compound. |
| clinical_max_phase | float | Highest clinical phase reached (ChEMBL ATC-N reference), if applicable. |
| atc | text | Anatomical Therapeutic Chemical code (reference set). |
| disease | text | Disease grouping from the ATC-N reference set. |
| is_flavonoid_core | bool | True if the structure matches a flavonoid/polyphenol core substructure (flavone/flavonol, isoflavone, flavanone, flavan-3-ol, chalcone). |
| mw | float | Molecular weight. |
| clogp | float | Crippen calculated logP. |
| tpsa | float | Topological polar surface area. |
| hbd | int | Hydrogen-bond donors. |
| hba | int | Hydrogen-bond acceptors. |
| rotatable_bonds | int | Rotatable bond count. |
| aromatic_rings | int | Aromatic ring count. |
| fraction_csp3 | float | Fraction of sp3-hybridised carbons. |
| qed | float | Quantitative estimate of drug-likeness. |

## Standardisation

SMILES are parsed with RDKit; the largest organic fragment is kept (salt/counter-ion removal), the structure is sanitised, and the canonical SMILES and InChIKey are computed. Deduplication is by full standard InChIKey.

## Flavonoid-core substructures

- flavone/flavonol: `O=c1cc(-c2ccccc2)oc2ccccc12`
- isoflavone: `O=c1c(-c2ccccc2)coc2ccccc12`
- flavanone: `O=C1CC(c2ccccc2)Oc2ccccc21`
- flavan-3-ol: `OC1Cc2ccccc2OC1c1ccccc1`
- chalcone: `O=C(C=Cc1ccccc1)c1ccccc1`