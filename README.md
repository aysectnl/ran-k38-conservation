# Evolutionary conservation of Ran K38 across eukaryotes

A small phylogenetics / sequence-analysis project by Ayşe Çetinel.

## Motivation

In a protein structure course project I modeled the Ran–XPO1 complex with AlphaFold3 and ChimeraX and introduced a Ran K38A point mutation. In the model, the Lys38–Glu840 distance increased from 3.55 Å to 6.50 Å. Structural change alone does not tell whether a residue matters, so this project asks an evolutionary question:

**Is Lys38 conserved across Ran orthologs, and do Ran sequences form a distinct clade separate from related small GTPases?**

The project also serves as a hands-on exercise in multiple sequence alignment, maximum-likelihood phylogenetics and profile-HMMs.

## Workflow

| Step | Task | Tool |
|---|---|---|
| 1 | Sequence collection | UniProtKB (Swiss-Prot) |
| 2 | Multiple sequence alignment | MAFFT |
| 3 | Maximum-likelihood tree | IQ-TREE (ModelFinder) |
| 4 | Profile-HMM | HMMER (hmmbuild, hmmsearch) |
| 5 | K38 conservation analysis | Python script |

## Step 1: Sequence collection

**Source:** UniProtKB, reviewed (Swiss-Prot) entries only, query `gene:RAN`.

**Why reviewed entries only?** Manually curated sequences have fewer annotation errors than unreviewed TrEMBL entries (808 results vs. 23 reviewed), which matters for alignment and tree quality.

**Why 9 species?** The search returned 23 reviewed entries. I excluded Ran-interacting proteins (RanGAP, RCC1, NTF2) and the paralog Ran-like (RANL_DROME), because they are not Ran orthologs. Among the remaining Ran entries, many are near-identical mammals or fish, which add redundancy without adding phylogenetic information. I kept one representative per major lineage to span a broad evolutionary range:

| Entry name | UniProt ID | Group |
|---|---|---|
| RAN_HUMAN | P62826 | Mammal |
| RAN_MOUSE | P62827 | Mammal |
| RAN_CHICK | P42558 | Bird |
| RAN_XENLA | P52301 | Amphibian |
| RAN_DANRE | P79735 | Fish |
| RAN_DROME | Q9VZ23 | Insect |
| RAN_CAEEL | O17915 | Nematode |
| RANA1_TOBAC | P41918 | Plant |
| GSP1_NEUCR | Q7RVL0 | Fungus |

**Outgroup:** HRAS (P01112), RAB1A (P62820) and ARF1 (P84077), all human, reviewed entries.

Why an outgroup? It roots the tree, tests whether Ran orthologs form their own clade, and later lets us check that a Ran profile-HMM separates Ran from non-Ran sequences.

Why these three? Ran belongs to the Ras superfamily of small GTPases. All members descend from a common ancestral gene that duplicated and diversified into subfamilies (Ras, Rab, Arf, Rho, Ran). I picked one well-characterized representative from three other major subfamilies. They are homologous to Ran, so they can be aligned reliably, but divergent enough to be expected to form separate clades. The selection criterion is shared ancestry (homology), not GTP binding, XPO1 interaction or conservation level. Conservation is a result of this project, so using it as a selection criterion would be circular. This choice is a practical one; other subfamily members (for example RHOA) would work equally well.

**Output:** `data/all.fasta`, 12 sequences (9 Ran + 3 outgroup), downloaded as a single UniProt FASTA (canonical) file. Lengths: Ran 215-216 aa; outgroup ARF1 181 aa, HRAS 189 aa, RAB1A 205 aa. Note that HRAS's UniProt entry name is RASH_HUMAN.

## Step 2: Alignment

**Tool:** MAFFT via the EMBL-EBI Job Dispatcher web server, default parameters, protein input, all 12 sequences (job run on 2026-10-28). Output saved as `data/all.aln.fasta`.

**Why align first?** Residue numbers differ between species because of N-terminal length differences and indels (for example zebrafish and *C. elegans* Ran are one residue shorter, so the same lysine is K37 there, not K38). To compare "the same" position across species, I use alignment columns rather than raw residue numbers.

**Sanity check:** The conserved GTPase motifs line up across all 12 sequences, including the outgroups: the P-loop (GxxxxGK[S/T]) and the DxxG motif (DTAGQ in Ran, HRAS and RAB1A). Outgroup N- and C-terminal tails align poorly (many gaps), which is expected since these regions diverged between subfamilies. This suggests the alignment of the core GTPase domain is biologically sensible.

**Note:** The MAFFT guide tree shown by the server is a fast heuristic used only to build the alignment. It is not used as a phylogenetic result; the tree in Step 3 is built separately with IQ-TREE.

## Step 3: Phylogenetic tree

**Tool:** IQ-TREE 2.0.7, run in Google Colab. Input: `data/all.aln.short.fasta`, the Step 2 alignment with FASTA headers shortened to plain names (some web servers reject special characters such as `|` and `=`).

```
iqtree2 -s data/all.aln.short.fasta -m MFP -B 1000 -T AUTO --prefix results/ran_mfp
```

ModelFinder (`-m MFP`) selected **LG+R2** as the best-fit model by BIC. Branch support is ultrafast bootstrap (UFBoot, 1000 replicates). A first pass on the Los Alamos web server with a hand-chosen LG+G4 model (that server has no automatic model selection) gave the same tree topology (`results/lanl_LG_G4.treefile`, `results/lanl_LG_G4_tree.pdf`). Support for the deep nodes was somewhat higher there (fungus + animals 74, nematode 72, versus 56 and 45 under LG+R2), so these nodes are sensitive to the model. The ModelFinder run is the main result.

**Rooting:** Rooting on the three outgroups (HRAS, RAB1A, ARF1) was not possible on the web server, because their last common ancestor is the root of the unrooted tree: the three outgroups do not form one clade. I therefore interpret the tree as unrooted.

**Results:**
- The 9 Ran sequences form one clade (UFBoot = 100), separate from HRAS, RAB1A and ARF1. The branch leading to this clade is long (about 4.1 substitutions per site), while branches inside the Ran clade are short (at most about 0.25). This supports Ran as a distinct subfamily.
- Deep splits inside Ran have weak support: fungus grouped with animals (56) and nematode placed outside insect + vertebrates (45). A species tree would place the nematode with the insect, so this node is probably unreliable. Insect + vertebrates (84) and the vertebrate/frog group (95) are better supported. I do not rely on the weakly supported nodes.
- Among vertebrates, branches are almost zero length (human and mouse Ran are identical) and zebrafish is placed inside the tetrapods (84), which differs from the species tree. The sequences are too similar to resolve their order, so I do not interpret this part.
- ARF1 sits next to the Ran clade (UFBoot = 80), but on very long branches (ARF1 about 5.3). With only three divergent outgroups, I draw no conclusion about which Ras-superfamily subfamily is closest to Ran.

**Files:** `results/ran_mfp.treefile` (Newick), `results/ran_mfp.iqtree` (IQ-TREE report).

## Step 4: Profile-HMM

**Tool:** HMMER 3.4, run in Google Colab.

**Build:** the 9 Ran sequences (without the outgroups) were aligned separately with MAFFT, so that the Ran-specific columns are not distorted by the divergent outgroups. A profile-HMM was then built from this alignment (`data/ran_only.aln.fasta`, 223 alignment columns, model length 216).

```
hmmbuild --amino --informat afa results/ran.hmm data/ran_only.aln.fasta
```

The effective sequence number is only 0.44, which reflects how similar the nine sequences are: this small set carries limited diversity, so the model describes Ran tightly.

**Search against all 12 sequences:**

```
hmmsearch --tblout results/ran_hmmsearch.tbl -E 10 results/ran.hmm data/all.targets.fasta
```

| Sequence | Bit score | E-value |
|---|---|---|
| RAN_HUMAN, RAN_MOUSE, RAN_XENLA | 483.6 | 3.7e-149 |
| RAN_CHICK | 483.3 | 4.7e-149 |
| RAN_DANRE | 475.0 | 1.6e-146 |
| RAN_DROME | 471.7 | 1.6e-145 |
| RAN_CAEEL | 467.8 | 2.6e-144 |
| GSP1_NEUCR (fungus) | 442.4 | 1.5e-136 |
| RANA1_TOBAC (plant) | 440.3 | 6.5e-136 |
| RAB1A_HUMAN (outgroup) | 95.8 | 1.6e-30 |
| RASH_HUMAN (outgroup) | 67.6 | 6.6e-22 |
| ARF1_HUMAN (outgroup) | 28.7 | 5.5e-10 |

**Interpretation:**
- Ran sequences score 440-484 bits, the three outgroups 29-96 bits, a gap of about 345 bits.
- The nine Ran sequences were used to build the model, so their high scores are expected and are not an independent test.
- All 12 sequences pass the E-value threshold used (E <= 10), and even at a stricter cutoff such as 1e-5 the three outgroups would still be reported. This fits their homology (same superfamily). Separating Ran from other small GTPases therefore needs a bit-score cutoff, not just an E-value. Any cutoff placed between 96 and 440 bits would be chosen after seeing the data and is not validated here.

**Held-out check:** three Ran sequences that were not used to build the model were searched: rat Ran (P62828), dog Ran (P62825) and tomato Ran1 (P38546). They scored 483.5, 483.5 and 434.3 bits (E-values 1e-149 or better, 1.1e-134 for tomato). Rat and dog are almost identical to human, so they mainly confirm the search works. Tomato is in the same plant family as tobacco, so it is only a modest test of generalisation.

**Files:** `results/ran.hmm`, `results/ran_hmmsearch.tbl`, `results/test_hmmsearch_new.tbl`.

**Notebook:** `notebooks/ran_analysis.ipynb` holds the HMMER and IQ-TREE commands and outputs as run in Google Colab. In the notebook the input files sit in the working directory (flat layout), while this repository groups them into `data/` and `results/`, so paths differ.

## Step 5: K38 conservation

**Method:** `scripts/residue_column.py` maps human Ran residue 38 to its alignment column (column 51 in `data/all.aln.fasta`) and reports the residue in that column for every sequence.

```
python scripts/residue_column.py data/all.aln.fasta RAN_HUMAN 38
```

**Result:**

| Sequence | Residue at column 51 | Own residue number |
|---|---|---|
| RAN_HUMAN, RAN_MOUSE, RAN_CHICK, RAN_XENLA, RAN_DROME | K | 38 |
| RAN_DANRE, RAN_CAEEL | K | 37 |
| GSP1_NEUCR (fungus) | K | 35 |
| RANA1_TOBAC (plant) | K | 41 |
| RASH_HUMAN (HRAS, outgroup) | E | 31 |
| RAB1A_HUMAN (outgroup) | S | 39 |
| ARF1_HUMAN (outgroup) | T | 45 |

**Interpretation:** Lys38 is conserved in all 9 Ran sequences, spanning animals, a fungus and a plant, and the lysine just before it (K37 in human) sits in the same `EFEKKY` motif. In the three outgroups the same column holds E, S and T, so this lysine looks specific to the Ran subfamily rather than a general feature of Ras-superfamily GTPases. Residue numbers differ between species (35 to 41), which is why columns, not raw numbers, were compared. This supports the choice of K38 as a functionally relevant residue for the K38A mutation modeled in the Ran–XPO1 structure project.

**Caveat:** This is based on 9 Ran sequences. The alignment around this region is less reliable for the outgroups, which diverged strongly from Ran, so the "Ran-specific" statement is a preliminary observation, to be checked against the phylogeny (Step 3) and a larger sequence set.

## Limitations

- Small dataset: 9 Ran sequences and 3 outgroups. The vertebrate Ran sequences are nearly identical, so they add little phylogenetic information.
- The tree is unrooted, the outgroups are few and divergent (long branches), and several deep nodes inside Ran have weak support.
- One alignment method with default settings, no trimming and no test of alignment uncertainty.
- The profile-HMM was built from, and first scored on, the same Ran sequences. The held-out test is small and uses close relatives, and no systematic set of negative sequences was tested.
- Conservation of Lys38 across Ran orthologs shows that the residue is conserved. It does not by itself show what the residue does in the Ran-XPO1 interface, and it does not validate the structural model.
- Not done: ancestral sequence reconstruction, formal tests of alternative trees, and experimental validation.

## Software

IQ-TREE 2.0.7, HMMER 3.4, MAFFT (default settings; version not recorded), Python 3 (see `scripts/`). UniProtKB sequences retrieved on 2026-10-28.
