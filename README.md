# buscophylo

A streamlined workflow for conducting phylogenetic analysis based on BUSCO orthologs (protein sequences, **not** nucleotide sequences).

## 🔁 Workflow
![Workflow](images/gpa_workflow.jpg)

---
## ⚙️ Step 0: Installation
```bash
gh repo clone ypchan/buscophylo
cd buscophylo
bash setup.sh
``` 

---
## 🚩 Step 1: Merge Genome Pieces (if needed)

For some genomes, contigs or chromosomes are stored in separate FASTA files. Use `merge_assembly_pieces.py` to merge these into a single FASTA file.

```bash
cat file_lst | merge_assembly_pieces.py - -o output_directory
```

**Input file format:**

| Column 1           | Column 2                                                         |
|--------------------|------------------------------------------------------------------|
| Genome prefix      | Path to each genome piece                                        |
|GCA_021398005.1 | ./GCA_021398005.1/ncbi_dataset/data/GCA_021398005.1/GCA_021398005.fna|
|GCA_021436885.1 | ./GCA_021436885.1/ncbi_dataset/data/GCA_021436885.1/chr1.fna|
...

---

## 🧮 Step 2: Obtain Basic Genome Statistics

Assess genome assembly quality using `genome_stat.py`.

```bash
genome_stat.py -h
```

**Usage:**
```bash
ls 00_genomes/*.fna | genome_stat.py -t 4 - > genome.statistics.tsv
find 00_genomes -name "*.fna" -type f | genome_stat.py -t 4 - > genome.statistics.tsv
```

**Output columns:**

| Genome | #Contigs | GC | #N_character | Contig_longest | Contig_minimum | N90 | L90 | N50 | L50 | N75 | L75 |
|--------|----------|----|--------------|----------------|----------------|-----|-----|-----|-----|-----|-----|
|GCA_000002515.1.fna|  6|  0.3642|  44 | 2,602,197|  1,062,590 | 1,320,834 | 5 | 1,753,957 | 3 | 1,715,506|  4 |

---

## 🧬 Step 3: Run BUSCO Analysis

Evaluate genome completeness and identify BUSCO genes.

```bash
cat genome_accession.list | while read genome_abbrev; do
  echo "busco -i ${genome_dir}/${genome_abbrev}.fna \
-o ${genome_abbrev}_busco \
--out_path ${output_dir} \
--offline --cpu 8 \
--mode geno \
-l ${lineage} &> ${output_dir}/${genome_abbrev}_busco.log"
done > run_busco.sh

# Run in parallel
nohup ParaFly -c run_busco.sh -failed_cmds run_busco.sh.failed &
```

---

## 📊 Step 4: Summarize BUSCO Results

Summarize BUSCO results and filter out poor assemblies using `summary_BUSCO_results.py`.

```bash
find . -name 'short_summary.specific*txt' | summary_BUSCO_results.py - | \
sed 's/genome_label/Assembly/;s/.fna//' > busco_statistics.tsv
```

**Example output:**

| Assembly         | C     | S     | D     | F     | M     | Complete BUSCOs | Single-copy BUSCOs | Duplicated BUSCOs | Fragmented BUSCOs | Missing BUSCOs | Total BUSCO groups |
|------------------|-------|-------|-------|-------|-------|-----------------|--------------------|-------------------|-------------------|----------------|--------------------|
| GCA_013390195.1  | 97.8% | 97.5% | 0.3%  | 0.2%  | 2.0%  | 1669            | 1664               | 5                 | 3                 | 34             | 1706               |
...

---

## 🗃️ Step 5: Construct Single-Copy BUSCO Datasets

Generate datasets of single-copy BUSCOs.

```bash
ls genome_labels.list | awk '{print $1"\t../01_busco/"$1"_busco/run_fungi_odb10/full_table.tsv"}' | \
single_copy_busco_datasets.py -d ~/database/busco/fungi_odb10/links_to_ODB10.txt \
-m class_matrix.tsv -c 50 -o 02_single_copy_busco_datasets -t 8 -

# Or
single_copy_busco_datasets.py -d ~/database/busco/fungi_odb10/links_to_ODB10.txt \
-m class_matrix.tsv -c 50 -o 02_single_copy_busco_datasets -t 8 label_full_table.path.list
```

**Input file example:**
```
GCA_902806535.1    ./GCA_902806535.1_HR_busco/run_ascomycota_odb10/full_table.tsv
GCA_002246955.1    ./GCA_002246955.1_ASM224695v1_busco/run_ascomycota_odb10/full_table.tsv
...
```

---

## 🧩 Step 6: Multiple Sequence Alignment with MAFFT

Align single-copy BUSCO protein sequences.

```bash
mkdir 03_mafft
ls 02_single_copy_busco_datasets/*.faa | sed "s/02_single_copy_busco_datasets\///" | \
while read file; do
  echo "mafft --thread 4 --auto 02_single_copy_busco_datasets/${file} > 03_mafft/${file%.faa}.mafft.faa 2> /dev/null"
done > run_mafft.sh

nohup ParaFly -c run_mafft.sh -CPU 8 --failed_cmds run_mafft.sh.failed &
```

---

## ✂️ Step 7: Trim Alignments with trimal

Remove poorly aligned regions.

```bash
mkdir 04_trimal
ls 03_mafft/ | sed 's/.mafft.faa//' | \
xargs -I {} trimal -in 03_mafft/{}.mafft.faa -out 04_trimal/{}.mafft.trimal.faa -gappyout
```

---

## 🧪 Step 8: Model Selection with IQ-TREE

Identify the best-fit evolutionary model.

```bash
mkdir 05_modelfinder
ls 04_trimal/*.faa | msa_length.py - | awk '$2>=300 {print $1}' | \
sed "s/.mafft.trimal.faa//;s/04_trimal\///" | \
while read file; do
  echo "iqtree -s 04_trimal/${file}.mafft.trimal.faa \
--seqtype AA \
-T 4 \
--prefix 05_modelfinder/${file} \
-m TESTONLY \
--msub nuclear \
--mset mrbayes"
done > run_modelfinder.sh

nohup ParaFly -c run_modelfinder.sh -CPU 8 --failed_cmds run_modelfinder.sh.failed &
```

---

## 🔗 Step 9: Concatenate MSAs

Concatenate alignments for phylogenetic inference.

```bash
mkdir 06_iqtree
ls 04_trimal/*.faa | msa_length.py - | awk '$2>=300 {print $1}' | \
concatenate_msa.py -o 06_iqtree/concatenated.faa -
```

---

## 🌳 Step 10: Maximum Likelihood Tree Construction with IQ-TREE

Build the phylogenetic tree using the best-fit model.

```bash
mkdir 06_iqtree
model=$(grep "Best-fit model: .* chosen according to BIC" 05_modelfinder/*log | \
sed 's/.*Best-fit model://;s/ chosen according to BIC//' | \
awk '{a[$0]+=1}END {for (i in a) print i, a[i]}' | sort -k 2 -nr)

nohup iqtree -s 06_iqtree/concatenated.faa \
--seqtype AA \
-o ${outgroup} \
--prefix 06_iqtree/iqtree.ml \
--mem 150G \
-T AUTO --threads-max 16 \
--ninit 50 \
--ufboot 1000 \
--abayes \
-m WAG+I+G4 &
```

---

## Citation

If you use BUSCOphylo in your research, please cite:

> Yanpeng Chen. BUSCOphylo: A workflow for phylogenetic analysis based on BUSCO orthologs. GitHub repository: [https://github.com/buscophylo](https://github.com/buscophylo)

---

