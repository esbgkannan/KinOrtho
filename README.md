## KinOrtho: a combination of full-length and domain-based orthology inference methods

## Requirement

Please ensure the following software is installed:

- `Python (v3.7.4 or later)` [https://www.python.org/downloads/](https://www.python.org/downloads/)
	- Instruction: [https://devguide.python.org/](https://devguide.python.org/)
- `NCBI BLAST+ (v2.9.0 or later)` [https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/](https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/)
	- Instruction: [https://www.ncbi.nlm.nih.gov/books/NBK279690/](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- `MCL (v14.137 or later)` [https://micans.org/mcl/src/](https://micans.org/mcl/src/)
	- Instruction: [https://micans.org/mcl/](https://micans.org/mcl/)
- `Biopython (v1.75 or later)` [https://biopython.org/wiki/Download](https://biopython.org/wiki/Download)
- `NumPy (v1.18.5 or later)` [https://numpy.org/install/](https://numpy.org/install/)

## Usage

* KinOrtho (Full-length): full-length orthology inference with query sequences
	* **`python kinortho.py <reference_proteomes> -f <full_length_query_seqs>`**
		* Example: python kinortho.py ./example/reference_proteomes/ -f ./example/HumanProteinKinase.fasta
		* Output: [https://github.com/esbgkannan/KinOrtho/tree/main/example/output/full-length/](https://github.com/esbgkannan/KinOrtho/tree/main/example/output/full-length/)
		
* KinOrtho (Domain-based): domain-based orthology inference with query sequences
	* **`python kinortho.py <reference_proteomes> -d <domain_based_query_seqs>`**
		* Example: `python kinortho.py ./example/reference_proteomes/ -d ./example/HumanKinaseDomain.fasta`
		* Output: [https://github.com/esbgkannan/KinOrtho/tree/main/example/output/domain-based/](https://github.com/esbgkannan/KinOrtho/tree/main/example/output/domain-based/)

* KinOrtho (Overlapping): combines the full-length and domain-based orthology inference results
	* **`python kinortho.py <reference_proteomes> -f <full_length_query_seqs> -d <domain_based_query_seqs>`**
		* Example: `python kinortho.py ./example/reference_proteomes/ -f ./example/HumanProteinKinase.fasta -d ./example/HumanKinaseDomain.fasta`
		* Output: [https://github.com/esbgkannan/KinOrtho/tree/main/example/output/overlapping/](https://github.com/esbgkannan/KinOrtho/tree/main/example/output/overlapping/)

## Options

* **-f <full_length_query_seqs>**
	* Full-lengt query sequences (FASTA format)
* **-d <domain_based_query_seqs>**
	* Domain-based query sequences (FASTA format)
* **-o <out_file>**
	* Output file (default: ./results.txt)
* **-E <e_value>**
	* E-value threshld (default: 1e-5)
* **-t <num_threads>**
	* Number of threads (default: 1)
* **-i <inflation_value>**
	* This value handles for affecting cluster granularity. (default: 1.5)
* **-e <min_ev>**
	* Minimal E-value. This value will replace the E-value '0' in BLAST output. (default: 1e-200)

## Output format

* KinOrtho (Full-length):
1. Species_1 - the proteome file name of protein 1
2. Protein_1 - the sequence ID of protein 1
3. Species_2 - the proteome file name of protein 2
4. Protein_2 - the sequence ID of protein 2
5. E-value_Full - the E-value of BLAST result (between protein 1 and protein 2)
6. Weight_Full - the weight between protein 1 and protein 2 in the graph
7. Relationship_Full - orthologous relationship between protein 1 and protein 2 (Ortholog/In-paralog/Co-ortholog)

* KinOrtho (Domain-based):
1. Species_1 - the proteome file name of protein 1
2. Protein_1 - the sequence ID of protein 1
3. Domain_1 - the domain region of protein 1 (dom[INDEX]-[START]-[END])
4. Species_2 - the proteome file name of protein 2
5. Protein_2 - the sequence ID of protein 2
6. Domain_2 - the domain region of protein 2 (dom[INDEX]-[START]-[END])
7. E-value_Domain - the E-value of BLAST result (between domain 1 and domain 2)
8. Weight_Domain - the weight between domain 1 and domain 2 in the graph
9. Relationship_Domain - orthologous relationship between domain 1 and domain 2 (Ortholog/In-paralog/Co-ortholog)

* KinOrtho (Overlapping):
1. Species_1 - the proteome file name of protein 1
2. Protein_1 - the sequence ID of protein 1
3. Domain_1 - the domain region of protein 1 (dom[INDEX]-[START]-[END])
4. Species_2 - the proteome file name of protein 2
5. Protein_2 - the sequence ID of protein 2
6. Domain_2 - the domain region of protein 2 (dom[INDEX]-[START]-[END])
7. E-value_Full - the E-value of BLAST result (between protein 1 and protein 2)
8. E-value_Domain - the E-value of BLAST result (between domain 1 and domain 2)
9. Weight_Full - the weight between protein 1 and protein 2 in the full-length graph
10. Weight_Domain - the weight between domain 1 and domain 2 in the domain-based graph
11. Relationship_Full - orthologous relationship between protein 1 and protein 2 (Ortholog/In-paralog/Co-ortholog)
12. Relationship_Domain - orthologous relationship between domain 1 and domain 2 (Ortholog/In-paralog/Co-ortholog)

## Example datasets

* [Reference proteomes](https://github.com/leon1003/KinOrtho/tree/master/example/reference_proteomes/) `(human, mouse, and fruit fly)`
* [Query sequences (full-length)](https://github.com/leon1003/KinOrtho/blob/master/example/HumanKinaseDomain.fasta) `(545 human protein kinases)`
* [Query sequences (domain-based)](https://github.com/leon1003/KinOrtho/blob/master/example/HumanProteinKinase.fasta) `(558 human kinase domains)`

## Updates

* [v1.0.0](https://github.com/esbgkannan/KinOrtho/tree/main/version/v1.0.0/)
	* The first version of KinOrtho

## Citation

To cite our work, please refer to:

> KinOrtho: a method for mapping human kinase orthologs across the tree of life and illuminating understudied kinases. Liang-Chin Huang, Rahil Taujale, Nathan Gravel, Aarya Venkat, Wayland Yeung, Dominic P Byrne, Patrick A Eyers, and Natarajan Kannan.
