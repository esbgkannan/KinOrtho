## KinOrtho: a combination of full-length and domain-based orthology inference methods

## Requirement

Please ensure the following software is installed:

- `Python (v3.8.3 or later)` [https://www.python.org/downloads/](https://www.python.org/downloads/)
- `Biopython (v1.75 or later)` [https://biopython.org/wiki/Download/](https://biopython.org/wiki/Download/)
- `NumPy (v1.18.5 or later)` [https://numpy.org/install/](https://numpy.org/install/)
- `NCBI BLAST+ (v2.6.0 or later)` [https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/)
- `MCL (v14.137 or later)` [https://micans.org/mcl/](https://micans.org/mcl/)

## Usage

* Identify all orthologous relationships without query sequences
	* **`python kinortho.py <reference_proteomes> [OPTIONS]`**
		* Example: python kinortho.py ./example/reference_proteomes/

* KinOrtho (Full-length): full-length orthology inference with query sequences
	* **`python kinortho.py <reference_proteomes> -f <full_length_query_seqs>`**
		* Example: python kinortho.py ./example/reference_proteomes/ -f ./example/HumanProteinKinase.fasta
		
* KinOrtho (Domain-based): domain-based orthology inference with query sequences
	* **`python kinortho.py <reference_proteomes> -d <domain_based_query_seqs>`**
		* Example: `python kinortho.py ./example/reference_proteomes/ -d ./example/HumanKinaseDomain.fasta`

* KinOrtho (Overlapping): combines the full-length and domain-based orthology inference results
	* **`python kinortho.py <reference_proteomes> -f <full_length_query_seqs> -d <domain_based_query_seqs>`**
		* Example: `python kinortho.py ./example/reference_proteomes/ -f ./example/HumanProteinKinase.fasta -d ./example/HumanKinaseDomain.fasta`

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
