# MUTaTe Pathogenetic SEQuences

The purpose of this project is to simulate pathogenic mutations given a sequence of DNA and its history of mutation rates at each loci in its DNA sequence to determine if it is possible to predict future mutations. While mutations are generally considered to be random, there are factors that influence mutation rates at specific loci. This can be seen when look at historical envolutionary records for dengue virus, where mutation rates at specific loci appear to have much high incidents of nucleotide changes. It is important to note that some mutations, which occur at the nucleotide level, may not propogate to the amino acid level, thus not impacting fitness. If it is possible to accurately predict mutation rates at specific loci, the findings in the project may be useful to pharmaceutical companies who design vaccines. It can be speculated that such companies would like to anticipate where the most likely mutations are to occur for pathogenic targets. 

This is a command line file. It will take in some required information and return mutated files.

How to use:

clone this git repository
```bash
git clone https://github.com/mikeaco/muttseq.git

cd muttseq
#activate virtual environment
./venv/Scripts/activate
#install requirements
pip install -r requirements.txt
```
mutater.py is the script you need to run the program.
* origin file - the file that correlates to the diversity file
* next strain diversity file - file that has the mutation rate weights
* mutation rate - how often this gene mutates per site per year (can be found on nextStrain website)
* mutation time - number of years in the future to mutate
* path - path to the genbank files
* save path - path to where the output files should be saved

This is how the parameters should look for the dengue virus:
```bash
python mutater.py DENV-1_IND_826883_1982.gb nextstrain_dengue_denv1_diversity.tsv 0.000621 1 genbank_files output
```

# Logic writeup
We used "get_numbers.py" to get the proper probabilities for mutation given the start nucelotide. The probability chart
is used to get the probability of what a nuceotide will change into given its current state. This is
using [Markov chain properties for amino acid sequences](https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter10.html). The get_numbers.py file takes hours (close to 20) to complete with one alignment so we ran it once where
we compared the original sequence to the sequences of ten other DENV1 files. We added a "pre_set" variable that has the counter for us 
so you do not have to go through the process. If you would like to run get_numbers.py then feel free you must change the files within the program since we did not
create it for users.

We first take the input data then parse through to get the useful information we need. For example, out of the genbank files we need 
to find the CDS (Coding Sequence) which is the only part that actually has
any inpact on the [sequence](https://en.wikipedia.org/wiki/Coding_region). We also isolated the amino acid sequence just for reference. Lastly, we took the date the sequence was stored as well as the start and stop locations.
We were able to then get the diversity file and use the bases and the frequency of mutation for each base as well. This served as weights
to our choices because some spots are more likely to mutate than others. We also took the mutation rate per year per base which can be found on the [NextStrain website](https://nextstrain.org/dengue/denv1?l=clock). For dengue we noticed it is 6.21x10^-4 subs per site per year. If we multiply this by the total number of sites and years we can get
the total mutations that the sequence will undergo. 

Once we have all of this information, we are able to start with the actual mutation process. 
Our program is set to be able to run on one file or multiple files if needed. 
We then take the sequences and mutate them in the areas chosen by our predictions and mutate it to our predicted nucleotide.
Since a point mutation is defined as a change in the nucleotide we did not allow the mutated points to change 
into themselves (A -> A).

Next, we aligned the sequences together which was the lengthiest process but a necessary step. After we have teh sligned sequence and
mutation algorithm we mutated the area found the codon it was representing and checked the blosom score for the change
and used that as a fitness function of how good the mutation is since the higher teh score the more likely
that it would happen in nature. In a future step, we will be isolating the top scoring sequences and rerunning them only for the 
future mutations to give us even more realistic mutations. 

Lastly, we changed the actual sequence file and the amino acid sequence and saved it to a genbank file so that it could be read again in the software
depending on the reads. 

#### Example:
The following picture shows the similarity between the mutated file and a random dengue_1 sequence from the same time period. 
The first line is the original file we used and mutated according to the mutation rate and set it to 5 years in the future of our 2004 dengue sample. 
The second line is our original 2004 file. 
The third line is our 2009 random dengue_1 sample file we found from NextStrain. 
The yellow color shows how it is similar to the original. The green color shows that our file mutated the sequences correctly in the correct locations. The blue color shows that our file mutated an area that did not align with the random 2009 file. The red color shows that our file did not mutate an area that was mutated in the random 2009 sample. We also had areas where we mutated the file but it mutated to an incorrect nucleotide which can be seen in grey.
![picture of how different the three sequences are](https://github.com/mikeaco/muttseq/blob/master/output/t2_seq_change.jpg?raw=true)
(The full alignment can be found the output file called 3_sequence_alignment)

### Possible TODO:
None of this has any impact on our current application but are nice things we could see in our program.
* Phylogenetic tree creation - creating a phylogenetic tree to see the divergence from the sequences
* Different mutation algorithm - our current set up is not as quick nor as actively changing as we had hoped so we can try different tasks.
* Faster counter writeup - our current file reader and mutator takes for the setup takes hours to complete. This needs to be piped into multiple processes.
