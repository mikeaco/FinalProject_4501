# FinalProject_4501

This is a command line file. It will take in some required information and return mutated files.

How to use:

clone this git repository
```bash
git clone https://github.com/mikeaco/FinalProject_4501.git

cd FinalProject_4501
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

