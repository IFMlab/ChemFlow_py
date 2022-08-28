# Speedy ChemFlow_py tutorial

## What you are going to do

The speedy tutorial shows how to build a python script to analyze the dataset of [CDK2](http://dude.docking.org/targets/cdk2) from the [DUD-E database](http://dude.docking.org/).

### Step 0 - Activate the anaconda environmet

Activate the conda environment
```
conda activate chemflowpy
```

## Step 1 - Download and prepare the files

Download and extract the [dataset](http://dude.docking.org/targets/cdk2/cdk2.tar.gz)
(untar the file cdk2 archive and ligand/decoys archives contained)
```
tar -xvf cdk2.tar
cd cdk2
tar -xvf actives_final.mol2.tar
tar -xvf decoys_final.mol2.tar
```
### Open the python interpreter
```
ipython
```
To create a file with 5 random active and 25 random decoys, paste:
```
import chemflow as cf
import random
# read the ligand and decoys as input
actives_dict = cf.read_ligands('actives_final.mol2')
decoys_dict = cf.read_ligands('decoys_final.mol2')
ligand_dict = dict(actives_dict, **decoys_dict)
# select 10 actives and 30 decoys randomly
random_act = random.sample(list(actives_dict.keys()), k=5)
random_dec = random.sample(list(decoys_dict.keys()), k=25)
random_ligand = random_act + random_dec
# create the compounds.mol2 file
compunds_text= '\n'.join([ligand_dict[key] for key in random_ligand])
with open('compounds.mol2', 'w') as file:
    file.write(compunds_text)
```

### Step 2 - Run Docking

To compute docking with all the docking program, use:
  
  ```
  import chemflow as cf
  for dock_progam in ('autodock', 'smina', 'plants'):
      docking = cf.Dock('receptor.pdb', 'compounds.mol2', program=dock_progam)
      docking.bounding_shape('crystal_ligand.mol2', shape='both')
      docking.setup(overwrite=False)
      docking.run(computed=True, keep_poses=3)
  ```

ChemFlow_py produce csv files containing the ranking and a 'computed.mol2' file containing the docked molecules, in each directory.

### Step 3 - Scoring function and free energy rescoring 
Rescoring relies on different scoring function to evaluate the protein-ligand interaction.
Use __smina__ to evaluate: dkoes_fast, dkoes_scoring, vina, vinardo
\
Use __plants__ to evaluate: chemplp, plp, plp95
\
Use __autodock__ to evaluate: autodock
\
Use __amber__ or __gromacs__ to evaluate: mmpbsa, mmgbsa 
  ```
  import chemflow as cf
  rescoring = cf.Rescore('receptor.pdb', 'smina/computed_ligands.mol2', program='plants', protocol='plp_rescoring')
  rescoring.bounding_shape('crystal_ligand.mol2', shape='both')
  rescoring.data['scoring'] = 'plp'
  rescoring.run(setup=True)
  ```

__*Remark:*__

Amber and Gromacs are very punctilious about the pdb structure of the protein. You may have to fix the missing atoms.

### Step 4 - Consensus ranking

Consensus approach is based on the combination of the different ranking obtained by docking.

```
import chemflow as cf
protocol_list = ['autodock', 'smina', 'plants']
cf.consensus(protocol_list, methods=['zscore', 'rbv', 'ass'])
```

If you have the list with all the decoys, you can analyse the results of docking, rescoring and consensus:
```
import chemflow as cf
decoy_list = cf.ligand_list('decoys_final.mol2')
protocol_list = ['autodock', 'smina', 'plants', 'plp_rescoring']
cf.analysis(protocol_list, decoy_list, consensus_deep=3, methods=['zscore', 'rbv'])
