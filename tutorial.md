# ChemFlow_py tutorial

## What you are going to do

ChemFlow_py is a toolkit for virtual screening, thus in this tutorial we will show some useful functions for molecular docking, 
consensus ranking and free energy rescoring. In this tutorial we will show you ChemFlow_py functionalities using a common dataset in virtual screening:
the protein [CDK2](http://dude.docking.org/targets/cdk2) from the [DUD-E database](http://dude.docking.org/).

## Step 0 - Activate the anaconda environmet

ChemFlow_py is designed to run in the anaconda environment provided, to do so type in the terminal:

```
conda activate chemflowpy
```

## Step 1 - Download and prepare the files

For this tutorial, you will need some files contained in this [DUD-E page](http://dude.docking.org/targets/cdk2)

Otherwise, all the files all already available from [this archive](http://dude.docking.org/targets/cdk2/cdk2.tar.gz)
(untar the file cdk2 archive and ligand/decoys archives contained)

The file we are going to use are:

| File                  | Description                                                    |
|-----------------------|----------------------------------------------------------------|
| `actives_final.mol2`  | File containing the known actives compounds in the mol2 format |
| `decoys_final.mol2`   | File containing the decoys compounds in the mol2 format        |
| `receptor.pdb`        | The pdb structure of the protein/receptor                      |
| `crystal_ligand.mol2` | The crystallographic structure of a ligand in the bindig syte  |

ChemFlow_py only takes only one compounds.mol2 file, thus we have to select the decoys and actives for the screening

To select the ligands you can:
- select manually the actives and decoys and insert them into a new file named 'compounds.mol2'
- Select 5 random actives and 25 random decoys using a python interpreter with the following code.

```
ipython # to enter the python interpreter
```

In the python interpreter paste:
```
import chemflow as cf
import random
# read the ligand and decoys as input
actives_dict = cf.read_ligands('actives_final.mol2')
decoys_dict = cf.read_ligands('decoys_final.mol2')
ligand_dict = {actives_dict, decoys_dict}
# select 10 actives and 30 decoys randomly
random_act = random.sample(list(actives_dict.keys()), k=5)
random_dec = random.sample(list(decoys_dict.keys()), k=25)
random_ligand = random_act + random_dec 
# create the compounds.mol2 file
compunds_text= '\n'.join([complete_dict[key] for key in random_ligand])
with open('compounds.mol2', 'w') as file:
    file.write(compunds_text)
```

## Step 2 - Run Docking

The docking programs currently supported are:

| Program    | characteristic                                                  |
|------------|-----------------------------------------------------------------|
 | `autodock` | Lamarckian Genetic Optimization, box grid with explicit spacing |
| `plants`   | Ant Colony Optimization, spherical grid                         |
| `qvina`    | Accelerated fork of Vina, box grid                              |
| `smina`    | Flexible fork of Vina, box grid                                 |
| `vina`     | Iterated Local Search Optimizer, box grid                       |


Most docking programs require the xyz coordinates of the center of the binding site, with 
some parameters to define the docking grid.
You can use ChemFlow tools and a reference ligand to compute these data.  
You may also provide the coordinates manually.

__*Using a reference ligand:*__

- command line:
  ```
    chemflow.py dock -r receptor.pdb -l compounds.mol2 -p smina --protocol smina_test --crystal crystal_ligand.mol2 --postprocess -k 3`
  ```
  - python interpreter/script
  
  ```
  import chemflow as cf
  docking = cf.Dock('receptor.pdb', 'compounds.mol2', program='smina', protocol='smina_test')
  docking.bounding_shape('crystal_ligand.mol2', shape='both')
  docking.setup(overwrite=False)
  docking.run(computed=True, keep_poses=3)
  ```
\
__*Manual coordinates:*__

To see the coordinates of the reference ligand center use (the program uses a default padding of 15 Å):
```
python $(which bound_shape.py) reference_ligand.mol2 --shape both 
```

Substitute 'x y z' with the  coordinates that you want
- command line:
  ```
  chemflow.py dock -r receptor.pdb -l compounds.mol2 -p smina --protocol smina_test --crystal crystal_ligand.mol2 --center x y z --postprocess -k 3`
  ```
  - python interpreter/script
  
  ```
  import chemflow as cf
  docking = cf.Dock('receptor.pdb', 'compounds.mol2', program='smina', protocol='smina_test')
  docking.data['center] = [x, y, z]
  docking.setup(overwrite=False)
  docking.run(computed=True, keep_poses=3)
  ```
  
__*Remark:*__

Each docking program has different settings, parameters and defualt values. To check them type:

```
chemflow.py dock -p program_name --help
```

__*Results:*__

You can try to run the code with different programs. 
ChemFlow_py organizes the file in a directory (the 'protocol' argument) and produce csv files 
containing the ranking and a 'computed.mol2' file containing the docked molecules

## Step 3 - Scoring function and free energy rescoring 

Rescoring relies on different scoring function to evaluate the protein-ligand interaction.
\
The scoring function currently supported are:
- autodock
- chemplp
- dkoes_fast
- dkoes_scoring
- plp
- plp95 
- vina
- vinardo

To do free energy rescoring:
- MMPB-SA
- MMGB-SA

Use __smina__ to evaluate: dkoes_fast, dkoes_scoring, vina, vinardo
\
Use __plants__ to evaluate: chemplp, plp, plp95
\
Use __autodock__ to evaluate: autodock
\
Use __amber__ or __gromacs__ to evaluate: mmpbsa, mmgbsa
\
\
__*The procedure is the same as for docking, we just use docked structures and specify the scoring function:*__

- command line:
  ```
  chemflow.py rescore -r receptor.pdb -l smina/computed_ligands.mol2 -p plants --scoring plp --protocol plp_rescoring --crystal crystal_ligand.mol2
  ```

  - python interpreter/script
  
  ```
  import chemflow as cf
  rescoring = cf.Rescore('receptor.pdb', 'smina/computed_ligands.mol2', program='plants', protocol='plp_rescoring')
  rescoring.bounding_shape('crystal_ligand.mol2', shape='both')
  rescoring.data['scoring'] = 'plp'
  rescoring.setup(overwrite=False)
  rescoring.run()
  ```
  
__*Advice:*__

Have a look at `chemflow.py rescore -p program_name --help` to see the full options

__*Remark:*__

Amber and Gromacs are very punctilious about the pdb structure of the protein. You may have to fix the missing atoms.

## Step 4 - Consensus ranking

Consensus approach is based on the combination of the different ranking obtained by docking.
Thus, to run consensus properly, run different docking program before

The consensus method supported are: 
\
ans, ass, ecr, rbn, rbr, rbv, z_score

Consensus is supported by __*the python interpreter only:*__
```
import chemfloww as cf
protocol_list = ['smina_test', 'plants_test', ...]
cf.consensus(protocol_list, methods=['zscore', 'rbv'])
```

If you have the list with all the decoys, you can analyse the results of docking, rescoring and consensus:
```
import chemflow as cf
decoy_list = cf.ligand_list('decoys_final.mol2')
protocol_list = ['smina_test', 'plants_test', ..., 'plp_rescoring', ...]
analysis(protocol_list, decoy_list, consensus_deep=3, methods=['zscore', 'rbv']):
```
