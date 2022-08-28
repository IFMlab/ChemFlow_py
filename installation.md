
# Installation instructions in 3 steps

## Step 1 - Download ChemFlow sources.

The sources code for ChemFlow_py can be downloaded from the [Github Repository](https://github.com/lmonari5/ChemFlow_py.git).    

To clone the repository:

```
git clone https://github.com/lmonari5/ChemFlow_py.git
```
    
Then add the ChemFlow_py scripts to your PATH and PYTHONPATH. To do so modify your ~/.bashr with your favourite text editor:

```
nano ~/.bashrc
```

Add:
```
export CHEMFLOW_PY_HOME="/your_path_to_ChemFlow_py"
export PYTHONPATH="${PYTHONPATH}:/${CHEMFLOW_PY_HOME}/modules"
export PATH="${PATH}:/${CHEMFLOW_PY_HOME}/modules"
```
Source the bashrc to update the PATH in the current terminal:
```
source ~/.bashrc
```
__*Remark*__ : make sure the script you want to run are executables:
```
chmod +x your_path_to_ChemFlow_py/modules/chemflow.py
chmod +x your_path_to_ChemFlow_py/modules/bound_shape.py
```
## Step 2 - Install Anaconda3 and set the environment

If you don't have Anaconda installed, follow the steps of the [Anaconda Installation Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

Then create the ChemFlow_py environment:
```
conda env create -n chemflowpy --file chemflowpy.yml
```

__*Remember to activate the environment before running the code*__

To activate the environment run:
```
conda activate chemflowpy
```

__*If the environment creation fails, tries these fix in order:*__
- set `conda config --set channel_priority false`
- you can find useful resources in the directory ChemFlow_py/templates/help/environment. Try again with:
   > `conda create -n chemflowpy --file chemflowpy_advanced.txt`
- In the directory ChemFlow_py/templates/help/environment you can find
'manual_install_environment.txt'. Just open the file and run every line in the terminal

## Step 3 - Install the docking/rescoring/MD programs supported

This step is quite simple: download the programs you want, make them executable and add them to the PATH.

__*Each program is optional*__

-  __*Vina*__ (free license):
    - Type:
        Docking
    - [Documentation](https://vina.scripps.edu/manual/):
        >https://vina.scripps.edu/manual/
    - [Download](https://vina.scripps.edu/downloads/#):
        >https://vina.scripps.edu/downloads/#
    - Make executable:
        > chmod +x your_path_to_vina_folder/bin/vina
    - Add to ~/.bashrc:
        > export PATH="your_path_to_vina_folder/bin:$PATH"
    - Command to run:
        > vina --help
    - Requirements:
        > None

-  __*QVina2.1*__ (free license):
    - Type:
        Docking
    - [Documentation](https://github.com/QVina/qvina):
        >https://github.com/QVina/qvina
    - [Download](https://github.com/QVina/qvina.git):
        >git clone https://github.com/QVina/qvina.git
    - Make executable:
        > chmod +x your_path_to_qvina_folder/bin/qvina2.1
    - Add to ~/.bashrc:
        > export PATH="your_path_to_qvina_folder/bin:$PATH"
    - Command to run:
        > qvina2.1 --help
    - Requirements:
        > None

-  __*Smina*__ (free license):
    - Type:
        Docking and Rescoring
    - [Documentation](https://sourceforge.net/p/smina/code/ci/master/tree/):
        >https://sourceforge.net/p/smina/code/ci/master/tree/
    - [Download](https://sourceforge.net/projects/smina/):
        >https://sourceforge.net/projects/smina/
        
        It downloads directly the executable, is adviced to put it in a specific directory
    - Make executable:
        > chmod +x your_path_to_smina_folder/smina.static
    - Add to ~/.bashrc:
        > export PATH="your_path_to_smina_folder:$PATH"
    - Command to run:
        > smina.static
    - Requirements:
        > None
        
-  __*Autodock4*__ (free license):
    - Type:
        Docking and Rescoring
    - [Documentation](https://autodock.scripps.edu/wp-content/uploads/sites/56/2022/04/AutoDock4.2.6_UserGuide.pdf):
        >https://autodock.scripps.edu/wp-content/uploads/sites/56/2022/04/AutoDock4.2.6_UserGuide.pdf
    - [Download](https://autodock.scripps.edu/download-autodock4/):
        >https://autodock.scripps.edu/download-autodock4/
    - Make executable:
        > chmod +x your_path_to_autodock/autodock4
    
        > chmod +x your_path_to_autodock/autogrid4
    - Add to ~/.bashrc:
        > export PATH="your_path_to_autodock_folder:$PATH"
    - Command to run:
        > autodock4

        > autogrid4
    - Requirements:
        > None
        
-  __*PLANTS*__ (free academic license):
    - Type:
        Docking and Rescoring
    - [Documentation](http://www.tcd.uni-konstanz.de/research/plants.php):

        Online: 

        >https://www.yumpu.com/en/document/view/27077942/protein-ligand-ant-system-user-manual-for-version-11-1-

        The full PDF is available at the download page

    - [Download](http://www.tcd.uni-konstanz.de/research/plants.php):
        >http://www.tcd.uni-konstanz.de/plants_download/

        It downloads directly the executable, is adviced to put it in a specific directory
    - Make executable:
        > chmod +x your_path_to_plants/PLANTS1.2_64bit
    - Add to ~/.bashrc:
        > export PATH="your_path_to_plants_folder/:$PATH"
    - Command to run:
        > PLANTS1.2_64bit
    - Requirements:
        > None
        
-  __*AMBER*__ (paid license, tested with Amber18):
    - Type:
        MD for Free Energy Rescoring
    - [Documentation](https://ambermd.org/):
        >https://ambermd.org/Manuals.php
    - [Download](https://ambermd.org/GetAmber.php):
        >https://ambermd.org/GetAmber.php
    - Make executable:
        follow installation according to your version
    - Add to ~/.bashrc:
        > source path_to_amber/amber.sh
    - Command to run (according to the installation):
        > pmemd
        > 
        > pmemd.cuda
        > 
        > pmemd.cuda_SPFP
    - Requirements:
        > [Ambertools](https://ambermd.org/AmberTools.php) (free license, already included in the environment)
        
        > [Gaussian_09](https://gaussian.com/glossary/g09/) (paid license) only for RESP charges:
        
         Add to ~/.bashrc:
         >>export g09root="path_to_gaussian"
         >>
         >>export GAUSS_SCRDIR=/tmp
         >>
         >>source $g09root/g09/bsd/g09.profile

        
-  __*GROMACS*__ (free license, tested with Gromacs 2021):
    - Type:
        MD for Free Energy Rescoring
    - [Documentation](https://www.gromacs.org/):
        >https://manual.gromacs.org/current/index.html
    - [Download](https://www.gromacs.org/):
        >https://manual.gromacs.org/current/download.html
    - Make executable:
        follow installation according to your version
    - Add to ~/.bashrc (in case you don't use "sudo make install"):
        > export PATH="your_path_to_gromacs_build/bin/:$PATH"
    - Command to run:
        > gmx
    - Requirements:
        > [Ambertools](https://ambermd.org/AmberTools.php) (free license, already included in the environment)
        
        > [CGenFF program and cgenff_charmm2gmx_py3_nx2.py](https://gaussian.com/glossary/g09/) (free academic license):
        
         >> Dowload and Intsallation:
         >>
         >>      http://kenno.org/pro/cgenff/program.php
         >>
         >> **replace the correspondent directories and files inside ChemFlow_py/templates/config/gromacs/cgenff**
                
