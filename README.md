<img src="/templates/help/logo.png" width=35% height=35%>

# ChemFlow_py 

ChemFlow_py is a software to simplify, collect and standardize the interface of 
virtual High Throughput Screening tools.

It includes:
 - molecular docking 
 - scoring function rescoring 
 - consensus ranking
 - free energy rescoring (MM/PB(GB)SA)

You can build your own custom workflow within a few steps.

It works both locally or on cluster (with slurm or pbs job manager).


## How to get started

You can follow the installation instruction (installation.md) to install the ChemFlow_py and
all the programs that you want to use.

Once you get it, you can test a little benchmark following the tutorial (tutorial.md), for both 
command line and python interpreter workflows.

If you are in rush, you may use the speedy tutorial (only for python interpreter).

We do not provide any of the licensed softwares used by ChemFlow_py.

## Build your own custom workflow (advanced).

ChemFlow_py is based on python modules, therefore it is completely customizable.

You can modify the default configuration templates inside ChemFlow_py/templates/config/ or you can even add
you a custom program to the workflow.

To add a program to run inside ChemFlow_py:
- use the template ChemFlow_py/templates/help/template_docking_module.py
- write the functions and variable defined in the template
- save your module inside ChemFlow_py/modules (i.e.: 'custom_module.py')
- simply run ChemFlow_py calling you module (i.e: 'program = custom_module')
