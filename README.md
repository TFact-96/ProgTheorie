# Protein Pow(d)er
## The Problem
- Proteins are long stains of amino-acids that control important proceses in the human body. One of the characteristics of proteins is that they fold themselves so they can be transported by cells in the human body. Proteins that have been folded in the wrong way cause some chronic illnesses such as Alzheimers and cancer. The protein folding problem consists of a sequence of amino acids, each labeled as either hydrophobic (H) or polarised (P). The sequence must be placed on a two-dimensional grid without overlapping, so that adjacent amino acids in the sequence remain horizontally or vertically adjacent in the grid. The goal is to minimize the energy, which in the simplest variation corresponds to maximizing the number of adjacent hydrophobic pairs.<br><br>

## Getting started<br>
- Install the requirements with `pip install -r requirements.txt`

- Format for running: `python main.py [protein_string] [optimalization_type]`<br>
- `protein_string`: A combination of C's, H's and P's. The desired protein chain. Example: `PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP`<br>
- `optimalization_type`:<br>
`R` = Just generate a random configuration of the protein string.<br>
`RHC` = Perform a Restart Hill Climb on the protein string.<br>
`SA` = Perform a (Bruteforce) Simulated Annealing on the protein string.<br>
`SL` = Start a self learning agent to find the best configuration for the given protein.<br>
- Further instructions appear when you run the program for each optimalization algorithm, algorithms are explained below.<br>
- At the end the user will be queried to plot a 3D plot of the chain with its bonds, and a stability per iteration plot if its `RHC` or `SA`.

## Algorithms

- #### Random Chain Generation<br>

- #### Pull Moves on generated chains<br>

- #### Restart Hill Climbing algorithm<br>

- #### (Bruteforce) Simulated Annealing algorithm

- #### Self Learning algorithm

These can be found in the folder `algorithms`. The `visualization` folder contains files to graph chains and stability over time.
In the `classes` folder the gridpoint, grid and node classes created.


