# Protein Pow(d)er
## The Problem
- Todo<br>

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

- Random Chain Generation<br>

- #### Pull Moves on generated chains<br>

- #### Restart Hill Climbing algorithm<br>

- #### (Bruteforce) Simulated Annealing algorithm

- #### Self Learning algorithm

These can be found in the folder `Algorithms`


