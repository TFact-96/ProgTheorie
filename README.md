# Protein Pow(d)er

## Getting started<br>
- Do `pip3 install -r requirements.txt` <br>

- Format: `python main.py [protein_string] [optimalization_type]`<br>
- `protein_string`: A combination of C's, H's and P's. The desired protein chain.<br>
- `optimalization_type`:<br>
R = Generate a random configuration of the protein string.<br>
RHC = Perform a Restart Hill Climb on the protein string.<br>
SA = Perform a (Bruteforce) Simulated Annealing on the protein string.<br>
- The random moves with both the Restart Hill Climb and the Simulated Annealing are Pull Moves on nodes.<br>
- Further instructions appear when you run the program.<br>