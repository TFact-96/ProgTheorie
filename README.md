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

## Algorithms

### Random Chain Generation<br>
- To be seen in the file `algorithms/RandomChain.py`
- The program starts with a N x N Grid, contained in the Grid class. First it places a Node object in the grid at (0, 0, 0), with the type as the first char in the protein string. Then for each char in the protein string, it calculates the next coordinate a random move away from the last Node's coordinates (could be 1 in z-direction, -1 in y, etc...). It then checks if the coordinates overlap with another Node, after which it just chooses another move. If all moves are overlapped, the chain is stuck and it just starts all over again. It does this until all Nodes from the protein string are placed in the Grid.

### Pull Moves on generated chains<br>
- TODO

### Restart Hill Climbing algorithm<br>
- TODO

### (Bruteforce) Simulated Annealing algorithm<br>
- To be seen in the file `algorithms/SimAnnealing.py`
- This algorithm uses the simulated annealing technique of finding a global minimum (best stability).<br>
1. It first creates a random chain
2. It does a pullmove on a random node in this chain calculates new score (stability). It then either accepts that move, or undoes that move by a formula:<br>
![Large accept = 2^{(oldScore - newScore) / temperature}]<br>
3. It does a random shot between 0 and 1, if that random shot is above the accept value, it undoes that move.
4. After this process it lowers the temperature by either of these formulas:<br>
Exponential: ![\Large Temperature = StartTemp \ast c^{iteration}]<br>
Linear: ![\Large Temperature = StartTemp - iteration * c]<br>
5. The process repeats back to 2, until `iteration` amount of times.<br>
- The `StartTemp` (between 0 and inf) and `c` (between 0 and 1) are arbitrarily chosen. It differs per chain and configuration, but usually it's good to start at a temperature around 2 and a coefficient `c` chosen such that after the `iteration` amount the Temperature will be extremely low (around 0.001). This can easily be solved by plugging in 0.001 as `Temperature`, `iteration` as amount of iterations, and 2 as `StartTemp`, then solving for `c`<br> 
It's easy to see that the lower the temperature, the higher the chance is that a better newScore will be accepted. Also a worse newScore would be rejected (plug numbers in with low temp and see for yourself). So at the beginning, the program will still accept bad moves (jump out of local minimum), but at the end only better scores are accepted (reach for global minimum). Which is why we want to choose the coefficient `c` like explained earlier.
- There is also function to bruteforce the Sim. Annealing. It's contained in the same file, and if SA is chosen as the argument in the program, it will ask immediately how many Sim. Annealings should be bruteforced. Don't want it? Choose 1 as the value.