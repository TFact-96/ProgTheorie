### Random Chain Generation<br>
- Can be found in the file: `algorithms/RandomChain.py`
- The program starts with a N x N Grid, contained in the Grid class. First it places a Node object in the grid at (0, 0, 0), with the type as the first char in the protein string. Then for each char in the protein string, it calculates the next coordinate a random move away from the last Node's coordinates (could be 1 in z-direction, -1 in y, etc...). It then checks if the coordinates overlap with another Node, after which it just chooses another move. If all moves are overlapped, the chain is stuck and it just starts all over again. It does this until all Nodes from the protein string are placed in the Grid.

### Pull Moves on generated chains<br>
- Can be found in the file: `algorithms/PullMove.py`
- Following the results of articles such as `A Complete and Effective Move Set for Simplified Protein Folding` by N. Lesh et al, the pull move has been chosen as an effictive move set for folding a chain. The pullmove is an iterative move which transfers a point `i` diagonally to the next point, which is next to `i + 1` or `i - 1`. After this move, the residue of the chain follows in its footsteps. Which part of the chain follows depends on which half the pull move was made.

### Restart Hill Climbing algorithm<br>
- Can be found in the file: `algorithms/RestartHillClimb.py`
- The hillclimbing method creates a randomly folded chain and performs a set amount of pull moves on random nodes. This set of pulls is repeated an arbitrary amount of times. If there are any stability improvements, the next set of pull moves will be performed on this very chain. If the next set does not result in an improvement, a new chain is created and this algorithm repeats.

### (Bruteforce) Simulated Annealing algorithm
- Can be found in the file: `algorithms/SimAnnealing.py`<br>
- This algorithm uses the simulated annealing technique of finding a global minimum (best stability).<br>

#### Algorithm
1. It first creates a random chain.
2. It deepcopies the old (this) chain and the old score (this stability).
3. It does a given amount of pullmoves on random nodes in this chain (amount is queried to the user when the program runs), then calculates new score (stability). It then either accepts that move, or undoes that move by a formula:<br>
![AcceptEquation](https://latex.codecogs.com/gif.latex?accept%20%3D%202%5E%7B%28oldScore%20-%20newScore%29%20/%20temperature%7D)<br>
4. It generates a random value between 0 and 1, if that random shot is above the `accept` value, it reverts to the old state. Otherwise, it goes on with the new state.
5. After this process it lowers the temperature by either of these formulas, chosen at the beginning:<br><br>
Exponential: ![ExpEquation](https://latex.codecogs.com/gif.latex?Temperature%20%3D%20StartTemp%20%5Cast%20c%5E%7Biteration%7D)<br><br>
Linear: ![LinearEquation](https://latex.codecogs.com/gif.latex?Temperature%20%3D%20StartTemp%20-%20iteration%20*%20c)<br><br>
6. The process repeats back to 2, until `iteration` amount of times.

#### Choosing the Temperature and Coefficient
- The `StartTemp` (between 0 and inf) and `c` (between 0 and 1 for exponential, 0 to inf for linear) are arbitrarily chosen. It differs per chain and configuration, but usually it's good to start at a temperature around 2 and a coefficient `c` chosen such that after the `iteration` amount the Temperature will be extremely low (around 0.001). This can easily be solved by plugging in 0.001 as `Temperature`, `iteration` as amount of iterations, and 2 as `StartTemp`, then solving for `c`.<br>
- It's easy to see that the lower the temperature, the higher the chance is that a better newScore will be accepted. Also a worse newScore would be rejected (plug numbers in with low temp and see for yourself). So at the beginning (`Temp > 1`), the program will still accept bad moves (jump out of local minimum), but at the end (`Temp << 1`) only better scores are accepted (reach for global minimum). Which is why we want to choose the coefficient `c` like explained earlier.

#### Bruteforcing
- There is also function to bruteforce the Sim. Annealing. It's contained in the same file, and if `SA` is chosen as the argument in the program, it will ask immediately how many Sim. Annealings should be bruteforced. After the bruteforce only the best result is returned, so more bruteforcing is always better, but at the cost of computing time.

### Self Learning algorithm
- Can be found in the folder `algorithms/SelfLearning`
- This algorithm uses a Q learning agent that finds the step sequence that finds the lowest stability.<br>

#### Algorithm
1. Get all possible directions to move in.
2. Take a step according to the highest Q value. Unless there are no more steps to be taken, in that case save Q and frequency tables and start over with the same chain.<br>
3. Every step calculate reward
4. With that reward the Q-table and frequency table are updated with the folowing formula:
![QValueUpdate](https://latex.codecogs.com/gif.latex?Q%5Bs%2Ca%5D%20%3D%20Q%5Bs%2Ca%5D&plus;%5Calpha%28N_%7Bsa%7D%5Bs%2Ca%5D%29%28r&plus;%5Cgamma%20max_%7Ba%27%7D%20Q%5Bs%2Ca%5D%20-%20Q%5Bs%27%2Ca%27%5D%29)
5. Repeat from step 1.<br>

#### Explanation
Q-table: indicates what kind of reward is possible if that step is taken, the higher the Q value the better the chance of a good end result.<br>
Frequency table: a table that tracks how many times a certain state is visited.<br>
Epsilon: randomness value. Every step there is a 10% (epsilon = 0.1) chance a random step will be taken to keep the program from overfitting<br>
First run: in the very first run a completely random chain is generated.<br>
Saving data: the Q and frequency tables are saved to .csv files.<br>
