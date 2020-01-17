from algorithms.chaingeneration.chaingenerate import generate_chain

def chain_pulling(amino, use_greedy, optimalization_tries, chain_generations, pull_times_per_chain):
    best_stability = 0
    best_chain_list = []

    for i in range(0, chain_generations):
        # generate_chain(amino, Greedymove=True, optimalization_tries bij het greedy genereren, 3D=True/2D=False)
        Chain = generate_chain(amino, use_greedy, optimalization_tries)

        if Chain:
            # set bonds and stability
            Chain.set_stability_and_bonds()
            
            if use_greedy:
                print(f"Generated greedy chain nr {i + 1}: Stability {Chain.stability} (before pulling)")
            else:
                print(f"Generated random chain nr {i + 1}: Stability {Chain.stability} (before pulling)")
            print("Pulling...")
            # pulling this chain
            Chain.random_pull(pull_times_per_chain)


            if Chain.stability < best_stability:
                best_chain = Chain
                best_chain_list.append(Chain)
                best_stability = Chain.stability
                print(f"After pulling this chain: Stability {Chain.stability} (Record)\n")
            else:
                print(f"After pulling this chain: Stability {Chain.stability}\n")
        else:
            # stuck chains will not be counted, so the count goes one lower
            i -= 1

    return best_chain_list