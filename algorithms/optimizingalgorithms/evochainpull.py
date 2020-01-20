from algorithms.chaingeneration.chaingenerate import generate_chain


######### Finds best pulled chain
def find_best_pulled_chain(amino, max_iteration, greedy_chain):
    best_chain = []
    
    # Current protein chain
    if greedy_chain:
        current_hilltop = generate_chain(amino, True, 20)
    else:
        current_hilltop = generate_chain(amino, False, 20)
        
    current_hilltop.set_stability_and_bonds()

    # Save as best chain
    best_c = current_hilltop
    best_stab_c = current_hilltop.stability    
    best_chain.append([best_c])

    for iteration in range(max_iteration):
        best_c_found = False
        print(iteration)

        for it in range(100):
            print(it)
            for index in range(1, len(current_hilltop.state) - 1):
                current_hilltop.pull_move(current_hilltop.state[index], current_hilltop.state)
                temp_stability = current_hilltop.get_stability_and_bonds(True)

                if temp_stability < best_stab_c:
                    # if this stability is better; set all bonds and stability in this object
                    current_hilltop.set_stability_and_bonds()
                    print(temp_stability)
                    best_c = current_hilltop
                    best_stab_c = temp_stability
                    best_c_found = True

        if best_c_found:
            current_hilltop = best_c
            current_stability = best_stab_c
            best_c_found = False

        else:
            best_chain.append([best_c, best_stab_c])
            
            if greedy_chain:
                current_hilltop = generate_chain(amino, True, 20)
            else:
                current_hilltop = generate_chain(amino, False, 20)
            
            current_hilltop.set_stability_and_bonds()

            best_c = current_hilltop
            best_stab_c = current_hilltop.stability

def find_best_c(self):
    for list_ in self.best_chain:
        print(list_[1])
