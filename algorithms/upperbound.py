# In 3D, each amino can maximally bond 3 times with others. This one calculates that each amino can only bond once with another amino.
# So kind of naive upperbound calculator

def calc_upperbound(protein):
    c_count = len([amino for amino in protein if amino == "C"])
    h_count = len([amino for amino in protein if amino == "H"])
    
    # no upperbound if no bonds possible
    if c_count == 0 and h_count == 0:
        return 0
        
    # check if there are an even amount of C's
    if (c_count % 2) == 0:
        c_even = True
    else:
        c_even = False
    
    # same with H's
    if (h_count % 2) == 0:
        h_even = True
    else:
        h_even = False
        
    # calc upperbound of stability by combinations of c's and h's
    if c_even and h_even:
        # first bond all C's --> S=-5*(c_count/2), 
        # then all H's --> S=-(h_count/2), Total:
        min_stability = -(5*c_count + h_count) / 2 # no unbound rest aminos
        
    elif c_even and not h_even:
        # first bond all C's --> S=-5*(c_count/2), 
        # then all H's but one --> S=-(h_count - 1)/2, Total:
        min_stability = -(5*c_count + (h_count - 1)) / 2 # always an unbound H
        
    elif h_even and not c_even:
        # first bond all C's but one --> S=-5*(c_count - 1)/2, 
        # then all H's --> S=-(h_count/2). Total:
        min_stability = -(5*(c_count - 1) + h_count) / 2 # always an unbound H or C
        
    elif not h_even and not c_even:
        # first bond all C's but one --> S=-5*(c_count - 1)/2, then all H's but one --> S=-(h_count - 1)/2
        # Then rest C and H makes a bond: S= -1. Total:
        min_stability = - 1 - ((5*(c_count - 1) + (h_count - 1)) / 2) # no unbound rest aminos
        
    return min_stability