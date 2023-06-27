
def simple_reaction(source_species, source_amount, reactants, products, reaction_rate_constants):
    reaction_prod = []
    branches_i = []
    sum_rate_constants = 0

    # find indices of each branch
    for branch in range(reactants.count(source_species)):
        working_branch = reactants.index(source_species)
        branches_i.append(working_branch)
        sum_rate_constants += reaction_rate_constants[working_branch]
        reactants[reactants.index(source_species)] = "-"
    
    # solve each branch
    for b in branches_i:
        p_amount = (reaction_rate_constants[b]/sum_rate_constants) * source_amount
        reaction_prod.append([products[b][0], p_amount])

    return reaction_prod, len(branches_i)

def sim(species, reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants):
    # create empty final list with same len as species
    final_species_amount = [1] * len(species)

    # find source species by setting all product species to 0 in final list
    for i in range(len(products)):
        for prod in products[i]:
            final_species_amount[species.index(prod)] = 0

    source_species_index = final_species_amount.index(1)
    source = [species[source_species_index]]

    r_check = 0

    for r, react in enumerate(reactants):
        if r_check <= r:
            next_source = []
            for s in source:
                reaction_prod, num_branches = simple_reaction([s], final_species_amount[species.index(s)], reactants, products, reaction_rate_constants)
                for p in reaction_prod:
                    final_species_amount[species.index(p[0][0])] += p[1]
                    next_source.append(p[0])
        
            # eliminate duplicate sources
            next_source = list(set(next_source))
            source = next_source
            r_check += num_branches
    
    final_species_amount[source_species_index] = 0
           
    return species, final_species_amount


print(sim(["A", "B", "C"], [["A"],["A"]], [[1],[1]], [["B"],["C"]], [[1],[1]], [1,3]))
print(sim(["A", "B", "C","D"], [["A"],["B"],["A"]], [[1],[1]], [["B"],["D"],["C"]], [[1],[1]], [1,8,3]))

print(sim(["A", "B", "C","D"], [["A"],["B"], ["C"],["A"]], [[1],[1]], [["B"], ["D"],["D"],["C"]], [[1],[1]], [1,3,4,2]))
print(sim(["A", "B", "C","D","E","F"], [["A"],["A"],["B"], ["C"], ["D"],["A"]], [[1],[1]], [["C"],["D"],["F"],["E"],["E"],["B"]], [[1],[1]], [2,1,3,3,3,1]))