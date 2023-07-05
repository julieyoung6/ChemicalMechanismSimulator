
def find_reversible_reaction(reactants, products):
    #NOT CURRENTLY USED IN THIS VERSION
    # finds the species that are involved in a reversible reaction

    collapsed_reactions = []
    # reversable_reactions = []
    rr_sources = []

    for i, reactant in enumerate(reactants):
        reaction_str = reactant[0] + products[i][0] #something to check when mult prods
        collapsed_reactions.append(reaction_str)

    for r in collapsed_reactions:
        if collapsed_reactions.count(r[::-1]) > 0:
            collapsed_reactions[collapsed_reactions.index(r[::-1])] = "-."
            # reversable_reactions.append(r)
            rr_sources.append([r[0]])
            rr_sources.append([r[1]])
    return rr_sources


def simple_reaction(source_species, source_amount, reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants):
    reactants_check = reactants.copy()
    reaction_prod = []
    branches_i = []
    sum_rate_constants = 0

    # find indices of each branch
    for branch in range(reactants_check.count(source_species)):
        working_branch = reactants_check.index(source_species)
        branches_i.append(working_branch)
        sum_rate_constants += reaction_rate_constants[working_branch]
        reactants_check[reactants_check.index(source_species)] = "-"
    
    # solve each branch
    for b in branches_i:
        for p in range(len(products[b])):
            p_amount = (reaction_rate_constants[b]/sum_rate_constants) * source_amount * product_coeffs[b][p]
            reaction_prod.append([products[b][p], p_amount])
    
    return reaction_prod, len(branches_i)

def sim(species, reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants, steps):
    # create empty lists with same len as species
    final_species_amount = [1] * len(species)  # cumulative amount
    current_species_amount = [0] * len(species)

    # find source species by setting all product species to 0 in final list
    for i in range(len(products)):
        for prod in products[i]:
            final_species_amount[species.index(prod)] = 0

    # find the source species
    source_species_index = final_species_amount.index(1)
    current_species_amount[source_species_index] = 1
    source = [species[source_species_index]]

    for r in range(steps):
        next_source = []
        for s in source:
            reaction_prod, num_branches = simple_reaction([s], current_species_amount[species.index(s)], reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants)
            current_species_amount[species.index(s)] = 0
            for p in reaction_prod:
                final_species_amount[species.index(p[0][0])] += p[1]
                current_species_amount[species.index(p[0][0])] = p[1]
                next_source.append(p[0])
    
        # eliminate duplicate sources
        next_source = list(set(next_source))
        source = next_source
           
    return species, final_species_amount


print(sim(["A", "B", "C","D", "E"],[["A"],["B"],["B"],["C"],["C"]],[[1],[1],[1],[1],[1]], [["B"],["C"],["D"],["B"],["E"]],[[1],[1],[1],[1],[1]], [1,1,1,1,1,], 5))
# print(find_reversible_reaction([["A"],["B"],["B"],["C"],["C"]], [["B"],["C"],["D"],["B"],["E"]]))