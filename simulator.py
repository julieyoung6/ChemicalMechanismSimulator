
def find_reversible_reaction(reactants, products):
    # finds the species that are involved in a reversible reaction

    collapsed_reactions = []
    rr_sources = []

    # create a string for each reaction
    for i, reactant in enumerate(reactants):
        reaction_str = reactant[0] + products[i][0] #something to check when mult prods
        collapsed_reactions.append(reaction_str)

    # find the reversible reaction
    for r in collapsed_reactions:
        if collapsed_reactions.count(r[::-1]) > 0:
            collapsed_reactions[collapsed_reactions.index(r[::-1])] = "-."
            rr_sources.append([r[0]])
            rr_sources.append([r[1]])

    return rr_sources

def external_species_extrapolation(starting_species, rr_sources, species, reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants, final_species_amount, current_species_amount):
    # extrapolates final yield for external species of a reversible reaction

    #rearrange rr_sources to have the starting species first
    if starting_species != rr_sources[0]:
        rr_sources = [starting_species] + rr_sources[:-1]
        

    # find external species
    external_species = []
    for i, reactant in enumerate(reactants):
        if reactant in rr_sources and products[i] not in rr_sources:
            external_species.append(products[i])

    ten_step_vals = []
    ten_step_source_out_of_cycle = 0
    twenty_step_vals = []
    twenty_step_source_out_of_cycle = 0
    
    for i in range(20):
        final_species_amount, current_species_amount, next_sources = simple_reaction(rr_sources[i%2], current_species_amount[species.index(rr_sources[i%2][0])], species, reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants, final_species_amount, current_species_amount)
        # collect data at 10 steps
        if i == 9:
            for es in external_species:
                ten_step_vals.append(final_species_amount[species.index(es[0])])
                ten_step_source_out_of_cycle += final_species_amount[species.index(es[0])]
    # collect data after 20 steps
    for es in external_species:
        twenty_step_vals.append(final_species_amount[species.index(es[0])])
        twenty_step_source_out_of_cycle += final_species_amount[species.index(es[0])]

    # extrapolate each external species value if 100% of the source left the cycle
    for i, es in enumerate(external_species):
        final_species_amount[species.index(es[0])] = (twenty_step_vals[i] - ten_step_vals[i])/(twenty_step_source_out_of_cycle - ten_step_source_out_of_cycle)


    return external_species, final_species_amount

def simple_reaction(source_species, source_amount, species, reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants, final_species_amount, current_species_amount):
    reactants_check = reactants.copy()
    reaction_prod = []
    next_sources = []
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

             # assign values to final and current amounts
            final_species_amount[species.index(products[b][p])] += p_amount
            current_species_amount[species.index(products[b][p])] = p_amount
            next_sources.append(products[b][p])


    return final_species_amount, current_species_amount, next_sources

def sim(species, reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants):
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

    rr_sources = find_reversible_reaction(reactants, products)
    r_check = 0

    for r in range(len(reactants) - len(rr_sources)):
        next_source = []
        for s in source:
            #checks if species is involved in the reversible reaction
            if [s] in rr_sources:
                external_species, final_species_amount = external_species_extrapolation([s], rr_sources, species, reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants, final_species_amount, current_species_amount)
                current_species_amount = final_species_amount.copy() 
                for rr_source in rr_sources:
                    current_species_amount[species.index(rr_source[0])] = 0
                for es in external_species:
                    next_source.append(es[0]) 
                break
            final_species_amount, current_species_amount, next_sources = simple_reaction([s], current_species_amount[species.index(s)], species, reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants, final_species_amount, current_species_amount)
            for ns in next_sources:
                next_source.append(ns)
            current_species_amount[species.index(s)] = 0
    
        # eliminate duplicate sources
        next_source = list(set(next_source))
        source = next_source

    return species, final_species_amount

print(sim(["A", "B", "C","D", "E", "F"],[["A"],["B"],["B"],["C"],["C"], ["E"]],[[1],[1],[1],[1],[1], [1]], [["B"],["C"],["D"],["B"],["E"], ["F"]],[[1],[1],[1],[1],[1],[1]], [1,1,1,1,1,1]))
print(sim(["A", "B", "C","D", "E", "F"],[["A"],["B"],["C"],["C"], ["E"],["B"]],[[1],[1],[1],[1],[1], [1]], [["B"],["D"],["B"],["E"], ["F"],["C"]],[[1],[1],[1],[1],[1],[1]], [1,1,20,1,1,20]))
