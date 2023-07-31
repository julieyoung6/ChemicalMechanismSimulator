def find_reaction_index(species, reactants):
    for i, reactant_list in enumerate(reactants):
        if species in reactant_list:
            return i
    
def find_reversible_reaction(reactants, products):
    # finds the species that are involved in a reversible reaction

    collapsed_reactions = []
    reversible_reactions = []

    # create a string for each reaction
    for i, reactant in enumerate(reactants):
        for j in range(len(products[i])):
            reaction_str = reactant[0] + products[i][j]
            collapsed_reactions.append(reaction_str)

    # find reversible reactions
    rr_species = []
    for r in collapsed_reactions:
        if collapsed_reactions.count(r[::-1]) > 0:
            collapsed_reactions[collapsed_reactions.index(r[::-1])] = "-."
            print(r)
            if r[0] not in rr_species and r[1] not in rr_species:
                reversible_reactions.append([[r[0]],[r[1]]])
                rr_species.append(r[0])
                rr_species.append(r[1])
            elif r[0] in rr_species:
                for rr in reversible_reactions:
                    if [r[0]] in rr:
                        rr.append([r[1]])
                        rr_species.append(r[1])
            elif r[1] in rr_species:
                for rr in reversible_reactions:
                    if [r[1]] in rr:
                        rr.append([r[0]])
                        rr_species.append(r[0])

    return reversible_reactions, rr_species

def external_species_extrapolation(starting_species, rr_sources, rr_species, species, reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants, final_species_amount, current_species_amount):
    # extrapolates final yield for external species of a reversible reaction

    cycle_starting_amount = current_species_amount[species.index(starting_species[0])]
    # find external species
    external_species = []
    for i, reactant in enumerate(reactants):
        if reactant in rr_sources and products[i] not in rr_sources:
            external_species.append(products[i])
    ten_step_vals = []
    ten_step_source_out_of_cycle = 0
    twenty_step_vals = []
    twenty_step_source_out_of_cycle = 0

    # run cycle for 20 steps
    source = starting_species
    for i in range(20):
        next_source = []
        validated_sources = [x for x in source if [x] in rr_sources]
       
        final_species_amount, current_species_amount, next_sources = cycle(validated_sources, species, reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants, final_species_amount, current_species_amount)
        for ns in next_sources:
            next_source.append(ns)
        next_source = list(set(next_source))
        source = next_source
 
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
        print(es, twenty_step_vals[i], ten_step_vals[i], twenty_step_source_out_of_cycle, ten_step_source_out_of_cycle)
        final_species_amount[species.index(es[0])] = ((twenty_step_vals[i] - ten_step_vals[i])/(twenty_step_source_out_of_cycle - ten_step_source_out_of_cycle) * (1-twenty_step_source_out_of_cycle) + twenty_step_vals[i])  * cycle_starting_amount 

    return external_species, final_species_amount

def cycle(source_species, species, reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants, final_species_amount, current_species_amount):
    next_sources = []

    # solve each source branch
    for source in source_species:
        branches_i = []
        sum_rate_constants = 0

        # constant multiplier in place of concentration of the second reactant (if there is one)
        x = .5

        # find indices of each branch
        for i, r in enumerate(reactants):
            if source in r:
                branches_i.append(i)
                if len(r) == 1:
                    sum_rate_constants += reaction_rate_constants[i]
                else: # assuming max of 2 reactants
                    sum_rate_constants += (reaction_rate_constants[i] * x)
                    reaction_rate_constants[i] *= x

        # solve each branch
        for b in branches_i:
            for p in range(len(products[b])):
                p_amount = (reaction_rate_constants[b]/sum_rate_constants) * current_species_amount[species.index(source)] * product_coeffs[b][p]

                # assign values to final and current amounts
                final_species_amount[species.index(products[b][p])] += p_amount
                current_species_amount[species.index(products[b][p])] += p_amount
                next_sources.append(products[b][p])

        current_species_amount[species.index(source)] = 0

    return final_species_amount, current_species_amount, next_sources


def simple_reaction(source_species, source_amount, species, reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants, final_species_amount, current_species_amount):
    reaction_prod = []
    next_sources = []
    branches_i = []
    sum_rate_constants = 0

    # constant multiplier in place of concentration of the second reactant (if there is one)
    x = .5

    # find indices of each branch
    for i, r in enumerate(reactants):
        if source_species[0] in r:
            branches_i.append(i)
            if len(r) == 1:
                sum_rate_constants += reaction_rate_constants[i]
            else: # assuming max of 2 reactants
                sum_rate_constants += (reaction_rate_constants[i] * x)
                reaction_rate_constants[i] *= x

    # solve each branch
    for b in branches_i:
        
        for p in range(len(products[b])):
            p_amount = (reaction_rate_constants[b]/sum_rate_constants) * source_amount * product_coeffs[b][p]
            reaction_prod.append([products[b][p], p_amount])

             # assign values to final and current amounts
            final_species_amount[species.index(products[b][p])] += p_amount
            current_species_amount[species.index(products[b][p])] += p_amount
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

    reversible_reactions, rr_species = find_reversible_reaction(reactants, products)

    for r in range(len(reactants) - len(rr_species)):

        if source == []:
            break

        next_source = []
        for s in source:
            # if species is involved in a reversible reaction use extrapolation method
            if s in rr_species:
                # finds the corresponding reversible reaction
                for i, rr in enumerate(reversible_reactions):
                    if [s] in rr:
                        rr_index = i
                external_species, final_species_amount = external_species_extrapolation([s], reversible_reactions[rr_index], rr_species, species, reactants, reactant_coeffs, products, product_coeffs, reaction_rate_constants, final_species_amount, current_species_amount)
                current_species_amount = final_species_amount.copy() 
                for rr_source in reversible_reactions[rr_index]:
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

    # prints each species and final value 
    for i, s in enumerate(species):
        print(s, final_species_amount[i])

    return species, final_species_amount


# print(sim(["A", "B", "C","D", "E", "F"],[["A"],["B"],["B"],["C"],["C"], ["E"]],[[1],[1],[1],[1],[1], [1]], [["B"],["C"],["D"],["B"],["E"], ["F"]],[[1],[1],[1],[1],[1],[1]], [1,1,1,1,1,1]))
# print(sim(["A", "B", "C","D", "E", "F"],[["A"],["B"],["C"],["C"], ["E"],["B"]],[[1],[1],[1],[1],[1], [1]], [["B"],["D"],["B"],["E"], ["F"],["C"]],[[1],[1],[1],[1],[1],[1]], [1,1,20,1,1,20]))


# # Cycle w/ more than 2 species B<-->C<-->E
# print(sim(["A", "B", "C","D", "E", "F"],[["A"],["B"],["B"],["C"],["C"], ["E"], ["E"]],[[1],[1],[1],[1],[1], [1],[1]], [["B"],["C"],["D"],["B"],["E"], ["C"], ["F"]],[[1],[1],[1],[1],[1],[1],[1]], [1,200,1,200,1,1,1]))

# # Disconnected Cycle Test B<-->C and E<-->F
# print(sim(["A", "B", "C","D", "E", "F", "G"],[["A"],["B"],["B"],["C"],["C"], ["E"], ["F"], ["F"]],[[1],[1],[1],[1],[1], [1],[1],[1]], [["B"],["C"],["D"],["B"],["E"], ["F"], ["E"], ["G"]],[[1],[1],[1],[1],[1],[1],[1],[1]], [1,1,1,1,1,1,1,1]))

# Multiple Reactants Step 1
# print(sim(["A", "B", "C", "D"], [["A","B"],["A"]], [[1,1],[1]], [["D"],["C"]], [[1],[1]],[1,1]))


# Test Mini Dag Reaction
# sim(["A", "B", "C","D"], [["A"],["A"],["A"],["B"],["C"]], [[1],[1],[1],[1],[1]], [["B"],["C"],["D"],["C"],["D"]],[[1],[1],[1],[1],[1]],[1,1,1,1,1] )
# create_f0am_file(["A", "B", "C","D"], [["A"],["A"],["A"],["B"],["C"]], [[1],[1],[1],[1],[1]], [["B"],["C"],["D"],["C"],["D"]],[[1],[1],[1],[1],[1]],[1,1,1,1,1], 'mini_test_dag_reaction')

import pandas as pd

# df = pd.read_csv('tree_reactions.csv', header=None)
# df = pd.read_csv('dag_reactions.csv', header=None)
df = pd.read_csv('cyclic_graph_reactions.csv', header=None)
df.columns = ["reactants", "products"]

reactants = []
products = []
rrconstants = []
coeffs = []
species = []

for i, r in enumerate(df["reactants"]):
    if r not in species:
        species.append(r)
    if df["products"][i] not in species:
        species.append(df["products"][i])
    reactants.append([r])
    products.append([df["products"][i]])
    rrconstants.append(1)
    coeffs.append([1])
species.sort()

sim(species, reactants, coeffs, products, coeffs.copy(), rrconstants)
# create_f0am_file(species, reactants, coeffs, products, coeffs.copy(), rrconstants, 'cyclic_reactions')