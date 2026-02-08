import random
import os

RNA_codon_table = {}

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.normpath(os.path.join(script_dir, '..', 'data'))
results_dir = os.path.normpath(os.path.join(script_dir, '..', 'results'))

#### functions ####

def Mutate_DNA(seq):
    if len(seq) == 0:
        return seq
    nuc_list = ["A", "T", "G", "C"]
    seq = seq.upper()

    rand_index = random.randrange(len(seq))
    rand_nuc = seq[rand_index]
    nuc_list.remove(rand_nuc)
    new_nuc = random.choice(nuc_list)

    return seq[:rand_index] + new_nuc + seq[rand_index + 1:]


def Insert_DNA(seq):
    if len(seq) == 0:
        return seq
    nuc_list = ["A", "T", "G", "C"]
    seq = seq.upper()

    rand_index = random.randrange(len(seq) + 1)
    new_nucs = random.choice(nuc_list) + random.choice(nuc_list)
    return seq[:rand_index] + new_nucs + seq[rand_index:]

def Delete_DNA(seq):
    if len(seq) == 0:
        return seq
    seq = seq.upper()

    k = random.randint(1, 3)
    rand_index = random.randrange(len(seq))
    return seq[:rand_index] + seq[rand_index + k:]

def Comp_seq(a, b):
    count = abs(len(a) - len(b))
    for x, y in zip(a, b):
        if x != y:
            count += 1
    return count

def Read_dict():
    global RNA_codon_table
    with open(os.path.join(data_dir, "codon_AA.txt")) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                RNA_codon_table[parts[0]] = parts[1]

def DNA_RNA_Cod(dna_seq):
    return dna_seq.upper().replace("T", "U")

def RNA_prot(rna_seq):
    protein = ""
    for i in range(0, len(rna_seq) - 2, 3):
        codon = rna_seq[i:i+3]
        if codon not in RNA_codon_table:
            break
        protein += RNA_codon_table[codon]
    return protein


#### main program ####

def simu53_BRCA():

    Read_dict()

    orgnl_seq = ""
    with open(os.path.join(data_dir, "human_p53_coding.txt"), "r") as f:
        for line in f:
            line = line.strip().upper()
            if not line.startswith(">"):
                orgnl_seq += line

    orig_prot = RNA_prot(DNA_RNA_Cod(orgnl_seq))

    user_input = input("Does the Female has a BRCA1,2 mutation? (Y=Yes, N=No) ").upper()

    has_BRCA = (user_input == "Y")

    runs = 1000
    total_generations = 0

    for run in range(runs):

        generations = 0
        mutated_seq = orgnl_seq
        mutation_count = 0

        while True:
            generations += 1

            if random.random() < 0.0001:

                mutations = random.random()

                if mutations < 0.98:
                    mutated_seq = Mutate_DNA(mutated_seq)
                    mutated_prot = RNA_prot(DNA_RNA_Cod(mutated_seq))

                    if Comp_seq(orig_prot, mutated_prot) > 0:
                        mutation_count += 1
                        orig_prot = mutated_prot

                elif mutations < 0.99:
                    mutated_seq = Insert_DNA(mutated_seq)
                    mutation_count += 1

                else:
                    mutated_seq = Delete_DNA(mutated_seq)
                    mutation_count += 1

            if has_BRCA and mutation_count >= 1:
                break

            if not has_BRCA and mutation_count >= 2:
                break

        total_generations += generations

    avg_generations = total_generations / runs

    avg_years = avg_generations / 365

    with open(os.path.join(results_dir, "results.txt"), "w") as f:
        if has_BRCA:
            f.write("For a female that DOES have BRCA1,2 Mutation:\n")
        else:
            f.write("For a female that does NOT have BRCA1,2 Mutation:\n")

        f.write("The mutation that will change the P53 protein will take in average " + str(round(avg_years, 2)) + " years.")


simu53_BRCA()