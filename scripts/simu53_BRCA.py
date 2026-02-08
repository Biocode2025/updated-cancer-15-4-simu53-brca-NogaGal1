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
