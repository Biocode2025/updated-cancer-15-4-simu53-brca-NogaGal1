import random
import os

RNA_codon_table = {}

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.normpath(os.path.join(script_dir, '..', 'data'))
results_dir = os.path.normpath(os.path.join(script_dir, '..', 'results'))

