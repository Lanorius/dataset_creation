import pandas as pd
import matplotlib.pyplot as plt

frame = pd.read_csv("p1_cleaned_frame.tsv", sep="\t")

prots = frame["BindingDB Target Chain  Sequence"].tolist()
comps = frame["Ligand SMILES"].tolist()

prots = list(set(prots))
comps = list(set(comps))

prot_lens = [len(x) for x in prots]
comp_lens = [len(x) for x in comps]

prot_lens.sort()
comp_lens.sort()

prot_cut = int(len(prot_lens)*0.9)+1
comp_cut = int(len(comp_lens)*0.9)+1

print(len(prot_lens))
print(len(comp_lens))

print(prot_lens[prot_cut])
print(comp_lens[comp_cut])

plt.hist(prot_lens)
plt.show()
