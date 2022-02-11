import pandas as pd
import matplotlib.pyplot as plt

frame = pd.read_csv("p1_cleaned_frame.tsv", sep="\t")

prots = frame["BindingDB Target Chain  Sequence"].tolist()
comps = frame["Ligand SMILES"].tolist()

prots = list(set(prots))
comps = list(set(comps))

prot_lens = [len(x) for x in prots]
comp_lens = [len(x) for x in comps]

prots.sort()
comps.sort()

print(prot_lens)

plt.hist(prot_lens)
plt.show()
