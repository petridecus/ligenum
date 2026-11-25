# ligenum

A hackathon project for enumerating small-molecule variants. The goal is to implement the chemistry kernels in Rust and let Python users compose them into their own ligand-preparation workflows.

I’m exploring whether a structure-of-arrays ligand graph can speed up tautomer enumeration for a single molecule. Gypsum-DL is the first performance baseline.

Today, `ligenum.parse_smiles("CC(=O)C")` parses acetone in Rust and returns a `LigandGraph` to Python. The parser supports a small SMILES subset. Next up is keto–enol tautomer enumeration, with protonation, stereoisomer, and conformer enumeration planned.
