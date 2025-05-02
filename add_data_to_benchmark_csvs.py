from rdkit.Chem import Descriptors

for n in tqdm([20, 30], desc="Benchmarking"):
    n_analogues = analogue_search(explorer=explorer, n=n)
    n_analogues = generate_analogue_mols(n_analogues, n=n)
    cse = ChemSpaceExplorerBenchmarkParams(
        n=n, analogues=n_analogues, method_params=ssmmcs_params
    )

    filename = f"{len(query_mols)}_results_for_{cse}.csv"
    # Sanitize the whole filename again just to be safe
    filename = sanitize_filename(filename)
    load_file = os.path.join(DATA_FOLDER, "MCS", filename)
    save_file = os.path.join(DATA_FOLDER, filename)

    df = pd.read_csv(load_file)
    print(f"Loaded {len(df)} {cse.method_params.method_name} results from {load_file}.")

    df.drop(
        columns=[
            "mcs_matched_query",
            "mss_matched_query",
            "failed",
            "mcs_size",
            "mcs_smarts",
            "mss_size",
            "mss_smarts",
            "query_match_mcs",
            "query_match_mss",
        ],
        inplace=True,
    )

    query_sizes = []
    query_masses = []
    analogue_smiles = []
    analogue_sizes = []
    mean_analogue_sizes = []
    analogue_masses = []
    mean_analogue_masses = []
    all_isf_scores = []
    all_distances = []
    for query_id, query_mol in enumerate(query_mols):
        query_analogues = n_analogues[n_analogues["query_spectrum_id"] == query_id]
        if len(query_analogues) <= 1:
            continue

        query_sizes.append(query_mol.GetNumAtoms())
        query_masses.append(Descriptors.MolWt(query_mol))
        analogue_smiles.append(query_analogues["smiles"].to_list())

        analogue_mols = query_analogues["mol"]

        sizes = [mol.GetNumAtoms() for mol in analogue_mols]
        analogue_sizes.append(sizes)
        mean_analogue_sizes.append(sum(sizes) / len(sizes) if len(sizes) > 0 else 0)

        masses = [Descriptors.MolWt(mol) for mol in analogue_mols]
        analogue_masses.append(masses)
        mean_analogue_masses.append(sum(masses) / len(masses) if len(masses) > 0 else 0)

        all_isf_scores.append(query_analogues["isf"].to_list())
        all_distances.append(query_analogues["predicted_distance"].to_list())

    df["query_size"] = query_sizes
    df["query_mass"] = query_masses
    df["analogue_smiles"] = analogue_smiles
    df["analogue_sizes"] = analogue_sizes
    df["mean_analogue_size"] = mean_analogue_sizes
    df["analogue_masses"] = analogue_masses
    df["mean_analogue_mass"] = mean_analogue_masses

    df["isf_scores"] = all_isf_scores
    df["mean_isf"] = df["isf_scores"].apply(
        lambda x: sum(x) / len(x) if len(x) > 0 else 0
    )
    df["distances"] = all_distances
    df["mean_distance"] = df["distances"].apply(
        lambda x: sum(x) / len(x) if len(x) > 0 else 0
    )

    # rearrange colums so that the columns "query_weight_diff", "analogue_weight_diffs"	and "stacked_weights" are at the end
    cols = df.columns.tolist()
    cols = [
        col
        for col in cols
        if col not in ["query_weight_diff", "analogue_weight_diffs", "stacked_weights"]
    ]
    cols += ["query_weight_diff", "analogue_weight_diffs", "stacked_weights"]
    df = df[cols]

    df.to_csv(save_file, index=False)
