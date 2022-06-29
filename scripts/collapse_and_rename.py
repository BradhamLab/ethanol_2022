import pandas as pd
import numpy as np

import os


def is_null(val):
    if isinstance(val, str) and val.lower() in ["none", "na", ""]:
        return True
    if val is None:
        return True
    if pd.isnull(val):
        return True
    return False


def get_gene_name(x, annos, name_columns):
    for col in name_columns:
        if not is_null(annos.at[x, col]):
            return annos.at[x, col]
    return x


if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        # is (gene x sample)
        counts = pd.read_csv(snakemake.input["counts"], index_col=0)
        annos = pd.read_csv(snakemake.input["annos"], index_col=0, delimiter="\t")
        annos["gene_id"] = annos.apply(lambda x: x.name.split("-RA")[0], axis=1)
        annos.set_index("gene_id", drop=True, inplace=True)
        names = snakemake.params["name_columns"]
        print("-" * 20 + " Getting gene names " + "-" * 20)
        samples = counts.columns.values
        annos["gname"] = annos.apply(
            lambda x: get_gene_name(x.name, annos, names), axis=1
        )
        counts["name"] = counts.apply(
            lambda x: get_gene_name(x.name, annos, names), axis=1
        )
        # counts = counts.set_index('name', drop=True)

        to_drop = []
        print("-" * 20 + " Collapsing duplicates " + "-" * 20)
        for dup in counts.index[counts["name"].duplicated()]:
            if dup not in to_drop:
                name = counts.at[dup, "name"]
                to_collapse = counts.index[counts["name"] == name]
                to_drop += [x for x in to_collapse if x != dup]
                counts.loc[dup, samples] = counts.loc[to_collapse, samples].sum(axis=0)
        print("-" * 20 + " Removing duplicates " + "-" * 20)
        out = counts.loc[~counts.index.isin(to_drop), :].copy()
        out.set_index("name", inplace=True)
        out.rename(columns={x: f"sample-{x}" for x in out.columns}, inplace=True)

        out.to_csv(snakemake.output["counts"])
        annos.to_csv(snakemake.output["annos"])
        # out[snakemake.params['uniprot']][ ,snakemake.params['uniprot']]
