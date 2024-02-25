import pathlib
import pandas as pd

hpfs = [15, 18, 21]
gsva_dir = pathlib.Path(
    "/projectnb/bradham/data/RNASeq/processed/ethanol-2020-11-20/uniprot/gsva"
)
enrichR_dir = pathlib.Path(
    "/projectnb/bradham/data/RNASeq/processed/ethanol-2020-11-20/uniprot/enrichr"
)
toxic_terms = pd.read_csv(
    "/projectnb/bradham/analysis/ethanol-2020-11/new_new_toxic_go_terms.txt",
    sep="\t",
    header=None,
    names=["term", "desc"],
)
gsva_dfs = []
enrichR_dfs = []
for hpf in hpfs:
    # gsva search
    gsva_df = (
        pd.read_csv(next(gsva_dir.glob(f"*{hpf}*{hpf}*_results.csv")), index_col=0)
        .assign(
            enr_rank=lambda df_: df_["Ethanol.mean.enrichment"].rank(ascending=False),
            ratio=lambda df_: df_.enr_rank.astype(int).astype(str) + f"/{len(df_)}",
            hpf=hpf,
        )
        .filter(toxic_terms.term, axis=0)
        .sort_values("enr_rank")[
            ["Ethanol.mean.enrichment", "desc", "enr_rank", "ratio", "hpf"]
        ]
    )
    gsva_dfs.append(gsva_df)
    # enrichr search
    enrichr_df = pd.read_csv(
        next(enrichR_dir.glob(f"*{hpf}*{hpf}*_significant.csv")), index_col=0
    ).filter(toxic_terms.term, axis=0)
    enrichR_dfs.append(enrichr_df)
pd.concat(gsva_dfs).to_csv("toxic_gsva_enrichment.csv")
