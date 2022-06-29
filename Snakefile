import os


configfile: "config.yaml"


comparisons = list(config["comparisons"].values())[0]
compare_strings = ["_vs_".join(x) for x in comparisons]


rule all:
    input:
        annos=os.path.join(config["workdir"], "collapsed", "annotations.csv"),
        deseq=expand(
            os.path.join(
                config["workdir"], "{collapse}", "deseq2", "normalized_counts.csv"
            ),
            collapse=["collapsed", "uniprot"],
        ),
        batch=expand(
            os.path.join(
                config["workdir"], "{collapse}", "combat", "corrected_counts.csv"
            ),
            collapse=["collapsed", "uniprot"],
        ),
        go=expand(
            os.path.join(config["workdir"], "uniprot", "gsva", "{comp}_results.csv"),
            comp=["_vs_".join(comp) for comp in comparisons],
        ),
        results=[
            os.path.join(
                config["workdir"],
                "uniprot",
                "plots",
                "gsva",
                "{}_top_enrichment.png".format("_vs_".join(comp)),
            )
            for comp in comparisons
        ],
        top=os.path.join(config["workdir"], "uniprot", "gsva", "top_pathways.csv"),


rule collapse_to_anno:
    input:
        counts=config["input"]["counts"],
        annos=config["input"]["annotations"],
    output:
        counts=os.path.join(config["workdir"], "collapsed", "counts.csv"),
        annos=os.path.join(config["workdir"], "collapsed", "annotations.csv"),
    params:
        name_columns=config["params"]["annotations"]["name_columns"],
    script:
        "scripts/collapse_and_rename.py"


rule collapse_to_uniprot:
    input:
        counts=config["input"]["counts"],
        annos=config["input"]["annotations"],
    output:
        counts=os.path.join(config["workdir"], "uniprot", "counts.csv"),
        annos=os.path.join(config["workdir"], "uniprot", "annotations.csv"),
    params:
        name_columns=["uniprot.hit"],
    script:
        "scripts/collapse_and_rename.py"


rule run_deseq2:
    input:
        counts=os.path.join(config["workdir"], "{collapse}", "counts.csv"),
        samples=config["input"]["samples"],
    params:
        control="Control",
        treatment="Treatment",
        batch="Sample",
        comparisons=comparisons,
        col=list(config["comparisons"].keys())[0],
        outdir=os.path.join(config["workdir"], "{collapse}", "deseq2"),
        plotdir=os.path.join(config["workdir"], "{collapse}", "plots", "deseq2"),
    output:
        res=expand(
            os.path.join(
                config["workdir"], "{collapse}", "deseq2", "{comparison}.csv"
            ),
            comparison=compare_strings,
            allow_missing=True,
        ),
        volcano=expand(
            os.path.join(
                config["workdir"], "{collapse}", "plots", "deseq2", "{comparison}.png"
            ),
            comparison=compare_strings,
            allow_missing=True,
        ),
        pca=os.path.join(
            config["workdir"],
            "{collapse}",
            "plots",
            "deseq2",
            "treatment-sample-pca.png",
        ),
        corr=os.path.join(
            config["workdir"],
            "{collapse}",
            "plots",
            "deseq2",
            "sample-diss-matrix.png",
        ),
        counts=os.path.join(
            config["workdir"], "{collapse}", "deseq2", "normalized_counts.csv"
        ),
    script:
        "scripts/deseq2.R"

DESC_COLUMN_NAME = {'go': 'desc', 'kegg': 'ko_names'}

rule batch_correct:
    input:
        samples=config["input"]["samples"],
        counts=os.path.join(
            config["workdir"], "{collapse}", "deseq2", "normalized_counts.csv"
        ),
    output:
        corrected=os.path.join(
            config["workdir"], "{collapse}", "combat", "corrected_counts.csv"
        ),
        pca=os.path.join(config["workdir"], "{collapse}", "plots", "combat", "pca.png"),
        dist=os.path.join(
            config["workdir"], "{collapse}", "plots", "combat", "distance.png"
        ),
    params:
        batch="Sample",
        group="Group",
        cofactors=["HPF", "Treatment", "Sample"],
    script:
        "scripts/batch_correct.R"


rule run_gsva:
    input:
        corrected=os.path.join(
            config["workdir"], "uniprot", "combat", "corrected_counts.csv"
        ),
        samples=config["input"]["samples"],
        annos=os.path.join(config["workdir"], "uniprot", "annotations.csv"),
        genesets=config["input"]['go']["genesets"],
        golookup=config["input"]['go']["lookup"],
    params:
        uniprot="uniprot.hit",
        comparisons=comparisons,
        group="Group",
        compare="Treatment",
        outdir=os.path.join(config["workdir"], "uniprot", "gsva"),
        plotdir=os.path.join(config["workdir"], "uniprot", "plots", "gsva"),
    output:
        X=expand(
            os.path.join(
                config["workdir"], "uniprot", "gsva", "{comparison}_gsva.csv"
            ),
            comparison=compare_strings,
        ),
        results=expand(
            os.path.join(
                config["workdir"], "uniprot", "gsva", "{comparison}_results.csv"
            ),
            comparison=compare_strings,
        ),
        png=expand(
            os.path.join(
                config["workdir"],
                "uniprot",
                "plots",
                "gsva",
                "{comparison}_top_enrichment.png",
            ),
            comparison=compare_strings,
        ),
        top=os.path.join(config["workdir"], "uniprot", "gsva", "top_pathways.csv"),
    script:
        "scripts/gsva_analysis.R"
