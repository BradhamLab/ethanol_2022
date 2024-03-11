import os


configfile: "config.yaml"


comparisons = {
    f"{values['reference']}_vs_{key}": dict(
        column=values["column"],
        reference=values["reference"],
        experimental=values["experimental"],
    )
    for key, values in config["comparisons"].items()
}

comparisons = config["comparisons"]


rule all:
    input:
        annos=os.path.join(config["workdir"], "collapsed", "annotations.csv"),
        deseq=expand(
            os.path.join(
                config["workdir"],
                "{collapse}",
                "{comparison}",
                "deseq2",
                "normalized_counts.csv",
            ),
            collapse=["collapsed", "uniprot"],
            comparison=comparisons.keys(),
        ),
        # go=expand(
        #     os.path.join(config["workdir"], "uniprot", "gsva", "{comp}_results.csv"),
        #     comp=["_vs_".join(comp) for comp in comparisons],
        # ),
        # results=[
        #     os.path.join(
        #         config["workdir"],
        #         "uniprot",
        #         "plots",
        #         "gsva",
        #         "{}_top_enrichment.png".format("_vs_".join(comp)),
        #     )
        #     for comp in comparisons
        # ],
        # top=os.path.join(config["workdir"], "uniprot", "gsva", "top_pathways.csv"),


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
        batch=config["params"]["deseq2"]["batch"],
        comparisons=comparisons,
        outdir=os.path.join(config["workdir"], "{collapse}"),
        plotdir=os.path.join(config["workdir"], "{collapse}"),
    output:
        res=expand(
            os.path.join(
                config["workdir"],
                "{{collapse}}",
                "{comparison}",
                "deseq2",
                "results.csv",
            ),
            comparison=comparisons.keys(),
            allow_missing=True,
        ),
        pca=expand(
            os.path.join(
                config["workdir"],
                "{{collapse}}",
                "{comparison}",
                "plots",
                "deseq2",
                "treatment-sample-pca.png",
            ),
            comparison=comparisons.keys(),
            allow_missing=True,
        ),
        counts=expand(
            os.path.join(
                config["workdir"],
                "{{collapse}}",
                "{comparison}",
                "deseq2",
                "normalized_counts.csv",
            ),
            comparison=comparisons.keys(),
            allow_missing=True,
        ),
    script:
        "scripts/deseq2.R"


# DESC_COLUMN_NAME = {"go": "desc", "kegg": "ko_names"}
# GSVA_INPUT = os.path.join(
#     config["workdir"], "uniprot", "{comparison}", "deseq2", "normalized_counts.csv"
# )
# # # batch correct if necessary
# # if config["params"]["deseq2"]["batch"] is not None:
# #     GSVA_INPUT = os.path.join(
# #         config["workdir"], "uniprot", "combat", "corrected_counts.csv"
# #     )
# #     rule batch_correct:
# #         input:
# #             samples=config["input"]["samples"],
# #             counts=os.path.join(
# #                 config["workdir"], "{collapse}", "deseq2", "normalized_counts.csv"
# #             ),
# #         output:
# #             corrected=os.path.join(
# #                 config["workdir"], "{collapse}", "combat", "corrected_counts.csv"
# #             ),
# #             pca=os.path.join(
# #                 config["workdir"], "{collapse}", "plots", "combat", "pca.png"
# #             ),
# #             dist=os.path.join(
# #                 config["workdir"], "{collapse}", "plots", "combat", "distance.png"
# #             ),
# #         params:
# #             batch="Sample",
# #             group="Group",
# #             cofactors=["HPF", "Treatment", "Sample"],
# #         script:
# #             "scripts/batch_correct.R"
# rule run_gsva:
#     input:
#         corrected=GSVA_INPUT,
#         samples=config["input"]["samples"],
#         annos=os.path.join(config["workdir"], "uniprot", "annotations.csv"),
#         genesets=config["input"]["go"]["genesets"],
#         golookup=config["input"]["go"]["lookup"],
#     params:
#         uniprot="uniprot.hit",
#         comparisons=comparisons,
#         group="Group",
#         compare="Treatment",
#         outdir=os.path.join(config["workdir"], "uniprot", "gsva"),
#         plotdir=os.path.join(config["workdir"], "uniprot", "plots", "gsva"),
#     output:
#         X=expand(
#             os.path.join(
#                 config["workdir"], "uniprot", "gsva", "{comparison}_gsva.csv"
#             ),
#             comparison=comparisons.keys(),
#         ),
#         results=expand(
#             os.path.join(
#                 config["workdir"], "uniprot", "gsva", "{comparison}_results.csv"
#             ),
#             comparison=comparisons.keys(),
#         ),
#         png=expand(
#             os.path.join(
#                 config["workdir"],
#                 "uniprot",
#                 "plots",
#                 "gsva",
#                 "{comparison}_top_enrichment.png",
#             ),
#             comparison=comparisons.keys(),
#         ),
#         top=os.path.join(config["workdir"], "uniprot", "gsva", "top_pathways.csv"),
#     script:
#         "scripts/gsva_analysis.R"
