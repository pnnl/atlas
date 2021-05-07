import sys
import pandas as pd
import numpy as np


with open(snakemake.log[0],"w") as log:
    
    sys.stdout= sys.stderr= log

    d= pd.read_csv(snakemake.input[0],index_col=0, squeeze=True, header=None,sep='\t')

    assert type(d) == pd.Series, "expect the input to be a two column file: {}".format(snakemake.input[0])

    old_cluster_ids = list(d.unique())
    if 0 in old_cluster_ids:
        old_cluster_ids.remove(0)

    map_cluster_ids = dict(zip(old_cluster_ids,
                               utils.gen_names_for_range(
                                   len(old_cluster_ids),
                                    prefix="{sample}_{binner}_".format(**snakemake.wildcards)
                                     )
                               ))

    new_d= d.map(map_cluster_ids)
    new_d.dropna(inplace=True)
    if new_d.shape[0]==0:
        print(f"No bins detected with binner {snakemake.wildcards.binner} in sample {snakemake.wildcards.sample}.\n"
                      "This will break the continuationof the pipeline. "
                      "Check what happened. Maybe the the assembly is too small. "
                      "You can either remove the binner (for all samples) from the config.yaml file or the sample from the sample.tsv"
                        )
        raise Exception(f"No bins detected with binner {snakemake.wildcards.binner} in sample {wildcards.sample}.")
    new_d.to_csv(snakemake.output[0],sep='\t')
