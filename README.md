## Signac paper

Code to reproduce analyses shown in [Stuart et al. 2021, Nature Methods](https://doi.org/10.1038/s41592-021-01282-5)

To run the workflow, first create a new conda environment containing the dependencies:

```
mamba env create -f environment.yaml
```

The entire workflow can be run by executing:

```
snakemake --cores 8
```

For information about the Signac package, see the [Signac repository](https://github.com/timoast/signac)

![](dag.svg)
