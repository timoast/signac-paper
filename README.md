## Signac paper

Code to reproduce analyses shown in [Stuart et al. 2020, bioRxiv](https://www.biorxiv.org/content/10.1101/2020.11.09.373613v1)

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
