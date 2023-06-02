# p2m
Map protein identifiers to metabolites

## Installation

Create and activate conda environment:

```bash
conda create -n p2m -c rdkit -c openbabel python rdkit openbabel pandas
conda activate p2m
```

Clone and setup p2m:

```bash
pip install git+https://github.com/pnnl/p2m
```

## Run

```bash
p2m run [OPTIONS]
```

See --help for more options:

```bash
p2m --help
p2m run --help
p2m download --help
```


