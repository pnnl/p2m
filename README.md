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
git clone https://github.com/pnnl/p2m
cd p2m
python setup.py install
```

## Download required local data

Modify DOWNLOAD_PATH to the folder you want to store mapping and ChEBI files:

```bash
p2m download --download_path [DOWNLOAD_PATH]
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


