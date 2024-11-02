# DNA-ESA

This repository contains code files for DNA-ESA. The main code files can be found under `src/dna2vec`, while the data is saved under the `data` directory.

## Installation

```
git clone {repo url}
cd {repo}
pip install -e .
```

Furthermore, if you wish to simulate reads (used for training and evaluation) you will have to install [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm). We do this using: 

```bash
apt install art-nextgen-simulation-tools
```


## Training the model

To train the model, simply run:

```
python train.py
```

The training script includes the config object. The config object contains all the hyperparameters for the model.

## Evaluate
Please refer to the `evaluate` sub-directory.


# Next Steps


- [ ] Training optimizations
    - [x] Add gradient clipping
    - [x] Add learning rate scheduler
    - [x] Add gradient accumulation
- [ ] Add more data
    - [x] Add all of the human genome
- [x] Issues

- [x] fix loss to mean loss? shouldn't change with batch size
- [ ] HPC implementation
