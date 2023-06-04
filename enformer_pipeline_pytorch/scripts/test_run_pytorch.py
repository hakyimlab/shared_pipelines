# %% [markdown]
# Testing running pytorch models using the common pipeline functions

# %%
import os, sys, json
import pandas as pd # for manipulating dataframes
import time
from datetime import date
import argparse
import torch


batch_utils_path = "/grand/gpu_hack/imlab/users/saideep/shared_folder/enformer_pipeline/scripts/modules"
sys.path.append(batch_utils_path)

import predictionUtils
import predictUtils_one
import predictUtils_two
import sequencesUtils
import enformer_pytorch
import attention

# %% [markdown]
# We can now start by building a test sequence

# %%
fasta_file = "/grand/TFXcan/imlab/data/hg_sequences/hg38/Homo_sapiens_assembly38.fasta"

test_region = "chr5_2000000_2000100"
seq_length_basenji = 131072
fasta_extractor = fasta_extractor = sequencesUtils.get_fastaExtractor(fasta_file)

test_seq = sequencesUtils.extract_reference_sequence(test_region, fasta_extractor, )

# %%
print(test_seq)

# %% [markdown]
# Great, we should load a test pytorch model now

# %%
torch.device("cpu")
model = enformer_pytorch.Enformer()
print(torch.cuda.is_available())
model.load_state_dict(torch.load("/grand/gpu_hack/imlab/users/saideep/test_train_OLD/models/model_0.pt"),map_location='cpu')
model.eval()


