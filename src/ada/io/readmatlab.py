import os
import re
import pandas as pd
import numpy  as np

from scipy.io import loadmat

from ..data.recording import RecordingSession
from ..registry       import registry

def file_info(f: str) -> dict:
    region = re.search(r"^zz_([a-zA-Z0-9]{2})", f).group(1)
    age = int(re.search(r"Textures2_([0-9]+)", f).group(1))
    n_subj = int(re.search(r"([0-9])subj\.mat$", f).group(1))

    return {
        "region": region, 
        "age": age, 
        "n_subj": n_subj
    }

def clean_stimuli(stim: np.ndarray) -> np.ndarray: 
    stim = pd.DataFrame(stim)
    stim.columns = ["coherence", "family", "id", "coherence_t3", "x5", "id_blank"]
    return (
        stim
        .drop(columns=["coherence_t3", "x5"])
        .assign(repetition=lambda df: df.groupby(["coherence", "family", "id", "id_blank"]).cumcount() + 1)
        .assign(famsamp=lambda df: df.groupby(["family", "id"]).ngroup() + 1)
        .to_numpy()
    )

@registry.register
def load_matlab_file(f: str) -> RecordingSession: 
    f_info   = file_info(os.path.basename(f))
    mat_data = loadmat(f, simplify_cells=True)
    return RecordingSession(
        bins=mat_data['bins'], 
        subj=mat_data['incl_subj'], 
        n_subj=f_info['n_subj'],
        stim=clean_stimuli(mat_data['save_label']),
        area=f_info['region'],
        wks_old=f_info['age'],
        stim_reps=mat_data['session_label'][3],
    )

@registry.register
def load_matlab_dir(dir: str, silent: bool = False) -> dict:
    """ 
    Loads entire directory of MATLAB files into a dictionary where file
    names are keys of RecordingSessions.
    """ 
    recordings = {}
    for f in os.listdir(dir):
        if re.search(r'\.mat$', f) is None:
            continue
        f_info = file_info(f)
        f_key = f_info['region'] + '_' + str(f_info['age']) + 'wk'
        try: 
            recordings[f_key] = load_matlab_file(os.path.join(dir, f))
            if not silent: 
                print(f"Loaded {f} as key {f_key}.")
        except Exception as e:
            print(f"Failed to load {f}: {e}")
    
    return recordings