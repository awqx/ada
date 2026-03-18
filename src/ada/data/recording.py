import numpy as np

from dataclasses import dataclass
from typing      import List, Union

from ..registry import registry

@registry.register
@dataclass
class RecordingSession: 
    """
    A single neural recording session from Lee et al. (2024)
    """
    # neural data: [stimuli x bins x channels]
    bins: np.ndarray

    # subject ID
    subj: Union[str, List[str]]
    n_subj: int

    # stimuli labels: [stimuli x characteristics]
    stim: np.ndarray

    # other session details
    area: Union[str, List[str]]
    wks_old: int 
    stim_reps: int