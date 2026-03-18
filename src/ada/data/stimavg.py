import numpy as np

from dataclasses import dataclass
from typing      import List, Union

from ..registry import registry


@registry.register
@dataclass
class SpikeAverages:
    """Stores averaged neural responses indexed by stimulus family."""
    # Dictionary mapping 'family_name' -> np.ndarray [Bins x Channels]
    avgs: dict[str, np.ndarray]
    area: str
    subj: Union[str, List[str]]
    
    def get_family(self, name: str):
        return self.avgs.get(name)

@registry.register
def avg_by_family(session) -> 'SpikeAverages':
    """
    Groups neural 'bins' by the 'family' column in session.stim 
    and returns a SpikeAverages object.
    """
    # Assume 'stim' is [Samples x Descriptors] 
    # and column 0 is the 'family' string or ID
    families = session.stim[:, 0] 
    unique_families = np.unique(families)
    
    family_map = {}
    
    for fam in unique_families:
        # 1. Find indices where the family matches
        indices = np.where(families == fam)[0]
        
        # 2. Extract those samples: [Matches x Bins x Channels]
        family_data = session.bins[indices, :, :]
        
        # 3. Average across the samples (axis 0)
        # Result is [Bins x Channels]
        family_map[str(fam)] = np.mean(family_data, axis=0)
        
    return SpikeAverages(
        avgs=family_map,
        area=session.area,
        subj=session.subj
    )