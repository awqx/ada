import numpy as np
import pandas as pd

from dataclasses import dataclass
from typing      import List, Union

from ..registry import registry


@registry.register
@dataclass
class SpikeAverages:
    """Stores averaged neural responses indexed by some grouping."""
    # Dictionary mapping 'family_name' -> np.ndarray [Bins x Channels]
    avgs: dict[str, np.ndarray]
    area: str
    subj: Union[str, List[str]]
    grouping: str
    
    def get_group(self, name: str):
        return self.avgs.get(name)

@registry.register
def avg_by_coherence(session) -> 'SpikeAverages':
    """
    Groups neural 'bins' by the 'coherence' column in session.stim 
    and returns a SpikeAverages object.
    """

    coherences = session.stim[:, 0] 
    unique_coherences = np.unique(coherences)
    
    coherence_map = {}
    
    for i in unique_coherences:
        # 1. Find indices where the family matches
        indices = np.where(coherences == i)[0]
        
        # 2. Extract those samples: [Matches x Bins x Channels]
        family_data = session.bins[indices, :, :]
        
        # 3. Average across the samples (axis 0)
        # Result is [Bins x Channels]
        coherence_map[str(i)] = np.mean(family_data, axis=0)
        
    return SpikeAverages(
        avgs=coherence_map,
        area=session.area,
        subj=session.subj,
        grouping="coherence",
    )

@registry.register
def avg_by_famsamp(session) -> 'SpikeAverages':
    """
    Groups neural 'bins' by the 'famsamp' column in session.stim 
    and returns a SpikeAverages object.
    """

    temp = pd.DataFrame(session.stim)
    temp.columns = ["coherence", "family", "id", "id_blank", "fam", "samp"]
    temp = temp.assign(famsamp=lambda df: df.groupby(["family", "id"]).ngroup())
    famsamps = temp[["famsamp"]].to_numpy().flatten()
    unique_famsamps = np.unique(famsamps)
    
    famsamp_map = {}
    
    for i in unique_famsamps:
        # 1. Find indices where the family matches
        indices = np.where(famsamps == i)[0]
        
        # 2. Extract those samples: [Matches x Bins x Channels]
        family_data = session.bins[indices, :, :]
        
        # 3. Average across the samples (axis 0)
        # Result is [Bins x Channels]
        famsamp_map[str(i)] = np.mean(family_data, axis=0)
        
    return SpikeAverages(
        avgs=famsamp_map,
        area=session.area,
        subj=session.subj,
        grouping="famsamp",
    )

@registry.register
def avg_by_family(session) -> 'SpikeAverages':
    """
    Groups neural 'bins' by the 'family' column in session.stim 
    and returns a SpikeAverages object.
    """
    # Assume 'stim' is [Samples x Descriptors] 
    # and column 0 is the 'family' string or ID
    families = session.stim[:, 1] 
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
        subj=session.subj,
        grouping="family",
    )

@registry.register
def avg_by_id(session) -> 'SpikeAverages':
    """
    Groups neural 'bins' by the 'id' column in session.stim 
    and returns a SpikeAverages object.
    """
    ids = session.stim[:, 2] 
    unique_ids = np.unique(ids)
    
    id_map = {}
    
    for i in unique_ids:
        # 1. Find indices where the family matches
        indices = np.where(unique_ids == i)[0]
        
        # 2. Extract those samples: [Matches x Bins x Channels]
        id_data = session.bins[indices, :, :]
        
        # 3. Average across the samples (axis 0)
        # Result is [Bins x Channels]
        id_map[str(i)] = np.mean(id_data, axis=0)
        
    return SpikeAverages(
        avgs=id_map,
        area=session.area,
        subj=session.subj,
        grouping="id",
    )
