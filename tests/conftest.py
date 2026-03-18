import pytest
import sys
from unittest.mock import MagicMock

# prevents failure of registry.register
mock_registry = MagicMock()
mock_registry.register = lambda x: x  
sys.modules["ada.registry"] = mock_registry 

import numpy as np
import scipy.io as sio

@pytest.fixture
def mock_mat_file(tmp_path):
    """
    Generates a temporary .mat file that matches the regex patterns
    expected by file_info and load_recording.
    """
    file_name = "zz_V1_Textures2_12wk_data_3subj.mat"
    file_path = tmp_path / file_name
    
    # Create dummy data structures
    data = {
        'bins': np.random.rand(3 * 100 * 20).reshape(3, 100, 20),
        'incl_subj': np.array([1, 2, 3]),
        'save_label': np.array([
            # coherence, family, id, coherence_t3, x5, id_blank
            [0.5, 1, 101, 0.5, 0, 0],
            [0.5, 1, 101, 0.5, 0, 0], 
            [0.8, 2, 102, 0.8, 0, 0]
        ]),
        'session_label': [None, None, None, 50] # stim_reps is at index 3
    }
    
    sio.savemat(file_path, data)
    return str(file_path)