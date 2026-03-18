import pytest
import numpy as np
from unittest.mock import patch
from ada import load_recording, RecordingSession

def test_make_recording():
    """
    Test that we can manually instantiate our registered dataclasses.
    """
    bins = np.zeros((4600, 100, 96))
    stim = np.zeros((4600, 5))
    
    session = RecordingSession(
        bins=bins, 
        subj='BB', 
        n_subj=1,
        stim=stim,
        area='V1',
        wks_old=10,
        stim_reps=4,
    )
    
    assert session.region == "V1"
    assert session.bins.shape == (4600, 100, 96)

@patch("ada.io.loadmatlab.loadmat")
def test_factory_logic(mock_loadmat):
    """
    Test the load_recording factory by mocking the scipy.io.loadmat return value.
    """
    # Mock data resembling your MATLAB structure
    mock_loadmat.return_value = {
        'bins': np.zeros((10, 10, 10)),
        'save_label': np.zeros((10, 5)),
        'session': {'label': 'mock_session'},
        'incl': {'subj': 'MonkeyB'}
    }
    
    # This calls the registered function via the top-level namespace
    session = load_recording("fake_path.mat")
    
    assert isinstance(session, RecordingSession)
    assert session.incl_subj == "MonkeyB"