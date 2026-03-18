import pytest
import ada
from ada.registry import registry

def test_namespace_population():
    """
    Verify that the top-level 'ada' namespace contains the registered 
    classes and functions without needing deep imports.
    """
    # Check for Data Classes
    assert hasattr(ada, "RecordingSession"), "RecordingSession failed to register at top level"
    assert hasattr(ada, "load_recording"), "load_recording failed to register at top level"

def test_registry_contents():
    """
    Verify the ObjectRegistry instance actually collected the objects.
    """
    assert "RecordingSession" in registry.objects
    assert "load_recording" in registry.objects
    
def test_all_declaration():
    """
    Ensure __all__ is populated for IDE support and 'from ada import *'.
    """
    assert "RecordingSession" in ada.__all__
    assert "load_recording" in ada.__all__

def test_module_locations():
    """
    Verify that the registered objects still know where they came from
    (Ensures the importlib logic didn't break their metadata).
    """
    from ada.data.recording import RecordingSession
    # The version at the top level should be the EXACT same object as the one in the sub-module
    assert ada.RecordingSession is RecordingSession