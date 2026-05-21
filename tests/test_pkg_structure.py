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
    assert hasattr(ada, "load_matlab_file"), "load_matlab_file failed to register at top level"

def test_registry_contents():
    """
    Verify the ObjectRegistry instance actually collected the objects.
    """
    ada.get_tools()
    assert "RecordingSession" in registry.objects
    assert "load_matlab_file" in registry.objects

def test_module_locations():
    """
    Verify that the registered objects still know where they came from
    (Ensures the importlib logic didn't break their metadata).
    """
    from ada.data.recording import RecordingSession
    # The version at the top level should be the EXACT same object as the one in the sub-module
    assert ada.RecordingSession is RecordingSession