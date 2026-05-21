import pytest
import ada
from ada.registry import registry

def test_registry_is_initially_empty():
    """Ensure the conftest fixture is working and we start clean."""
    assert len(registry.objects) == 0

def test_discovery_fills_registry():
    """Verify that calling discover actually finds the methods."""
    # Trigger discovery manually
    tools = registry.discover("ada", ['data', 'io', 'processing'])
    
    assert len(tools) > 0
    assert any(callable(obj) for obj in tools.values())

def test_lazy_attribute_access():
    """Verify that accessing 'ada.tool' triggers the registry."""
    # 'load_recordings' is decorated with @registry.register
    # trigger __getattr__ in ada/__init__.py
    tool = getattr(ada, "RecordingSession", None)
    
    assert tool is not None
    assert registry._is_initialized is True
    assert "RecordingSession" in registry.objects

def test_import_from_ada():
    """Verify that 'from ada import ...' works with our new setup."""
    # This is the 'acid test' for the __getattr__ logic
    try:
        from ada import load_matlab_file
        assert callable(load_matlab_file)
    except ImportError:
        pytest.fail("Failed to import registered tool directly from package")