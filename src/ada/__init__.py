from .registry import registry
import sys

_SUB_PACKAGES = ['data', 'io', 'processing']

def __getattr__(name):
    """
    Triggered when an attribute (like a registered tool) 
    is accessed but not yet found in the module namespace.
    """
    # Trigger discovery only once
    registry.discover(__name__, _SUB_PACKAGES)
    
    if name in registry.objects:
        return registry.objects[name]
    
    raise AttributeError(f"module {__name__} has no attribute {name}")

def __dir__():
    """
    Ensures that registered tools show up in dir(ada) 
    and tab-completion in IDEs.
    """
    registry.discover(__name__, _SUB_PACKAGES)
    return list(globals().keys()) + list(registry.objects.keys())

def get_tools():
    """
    Public API to get registered objects. 
    """
    # This triggers the import of submodules safely
    tools_dict = registry.discover(__name__, _SUB_PACKAGES)
    
    # Optional: If you want them available as ada.tool_name
    curr = sys.modules[__name__]
    for name, obj in tools_dict.items():
        if not hasattr(curr, name):
            setattr(curr, name, obj)
            
    return list(tools_dict.keys())

if not hasattr(sys, "_ada_importing"):
    sys._ada_importing = True
    try:
        # Only run this if we aren't already in the middle of a nested import
        registry.discover(__name__, _SUB_PACKAGES)
    finally:
        del sys._ada_importing