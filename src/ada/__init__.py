import pkgutil
import importlib
import sys
from .registry import registry

# The sub-directories to scan
_SUB_PACKAGES = ['data', 'io', 'processing']
__all__ = []

def _automate_imports():
    """
    Crawls modules and detected @registry decorators to populate __all__.
    """
    curr  = sys.modules[__name__]
    tools = []
    
    for sub_pkg in _SUB_PACKAGES:
        try:
            pkg_path = f"{__name__}.{sub_pkg}"
            sub_pkg = importlib.import_module(pkg_path)
            
            # Walk through all .py files in the sub-folder and execute decorators
            for _, mod_nm, _ in pkgutil.walk_packages(sub_pkg.__path__, sub_pkg.__name__ + "."):
                importlib.import_module(mod_nm)

                
            for nm, obj in registry.objects.items():
                # safety check: don't overwrite if already imported
                setattr(curr, nm, obj)
                tools.append(nm)     
                      
        except Exception as e:
            print(f"Warning: Registration failed for {sub_pkg}: {e}")

    return tools

__all__ = _automate_imports()