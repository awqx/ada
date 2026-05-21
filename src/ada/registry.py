import importlib
import pkgutil

class ObjectRegistry:
    def __init__(self):
        self.objects = {}
        self._is_initialized = False

    def register(self, obj):
        self.objects[obj.__name__] = obj
        return obj

    def discover(self, package_name, sub_packages):
        if self._is_initialized:
            return self.objects

        print(f"DEBUG: Starting discovery for {package_name}...")
        for sub in sub_packages:
            full_path = f"{package_name}.{sub}"
            try:
                pkg = importlib.import_module(full_path)
                for _, mod_nm, _ in pkgutil.walk_packages(pkg.__path__, pkg.__name__ + "."):
                    print(f"  > Importing module: {mod_nm}")
                    importlib.import_module(mod_nm)
            except Exception as e:
                print(f"  ! Error in {sub}: {e}")

        self._is_initialized = True
        return self.objects

    def _reset(self):
        self.objects.clear()
        self._is_initialized = False

    def debug_status(self):
        """Prints the current state of the registry."""
        print(f"--- Registry Debug ---")
        print(f"Initialized: {self._is_initialized}")
        print(f"Objects found ({len(self.objects)}): {list(self.objects.keys())}")
        return self.objects

registry = ObjectRegistry()