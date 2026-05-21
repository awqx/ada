import ada
from ada.registry import registry

# 1. Check if the registry is empty initially (it should be)
print("Initial objects:", registry.objects.keys())

# 2. Force discovery
print("\nTriggering discovery...")
# registry.discover('ada', ['data', 'io', 'processing'])
ada.get_tools()

# 3. Check again
registry.debug_status()

# 4. Check attribute access
try:
    print(f"\nTesting access: {ada.load_recording}")
except AttributeError:
    print("\nFAILED: Tool not found in ada namespace.")