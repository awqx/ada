class ObjectRegistry:
    def __init__(self):
        self.objects = {}

    def register(self, obj):
        self.objects[obj.__name__] = obj
        return obj

# Global instance
registry = ObjectRegistry()