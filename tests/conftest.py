def pytest_collection_modifyitems(items):
    """Runs the fast offline unit tests before the end-to-end tests that download data"""
    items.sort(key=lambda item: item.path.name != "test_units.py")
