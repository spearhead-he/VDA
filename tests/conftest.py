def pytest_collection_modifyitems(items):
    """Runs the fast offline tests before the network tests, which download data"""
    items.sort(key=lambda item: item.get_closest_marker("network") is not None)
