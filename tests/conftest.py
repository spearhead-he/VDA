def pytest_collection_modifyitems(items):
    """Runs the fast offline tests before the end-to-end test, which downloads data"""
    items.sort(key=lambda item: item.path.name == "test_end_to_end.py")
