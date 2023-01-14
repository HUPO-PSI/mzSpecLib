import os

data_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "test_data"))


def datafile(name):
    path = os.path.join(data_path, name)
    if not os.path.exists(path):
        raise Exception(
            ("The path %r could not be located in the test suite's data package." % (path, )) +
            "If you are NOT running the test suite, you should not"
            " be using this function.")
    return path
