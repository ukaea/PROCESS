def pytest_addoption(parser):
    parser.addoption(
        "--local-repository",
        action="store",
        default=None,
        help="Location of the local repository of reference MFiles.",
    )
