"""The notebooks of the VDA tool and their example files, for an installation of the package (e.g. from PyPI).

The notebooks are not part of the package: copy_notebooks() downloads them, with the example files, from the GitHub
release of the installed version. It is also the command spearhead-vda-notebooks:

    spearhead-vda-notebooks              # into the current folder
    spearhead-vda-notebooks my_folder    # into my_folder
"""
import argparse
import os
from urllib.error import HTTPError
from urllib.request import urlopen

from . import __version__

REPOSITORY_RAW_URL = "https://raw.githubusercontent.com/spearhead-he/VDA"

# The notebooks and their example files, as paths in the repository
FILES = (
    "vda_multi_event.ipynb",
    "vda_multi_spacecraft.ipynb",
    "vda_custom_points.ipynb",
    "examples/datetime_range_example.csv",
    "examples/datetime_range_bg_example.csv",
    "examples/datetime_range_only_example.csv",
    "examples/reference_times_example.csv",
    "examples/custom_points_example.csv",
)


def files_url(version) -> str:
    """URL of the files of the repository at the release of the version"""
    return f"{REPOSITORY_RAW_URL}/v{version}"


def copy_notebooks(folder=".", overwrite=False, version=__version__) -> dict:
    """Downloads the notebooks and their example files of the version (by default the installed one) into the folder.

    Existing files are kept, unless overwrite is True. Returns {path: "downloaded" or "kept"}.
    """
    url = files_url(version)
    copied = {}
    for name in FILES:
        path = os.path.join(folder, *name.split("/"))
        if os.path.exists(path) and not overwrite:
            copied[path] = "kept"
            continue
        try:
            with urlopen(f"{url}/{name}", timeout=60) as response:
                data = response.read()
        except HTTPError as e:
            if e.code == 404:
                raise ValueError(f"{name} of version {version} not found on GitHub ({url}/{name}): "
                                 f"is v{version} a release of the VDA tool?") from None
            raise
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open(path, "wb") as f:
            f.write(data)
        copied[path] = "downloaded"
    return copied


def main(argv=None) -> None:
    """The command spearhead-vda-notebooks"""
    parser = argparse.ArgumentParser(
        prog="spearhead-vda-notebooks",
        description="Downloads the notebooks of the SPEARHEAD VDA tool and their example files from GitHub, "
                    "for the installed version of the tool.",
    )
    parser.add_argument("folder", nargs="?", default=".", help="folder to download them into (default: the current one)")
    parser.add_argument("--overwrite", action="store_true", help="replace the files that already exist")
    parser.add_argument("--version", default=__version__,
                        help=f"version of the tool whose notebooks are downloaded (default: the installed one, {__version__})")
    args = parser.parse_args(argv)
    try:
        copied = copy_notebooks(args.folder, args.overwrite, args.version)
    except ValueError as e:
        parser.exit(1, f"Error: {e}\n")
    for path, action in copied.items():
        print(f"{'Downloaded' if action == 'downloaded' else 'Kept the existing'} {path}")
    if "kept" in copied.values():
        print("Use --overwrite to replace the existing files.")
    print(f"Open a notebook with e.g.: jupyter-lab {os.path.join(args.folder, FILES[0])}")
