# Development

- [Installation](#installation)
- [Tests](#tests)
- [Baseline image](#baseline-image)
- [Continuous integration](#continuous-integration)
- [Publishing to PyPI](#publishing-to-pypi)

## Installation

From the base directory of the repository, install the package in editable mode with the dependencies of the notebook and of the tests:

```bash
pip install -e ".[notebook,test]"
```

## Tests

The tests are in `tests/`, one file per part of the tool:

| File | Tests |
|---|---|
| `test_onsets.py` | Onset determination: Sigma threshold and Poisson-CUSUM, and the number of channels with onsets |
| `test_events.py` | Events input (files and datetime ranges, and their errors) and background windows |
| `test_data.py` | Grouping of the energy channels, the checks before the download, saved data, the folder of the downloaded data, and when the energy ranges are read again |
| `test_observers.py` | Observers (spacecraft): data loading, viewings, default grouped channels, and the download of the SPICE kernels |
| `test_onset_selection.py` | Onset selection |
| `test_vda_fit.py` | VDA fit and plot: release times at the Sun and at 1 AU, errors, warnings, the printed results and the order of the steps |
| `test_views.py` | Plots of the background window, of the detected onsets (also with a zero background) and of the VDA fit, and their files |
| `test_parameters.py` | Parameters |
| `test_notebook.py` | Notebook interface: the widgets of the parameters form, the background window and onset selection, and the steps of the analysis |
| `test_multi_spacecraft.py` | Multi-spacecraft notebook: the parameters form, the analysis of each spacecraft and the combined VDA plot |
| `test_notebooks_download.py` | Download of the notebooks and their example files (`spearhead-vda-notebooks`), without the network |
| `test_custom_points.py` | Custom points: their VDA fit, energies and particles, the points file, the plot and the notebook |
| `test_end_to_end.py` | The default event, from the download of the data to the VDA plot, compared with its baseline image |

All the tests except the end-to-end test run offline, with synthetic data (the helpers are in `tests/helpers.py`). The tests that download data are marked `network` and run after the offline tests. To run all the tests, from the base directory of the repository:

```bash
pytest tests/
```

Only the offline tests (e.g. without an internet connection):

```bash
pytest -m "not network" tests/
```

With the comparison of the VDA plot with its baseline image (matplotlib 3.11 or newer; with older versions the comparison is skipped):

```bash
pytest -rP --mpl --mpl-baseline-path=baseline --mpl-baseline-relative --mpl-generate-summary=html tests/
```

The notebooks run headless with nbmake (they download data), all of them or one:

```bash
pytest --nbmake *.ipynb
```

```bash
pytest --nbmake vda_multi_event.ipynb
```

## Baseline image

The baseline image of the end-to-end test is `tests/baseline/test_vda_default.png`. When a change of the VDA plot is intended, update it (with matplotlib 3.11 or newer):

```bash
pytest --mpl-generate-path=tests/baseline tests/test_end_to_end.py
```

## Continuous integration

GitHub Actions (`.github/workflows/pytest.yml`) runs on every push and pull request, with Python 3.11 to 3.14: flake8, all the tests with the image comparison and coverage, and the notebook with nbmake. The coverage is uploaded to Codecov.

The workflow `.github/workflows/publish.yml` builds the package (sdist and wheel) and checks it with `twine check` on every push and pull request.

## Publishing to PyPI

The package `spearhead-vda` is published to [PyPI](https://pypi.org/project/spearhead-vda/) by `.github/workflows/publish.yml` when a GitHub release is published, with trusted publishing: the PyPI project trusts this workflow and its GitHub environment `pypi`, so no tokens are needed. Publishing waits for the approval of the environment's reviewers. The version is `spearhead.vda.__version__`, which must be new on PyPI.

Running the workflow manually (Actions > publish > Run workflow) publishes the package to [TestPyPI](https://test.pypi.org/project/spearhead-vda/) (environment `testpypi`), e.g. to check a release before publishing it. A version can be uploaded only once to each index. To install from TestPyPI, download the package without its dependencies (TestPyPI has broken test versions of some of them), then install it with its dependencies from PyPI, e.g. for version 1.0.0:

```bash
pip download --no-deps --index-url https://test.pypi.org/simple/ spearhead-vda==1.0.0 -d dist_test
```

```bash
pip install "dist_test/spearhead_vda-1.0.0-py3-none-any.whl[notebook]"
```

To build and check the package locally:

```bash
pip install build twine
python -m build
twine check --strict dist/*
```
