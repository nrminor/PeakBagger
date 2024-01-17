# PeakBagger
[![Open Source Files](https://github.com/nrminor/PeakBagger/actions/workflows/open-source.yaml/badge.svg)](https://github.com/nrminor/PeakBagger/actions/workflows/open-source.yaml)

### Setup Instructions
This project uses [Poetry](https://python-poetry.org) to manage dependencies and make its software environment reproducible. To make sure it runs on your machine, first install poetry, clone this repo, and then in the repo directory, run `poetry install`. This will make sure all the packages with all the right versions are on hand for the module to run. Once dependencies are installed, you're ready to run the module in one of two ways: with a Poetry shell or through `poetry run`.

Poetry shell looks like this:
```
poetry shell # this will activate a shell within the project's virtual environment
python3 peakbagger/main.py -d /path/to/results
```

Without activating the environment, you may also simply use `poetry run`:
```
poetry run python3 peakbagger/main.py -d /path/to/results
```
