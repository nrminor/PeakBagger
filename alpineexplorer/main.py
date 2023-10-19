#!/usr/bin/env python3

"""
ALPINE Explorer is a Python module specifying functions and classes for summarizing
results from one or many [ALPINE](https://github.com/nrminor/ALPINE) runs. This functions
can be run via the script itself, like so:

```
poetry shell
python3 alpineexplorer/main.py ~/path/to/results
```

They can also be accessed in a Jupyter or Quarto notebook with the following:
```
import alpineexplorer
```
"""

import os
import sys
import glob
import argparse
from typing import Optional, Iterable
from dataclasses import dataclass
from result import Result, Ok, Err
import polars as pl


@dataclass
class SearchBranch:
    """
    Dataclass `SearchBranch` contains the necessary file paths for each geography
    within the broader search tree.
    """

    parent_dir: str
    geography: str
    double: Optional[str]
    anachron: Optional[str]
    highdist: Optional[str]
    early_stats: Optional[str]
    late_stats: Optional[str]


@dataclass
class StarterPaths:
    """
    Dataclass `StarterPaths` contains iterables of geography names
    and the corresponding results paths. These iterables are used to
    construct a multi-level search tree of directories containing ALPINE
    results.
    """

    geo: Iterable[str]
    path: Iterable[str]


def parse_command_line_args() -> Result[argparse.Namespace, str]:
    """
        Parses command line arguments, defaulting to the current working directory.

    Args:
        `None`

    Returns:
        `Result[argparse.Namespace, str]`: Returns `Ok(argparse.Namespace)` if args could
        be parsed, else returns `Err(str)`.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results_dir",
        "-d",
        type=str,
        required=False,
        default=".",
        help="The directory to search within.",
    )
    args = parser.parse_args()

    return Ok(args)


def _clean_string(subdir_name: str) -> str:
    """
        Helper function that takes an ALPINE result subdirectory, removes 
        dataset-specific prefixes, and replaces underscores with spaces to
        get a clean geography name.

    Parameters:
        - `subdir_name (str)`: The subdirectory to clean.

    Returns:
        - `str`: The cleaned subdirectory.
    """

    return (
        subdir_name.replace("LocalDataset_", "")
        .replace("GISAID_", "")
        .replace("GenBank_", "")
        .replace("_", " ")
    )


def construct_file_paths(result_root: str) -> Result[StarterPaths, str]:
    """
    The function `construct_file_paths()` looks in the results directory
    provided by the user and parses out the results directories for each
    geography. It also parses the names the geographies, and uses them
    as keys in a dictionary to be iterated through downstream.
    """

    paths = [
        dir
        for dir in os.listdir(result_root)
        if os.path.isdir(os.path.join(result_root, dir))
    ]
    if len(paths) == 0:
        return Err("No subdirectories found in provided results directory.")

    geos = [_clean_string(path) for path in paths]

    starter_paths = StarterPaths(geo=geos, path=paths)

    return Ok(starter_paths)


def define_search_tree(
    start_paths: StarterPaths,
) -> Result[dict[str, SearchBranch], str]:
    """
        The function `define_search_tree()` iterates through the initial dictionary
        of geographies provided by `construct_file_paths()` and constructs each
        "branch" of the search tree. These branches contain all of the available
        results files to be assessed or None if a particular file was not
        generated.

    Args:
        `path_dict: dict[str, str]`: A dictionary where the keys are each geography
        searched by ALPINE, and the values are the file path corresponding each
        geography's results.

    Returns:
        `Result[dict[str, SearchBranch], str]`: A Result type containing either a
        search tree object (`Ok(dict[str, SearchBranch)`) or an error string
        (`Err(str)`)
    """

    search_tree = {}

    for geo, path in zip(start_paths.geo, start_paths.path):
        double = glob.glob(f"{path}/*double_candidates")
        anachron = glob.glob(f"{path}/*metadata_candidates")
        highdist = glob.glob(f"{path}/*high_distance_clusters")
        early_stats = glob.glob(f"{path}/*early_stats.tsv")
        late_stats = glob.glob(f"{path}/*late_stats.tsv")

        search_tree[geo] = SearchBranch(
            parent_dir=path,
            geography=geo,
            double=double[0] if double else None,
            anachron=anachron[0] if anachron else None,
            highdist=highdist[0] if highdist else None,
            early_stats=early_stats[0] if early_stats else None,
            late_stats=late_stats[0] if late_stats else None,
        )

    return Ok(search_tree)


def get_early_count(path: Optional[str]) -> Optional[int]:
    """
        The function `get_early_count()` reads a file of statistics from early in
        the pipeline meant to describe the number of input sequences. `get_early_count()`
        does so for each geography independently so that a dataframe of statistics
        can be generated in a vectorized manner.

    Args:
        `path: Optional[str]`: Either a string specifying a path to read from or `None`.

    Returns:
        `Optional[int]`: Either an integer specifying the number of input sequences for
        a geography, or `None`.
    """

    if path is None:
        return None

    # Read the TSV file into a DataFrame
    stats_df = pl.read_csv(path, separator="\t")

    # Get the count from the first row of the 'num_seqs' column
    count = stats_df.select(pl.col(["num_seqs"])).to_series().to_list()[0]

    return count


def get_late_count(path: Optional[str]) -> Optional[int]:
    """
        The function `get_late_count` reads a file of statistics from the end of
        the pipeline meant to describe the number of double candidate sequences,
        which is to say sequences from viral lineages that are both highly evolved
        and anachronistic. `get_late_count` does so for each geography independently
        so that a dataframe of statistics can be generated in a vectorized manner.

    Args:
        `path: Optional[str]`: Either a string specifying a path to read from or `None`.

    Returns:
        `Optional[int]`: Either an integer specifying the number of input sequences for
        a geography, or `None`.
    """

    # propagate None if that is the input
    if path is None:
        return None

    # Read the TSV file into a DataFrame
    stats_df = pl.read_csv(path, separator="\t")

    # Get the count from the first row of the 'num_seqs' column
    count = stats_df.select(pl.col(["num_seqs"])).to_series().to_list()[0]

    return count


def summarize_anachrons(path: Optional[str]) -> Optional[int]:
    """
        The function `summarize_anachrons` finds the metadata for anachronistic
        sequences and counts the number of entries, if any, for the provided
        geography and associated filepath.

    Args:
        `path: Optional[str]`: Either a string specifying a path to read from or `None`.

    Returns:
        `Optional[int]`: Either an integer specifying the number of input sequences for
        a geography, or `None`.
    """

    # propagate None if that is the input
    if path is None:
        return None

    # Construct the file path for the TSV file
    file_path = f"{path}/anachronistic_metadata_only_candidates.tsv"

    # Read the TSV file into a DataFrame
    anachron_df = pl.read_csv(file_path, separator="\t")

    # Return the number of rows in the DataFrame
    return anachron_df.shape[0]


def summarize_highdist(path: Optional[str]) -> Optional[int]:
    """
        The function `summarize_highdist` finds the metadata for high distance
        sequences and counts the number of entries, if any, for the provided
        geography and associated filepath.

    Args:
        `path: Optional[str]`: Either a string specifying a path to read from or `None`.

    Returns:
        `Optional[int]`: Either an integer specifying the number of input sequences for
        a geography, or `None`.
    """

    # propagate None if that is the input
    if path is None:
        return None

    # Construct the file path for the TSV file
    file_path = f"{path}/high_distance_candidates.tsv"

    # Read the TSV file into a DataFrame
    anachron_df = pl.read_csv(file_path, separator="\t")

    # Return the number of rows in the DataFrame
    return anachron_df.shape[0]


def main() -> None:
    """
        Main daisy-chains the above functions and controls the flow of
        data through them if this script is run on its own in the command
        line.

    Args:
        None

    Returns:
        None
    """

    # collect a result directory from the command line, defaulting to the
    # current working directory
    args = parse_command_line_args().expect("Unable to access filesystem.")

    # make sure the provided directory exists, and if not, exit the program
    assert os.path.exists(args.result_dir) and os.path.isdir(
        args.result_dir
    ), "Provided file path does not exist or is not a directory."

    # construct an initial dictionary of the geographies and top-level
    # per-geography directories to be searched while handling any errors.
    message_result = construct_file_paths(args.result_dir)
    match message_result:
        case Ok(result):
            search_tree = result
        case Err(message):
            sys.exit(
                f"Error originated while constructing file paths:\n\
                     {message}"
            )

    print(f"{search_tree}")


if __name__ == "__main__":
    main()
