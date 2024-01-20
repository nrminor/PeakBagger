#!/usr/bin/env python3

"""
ALPINE Explorer is a Python module specifying functions and classes for summarizing
results from one or many [ALPINE](https://github.com/nrminor/ALPINE) runs. This functions
can be run via the script itself, like so:

```
poetry install
poetry shell
python3 alpineexplorer/main.py ~/path/to/results
```

They can also be accessed in a Jupyter or Quarto notebook with the following:
```
import alpineexplorer
```

We recommend users access these functions via the Quarto notebook `ALPINEExplorer.qmd`
provided in this repository.
"""

import os
import sys
import glob
import argparse
import tempfile
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


@dataclass
class CompiledMeta:
    """
    This simple dataclass ensures that compiled metadata dataframes are
    stored in the correct order and unpacked properly.
    """

    anachron: pl.DataFrame | pl.LazyFrame
    highdist: pl.DataFrame | pl.LazyFrame
    double: pl.DataFrame | pl.LazyFrame


def parse_command_line_args() -> Result[argparse.Namespace, str]:
    """
        Parses command line arguments, defaulting to the current working directory.

    Args:
        - `None`

    Returns:
        - `Result[argparse.Namespace, str]`: Returns `Ok(argparse.Namespace)` if args could
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

    Args:
        - `result_root: str`: the "results root" directory to search within.
        This is typically a directory named with a YYYYMMDD-formatted date.
    Returns:
        - `Result[StarterPaths, str]`: a result type that is either an instance
        of class `StarterPaths` or an error message string.
    """

    subdirs = [
        dir
        for dir in os.listdir(result_root)
        if os.path.isdir(os.path.join(result_root, dir))
    ]
    if len(subdirs) == 0:
        return Err("No subdirectories found in provided results directory.")

    geos = [_clean_string(path) for path in subdirs]
    paths = [os.path.join(result_root, dir) for dir in subdirs]

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
        - `path_dict: dict[str, str]`: A dictionary where the keys are each geography
        searched by ALPINE, and the values are the file path corresponding each
        geography's results.

    Returns:
        - `Result[dict[str, SearchBranch], str]`: A Result type containing either a
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


def _get_early_count(path: Optional[str]) -> Optional[int]:
    """
        The function `get_early_count()` reads a file of statistics from early in
        the pipeline meant to describe the number of input sequences. `get_early_count()`
        does so for each geography independently so that a dataframe of statistics
        can be generated in a vectorized manner.

    Args:
        - `path: Optional[str]`: Either a string specifying a path to read from or `None`.

    Returns:
        - `Optional[int]`: Either an integer specifying the number of input sequences for
        a geography, or `None`.
    """

    if path is None:
        return None

    # Read the TSV file into a DataFrame
    stats_df = pl.read_csv(path, separator="\t")

    # Get the count from the first row of the 'num_seqs' column
    count = stats_df.select(pl.col(["num_seqs"])).to_series().to_list()[0]

    return count


def _get_late_count(path: Optional[str]) -> Optional[int]:
    """
        The function `get_late_count` reads a file of statistics from the end of
        the pipeline meant to describe the number of double candidate sequences,
        which is to say sequences from viral lineages that are both highly evolved
        and anachronistic. `get_late_count` does so for each geography independently
        so that a dataframe of statistics can be generated in a vectorized manner.

    Args:
        - `path: Optional[str]`: Either a string specifying a path to read from or `None`.

    Returns:
        - `Optional[int]`: Either an integer specifying the number of input sequences for
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


def _summarize_anachrons(path: Optional[str]) -> Optional[int]:
    """
        The function `summarize_anachrons` finds the metadata for anachronistic
        sequences and counts the number of entries, if any, for the provided
        geography and associated filepath.

    Args:
        - `path: Optional[str]`: Either a string specifying a path to read from or `None`.

    Returns:
        - `Optional[int]`: Either an integer specifying the number of input sequences for
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


def _summarize_highdist(path: Optional[str]) -> Optional[int]:
    """
        The function `summarize_highdist` finds the metadata for high distance
        sequences and counts the number of entries, if any, for the provided
        geography and associated filepath.

    Args:
        - `path: Optional[str]`: Either a string specifying a path to read from or `None`.

    Returns:
        - `Optional[int]`: Either an integer specifying the number of input sequences for
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


def stats_pipeline(search_tree: dict[str, SearchBranch]) -> Result[pl.DataFrame, str]:
    """
        Function `stats_pipeline` constructs a Polars dataframe progressively by column,
        where each column displays information about the various categories of results
        from ALPINE, e.g., anachronistic sequences, high-distance sequences, and
        sequences that match both criteria. It uses a series of helper functions to read
        files for each category.

    Args:
        - `search_tree: dict[str, SearchBranch]`: A dictionary containing each geography
        as its keys, and a `SearchBranch` dataclass instance specifying the paths to
        traverse as its values.

    Returns:
        - `Result[pl.DataFrame, str]`: A result type that is either a Polars data frame
        or an error message string.
    """

    # initialize the dataframe
    stats_df = pl.DataFrame({"Geography": search_tree.keys()})

    # bring in input sequence counts
    stats_df = stats_df.with_columns(
        pl.Series(
            "Input Sequence Count",
            [_get_early_count(branch.early_stats) for _, branch in search_tree.items()],
        )
    )

    # bring in double candidate counts
    stats_df = stats_df.with_columns(
        pl.Series(
            "Double Candidate Count",
            [_get_late_count(branch.late_stats) for _, branch in search_tree.items()],
        )
    )

    # compute a percentage of the total inputs that were flagged as double candidates
    stats_df = stats_df.with_columns(
        (
            (pl.col("Double Candidate Count") / pl.col("Input Sequence Count")) * 100
        ).alias("Double Candidate Prevalence (%)")
    )

    # Compute the rate of double candidates in a 1 in X format
    stats_df = stats_df.with_columns(
        pl.concat_str(
            [
                pl.lit("1 in "),
                (1 / (pl.col("Double Candidate Prevalence (%)") / 100))
                .floor()
                .cast(pl.Utf8),
            ]
        ).alias("Double Candidate Rate")
    )

    # count the number of anachronistic lineages
    stats_df = stats_df.with_columns(
        pl.Series(
            "Anachronistic Count",
            [
                _summarize_anachrons(branch.anachron)
                for _, branch in search_tree.items()
            ],
        )
    )

    # compute a percentage of the total inputs that were flagged as anachronistics
    stats_df = stats_df.with_columns(
        ((pl.col("Anachronistic Count") / pl.col("Input Sequence Count")) * 100).alias(
            "Anachronistic Prevalence (%)"
        )
    )

    # Compute the rate of anachronistic candidates in a 1 in X format
    stats_df = stats_df.with_columns(
        pl.concat_str(
            [
                pl.lit("1 in "),
                (1 / (pl.col("Anachronistic Prevalence (%)") / 100))
                .floor()
                .cast(pl.Utf8),
            ]
        ).alias("Anachronistic Rate")
    )

    # get the number of high distance candidates
    stats_df = stats_df.with_columns(
        pl.Series(
            "High Distance Count",
            [_summarize_highdist(branch.highdist) for _, branch in search_tree.items()],
        )
    )

    # compute the percentage of inputs that were high distance
    stats_df = stats_df.with_columns(
        ((pl.col("High Distance Count") / pl.col("Input Sequence Count")) * 100).alias(
            "High Distance Prevalence (%)"
        )
    )

    # Compute the rate of high distance candidates in a 1 in X format
    stats_df = stats_df.with_columns(
        pl.concat_str(
            [
                pl.lit("1 in "),
                (1 / (pl.col("High Distance Prevalence (%)") / 100))
                .floor()
                .cast(pl.Utf8),
            ]
        ).alias("High Distance Rate")
    )

    if stats_df.shape[0] == 0:
        return Err("No results could be compiled.")

    return Ok(stats_df)


def _search_tree_meta(
    search_tree: dict[str, SearchBranch],
    branch_name: str,
    filename: str,
    output_name: str,
) -> None:
    """
        `_search_tree_meta` is a private helper function that iterates through the search
        tree, reads, and converts metadata to produce one Arrow-formatted table.

    Args:
        - `search_tree: dict[str, SearchBranch]`: A dictionary containing each geography
        as its keys, and a `SearchBranch` dataclass instance specifying the paths to
        traverse as its values.
        - `branch_name: str`: The name of the branch to select at each `SearchBranch`.
        - `filename: str`: The name of the file expected in the search directories.
        - `output_name: str`: An output name for the Arrow file.

    Returns:
        None
    """

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_file_path = os.path.join(temp_dir, "tmp.tsv")

        ticker = 0
        for geo, branch in search_tree.items():
            meta_path = f"{getattr(branch, branch_name)}/{filename}"
            if not os.path.isfile(meta_path):
                continue

            ticker += 1
            write_mode = "a"
            header = bool(ticker == 1)
            if header:
                write_mode = "w"

            with open(temp_file_path, write_mode, encoding="utf-8") as tmp:
                pl.read_csv(
                    meta_path,
                    separator="\t",
                ).with_columns(
                    pl.lit(geo).alias("Geography")
                ).write_csv(tmp, separator="\t", has_header=header)

        if ticker > 0:
            pl.scan_csv(
                temp_file_path, separator="\t", infer_schema_length=250
            ).sink_ipc(output_name, compression="zstd")


def compile_metadata(search_tree: dict[str, SearchBranch]) -> Result[CompiledMeta, str]:
    """
        The function `compile_metadata` traverses the search tree and creates a queryable
        metadata in Polars LazyFrame format. It also saves these databases as compressed
        Apache Arrow representations, which can be referred to via Polars later.

    Args:
        - `search_tree: dict[str, SearchBranch]`: The search tree object, which is a
        dictionary where the keys are the geographies searched, and the values are
        instances of the `SearchBranch` dataclass containing any relevant results
        file paths.

    Returns:
        - `Result[CompiledMeta, str]`: A result type containing either an instance of
        the `CompiledMeta` dataclass, which stores each data/lazy frame with labels to
        control unpacking downstream, or an error message string (`Err(str)`).
    """

    # start with anachronistic candidates
    _search_tree_meta(
        search_tree,
        "anachron",
        "anachronistic_metadata_only_candidates.tsv",
        "anachronistics-meta.arrow",
    )

    # next, do high-distance candidates
    _search_tree_meta(
        search_tree, "highdist", "high_distance_candidates.tsv", "highdist-meta.arrow"
    )

    # and finally, double candidagtes
    _search_tree_meta(
        search_tree, "double", "double_candidate_metadata.tsv", "double-meta.arrow"
    )

    return Ok(
        CompiledMeta(
            anachron=pl.scan_ipc("anachronistics-meta.arrow"),
            highdist=pl.scan_ipc("highdist-meta.arrow"),
            double=pl.scan_ipc("double-meta.arrow"),
        )
    )


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
    assert os.path.exists(args.results_dir) and os.path.isdir(
        args.results_dir
    ), "Provided file path does not exist or is not a directory."

    # construct an initial dataclass of the geographies and top-level
    # per-geography directories to be searched while handling any errors.
    message_result = construct_file_paths(args.results_dir)
    match message_result:
        case Ok(starter_paths):
            pass
        case Err(message):
            sys.exit(
                f"Error originated while constructing file paths:\n\
                    {message}"
            )

    # use the starter paths and recur into subdirectories to map out the
    # "search tree" of results files
    tree_result = define_search_tree(starter_paths)
    match tree_result:
        case Ok(search_tree):
            pass
        case Err(message):
            sys.exit(
                f"Unable to search through provided results directories:\n\
                    {message}"
            )

    # traverse the search tree, reading files as they're found, and summarize
    # the results in a polars dataframe
    stats_result = stats_pipeline(search_tree)
    match stats_result:
        case Ok(stats_df):
            pass
        case Err(message):
            sys.exit(
                f"Compiling statistics failed with the following message:\n\
                    {message}"
            )
    stats_df.write_excel("alpine_run_statistics.xlsx", autofit=True)

    # compile all metadata for all geographies so that it is queryable downstream
    meta_result = compile_metadata(search_tree)
    match meta_result:
        case Ok(compiled):
            _anachron_meta = compiled.anachron
            _highdist_meta = compiled.highdist
            _double_meta = compiled.double
        case Err(message):
            sys.exit(
                f"Error encountered while compiling metadata for candidates:\n\
                    {message}"
            )
    
    # Count the number of candidates for the full date range of the dataset
    


if __name__ == "__main__":
    main()
