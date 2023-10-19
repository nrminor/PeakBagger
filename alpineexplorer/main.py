#!/usr/bin/env python3

"""
Placeholder script
"""

import sys
import glob
from typing import Optional
from dataclasses import dataclass
from result import Result, Ok, Err
import polars as pl

@dataclass
class SearchBranch:
    """
    Search branch contains the necessary file paths for each geography
    within the broader search tree.
    """
    parent_dir: str
    geography: str
    double: Optional[str]
    anachron: Optional[str]
    highdist: Optional[str]
    early_stats: Optional[str]
    late_stats: Optional[str]


def construct_file_paths(result_root: str) -> Result[dict[str, str], str]:
    """
    More to come soon. This function is just a placeholder.
    """

    tmp_dict = dict({f"{result_root}": "Hello, world!"})

    if len(tmp_dict) == 0:
        return Err("No subdirectories found in provided results directory.")

    return Ok(tmp_dict)

def define_search_tree(path_dict: dict[str, str]) -> Result[dict[str, SearchBranch], str]:
    """
    More to come soon. This function is just a placeholder.
    """

    search_tree = {}

    for geo, path in path_dict.items():
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
            late_stats=late_stats[0] if late_stats else None
        )

    return Ok(search_tree)

def get_early_count(path: Optional[str]) -> Optional[int]:
    """
    The function `get_early_count` reads a file of statistics from early in
    the pipeline meant to describe the number of input sequences. `get_early_count`
    does so for each geography independently so that a dataframe of statistics
    can be generated in a vectorized manner.
    """

    if path is None:
        return None

    # Read the TSV file into a DataFrame
    stats_df = pl.read_csv(path, separator='\t')

    # Get the count from the first row of the 'num_seqs' column
    count = stats_df.select(pl.col("num_seqs")).to_series().to_list()[0]

    return count

def get_late_count(path: Optional[str]) -> Optional[int]:
    """
    The function `get_late_count` reads a file of statistics from the end of
    the pipeline meant to describe the number of double candidate sequences,
    which is to say sequences from viral lineages that are both highly evolved 
    and anachronistic. `get_late_count` does so for each geography independently
    so that a dataframe of statistics can be generated in a vectorized manner.
    """

    # propagate None if that is the input
    if path is None:
        return None

    # Read the TSV file into a DataFrame
    stats_df = pl.read_csv(path, separator='\t')

    # Get the count from the first row of the 'num_seqs' column
    count = stats_df.select(pl.col("num_seqs")).to_series().to_list()[0]

    return count

def summarize_anachrons(path: Optional[str]) -> Optional[int]:
    """
    The function `summarize_anachrons` finds the metadata for anachronistic
    sequences and counts the number of entries, if any, for the provided
    geography and associated filepath.
    """

    # propagate None if that is the input
    if path is None:
        return None

    # Construct the file path for the TSV file
    file_path = f"{path}/anachronistic_metadata_only_candidates.tsv"

    # Read the TSV file into a DataFrame
    anachron_df = pl.read_csv(file_path, separator='\t')

    # Return the number of rows in the DataFrame
    return anachron_df.shape[0]

def summarize_highdist(path: Optional[str]) -> Optional[int]:
    """
    The function `summarize_highdist` finds the metadata for high distance
    sequences and counts the number of entries, if any, for the provided
    geography and associated filepath.
    """

    # propagate None if that is the input
    if path is None:
        return None

    # Construct the file path for the TSV file
    file_path = f"{path}/high_distance_candidates.tsv"

    # Read the TSV file into a DataFrame
    anachron_df = pl.read_csv(file_path, separator='\t')

    # Return the number of rows in the DataFrame
    return anachron_df.shape[0]

def main() -> None:
    """
    More to come soon. This script is just a placeholder.
    """

    message_result = construct_file_paths(".")
    match message_result:
        case Ok(result):
            message = result["."]
        case Err(message):
            sys.exit(f"Error originated while constructing file paths:\n\
                     {message}")

    print(f"{result}")

if __name__ == "__main__":
    main()
