module ALPINEExplorer

using DataFrames, CSV, Pipe, ErrorTypes, Glob, ExcelFiles, Arrow

export construct_file_paths, define_search_tree, stats_pipeline, make_arrow_databases

"""
The function `construct_file_paths` takes what can be thought of as the Results root directory,
checks for any subdirectories, and returns a dictionary of geography-to-filepath. If checks
were not successful, the function returns a helpful error message to be handled elsewhere.

Arguments:
 - `root_dir::String`: the file path to the results directory to be searched.

Returns:
 - `Result{Dict{String, String}, ErrorException}`: a Rust-like enum that is either
 the dictionary in the successful case, or an exception in the unsuccessful case.
"""
function construct_file_paths(root_dir::String)::Result{Dict{String, String}, ErrorException}

    paths = filter(isdir, readdir(root_dir; join=true))

    if length(paths) == 0
        return Err(ErrorException("No subdirectories found in the provided directory."))
    end

    geographies = @pipe paths |>
        split.(_, "/"; keepempty=false) |>
        last.(_) |>
        replace.(_, "GISAID_" => "") |>
        replace.(_, "LocalDataset_" => "") |>
        replace.( _, "GenBank_" => "") |>
        replace.(_, "_" => " ")
    
    if length(geographies) == 0
        return Err(ErrorException("No subdirectories present with the expected prefixes,\
        i.e., 'GISAID_', 'GenBank_', or 'LocalDataset_."))
    end

    geo_paths = Dict(geo => path for (geo, path) in zip(geographies, paths))

    return Ok(geo_paths)

end

"""
Here we define the immutable struct `SearchBranch`, which will store the necessary
file paths for each geography within the broader search tree. This struct may be
made mutable so that it can be filled recursively in the case of there being more
than two levels of results directories to traverse.
"""
struct SearchBranch
    parent_dir::String
    geography::String
    double::Union{String, Missing}
    anachron::Union{String, Missing}
    highdist::Union{String, Missing}
    early_stats::Union{String, Missing}
    late_stats::Union{String, Missing}
end

"""
The function `define_search_tree` creates a dictionary to guide the compilation of
results. This dictionary assumes only two levels of results directories at this time,
though this may change. The dictionary keys in this case are the geographies, and their
values are all available information stored in `SearchBranch` structs. Downstream
functions will iterate through this search tree and process each branch. Note that this
function is written such that it will return empty search branches, which must be
handled downstream.

Arguments:
 - `path_dict::Dict{String, String}`: The dictionary of geographies and the associated
 directories.
 - `root_dir::String`: a string specifying the results "root" directory within which to
 search subdirectories.

Returns:
 - `Dict{Symbol, SearchBranch}`: A dictionary specifying the search tree, where each
 geography has its own associated `SearchBranch` of files and directories to search.
"""
function define_search_tree(path_dict::Dict{String, String})::Dict{Symbol, SearchBranch}

    search_tree = Dict{Symbol,SearchBranch}()
    search_tree = Dict{Symbol, SearchBranch}()
    for (geo, path) in path_dict

        # check for paths
        double =
            if length(readdir(glob"*double_candidates", path)) >
                0
                readdir(glob"*double_candidates", path)[1]
            else
                missing
            end
        anachron =
            if length(readdir(glob"*metadata_candidates", path)) >
            0
                readdir(glob"*metadata_candidates", path)[1]
            else
                missing
            end
        highdist =
            if length(readdir(glob"*high_distance_clusters", path)) >
            0
                readdir(glob"*high_distance_clusters", path)[1]
            else
                missing
            end
        early_stats =
            if length(readdir(glob"*early_stats.tsv", path)) >
            0
                readdir(glob"*early_stats.tsv", path)[1]
            else
                missing
            end
        late_stats =
            if length(readdir(glob"*late_stats.tsv", path)) >
            0
                readdir(glob"*late_stats.tsv", path)[1]
            else
                missing
            end
        search_tree[Symbol(geo)] = SearchBranch(
            path,
            geo,
            double,
            anachron,
            highdist,
            early_stats,
            late_stats,
        )
    end

    return search_tree

end

"""
The function `get_early_count` reads a file of statistics from early in
the pipeline meant to describe the number of input sequences. `get_early_count`
does so for each geography independently so that a dataframe of statistics
can be generated in a vectorized manner.
"""
function get_early_count(path::Union{String, Missing})

    if path === missing
        return missing
    end

    stats_df = CSV.read(path, delim='\t', DataFrame)

    count = stats_df.num_seqs[1]

    return count

end

"""
The function `get_late_count` reads a file of statistics from early in
the pipeline meant to describe the number of double candidate sequences,
which is to say sequences from viral lineages that are both highly evolved 
and anachronistic. `get_late_count` does so for each geography independently
so that a dataframe of statistics can be generated in a vectorized manner.
"""
function get_late_count(path::Union{String, Missing})

    if path === missing
        return missing
    end

    stats_df = CSV.read(path, delim='\t', DataFrame)

    count = stats_df.num_seqs[1]

    return count

end

"""
The function `summarize_anachrons` finds the metadata for anachronistic
sequences and counts the number of entries, if any, for the provided
geography and associated filepath.
"""
function summarize_anachrons(path::Union{String, Missing})

    if path === missing
        return missing
    end

    anachron_df = CSV.read("$path/anachronistic_metadata_only_candidates.tsv",
                            delim='\t', DataFrame)

    return nrow(anachron_df)

end

"""
The function `summarize_highdist` finds the metadata for high distance
sequences and counts the number of entries, if any, for the provided
geography and associated filepath.
"""
function summarize_highdist(path::Union{String, Missing})

    if path === missing
        return missing
    end

    high_dist_df = CSV.read("$path/high_distance_candidates.tsv",
                            delim='\t', DataFrame)

    return nrow(high_dist_df)

end

"""
The function `stats_pipeline` makes heavy use of `Pipe.jl` to daisy-chain
a number of functions that help compile a data frame of statistics for all
ALPINE runs in the provided results directory.
"""
function stats_pipeline(search_tree::Dict{Symbol, SearchBranch})::Result{DataFrame, ErrorException}

    stats_df = @pipe DataFrame("Geography" => String.(keys(search_tree))) |>

        # Get the number of input sequences
        hcat(_, DataFrame("Input Sequence Count" => 
            [get_early_count(branch.early_stats) for (geo, branch) in search_tree])) |>
        
        # add a column for the number of double candidates
        hcat(_, DataFrame("Double Candidate Count" => 
            [get_late_count(branch.late_stats) for (geo, branch) in search_tree])) |>
        
        # generate a percentage of the total inputs that were flagged as double candidates
        transform(_, [:("Double Candidate Count"), :("Input Sequence Count")] =>
            ByRow((x, y) -> ifelse(ismissing(x) || ismissing(y), missing, (x / y) * 100)) => 
            :("Double Candidate Prevalence (%)")) |>
        
        # Compute the rate of double candidates in a 1 in X format
        transform(_, :"Double Candidate Prevalence (%)" => 
            ByRow(x -> ifelse(ismissing(x) || x == 0, "None", "1 in $(floor(1 / (x / 100)))")) => 
            :"Double Candidate Rate") |>

        # count the number of anachronistic lineages
        hcat(_, DataFrame("Anachronistic Count" =>
            [summarize_anachrons(branch.anachron) for (geo, branch) in search_tree])) |>

        # generate a percentage of the total inputs that were flagged as anachronistics
        transform(_, [:("Anachronistic Count"), :("Input Sequence Count")] =>
            ByRow((x, y) -> ifelse(ismissing(x) || ismissing(y), missing, (x / y) * 100)) => 
            :("Anachronistic Prevalence (%)")) |>
        
        # Compute the rate of anachronistic candidates in a 1 in X format
        transform(_, :"Anachronistic Prevalence (%)" => 
            ByRow(x -> ifelse(ismissing(x) || x == 0, "None", "1 in $(floor(1 / (x / 100)))")) => 
            :"Anachronistic Rate") |>
        
        # get the number of high distance candidates
        hcat(_, DataFrame("High Distance Count" =>
            [summarize_highdist(branch.highdist) for (geo, branch) in search_tree])) |>
        
        # compute the percentage of inputs that were high distance
        transform(_, [:("High Distance Count"), :("Input Sequence Count")] =>
            ByRow((x, y) -> (x / y) * 100) => :("High Distance Prevalence (%)")) |>
        
        # Compute the rate of high distance candidates in a 1 in X format
        transform(_, :"High Distance Prevalence (%)" => 
            ByRow(x -> ifelse(ismissing(x) || x == 0, "None", "1 in $(floor(1 / (x / 100)))")) => 
            :"High Distance Rate")
        
    # save as an excel file
    save("alpine_run_statistics.xlsx", stats_df)
    
    if nrow(stats_df) == 0
        return(Err(ErrorException(
            "No results could be compiled."
        )))
    end
    
    return Ok(stats_df)

end

"""
"""
function make_arrow_databases(search_tree::Dict{Symbol, SearchBranch})

    # First, we'll handle anachronistics
    anachron_arrow = "anachronistics-meta.arrow"
    mktempdir() do temp_path
        ticker = 0
        for (geo, branch) in search_tree
    
            if !isfile("$(branch.anachron)/anachronistic_metadata_only_candidates.tsv")
                continue
            end
    
            ticker += 1
    
            if ticker == 1
                @pipe CSV.read("$(branch.anachron)/anachronistic_metadata_only_candidates.tsv",
                            delim='\t', DataFrame) |>
                    transform(_, :Accession => ByRow(n -> String(geo)) => :Geography) |>
                    CSV.write("$temp_path/anachron.tsv", _, delim = '\t')
                continue
            end
    
            @pipe CSV.read("$(branch.anachron)/anachronistic_metadata_only_candidates.tsv",
                            delim='\t', DataFrame) |>
                    transform(_, :Accession => ByRow(n -> String(geo)) => :Geography) |>
                    CSV.write("$temp_path/anachron.tsv", _, delim = '\t', append=true)
    
        end
        Arrow.write(anachron_arrow, CSV.File("$temp_path/anachron.tsv"; delim='\t'))
    end
    
    # Next, the highly evolved/high-distance
    highdist_arrow = "highdist-meta.arrow"
    mktempdir() do temp_path
        ticker = 0
        for (geo, branch) in search_tree
    
            if !isfile("$(branch.highdist)/high_distance_candidates.tsv")
                continue
            end
    
            ticker += 1
    
            if ticker == 1
                @pipe CSV.read("$(branch.highdist)/high_distance_candidates.tsv",
                            delim='\t', DataFrame) |>
                    transform(_, :Accession => ByRow(n -> String(geo)) => :Geography) |>
                    CSV.write("$temp_path/highdist.tsv", _, delim = '\t')
                continue
            end
    
            @pipe CSV.read("$(branch.highdist)/high_distance_candidates.tsv",
                        delim='\t', DataFrame) |>
                transform(_, :Accession => ByRow(n -> String(geo)) => :Geography) |>
                CSV.write("$temp_path/highdist.tsv", _, delim = '\t', append=true)
    
        end
        Arrow.write(highdist_arrow, CSV.File("$temp_path/highdist.tsv"; delim='\t'))
    end
    
    # And finally, the double candidates
    double_arrow = "double-meta.arrow"
    mktempdir() do temp_path
        ticker = 0
        for (geo, branch) in search_tree
    
            if !isfile("$(branch.double)/double_candidate_metadata.tsv")
                continue
            end
    
            ticker += 1
    
            if ticker == 1
                @pipe CSV.read("$(branch.double)/double_candidate_metadata.tsv",
                            delim='\t', DataFrame) |>
                    transform(_, :Accession => ByRow(n -> String(geo)) => :Geography) |>
                    CSV.write("$temp_path/double.tsv", _, delim = '\t')
                continue
            end
    
            @pipe CSV.read("$(branch.double)/double_candidate_metadata.tsv",
                        delim='\t', DataFrame) |>
                transform(_, :Accession => ByRow(n -> String(geo)) => :Geography) |>
                CSV.write("$temp_path/double.tsv", _, delim = '\t', append=true)
    
        end
        Arrow.write(double_arrow, CSV.File("$temp_path/double.tsv"; delim='\t'))
    end

    return (anachron_arrow, highdist_arrow, double_arrow)
    
end

end
