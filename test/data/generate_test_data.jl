# Generate data for testing.
# These are a subset of Karney's dataset:
# - Karney, Charles F. F. (2010). Test set for geodesics [Data set].
#   Zenodo. https://doi.org/10.5281/zenodo.32156

import Downloads

function save_lines(filename, rows, row_indices; outdir=@__DIR__)
    open(joinpath(outdir, filename), "w") do io
        for irow in row_indices
            println(io, rows[irow])
        end
    end
end

if ! @isdefined all_data
    all_raw_data = IOBuffer()
    Downloads.download(
        "https://zenodo.org/record/32156/files/GeodTest.dat?download=1",
        all_raw_data
    )

    all_data = String(take!(all_raw_data))
    rows = split(all_data, '\n')
end

# Get subset of each type of test
lengths = [100_000; fill(50000, 8)]
names = ["random", "antipodal", "short_distance", "one_end_near_pole", "both_ends_near_poles", "nearly_meridional", "nearly_equatorial", "between_vertices", "close_to_vertices"]
i = 1
for (n, name) in zip(lengths, names)
    save_lines("$(name).txt", rows, i:100:(i+n-1))
    global i += n
end
