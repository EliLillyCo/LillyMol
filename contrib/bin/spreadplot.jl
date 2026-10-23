# Plot the distance trajectory produced by the -S option of nplotnn.

using ArgMacros
using CSV
using Dierckx
import Humanize: digitsep
using Plots

"""Parse one or more comma-separated distance thresholds."""
function parse_thresholds(value::AbstractString)
    fields = strip.(split(value, ','))
    if isempty(fields) || any(isempty, fields)
        throw(ArgumentError("--dist requires one or more comma-separated numbers"))
    end

    thresholds = Float64[]
    for field in fields
        threshold = tryparse(Float64, field)
        if isnothing(threshold)
            throw(ArgumentError("Invalid --dist value '$(field)'"))
        end
        if !isfinite(threshold) || !(0.0 <= threshold <= 1.0)
            throw(ArgumentError("--dist values must be between 0.0 and 1.0"))
        end
        push!(thresholds, threshold)
    end

    return thresholds
end

"""Return the first index at which `values` reaches or falls below `cutoff`."""
function crossing_point(values::AbstractVector{<:Real}, cutoff::Real)
    return findfirst(value -> value <= cutoff, values)
end

"""Return trajectory indices crossing regular distance intervals."""
function crossing_points(values::AbstractVector{<:Real}, delta::Real)
    if !isfinite(delta) || !(0.0 < delta <= 1.0)
        throw(ArgumentError(
            "--point must be greater than 0.0 and no greater than 1.0"))
    end

    nthresholds = ceil(Int, 1.0 / delta) - 1
    thresholds = (1.0 - i * delta for i in 1:nthresholds)
    indices = Int[]
    for threshold in thresholds
        threshold <= 0.0 && continue
        ndx = crossing_point(values, threshold)
        isnothing(ndx) || push!(indices, ndx)
    end

    # A sharp drop may cross several thresholds at the same trajectory point.
    return unique(indices)
end

"""Read and validate an nplotnn distance trajectory."""
function read_trajectory(input_file::AbstractString)
    data = CSV.File(input_file; header=true, delim=',')
    columns = propertynames(data)
    :sel in columns || throw(ArgumentError("Input must contain a 'sel' column"))
    :distance in columns || throw(ArgumentError("Input must contain a 'distance' column"))

    selected_raw = collect(data.sel)
    distance_raw = collect(data.distance)
    length(selected_raw) == length(distance_raw) ||
        throw(ArgumentError("The 'sel' and 'distance' columns have different lengths"))
    length(distance_raw) >= 2 ||
        throw(ArgumentError("Input must contain at least two trajectory rows"))
    any(ismissing, selected_raw) &&
        throw(ArgumentError("The 'sel' column contains missing values"))
    any(ismissing, distance_raw) &&
        throw(ArgumentError("The 'distance' column contains missing values"))

    all(value -> value isa Integer, selected_raw) ||
        throw(ArgumentError("The 'sel' column must contain integers"))
    all(value -> value isa Real, distance_raw) ||
        throw(ArgumentError("The 'distance' column must contain numbers"))

    selected = Int.(selected_raw)
    distance = Float64.(distance_raw)

    all(isfinite, distance) || throw(ArgumentError("Distances must be finite"))
    all(value -> 0.0 <= value <= 1.0, distance) ||
        throw(ArgumentError("Distances must be between 0.0 and 1.0"))
    all(value -> value >= 0, selected) ||
        throw(ArgumentError("The 'sel' column cannot contain negative values"))
    all(selected[i] > selected[i - 1] for i in 2:length(selected)) ||
        throw(ArgumentError("The 'sel' column must be strictly increasing"))
    all(distance[i] <= distance[i - 1] for i in 2:length(distance)) ||
        throw(ArgumentError("The distance trajectory must be monotonically non-increasing"))

    return selected, distance
end

function main()
    @inlinearguments begin
        @helpusage "nplotnn -S <file> file.spr > file.spr.smi; spreadplot.jl -S spread <file>    generates spread.png"
        @helpdescription """Consumes the output of the -S option of nplotnn, generating a plot of distance vs number selected.
        -S : specifies the output filename stem; the format extension will be added.
        --title : title of the plot, --title 'Title'
        --label : legend label for the curve, --label 'foo'
        --lwd : width of the trajectory line, --lwd 8
        --color : colour of the trajectory line, --color red
        --dist : draw horizontal lines where distances are crossed, --dist 0.10,0.15,0.20
        --spline : use spline interpolation, useful for large datasets and SVG output, --spline <npoints>
        --point : draw markers at regular distance intervals, --point 0.1
        """
        @argumentrequired String stem "-S"
        @argumentoptional String thresholds "--dist"
        @argumentoptional String label "--label"
        @argumentoptional String title "--title"
        @argumentoptional Int32 spline "--spline"
        @argumentoptional Float32 point "--point"
        @argumentdefault String "png" format "--format"
        @argumentdefault Int32 4 lwd "--lwd"
        @argumentdefault String "green" color "--color"

        @positionalrequired String input_file "input_file"
    end

    selected, distance = read_trajectory(input_file)
    curve_label = isnothing(label) ? "Selected" : label
    parsed_thresholds = isnothing(thresholds) ? Float64[] : parse_thresholds(thresholds)
    ymaximum = isempty(parsed_thresholds) ? distance[2] :
        max(distance[2], maximum(parsed_thresholds))
    ymaximum = max(ymaximum, eps(Float64))
    lwd > 0 || throw(ArgumentError("--lwd must be greater than zero"))

    if isnothing(spline)
        plt = plot(selected, distance;
                   ylim=(0.0, ymaximum), xlabel="Number Selected", ylabel="Distance",
                   linewidth=lwd, color=Symbol(color), label=curve_label)
    else
        spline >= 2 || throw(ArgumentError("--spline must be at least 2"))
        x = range(first(selected), last(selected); length=spline)
        interpolation = Spline1D(selected, distance; k=min(3, length(distance) - 1))
        interpolated_distance = interpolation.(x)
        plt = plot(x, interpolated_distance;
                   ylim=(0.0, ymaximum), xlabel="Number Selected", ylabel="Distance",
                   linewidth=lwd, color=Symbol(color), label=curve_label)
    end

    if !isnothing(title)
        plot!(plt; title=title)
    end

    if !isnothing(point)
        indices = crossing_points(distance, point)
        scatter!(plt, selected[indices], distance[indices]; label=false)
    end

    if !isempty(parsed_thresholds)
        for threshold in parsed_thresholds
            ndx = crossing_point(distance, threshold)
            threshold_label = if isnothing(ndx)
                "Dist $(threshold) not reached"
            else
                "Dist $(threshold) count $(digitsep(selected[ndx]))"
            end
            plot!(plt, [first(selected), last(selected)], [threshold, threshold];
                  linestyle=:dash, color=:grey, label=threshold_label)
        end
    end

    savefig(plt, "$(stem).$(format)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
