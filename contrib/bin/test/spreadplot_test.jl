using Test

const SPREADPLOT = normpath(joinpath(@__DIR__, "..", "spreadplot.jl"))
include(SPREADPLOT)

@testset "parse_thresholds" begin
    @test parse_thresholds("0.15") == [0.15]
    @test parse_thresholds("0.1, 0.25,1") == [0.1, 0.25, 1.0]

    for value in ("", "0.1,", ",0.1", "abc", "-0.1", "1.1", "NaN", "Inf")
        @test_throws ArgumentError parse_thresholds(value)
    end
end

@testset "crossing_point" begin
    values = [1.0, 0.7, 0.4]
    @test crossing_point(values, 1.0) == 1
    @test crossing_point(values, 0.7) == 2
    @test crossing_point(values, 0.5) == 3
    @test isnothing(crossing_point(values, 0.2))
end

@testset "crossing_points" begin
    values = [1.0, 0.8, 0.5, 0.2, 0.0]
    @test crossing_points(values, 0.1) == [2, 3, 4, 5]
    @test crossing_points(values, 0.25) == [3, 4]
    @test isempty(crossing_points(values, 1.0))

    for delta in (0.0, -0.1, 1.1, Inf, NaN)
        @test_throws ArgumentError crossing_points(values, delta)
    end
end

function write_csv(directory::AbstractString, name::AbstractString,
                   contents::AbstractString)
    fname = joinpath(directory, name)
    write(fname, contents)
    return fname
end

@testset "read_trajectory" begin
    mktempdir() do directory
        valid = write_csv(directory, "valid.csv",
                          "sel,distance\n0,1.0\n10,0.8\n20,0.5\n30,0.2\n")
        selected, distance = read_trajectory(valid)
        @test selected == [0, 10, 20, 30]
        @test distance == [1.0, 0.8, 0.5, 0.2]

        invalid_inputs = [
            ("missing_sel.csv", "count,distance\n0,1.0\n1,0.5\n"),
            ("missing_distance.csv", "sel,value\n0,1.0\n1,0.5\n"),
            ("one_row.csv", "sel,distance\n0,1.0\n"),
            ("missing_value.csv", "sel,distance\n0,1.0\n1,\n"),
            ("non_integer_sel.csv", "sel,distance\n0,1.0\n1.5,0.5\n"),
            ("negative_sel.csv", "sel,distance\n-1,1.0\n0,0.5\n"),
            ("repeated_sel.csv", "sel,distance\n0,1.0\n0,0.5\n"),
            ("decreasing_sel.csv", "sel,distance\n1,1.0\n0,0.5\n"),
            ("distance_too_large.csv", "sel,distance\n0,1.0\n1,1.1\n"),
            ("negative_distance.csv", "sel,distance\n0,1.0\n1,-0.1\n"),
            ("increasing_distance.csv", "sel,distance\n0,0.5\n1,0.6\n"),
        ]

        for (name, contents) in invalid_inputs
            fname = write_csv(directory, name, contents)
            @test_throws ArgumentError read_trajectory(fname)
        end
    end
end

@testset "command line rendering" begin
    mktempdir() do directory
        input_file = write_csv(directory, "trajectory.csv",
                               "sel,distance\n0,1.0\n10,0.8\n20,0.5\n30,0.2\n40,0.1\n")

        direct_stem = joinpath(directory, "direct")
        direct_command = `$(Base.julia_cmd()) $SPREADPLOT -S $direct_stem --dist 0.9,0.5,0.05 --point 0.1 $input_file`
        run(addenv(direct_command, "GKSwstype" => "100"))
        direct_output = direct_stem * ".png"
        @test isfile(direct_output)
        @test filesize(direct_output) > 0

        spline_stem = joinpath(directory, "spline")
        spline_command = `$(Base.julia_cmd()) $SPREADPLOT -S $spline_stem --dist 0.5 --spline 50 --format svg $input_file`
        run(addenv(spline_command, "GKSwstype" => "100"))
        spline_output = spline_stem * ".svg"
        @test isfile(spline_output)
        @test filesize(spline_output) > 0
    end
end
