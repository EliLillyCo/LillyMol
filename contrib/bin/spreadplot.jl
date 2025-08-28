# Consume the output of the -S option of nplotnn.

using ArgMacros

using CSV
using Plots

function crossing_point(values::Array, cutoff::Real)::Int
  for i in eachindex(values)
    values[i] <= cutoff && return i
  end

  return 1;
end

function main()
  @inlinearguments begin
    @helpusage "nplotnn -S <file> file.spr > file.spr.smi; spreadplot.jl -S spread <file>    generates spread.png"
    @helpdescription """Consumes the output of the -S option of nplotnn, generating a plot of distance vs number selected.
    -S : specifies the name of the output file - .png will be added.
    --title : option is the title of the plot, --title 'Title'
    --label : legend of the curve, --legend 'foo'
    --lwd : controls the width of the line, --lwd 8
    --color : the color of the line, --color red
    --dist : draw a horizontal line at the point where a distance is crossed, --dist 0.15
    The 
    """
    @argumentrequired String stem "-S"
    @argumentoptional Float32 threshold "--dist"
    @argumentoptional String label "--label"
    @argumentoptional String title "--title"
    @argumentdefault Int32 4 lwd "--lwd"
    @argumentdefault String "green" colour "--color"

    @positionalrequired String input_file "input_file"
  end

  data = CSV.File(input_file, header=true, delim= ",")
  distance = data.distance
  if isnothing(label)
    label = "Selected"
  end

  plot(range(1,length(distance)), distance, ylim=(0.0, distance[2]), xlabel="Number Selected",
       ylabel="Distance", lwd=lwd, color=Symbol(colour), label=label)

   if ! isnothing(title)
     plot!(title=title)
   end

  if ! isnothing(threshold)
    ndx = crossing_point(distance, threshold)
    y = distance[ndx]
    plot!([0, length(distance)], [y, y], linestyle=:dash, color=:grey, label="Dist $(threshold) count $(ndx)")
  end

  savefig("$(stem).png")
end

main()

