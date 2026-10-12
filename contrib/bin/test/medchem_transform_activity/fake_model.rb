#!/usr/bin/env ruby

input = ARGV.fetch(0)
output = ARGV.fetch(1)

predictions = {
  'add_oxygen' => 3.0,
  'add_nitrogen' => 1.5
}

File.open(output, 'w') do |out|
  out.puts 'id prediction'
  File.foreach(input) do |line|
    tokens = line.split
    next if tokens.empty?

    _smiles, name = tokens
    transformation = name.split('__MCTA__').fetch(1)
    prediction = predictions.fetch(transformation)
    out.puts [name, prediction].join(' ')
  end
end
