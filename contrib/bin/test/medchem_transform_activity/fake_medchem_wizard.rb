#!/usr/bin/env ruby

smiles = ARGV.fetch(-1)
separator = ARGV.fetch(ARGV.index('-W') + 1)
raise 'missing -V' unless ARGV.include?('-V')

File.foreach(smiles) do |line|
  tokens = line.split
  next if tokens.empty?

  _smiles, id = tokens
  case id
  when 'mol1'
    puts "CC mol1#{separator}add_carbon"
    puts "CO mol1#{separator}add_oxygen"
  when 'mol2'
    puts "CCC mol2#{separator}add_carbon"
    puts "CN mol2#{separator}add_nitrogen"
  end
end
