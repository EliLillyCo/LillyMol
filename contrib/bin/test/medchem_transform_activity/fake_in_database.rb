#!/usr/bin/env ruby

found = nil
notfound = nil
args = ARGV.dup

until args.empty?
  arg = args.shift
  found = args.shift if arg == '-F'
  notfound = args.shift if arg == '-U'
end

raise 'missing -F' unless found
raise 'missing -U' unless notfound

File.open(found, 'w') do |matched|
  File.open(notfound, 'w') do |unmatched|
    STDIN.each_line do |line|
      smiles, name = line.split
      if name.end_with?('__MCTA__add_carbon')
        target = name.start_with?('mol1__MCTA__') ? 'mol2' : 'mol3'
        matched.puts "#{smiles} #{name} #{target}"
      else
        unmatched.puts "#{smiles} #{name}"
      end
    end
  end
end
