#!/usr/bin/env ruby
raise 'missing column' unless ARGV[0, 2] == ['-c', '2']
raise 'missing limit' unless ARGV[2] == '-r'
raise 'wrong subset' unless ARGV[4, 2] == ['-subset', '.*__MCTA__(.*)']
raise 'missing error-exit flag' unless ARGV[6] == '-whash'
limit = Integer(ARGV.fetch(3))
raise 'invalid limit' unless limit.positive?
counts = Hash.new(0)
File.foreach(ARGV.fetch(7)) do |line|
  reaction = line.split.fetch(1).split('__MCTA__').fetch(1)
  next if counts[reaction] >= limit

  counts[reaction] += 1
  print line
  exit 1 if ENV['MCTA_SAMPLE_FAIL']
end
