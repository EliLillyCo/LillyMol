#!/usr/bin/env ruby

# Given a set of fingerprints in a file, construct canonical forms.

require_relative('lib/iwcmdline')

def usage(rc)
  $stderr << "Reads a file of gfp fingerprint specifications and writes them in a canonical order\n"
  $stderr << " -summarise       output is from svmfp_summarise_results\n"
  $stderr << " -skip <n>        skip the first <n> records, same as 'sed 1,<n>d | gfp_caninical_fingerprint -'\n"
  $stderr << " -v               verbose output - provides statistics on number times fingerprint found\n"
  exit rc
end

def parse_from_summarise(line)
  fingerprints = line.split[0]
  fps =  fingerprints.split('-').select{ |fp| fp.length > 0}.map { |fp|  "-#{fp}" }.sort.uniq

  return fps
end

def main
  cl = IWCmdline.new("-v-summarise-skip=ipos")
  if cl.unrecognised_options_encountered
    $stderr << "unrecognised_options_encountered\n"
    usage(1)
  end

  verbose = cl.option_present('v')

  from_summarise = cl.option_present('summarise')

  if ARGV.empty?
    $stderr << "Must specify file of fingerprints as an argument\n"
    usage(1)
  end

  skip = if cl.option_present('skip')
           cl.value('skip')
         else
           0
         end

  seen_fp = Hash.new(0)
  seen_line = Hash.new

  fingerprints_read = 0
  duplicates_suppressed = 0

  ARGF.drop(skip).each do |line|
    fingerprints_read += 1
    line = line.chomp
    if from_summarise
      f = parse_from_summarise(line)
    else
      f = line.split.sort.uniq
    end
    f.each do |fp|
      seen_fp[fp] += 1
    end
    if seen_line[f]
      duplicates_suppressed += 1
      next
    end
    seen_line[f] = true

    $stdout << "#{f.join(' ')}\n"
  end

  if verbose
    $stderr << "Read #{fingerprints_read} suppressed #{duplicates_suppressed} duplicates\n"
    seen_fp.sort_by { |k, v| -v}.each do |fp, count|
      $stderr << "#{fp} #{count}\n"
    end
  end
end

main
