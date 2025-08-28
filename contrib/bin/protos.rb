#!/usr/bin/env ruby

require 'find'

# Identify all the protos in a LillyMol distribution and report their information
# and possible links


# Look for `needle` in `fname`.
# If found, return the following line.
# grep -A 1 needle fname
def find_in_file(needle, fname)
  gathering_lines = false
  result = []
  File.foreach(fname) do |line|
    if gathering_lines
      return result if line.chomp.empty?
      result << line.gsub(/^\/\/ /, "").chomp
      next
    end

    next unless line.include?(needle)

    gathering_lines = true
  end

  return ""
end

def main
  repos = [
    'https://github.com/EliLillyCo/LillyMol',
    'https://github.com/EliLillyCo/LillyMolPrivate',
    'https://github.com/IanAWatson/LillyMol'
  ]

  $stdout << "Online documentation at\n\n"
  $stdout << repos.join("\n")
  $stdout << "\n\n"

  lillymol_home = if ENV.key?('LILLYMOL_HOME')
                    ENV['LILLYMOL_HOME']
                  else
                    # script is assumed to be in LILLYMOL_HOME/contrib/bon
                    File.realpath(File.dirname(File.dirname(File.dirname(__FILE__))))
                  end
  src = File.join(lillymol_home, 'src')
  unless File.directory?(src)
    $stderr << "Where is #{src}, cannot continue\n"
    return 1
  end

  proto_rx = Regexp.new('\.proto$')

  Find.find(src) do |fname|
    next unless proto_rx.match(fname)
    msg = find_in_file('LillyMol:UserDocumentation', fname)
    next if msg.empty?

    $stdout << "#{fname}\n"
    msg.each do |line|
      $stdout << "  #{line}\n"
    end
    $stdout << "\n"
  end
end

main
