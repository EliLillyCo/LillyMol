#!/usr/bin/env ruby

ianhome = "/home/rx87851"

require 'json/ext'

require "#{ianhome}/ruby/lib/iwcmdline.rb"
#require "#{ianhome}/ruby/lib/iwcmdline_v2.rb"
require "#{ianhome}/ruby/lib/bindir.rb"

$expert = false

def usage (rc)
  $stderr.print "Converts a query file to JSON format\n"
  $stderr.print " -nq            strip off quotes\n"
  $stderr.print " -al            arrays on one line\n"
  $stderr.print " -expert        more options\n" unless ($expert)
  $stderr.print " -v             verbose output\n"
  exit(rc)
end

# IWCmdlineV2  if using v2

cl = IWCmdline.new("-v-expert-bindir=dir-nq-al-sq")

$expert = cl.option_present('expert')

if cl.unrecognised_options_encountered()
  $stderr.print "Unrecognised options encountered\n"
  usage(1)
end

verbose = cl.option_present('v')

bindir = Bindir.new(ianhome)

if (cl.option_present('bindir'))
  cl.values('bindir').each do |d|
    bindir.add_dir(d)
  end
end

if (0 == ARGV.size)
  $stderr.print "Insufficient arguments\n"
  usage(2)
end

$quoted_string_rx = Regexp.new('^"(.+)"$')
$array_rx = Regexp.new('^\((.*\s.*)\)$')


class MSI_Attribute
  attr_accessor :name, :val
  def initialize(n, v)
    @name = n
    @val = v

    m = $array_rx.match(v)
    if (m)
      @val = Array.new
      m[1].split.each do |a|
        @val.push(a)
      end
      return
    end

    m = $quoted_string_rx.match(v)
    if (m)
      @val = m[1]
      return 
    end

  end
end

$start_of_object_rx = Regexp.new('^ *\((\d+) (\S+)')
$end_of_object_rx = Regexp.new('^ *\)')
$attribute_rx = Regexp.new(' *\(A [A-Z] (\S+) (.*)\)')
$comment_rx = Regexp.new('^#')

class MSI_Object
  def initialize(i, s)
    @object_number = i
    if ("Query" == s)
      @name = 'SubstructureQuery'
    else
      @name = s
    end
    @attr = Array.new
    @children = Array.new
  end

  def build(s)
    s.each do |line|
      return true if ($end_of_object_rx.match(line))

      next if ($comment_rx.match(line))

      next if (0 == line.chomp.length)

      m = $attribute_rx.match(line)
      if (m)
        next if (m[1].downcase == "version")
        a = MSI_Attribute.new(m[1], m[2])
        @attr.push(a)
        next
      end

      m = $start_of_object_rx.match(line)
      if (m)
        o = MSI_Object.new(m[1].to_i, m[2])
        startline = s.lineno
        if (! o.build(s))
          $stderr << "Cannot parse object '#{m[1]} #{m[2]}' at line #{startline}\n"
          return false
        end

        @children.push(o)
        next
      end

      $stderr << "Skipping unrecognised line '#{line.chomp}' at line #{s.lineno}\n"
    end

    return true
  end

  def to_hash
    rc = Hash.new

    rc[:object_type] = @name

    if ("Query_Atom" == @name)
      rc[:atom_number] = @object_number
    elsif ("Query" == @name || 'SubstructureQuery' == @name)
      true
    else
      true
#     rc[:object_number] = @object_number
    end

    @attr.each do |at|
      rc[at.name] = at.val
    end

    if (@children.size > 0)
      rc[:children]  = Array.new

      @children.each do |c|
        rc[:children].push(c.to_hash)
      end
    end

    return rc;
  end

  def to_json(*a)
    rc = to_hash

    rc.to_json(*a)
  end
end

o = false

ARGF.each do |line|
# $stderr << "Examining '#{line.chomp}'\n"
  m = $start_of_object_rx.match(line)
  if (m)
    o = MSI_Object.new(m[1].to_i, m[2])

    if (! o.build(ARGF))
      $stderr << "Cannot parse top level object\n";
      exit(2)
    end
  end
end

if (! o)
  $stderr << "Huh, no object found\n"
  exit 2
end

json = JSON.pretty_generate(o)

json.gsub!('"', "") if (cl.option_present('nq'))

if (cl.option_present('al'))
  start_of_array = Regexp.new(' \[$')
  end_of_array_comma = Regexp.new('\],$')
  end_of_array = Regexp.new('\]$')
  children_rx = Regexp.new('children: \[')
  f = json.split("\n")
  newf = Array.new

  i = 0;
  while (i < f.size)
    line = f[i]
#   $stderr << "Checking for array '#{line}'\n"
    i += 1

    if (! start_of_array.match(line) || children_rx.match(line))
      newf.push(line)
      next
    end

#   $stderr << "Detected start of array in '#{line}'\n"
    while (i < f.size)
      l = f[i]
      if (end_of_array_comma.match(l))
        line << " ],"
        newf.push(line)
        i += 1
        break
      elsif (end_of_array.match(l))
        line << " ]"
        newf.push(line)
#       $stderr << "Found end of array '#{line}'\n"
        i += 1
        break
      else
        line << ' ' << l.gsub(/^\s+/, ' ')
#       newf.push(line.gsub(/^\s*/, ' '))
      end

      i += 1
    end
  end

  json = newf.join("\n") if (newf.size < f.size)
end

if (cl.option_present('sq'))
  f = json.split("\n")

  to_requote = Regexp.new('(Comment|smarts): (.+)')
  to_requote_comma = Regexp.new('(Comment|smarts): (.+),$')

  f.each_index do |i|
    line = f[i]

#   $stderr << "Examining '#{line}' for quoting\n"

    m = to_requote_comma.match(line)
    if (m)
      f[i] = line.gsub(/: (.+),$/, ": \"\\1\",")
      next
    end

    m = to_requote.match(line)
    if (m)
      f[i] = line.gsub(/: (.+)$/, ": \"\\1\"")
      next
    end
  end

  json = f.join("\n")
end

$stdout << json << "\n"


