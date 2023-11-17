#!/usr/bin/env ruby
# $Id$
ianhome = ENV['C3TK_BIN']
# ianhome = "/home/rx87851"

require "#{ianhome}/ruby/lib/iwcmdline.rb"
require "#{ianhome}/ruby/lib/bindir.rb"

$expert = false

# xmg_template = "#{ianhome}/lib/spreadplot.agr"

def usage (rc)
  $stderr.print "Takes the output from gc3tk's gfp_spread and produces a distance plot\n"
  $stderr.print "gfp_spread.sh file.gfp > file.spr\n"
  $stderr.print "#{$0} file.spr > results\n"
  # $stderr.print " -xmg <stem>   create plot with xmgrace\n"
  # $stderr.print " -R <fname>    create an R script\n"
  $stderr.print " -Rplot <fname>  create plot using R (pdf or png) - must have suffix\n";
  $stderr.print " -np ... -np     options passed directly to nplotnn (-T for example)\n" if ($expert)
  $stderr.print " -bindir <dir>   directory for executables\n" if ($expert)
  $stderr.print " -expert         more options\n" unless ($expert)
  $stderr.print " -v              verbose output\n"
  exit(rc)
end

cl = IWCmdline.new("-v-expert-bindir=dir-Rplot=s-np=close")

if (cl.option_present('expert'))
  $expert = true
end

if cl.unrecognised_options_encountered()
  $stderr.print "Unrecognised options encountered\n"
  usage(1)
end

verbose = cl.option_present('v')

unless cl.option_present('Rplot')
  $stderr << "Must specify the name for the output file via the -Rplot option\n"
  usage(1)
end

img_out = cl.value('Rplot')
rfile = "#{img_out}.R"
if img_out.end_with?(".pdf")
  pdf_out = true
elsif img_out.end_with?(".png")
  png_out = true
else
  pdf_out = true # output to pdf by default
  img_out << ".pdf"
end
$stderr.print "Plot will be saved as #{img_out}"

# rstem = false
# rfile = false
# if cl.option_present('R')
#   rstem = cl.value('R')
#
#   rfile = rstem.dup
#
#   if ! rfile.end_with?(".r")
#     rfile << ".r"
#   end
# end

bindir = Bindir.new(ianhome)

if (cl.option_present('bindir'))
  cl.value('bindir').each do |d|
    bindir.add_dir(d)
  end
end

nplotnn = bindir.find_executable('nplotnn.sh')

if (0 == ARGV.size)
  $stderr.print "Insufficient arguments\n"
  usage(2)
end

if (ARGV.size > 1)
  $stderr << "Only takes one argument\n"
  usage(1)
end

fname = ARGV.shift

if (! File.size?(fname))
  $stderr << "Missing or empty input file '#{fname}'\n"
  exit 2
end

if (! /\.spr$/.match(fname))
  $stderr << "Warning, file name '#{fname}' does not match .spr\n"
end

nplotnn_options = "-n 0 -s"

nplotnn_options << ' -v' if (verbose)

if (cl.option_present('np'))
  nplotnn_options << ' ' << cl.value('np')
end

p = IO.popen("#{nplotnn} #{nplotnn_options} #{fname}", mode='r')
raise "Cannot open pipe" unless (p)

n = 0
dist = false

zdata = Array.new

nsel = Array.new
dists = Array.new

nmax = 0
dmax = 0.0

last_record = false

p.each do |line|
  f = line.split

  n += 1

  d = f.pop

  if (d != dist)
    $stdout.print "#{n} #{d}\n"
    zdata.push("#{n} #{d}\n") if (img_out)
    if rfile
      nsel.push(n)
      dists.push(d)
    end
    dist = d
    nmax = n
    dmax = d.to_f if (d.to_f > dmax)
    last_record = false
  else
    last_record = "#{n} #{d}\n"
  end
end

p.close

if (last_record)
  $stdout.print last_record 
  zdata.push(last_record)
  if rfile
    f = last_record.split
    nsel.push(f[0])
    dists.push(f[1])
  end
end

if img_out
  fname = "#{img_out}.R"
  outp = File.open(fname, mode='w')
  raise "Cannot open '#{fname}'" unless outp
  outp << 'dist=c('
  dists.shift
  ymax = dists[0]
  nsel.shift
  dists.each_with_index do |d, i|
    outp << ',' if i > 0
    outp << d
  end
  outp << ")\n"
  outp << "nsel=c("
  nsel.each_with_index do |n, i|
    outp << ',' if i > 0
    outp << n
  end
  outp << ")\n"
  if png_out
    outp << "png('#{img_out}',width=800,height=700)\n"
  elsif pdf_out
    outp << "pdf('#{img_out}',width=8,height=7)\n"
  end
  outp << "plot(nsel,dist,lwd=3,col='red',ylim=c(0,#{ymax}),xlab='Number Selected',ylab='Distance',main='Spread',type='l')\n"
  outp << "dev.off()\n"
  outp.close
  if ! system("Rscript #{fname}")
    $stderr.print "Smething went wrong with executing Rscript on #{fname} - is R in the path?\n"
  end
end

exit 0
