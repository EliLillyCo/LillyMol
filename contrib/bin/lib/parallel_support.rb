# Common functions used by the parallel tools

module PSupport

# Return true if zero of the tokens in `argv` are files.
# Return true if there are more than 1 files in `argv`.
def looks_like_multiple_files(argv)
  file_count = 0
  argv.each do |fname|
    file_count += 1 if File.exist?(fname)
  end

  if file_count == 0
    $stderr << "No files on command line\n"
    return true
  end

  if file_count > 1
    $stderr << "Cannot process multiple files\n";
    return true
  end

  return false
end

# A set of intermediate files has been generated.
# All such files must exist.
def all_files_exist(files)
  rc = true
  files.each do |fname|
    next if File.exist?(fname)
    $stderr << "Intermediate file #{fname} not found\n"
    rc = false
  end

  return rc
end

# Remove an array of files
def remove_files(files)
  files.each do |fname|
    File.unlink(fname) if File.exist?(fname)
  end
end

def quote_special_characters(argv)
  special_chars = Regexp.new('[$\[!><()@~]')
  
  result = ""
  
  argv.each do |a|
    if special_chars.match(a)
      result << " '#{a}'"
    else
      result << " #{a}"
    end
  end

  return result
end

def get_offsets(input_file, nthreads, verbose)

  file_size = File.size?(input_file)
  
  inp = File.open(input_file, mode='r')
  raise "Cannot open '#{input_file}'" unless inp
  
  chunk_size = file_size / nthreads
  
  if chunk_size < 10
    $stderr << "File '#{input_file}' too small for #{nthreads} threads\n"
    exit 2
  end
  
  
  offset = Array.new
  offset.push(0)
  
  (1...nthreads).each do |i|
    o = i * chunk_size
  
    if ( ! inp.seek(o, IO::SEEK_SET))
      raise "Cannot seek to #{o} in '#{input_file}'\n"
    end
    line = inp.gets
    offset.push(inp.pos)
  end
  
  if verbose
    $stderr << "'#{input_file}', size #{file_size} divided into #{nthreads} chunks, each #{chunk_size} bytes\n"
    offset.each do |o|
      $stderr << o << "\n"
    end
  end

  return offset
end

module_function :looks_like_multiple_files, :all_files_exist, :remove_files, :get_offsets, :quote_special_characters

end  # PSupport
