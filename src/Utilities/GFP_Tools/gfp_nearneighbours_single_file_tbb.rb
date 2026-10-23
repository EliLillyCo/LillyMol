#!/usr/bin/env ruby

# Write code for gfp_nearneighbours_single_file_tbb.

nthreads = 32
# We need to do all pairs.
todo = []
(0...nthreads).each do |i|
  (i + 1...nthreads).each do |j|
    todo << [i,j]
  end
end

# Write a group of tasks that are independent.
def write_group(group, this_group, active_this_group, diagonal_done)
  $stderr << this_group << "\n"
  grp = "g#{group}"
  $stdout << "  tbb::task_group #{grp};\n"
  this_group.each do |g|
    i = g[0]
    j = g[1]
    $stdout << "  #{grp}.run([&] { gfp_nearneighbours(pool, p[#{i}], p[#{i+1}], p[#{j}], p[#{j+1}]); });\n"
  end

  diagonal_done.each_with_index do |d, ndx|
    next if d
    next if active_this_group[ndx]
    $stdout << "  gfp_nearneighbours_diagonal(pool, p[#{ndx}], p[#{ndx+1}]);\n"
    diagonal_done[ndx] = true
    break
  end
  $stdout << "  #{grp}.wait();\n"
  $stdout << "  if (verbose) {\n"
  $stdout << "    cerr << \" end #{grp}\\n\";\n"
  $stdout << "  }\n"
end

# Diagonals actually get done all together at the end.
diagonal_done = Array.new(nthreads, false)
pair_done = Hash.new
group = 0
while todo.size > 0
  active_this_group = Array.new(nthreads, false)
  this_group = []
  this_group << todo.shift
  active_this_group[this_group[0][0]] = true
  active_this_group[this_group[0][1]] = true
  ndx = 0
  while ndx < todo.size
    i = todo[ndx][0]
    j = todo[ndx][1]
    if active_this_group[i] || active_this_group[j]
      ndx += 1
      next
    end
    this_group << todo.delete_at(ndx)
    active_this_group[i] = true
    active_this_group[j] = true
    break if this_group.size == 16   # seems to never happen
  end
  $stderr << "Found group of size #{this_group.size} items, active #{active_this_group.count(true)}\n"
  write_group(group, this_group, active_this_group, diagonal_done)
  group += 1
end

# Diagonals will be left over.
group_size = 0
grp = ""
diagonal_done.each_with_index do |d, ndx|
  next if d
  if group_size == 0
    group += 1
    grp = "g#{group}"
    $stdout << "  tbb::task_group #{grp};\n"
  end
  $stdout << "  #{grp}.run([&] { gfp_nearneighbours_diagonal(pool, p[#{ndx}], p[#{ndx+1}]); });\n"

  group_size += 1
  if group_size == nthreads / 2
    $stdout << "  #{grp}.wait();\n"
    $stdout << "  if (verbose) {\n"
    $stdout << "    cerr << \" end #{grp}\\n\";"
    $stdout << "  }\n"
    group_size = 0
  end
end

exit
