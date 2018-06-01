# Mimic many of the functions of the C++ cmdline object
# $Id$
class IWCmdline
  def initialize (opts)
    @_count = Hash.new{0}    # count of how many times each argument encountered

    opts.freeze

    if 0 == opts.length
      return
    end

    @_unrecognised_options = Array.new

    okargs = Hash.new

    opts.split('-').each { |opt|
      next unless opt.length > 0

      opt.gsub!(/ /, "")   # spaces cannot be part of this

      if opt =~ /\S+=\S+/      # of the form option=requirement
        a = opt.split('=')
        okargs[a[0]] = a[1]
      else
        okargs[opt] = ""
      end
    }

#   okargs.each { |key, value|
#     $stderr.print "Qualifiers for option '#{key}' '#{value}'\n"
#   }

    starts_with_dash = Regexp.new("^-+")

    @_option_value = Hash.new

    iptr = 0
    while (iptr < ARGV.size)
      opt = "#{ARGV[iptr]}"

# important design decision:  should it stop grabbing arguments once
# it encounters a non-option?

      unless starts_with_dash.match(opt)
        iptr += 1
        break
#       next
      end

      if ("-" == opt && iptr == ARGV.size - 1)
        break
      end

      opt.gsub!(/^-+/, "")

#     $stderr.print "Examining argument '#{opt}'\n"

      unless okargs.has_key?(opt)
#       $stderr.print "Unrecognised option '#{opt}'\n"
        @_unrecognised_options.push(opt)
        ARGV.delete_at(iptr)
        next
      end

      if @_count.has_key?(opt)
        @_count[opt] += 1
      else
        @_count[opt] = 1
        @_option_value[opt] = Array.new
      end

      qualifiers = okargs[opt]

      ARGV.delete_at(iptr)

      next unless qualifiers.length > 0

      if ("s" == qualifiers)
        tmp = ARGV[iptr]
        @_option_value[opt].push(tmp)
        ARGV.delete_at(iptr)
        next
      end

      if ("i" == qualifiers || "int" == qualifiers)
        tmp = ARGV[iptr]
        raise "Invalid integer qualifier for -#{opt} option '#{tmp}'" unless tmp =~ /^-*\d+$/
        @_option_value[opt].push(tmp.to_i)
        ARGV.delete_at(iptr)
        next
      end

      if ("u" == qualifiers || "uint" == qualifiers)
        tmp = ARGV[iptr]
        raise "Invalid unsigned integer qualifier for -#{opt} option '#{tmp}'" unless tmp =~ /^\d+$/
        @_option_value[opt].push(tmp.to_i)
        ARGV.delete_at(iptr)
        next
      end

      if ("ipos" == qualifiers)
        tmp = ARGV[iptr]
        raise "Invalid positive integer qualifier for -#{opt} option '#{tmp}'" unless tmp =~ /^\d+$/ && tmp.to_i > 0
        @_option_value[opt].push(tmp.to_i)
        ARGV.delete_at(iptr)
        next
      end


      if ('f' == qualifiers || 'float' == qualifiers)
        tmp = ARGV[iptr]
        begin
          f = Float(tmp)
        rescue
          raise "Invalid float for '-#{opt}', '#{tmp}'\n"
        end
        @_option_value[opt].push(tmp.to_f)
        ARGV.delete_at(iptr)
        next
      end

      if ('fraction' == qualifiers)
        tmp = ARGV[iptr];       # Must be some way of getting this to work. Want an exception if an invalid number
        begin
          f = Float(tmp);
        rescue
          raise "Invalid fraction for '-#{opt}', '#{tmp}'\n"
        end
        raise "Invalid fraction '#{f}'" unless (f >= 0.0 && f <= 1.0)
        @_option_value[opt].push(f)
        ARGV.delete_at(iptr)
        next
      end

      if ("sfile" == qualifiers)
        tmp = ARGV[iptr];
        raise "Must specify file name for option '-#{opt}'" unless (tmp)
        raise "Missing or empty file '#{tmp}'" unless (FileTest.size?(tmp))
        @_option_value[opt].push(tmp)
        ARGV.delete_at(iptr)
        next
      end

      if ("xfile" == qualifiers)
        tmp = ARGV[iptr];
        raise "Must specify file name for option '-#{opt}'" unless (tmp)
        raise "Missing or empty file '#{tmp}'" unless (FileTest.executable_real?(tmp))
        @_option_value[opt].push(tmp)
        ARGV.delete_at(iptr)
        next
      end

      if ("dir" == qualifiers)
        tmp = ARGV[iptr];
        raise "Must specify file name for option '-#{opt}'" unless (tmp)
        raise "Missing or invalid directory file '#{tmp}'" unless (FileTest.directory?(tmp))
        @_option_value[opt].push(tmp)
        ARGV.delete_at(iptr)
        next
      end

      if ("close" == qualifiers)
        gotclose = false
        closing_option = Regexp.new("^-+#{opt}$")
        tmp1 = String.new("")
        while (iptr < ARGV.size)
          tmp2 = ARGV[iptr]
          ARGV.delete_at(iptr)
          if (closing_option.match(tmp2))
            gotclose = true
            break
          else
            if (tmp1.length > 0)
              tmp1 << ' ' << tmp2
            else
              tmp1 << tmp2
            end
          end
        end
        raise "No closing -#{opt}" unless gotclose
        @_option_value[opt].push(tmp1)
        next
      end

      $stderr.print "IWCmdline:initialize:invalid qualifier '#{qualifiers}' for option '#{opt}', value discarded\n"
      ARGV.delete_at(iptr)    # just discard it???!
    end
  end

  def debug_print
  end

  def unrecognised_options_encountered
    return @_unrecognised_options.size > 0
  end

  def unrecognised_options
    return @_unrecognised_options
  end

  def option_present (opt)
    return (@_count.has_key?(opt))
  end

  def option_count (opt)
    return 0 unless @_count.has_key?(opt)

    return @_count[opt]
  end

  def value(opt, ndx = 0)
    return false unless @_option_value.has_key?(opt)

    tmp = @_option_value[opt]

    return tmp[ndx]
  end

  def value_or_empty_string (opt)
    return "" unless @_option_value.has_key?(opt)

    tmp = @_option_value[opt]

    return tmp[0]
  end

  def values (opt)
#   $stderr.print "Do we have '#{opt}' " << (@_option_value.has_key? (opt)).to_s << "\n"

    return [] unless @_option_value.has_key?(opt)

    return @_option_value[opt]
  end

  def number_of_these_options_set (opt)
    rc = 0
    opt.each { |o|
      rc += 1 if @_option_value.has_key?(o)
#     $stderr.print "After checking '#{o}' rc = #{rc}\n"
    }

    rc
  end

# We may decide to combine two options. 

  def combine_options (o1, o2)
#   $stderr.print "Combining '#{o1}' with '#{o2}'\n"
    if (! @_option_value.has_key?(o1) && ! @_option_value.has_key?(o2))
      return
    elsif (! @_option_value.has_key?(o1))
      @_option_value[o1] = Array.new
    elsif (! @_option_value.has_key?(o2))
      @_option_value[o2] = Array.new
    end
     
    tmp = @_option_value[o1] | @_option_value[o2]

    @_option_value[o1] = tmp
    @_option_value[o2] = tmp
    @_count[o1] = tmp.size
    @_count[o2] = tmp.size
  end

  def values_as_array(opt, separator)

    return [] unless @_option_value.has_key?(opt)

    rc = Array.new

    @_option_value[opt].each do |o|
      f = o.split(separator)
      f.each do |x|
        rc.push(x)
      end
    end

    return rc
  end

end
