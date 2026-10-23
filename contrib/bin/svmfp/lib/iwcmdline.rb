# frozen_string_literal: true

# Mimic many of the functions of the C++ cmdline object
class IWCmdline # rubocop:disable Metrics/ClassLength
  def initialize(opts)
    @_count = Hash.new(0) # count of how many times each argument encountered

    opts.freeze

    return unless opts.length.positive?

    @_unrecognised_options = []

    okargs = {}

    opts.split('-').each do |opt|
      next if opt.empty?

      opt.gsub!(/ /, '') # spaces cannot be part of this

      if opt =~ /\S+=\S+/ # of the form option=requirement
        a = opt.split('=')
        okargs[a[0]] = a[1]
      else
        okargs[opt] = ''
      end
    end

    starts_with_dash = Regexp.new('^-+')

    @_option_value = {}

    iptr = 0
    while iptr < ARGV.size
      opt = ARGV[iptr].to_s

      # important design decision:  should it stop grabbing arguments once
      # it encounters a non-option?

      unless starts_with_dash.match(opt)
        iptr += 1
        break
      end

      break if opt == '-' && iptr == ARGV.size - 1

      opt = opt.gsub(/^-+/, '')

      unless okargs.key?(opt)
        @_unrecognised_options.push(opt)
        ARGV.delete_at(iptr)
        next
      end

      if @_count.key?(opt)
        @_count[opt] += 1
      else
        @_count[opt] = 1
        @_option_value[opt] = []
      end

      qualifiers = okargs[opt]

      ARGV.delete_at(iptr)

      next unless qualifiers.length.positive?

      if qualifiers == 's'
        tmp = ARGV.delete_at(iptr)
        @_option_value[opt].push(tmp)
        next
      end

      if %w[i int].include?(qualifiers)
        tmp = ARGV.delete_at(iptr)
        raise "Invalid integer qualifier for -#{opt} option '#{tmp}'" unless tmp =~ /^-*\d+$/

        @_option_value[opt].push(tmp.to_i)
        next
      end

      if %w[u uint].include?(qualifiers)
        tmp = ARGV.delete_at(iptr)
        raise "Invalid unsigned integer qualifier for -#{opt} option '#{tmp}'" unless tmp =~ /^\d+$/

        @_option_value[opt].push(tmp.to_i)
        next
      end

      if qualifiers == 'ipos'
        tmp = ARGV.delete_at(iptr)
        unless tmp =~ /^\d+$/ && tmp.to_i.positive?
          raise "Invalid positive integer qualifier for -#{opt} option '#{tmp}'"
        end

        @_option_value[opt].push(tmp.to_i)
        next
      end

      if %w[f float].include?(qualifiers)
        tmp = ARGV.delete_at(iptr)
        begin
          f = Float(tmp)
        rescue ArgumentError
          raise "Invalid float for '-#{opt}', '#{tmp}'\n"
        end
        @_option_value[opt].push(tmp.to_f)
        next
      end

      if qualifiers == 'fraction'
        tmp = ARGV.delete_at(iptr)
        begin
          f = Float(tmp)
        rescue ArgumentError
          raise "Invalid fraction for '-#{opt}', '#{tmp}'\n"
        end
        raise "Invalid fraction '#{f}'" unless f >= 0.0 && f <= 1.0

        @_option_value[opt].push(f)
        next
      end

      if qualifiers == 'sfile'
        tmp = ARGV.delete_at(iptr)
        raise "Must specify file name for option '-#{opt}'" unless tmp
        raise "Missing or empty file '#{tmp}'" unless FileTest.size?(tmp)

        @_option_value[opt].push(tmp)
        next
      end

      if qualifiers == 'xfile'
        tmp = ARGV.delete_at(iptr)
        raise "Must specify file name for option '-#{opt}'" unless tmp
        raise "Missing or empty file '#{tmp}'" unless FileTest.executable_real?(tmp)

        @_option_value[opt].push(tmp)
        next
      end

      if qualifiers == 'dir'
        tmp = ARGV.delete_at(iptr)
        raise "Must specify file name for option '-#{opt}'" unless tmp
        raise "Missing or invalid directory file '#{tmp}'" unless FileTest.directory?(tmp)

        @_option_value[opt].push(tmp)
        next
      end

      if qualifiers == 'close'
        gotclose = false
        closing_option = Regexp.new("^-+#{opt}$")
        tokens = []
        while iptr < ARGV.size
          tmp2 = ARGV.delete_at(iptr)
          if closing_option.match(tmp2)
            gotclose = true
            break
          end
          tokens.push(tmp2)
        end
        raise "No closing -#{opt}" unless gotclose

        @_option_value[opt].push(tokens.join(' '))
        next
      end

      $stderr.print "IWCmdline:initialize:invalid qualifier '#{qualifiers}' for option '#{opt}', value discarded\n"
      ARGV.delete_at(iptr) # just discard it???!
    end
  end

  def debug_print; end

  def unrecognised_options_encountered
    !@_unrecognised_options.empty?
  end

  def unrecognised_options
    @_unrecognised_options
  end

  def option_present(opt)
    @_count.key?(opt)
  end

  def option_count(opt)
    return 0 unless @_count.key?(opt)

    @_count[opt]
  end

  def value(opt, ndx = 0)
    return false unless @_option_value.key?(opt)

    tmp = @_option_value[opt]

    tmp[ndx]
  end

  def value_or_empty_string(opt)
    return '' unless @_option_value.key?(opt)

    tmp = @_option_value[opt]

    tmp[0]
  end

  def values(opt)
    # $stderr.print "Do we have '#{opt}' " << (@_option_value.key?(opt)).to_s << "\n"

    return [] unless @_option_value.key?(opt)

    @_option_value[opt]
  end

  def number_of_these_options_set(opt)
    rc = 0
    opt.each do |o|
      rc += 1 if @_option_value.key?(o)
    end

    rc
  end

  # We may decide to combine two options.

  def combine_options(o1, o2) # rubocop:disable Naming/MethodParameterName
    # $stderr.print "Combining '#{o1}' with '#{o2}'\n"
    return unless @_option_value.key?(o1) && @_option_value.key?(o2)

    if !@_option_value.key?(o1)
      @_option_value[o1] = []
    elsif !@_option_value.key?(o2)
      @_option_value[o2] = []
    end

    tmp = @_option_value[o1] | @_option_value[o2]

    @_option_value[o1] = tmp
    @_option_value[o2] = tmp
    @_count[o1] = tmp.size
    @_count[o2] = tmp.size
  end

  def values_as_array(opt, separator)
    return [] unless @_option_value.key?(opt)

    rc = []

    @_option_value[opt].each do |o|
      f = o.split(separator)
      f.each do |x|
        rc.push(x)
      end
    end

    rc
  end
end
