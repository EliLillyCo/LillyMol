# frozen_string_literal: true

OptionValue = Struct.new(:option, :value) do
end

# Some other ways of doing the same thing.
# class OptionValue < Struct.new(:option, :value)
# end
# Struct.new("OptionValue", :option, :value)

# Support functions for processing the arguments to gfp_make.
module GfpMakeSupport
  # Argument parsing for gfp_make.
  # args is an array of options.
  # Some will be single options, -MK, others will be of the
  # form -oMK ... -cMK
  # Return an Array of Structs, where each struct constists of
  # the key and value.
  # -MK -IW -oMK x y z -cMK -QQ is returned as
  # [ [-MK, nil],
  #   [-IW, nil],
  #   [-MK, [x, y, z]],
  #   [-QQ, nil]]

  def self.group_args(args) # rubocop:disable Metrics/MethodLength, Metrics/AbcSize
    return nil if args.empty?

    to_be_returned = []
    argptr = 0
    while argptr < args.size
      arg = args[argptr]
      argptr += 1
      unless arg.start_with?('-')
        $stderr << "group_args:does not look like option #{arg}\n"
        return nil
      end
      if arg.start_with?('-c')
        $stderr << "group_args:unpaired closing option #{arg}\n"
        return nil
      end
      unless arg.start_with?('-o')
        to_be_returned.append(OptionValue.new(arg[1..-1], nil))
        next
      end
      # Find the corresponding closing -c flag.
      expected = arg.sub(/^-o/, '-c')
      found = args[argptr..args.size].find_index(expected)
      unless found
        $stderr << "group_args: no #{expected} found\n"
        return nil
      end
      to_be_returned.append(OptionValue.new(arg[2..-1], args[argptr...(argptr + found)].join(' ')))
      argptr += found + 1
    end

    to_be_returned
  end

  # Return a Hash of the unique options in `fp_args`.
  # Return nil if there are duplicate options present
  def self.all_options_unique(fp_args)
    options = {}
    fp_args.each do |fp|
      if options.key?(fp.option)
        $stderr << "Duplicate fingerprint #{fp}\n"
        return nil
      end
      options[fp.option] = fp.value
    end
    options
  end
end
