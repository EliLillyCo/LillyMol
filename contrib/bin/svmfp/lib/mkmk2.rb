# frozen_string_literal: true

# Part of gfp_make.
# If both -MK and -MK2 are specified, as defaults, consoldiate those entries.

def consoliate_mkmk2(args) # rubocop:disable Metrics/CyclomaticComplexity, Metrics/MethodLength
  mk_index = -1
  mk2_index = -1
  args.each_with_index do |fp, i|
    case fp.option
    when 'MK'
      return args if fp.value

      mk_index = i
    when 'MK2'
      return args if fp.value

      mk2_index = i
    end
  end

  # Nothing specified.
  return args if mk_index.negative? && mk2_index.negative?

  # Just -MK
  return args if mk2_index.negative?

  # If just MK2, must get MK as well.
  if mk_index.negative? && mk2_index >= 0
    args[mk2_index] = OptionValue.new('MK', '-J FPMK -J LEVEL2=FPMK2')
    return args
  end

  # Both specified

  args[mk_index] = OptionValue.new('MK', '-J FPMK -J LEVEL2=FPMK2')
  args.delete_at(mk2_index)
  args
end
