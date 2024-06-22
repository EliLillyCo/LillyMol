# SCAF class for gfp_make

require_relative 'lib/fp_common'

# Scaffold class.
class SCAF
  attr_reader :description

  def initialize
    @rx = Regexp.new('^SCAFI*')
    @description = 'Scaffold fingerprints'
    @executable = 'molecular_abstractions'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^SCAF/.match(fp)
    raise "Unrecognized SCAF fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)

    if fp.include?('I')
      path_length, atype, fixed = FpCommon.parse_fp_token(fp[5..])
      directive = 'ISO FP'
    else
      path_length, atype, fixed = FpCommon.parse_fp_token(fp[4..])
      directive = 'FP'
      cmd << " -a 'scaffold(FP)'"
    end

    if atype
      directive << " ATYPE=#{atype}"
    end

    cmd << " -a 'scaffold(#{directive})'"

    cmd
  end
end
