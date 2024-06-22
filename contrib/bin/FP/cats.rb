# CATS class for gfp_make

require_relative 'lib/fp_common'

# Class for CATS fingerprints.
class CATS
  attr_reader :description

  def initialize
    @rx = Regexp.new('^CATSP*')
    @description = 'CATS pharmacaphore fingerprint'
    @executable = 'jwcats.sh'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^CATSP*(\d+)*/.match(fp)
    raise "Unrecognized CATS fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)
    tag = "NCCATS"
    if /^CATSP(\d+)*/.match(fp)
      path_length, atype, fixed = FpCommon.parse_fp_token(fp[5..])
      tag << 'P'
      cmd << ' -p'
    else
      path_length, atype, fixed = FpCommon.parse_fp_token(fp[4..])
    end

    if path_length
      tag << path_length.to_s
    end

    cmd << " -J #{tag}"
    cmd << " -m #{path_length}" if path_length

    cmd
  end
end
