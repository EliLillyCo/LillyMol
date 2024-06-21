# frozen_string_literal: true

# IWFP class for gfp_make

require_relative 'lib/fp_common'

# Molecular properties fingerprint.
class MPR
  attr_reader :description

  def initialize
    @rx = Regexp.new('^MPR')
    @description = 'Molecular property fingerprint'
    @executable = 'temperature'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^MPR(\d+)*/.match(fp)
    raise "Unrecognized IW fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)
    # path_length, atype = FpCommon.parse_fp_token(fp[3..])

    cmd << ' -J MPR'
    cmd
  end
end
