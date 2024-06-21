# frozen_string_literal: true

# Xlogp class for gfp_make

require_relative 'lib/fp_common'

# Class for Xlogp fingerprints.
class XLOGP
  attr_reader :description

  def initialize
    @rx = Regexp.new('^XLOGP')
    @description = 'XlogP fingerprint'
    @executable = 'xlogp'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^XLOGP(\d+)*/.match(fp)
    raise "Unrecognized XLOGP fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)
    replicates, atype, fixed = FpCommon.parse_fp_token(fp[5..])

    # fixed is the default, and non colliding does not work with iwfp.
    cmd << ' -J NCXLOGP'
    cmd << " -p #{replicates}" if replicates
    cmd
  end
end
