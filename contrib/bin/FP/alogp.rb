# frozen_string_literal: true

# Alogp class for gfp_make

require_relative 'lib/fp_common'

# Class for Alogp fingerprints.
class ALOGP
  attr_reader :description

  def initialize
    @rx = Regexp.new('^ALOGP')
    @description = 'AlogP fingerprint'
    @executable = 'alogp'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^ALOGP(\d+)*/.match(fp)
    raise "Unrecognized ALOGP fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)
    replicates, atype, fixed = FpCommon.parse_fp_token(fp[5..])

    # fixed is the default, and non colliding does not work with iwfp.
    cmd << ' -J NCALOGP -Y alcacid -Y RDKIT.N+ -Y quiet'
    cmd << " -p #{replicates}" if replicates
    cmd
  end
end
