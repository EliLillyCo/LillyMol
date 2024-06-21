# frozen_string_literal: true

# SKE class for gfp_make

require_relative 'lib/fp_common'

# Skeleton class.
class SKE
  attr_reader :description

  def initialize
    @rx = Regexp.new('^SKE')
    @description = 'Skeleton fingerprints'
    @executable = 'molecular_abstractions'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^SKE/.match(fp)
    raise "Unrecognized SKE fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)

    cmd << " -c -a 'allbonds.allatoms(S FP=FPSKE)'"
    cmd
  end
end
