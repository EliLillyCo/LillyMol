# frozen_string_literal: true

# RS class for gfp_make

require_relative 'lib/fp_common'

# Atom pair class.
class RS
  attr_reader :description

  def initialize
    @rx = Regexp.new('^RS')
    @description = 'Ring substitution fingerprints'
    @executable = 'ring_substitution'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^RS/.match(fp)
    raise "Unrecognized RS fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)
    cmd
  end
end
