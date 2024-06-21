# frozen_string_literal: true

# MaccsKeys class for gfp_make

require_relative 'lib/fp_common'

# Class for MACCS keys.
class MK
  attr_reader :description

  def initialize
    @rx = Regexp.new('^MK')
    @description = 'MACCS keys fingerprint'
    @executable = 'maccskeys'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^MK(\d+)*/.match(fp)
    raise "Unrecognized MK fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                        extra_qualifiers: extra_qualifiers)
    if fp.match(/:sparse/) 
      cmd << " -J NC=NCMK"
    else
      cmd << " -J FPMK"
    end
    cmd
  end
end
