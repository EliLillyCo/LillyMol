# frozen_string_literal: true

# tsubstructure derived fingerprint class for gfp_make

require_relative 'lib/fp_common'

# Class for tsubstructure based fingerprints.
class TSUB
  attr_reader :description

  def initialize
    @rx = Regexp.new('^TSUB')
    @description = 'tsubstructure based fingerprint'
    @executable = 'tsubstructure'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^TSUB(\d+)*/.match(fp)
    raise "Unrecognized TSUB fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)
    cmd << " -M bitrep=#{m[2]}" if m[2]

    # If of the form TSUB:xxx, 'xxx' is assumed to be a file of smarts.
    tokens = fp.split(/[:=]/)
    raise 'No queries' unless tokens.length == 2 || extra_qualifiers

    cmd << " -q S:#{tokens[1]}" if tokens.length == 2
    cmd << ' -J FPSUB'

    cmd
  end
end
