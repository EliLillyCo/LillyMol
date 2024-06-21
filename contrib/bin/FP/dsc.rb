# frozen_string_literal: true

# iwdescr fingerprint class for gfp_make

require_relative 'lib/fp_common'

# Class for descriptor based fingerprints.
class DSC
  attr_reader :description

  def initialize
    @rx = Regexp.new('^DSC')
    @description = 'Descriptor based fingerprint'
    @executable = 'iwdescr'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^DSC(\d+)*/.match(fp)
    raise "Unrecognized DSC fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)
    # iwdescr already has a -f option that does something different.
    cmd = cmd.gsub(/-f/, '-G FILTER') unless first_in_pipeline
    cmd << ' -O dm -O complex'
    # TODO: ianwatson Figure out donor/acceptor things...
    cmd
  end
end
