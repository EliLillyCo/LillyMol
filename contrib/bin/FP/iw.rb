# frozen_string_literal: true

# IWFP class for gfp_make

require_relative 'lib/fp_common'

# Class for IW Linear fingerprints.
class IW
  attr_reader :description

  def initialize
    @rx = Regexp.new('^IW')
    @description = 'Linear path fingerprint'
    @executable = 'iwfp'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^IW(\d+)*/.match(fp)
    raise "Unrecognized IW fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)
    path_length, atype, fixed = FpCommon.parse_fp_token(fp[2..])

    # fixed is the default, and non colliding does not work with iwfp.
    cmd << ' -J FPIW'
    cmd << "#{path_length} -R #{path_length}" if path_length
    cmd << " -P #{atype}" if atype
    cmd
  end
end
