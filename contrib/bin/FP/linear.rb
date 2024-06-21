# Linear_fingerprint
# frozen_string_literal: true

# LINEAR class for gfp_make

require_relative 'lib/fp_common'

# Class for Linear fingerprints.
class LINEAR
  attr_reader :description

  def initialize
    @rx = Regexp.new('^LINEAR')
    @description = 'Linear fingerprint (faster than IW)'
    @executable = 'linear_fingerprint'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^LINEAR(\d+)*/.match(fp)
    raise "Unrecognized LINEAR fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)
    path_length, atype, fixed = FpCommon.parse_fp_token(fp[6..])

    if fixed
      cmd << ' -J fixed -J FPLN'
    else
      cmd << ' -J NCLN'
    end
    cmd << "#{path_length} -R #{path_length}" if path_length
    cmd << " -P #{atype}" if atype
    cmd
  end
end
