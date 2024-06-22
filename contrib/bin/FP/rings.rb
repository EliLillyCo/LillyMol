# frozen_string_literal: true

# RINGS class for gfp_make

require_relative 'lib/fp_common'

# Ring class.
class RINGS
  attr_reader :description

  def initialize
    @rx = Regexp.new('^RINGSI*')
    @description = 'Rings fingerprint'
    @executable = 'molecular_abstractions'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^RING(I)*/.match(fp)
    raise "Unrecognized RING fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)

    cmd << ' -c -Y nbits=1024'
    if m[1] == 'I'
      cmd << " -a 'rings(ISO FP)'"
    else
      cmd << " -a 'rings(FP)'"
    end
    cmd
  end
end
