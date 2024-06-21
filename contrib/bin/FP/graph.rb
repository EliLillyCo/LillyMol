# GRAPH class for gfp_make

require_relative 'lib/fp_common'

# Scaffold class.
class GRAPH
  attr_reader :description

  def initialize
    @rx = Regexp.new('^GRAPH$')
    @description = 'Molecular graph fingerprints'
    @executable = 'molecular_abstractions'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^GRAPH$/.match(fp)
    raise "Unrecognized GRAPH fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)

    cmd << " -Y nbits=1024 -a 'graph(FP)'"

    cmd
  end
end
