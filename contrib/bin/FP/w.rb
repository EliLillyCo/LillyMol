# frozen_string_literal: true

# W class for gfp_make

require_relative 'lib/fp_common'

# Class for W topological torsion fingerprints.
class W
  attr_reader :description

  def initialize
    @rx = Regexp.new('^W')
    @description = 'iwdescr fingerprints'
    @executable = 'iwdescr'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  # The command can be
  # W\d+ or W\d+,f1,f2,f3,... where the 'f' are feature names, natoms...
  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    replicates = 1
    m = /^W(\d+)*/.match(fp)
    raise "Unrecognized TT fp form '#{fp}'" unless m
    fp.gsub!(/^W(\d+)*/, "")

    cmd = @executable.dup
    cmd << " -G FILTER" unless first_in_pipeline
    cmd << " -G #{m[1]}" if m[1].to_i > 1
    if fp.empty?
      cmd << " -G ALL"
    else
      cmd << ' -G ' << fp[1..-1]
    end
    cmd << ' -G R=50'
    # path_length, atype, fixed = FpCommon.parse_fp_token(fp[2..])

    cmd
  end
end
