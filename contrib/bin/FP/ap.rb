# frozen_string_literal: true

# AP class for gfp_make

require_relative 'lib/fp_common'

# Atom pair class.
class AP
  attr_reader :description

  def initialize
    @rx = Regexp.new('^M*AP')
    @description = 'Atom pair fingerprints'
    @executable = 'atom_pair_fingerprint'
  end

  def match?(fp) # rubocop:disable Naming/MethodParameterName
    @rx.match?(fp)
  end

  def expand(fp, first_in_pipeline:, extra_qualifiers:) # rubocop:disable Naming/MethodParameterName
    m = /^M*AP(\d+)*/.match(fp)
    raise "Unrecognized AP fp form '#{fp}'" unless m

    cmd = FpCommon.initial_command_stem(@executable, first_in_pipeline: first_in_pipeline,
                                                     extra_qualifiers: extra_qualifiers)
    path_length, atype, fixed = FpCommon.parse_fp_token(fp.gsub(/^M*AP/, ''))
    $stderr << "path_length #{path_length} atype #{atype} fixed #{fixed}\n"

    if fp.match(/:fixed/)
      cmd << ' -J fixed -J FPAP'
    else
      cmd << ' -J NCAP'
    end
    cmd << "#{path_length} -R #{path_length}" if path_length
    cmd << " -P #{atype}" if atype
    cmd
  end
end
