# frozen_string_literal: true

# Part of gfp_make

# Report to stderr the recognized fingerprints in `fps`.
def report_recognized_fingerprints(fps)
  $stderr << "#{fps.length} fingerprints recognized\n"
  fps.each do |k, v|
    $stderr << "#{k} #{v} #{v.description}\n"
  end
end

# Return a map of fingerprint names to object that can process
# that kind of fingerprint.
# By convention, the name of the class is the uppercase of
# the config file name.
# config_dirs is an array of directories in which to look.
# If multiple items are found, the last will be silently used.

def config_fingerprints(config_dirs, verbose) # rubocop:disable Metrics/AbcSize, Metrics/CyclomaticComplexity, Metrics/MethodLength, Metrics/PerceivedComplexity
  fps = {} # To be returned.
  config_dirs.each do |dir|
    $stderr << "Fetching configs from #{dir}\n" if verbose.positive?
    Dir.entries(dir).each do |fname|
      next if /^\./.match(fname)
      next if fname == 'fp_common.rb'
      next if fname == 'lib'

      m = /^(\S+)\.rb/.match(fname)
      unless m
        $stderr << "No match file name #{fname}, ignored\n"
        next
      end
      require("#{dir}/#{fname}")
      class_name = m[1].upcase # By convention.
      fps[class_name] = eval("#{class_name}.new", binding, __FILE__, __LINE__)  # rubocop:disable Security/Eval
    end
  end

  report_recognized_fingerprints(fps) if verbose.positive?

  fps
end
