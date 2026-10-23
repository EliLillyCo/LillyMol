# Ruby
Most ruby scripts included in LillyMol will run under all recent versions of ruby,
including 4.x.x versions. The most significant problems encountered are when ruby
tools use Protocol Buffers. Some years ago there was a breaking change
in how the google-protbuf package functioned, and that has caused problems for
LillyMol releases.

The solution is to pin things to more up to date versions of google-protobuf. In
the top level directory, LILLYMOL_HOME there is Gemfile that pins a functioning
version of google-protbuf.

Here is a more detailed explanation from Claude as it guided me through this.

1. The Gemfile
If LillyMol doesn't already have a Gemfile, add one at the project root:

```ruby
# Gemfile
source "https://rubygems.org"

gem "google-protobuf", ">= 3.21"
```

>= 3.21 is a safe floor — that's comfortably within the era that emits/reads the serialized-descriptor format,
well past the 3.18 compatibility line.

Then generate a lockfile once, from a machine with a good gem environment:

```bash
bundle install
```

I have done this with LillyMol.

This creates a Gemfile.lock that pins the exact resolved version. Commit Gemfile.lock to
the repo. This is the part that actually solves your cross-machine drift problem — anyone who runs
bundle install afterward gets the identical google-protobuf version you tested with,
not whatever happens to be latest on PyPI-equivalent (RubyGems) that day.

2. Running the program with Bundler

Instead of just ruby svmfp_make.rb, it becomes:

```bash
bundle exec ruby svmfp_make.rb
bundle exec makes sure the script sees the gem versions from Gemfile.lock, not whatever's globally
```
installed under that Ruby version. This is the mechanism that makes the pin actually enforceable
— without it, Bundler's pinning is just a suggestion.

3. Regenerating _pb.rb in the modern (serialized) format
On your Ubuntu machine, confirm your protoc is new enough:

```bash
protoc --version
```
Anything 23.x or newer will emit the serialized form by default. Then, from
the `${LILLYMOL_HOME}/src` directory, regenerate:

```bash
protoc --ruby_out=. Utilities/GFP_Tools/gfp_model.proto
cp Utilities/GFP_Tools/gfp_model_pb.rb ../contrib/bin/svmfp/lib/
```
Quick sanity check that you got the new format and not the old DSL — open the generated
file and look for add_serialized_file rather than add_file(...) do:

```bash
head -20 Utilities/GFP_Tools/gfp_model_pb.rb
```

4. Documentation note (for your build docs)

Something like:

Ruby protobuf bindings: _pb.rb files are pre-generated (via protoc --ruby_out) and committed to the repo,
since not all build hosts have protoc installed. The google-protobuf gem version is pinned
in Gemfile.lock — always run via bundle exec rather than bare ruby, so the correct gem version
is loaded. If you must regenerate _pb.rb files, use a protoc version ≥ 23.0 so the output
uses the serialized-descriptor format (add_serialized_file), which is compatible
with google-protobuf gem versions back to 3.18 — do not use an older protoc that emits the
legacy DSL (build do...end blocks), as this breaks under current gem versions with an undefined method 'build' error.

One gotcha worth flagging now

If rbenv + bundler isn't already part of your workflow, you'll need gem install bundler once per Ruby version/machine (it's not bundled by default in all Ruby installs). Worth checking:

bash
bundle -v
If that's missing anywhere in your machine fleet, that's a one-time gem install bundler before any of this works.

