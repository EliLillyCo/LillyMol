# Copied from https://www.grahambrooks.com/software-development/2021/05/31/bazel-local-install.html


def _local_install_impl(ctx):
    target = ctx.attr.target
    shell_commands = ""
    env = ctx.configuration.default_shell_env

    # Enable via --action_env=BINDIR=/path/to/somewhere
    # but beware, changing it forces a recompile.
    if 'BINDIR' in env:
      target = env['BINDIR']

    for s in ctx.files.srcs:
        shell_commands += "dest=%s/$(basename %s)\n" % (target, s.short_path)
        shell_commands += "if [[ ! -s ${dest} || ${dest} -ot %s ]] ; then\n" % (s.short_path)

        shell_commands += "  echo 'Copy %s to %s'\n" % (s.short_path, target)
        shell_commands += "  cp -f %s %s\nfi\n" % (s.short_path, target)

    ctx.actions.write( 
        output = ctx.outputs.executable,
        is_executable = True,
        content = shell_commands,
    )
    runfiles = ctx.runfiles(files = ctx.files.srcs)
    return DefaultInfo(
        executable = ctx.outputs.executable,
        runfiles = runfiles,
    )

local_install = rule(
    executable = True,
    implementation = _local_install_impl,
    attrs = {
        "srcs": attr.label_list(allow_files = True),
        "target": attr.string(default = "/lrlhps/users/rx87690/LillyMolPrivate/bin/Linux", doc = "Installation target directory"),
    },
)
