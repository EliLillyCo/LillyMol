# Local git-hooks

## Compile Time
We wish to include the date compiled with executables. Historically
this could be done by referencing the preprocessor symbol __DATE__
and __TIME__. But bazel deliberately redacts these - in order to avoid
needless recompliation.

Theoretically it is possible to suppress this, perhaps by adding
```
--features=-determinism
```
but I could not get that to work. Probably a misunderstanding in my
behalf...

As an alternative we create a file, `Foundational/iwmisc/compile_time.cc`
which has today's date built in. The `pre-commit` hook in this directory
checks the date currently in that file, and if it is today, does nothing.

If it is out of date, it will update it to today's date. You will then
need to recompile, in order to get today's date included in executables.
That will naturally happen on the next compilation.

But note that if it updates the file to today's date, the repo is no
longer clean, and you should add/commit the updated compile_time.cc file.

To get started, 
**YOU NEED TO COPY `pre-commit` TO  .git.hooks** in order to enable
this in your repo. I think there might be ways of automating this
but not sure - open to ideas.

## Other approaches.
Previously the git hash was passed into the executable via a shell
variable. But then even a minor commit triggers a recompilation of
everything - because a shell variable has changed.

Now perhaps a fresh compile at the start of every day is not a bad
thing, if there were a TODAY variable.  But recompiling everything
after a minor commit is burdensome.

## Summary
It is an open question how to best accomplish the task of having
the git hash and the compile time available in executables.

In all cases, we face a `lag` problem. For example once you have
got code working, those binaries have the old date and git hash
compiled in. You
then need to commit the current work in progress, which generates a
new git hash and date, and then recompile.  But it does mean that
`compile_time.cc` is always 'dirty' - because the previous commit
changed it.  And checking a file with an embedded date into the repo
seems wrong.

Maybe something like
[stack overflow](https://stackoverflow.com/questions/1792838/how-do-i-enable-the-ident-string-for-a-git-repository/1792913#1792913)
but I also want the date of compilation. That may offer some options.

There might also be bazel based solutions
[stack overflow](https://stackoverflow.com/questions/66647916/how-to-get-git-commit-sha-inside-build-bazel)
