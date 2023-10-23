# Licensed Software

This directory contains LillyMol interfaces to some
licensed software.

If you have a license, adjust WORKSPACE to point to your
installation, and then you should be able to do
```
bazelisk build -c opt <other bazel options> Vendor:all
```
