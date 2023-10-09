# Why Protocol Buffers

Protocol Buffers are widely used at Google, and offer many advantages

* Cross language API (C++, Python, Ruby, Julia, Go ...)
* Can be serialized to binary form very efficiently

While LillyMol seldom uses serialized protos, the language bindings
are a compelling advantage for us.

Other config languages could have been used for LillyMol, JSon, XML,
yaml, etc, but protos offer significant benefits over those others.
Note that protos do inter-operate with JSon, and it would be possible
to do most of what we have done using textproto format with JSon.

For a general introduction to Protocol Buffers refer to the documentation
published by Google [Protocol Buffers](https://protobuf.dev/) and related
sites.

