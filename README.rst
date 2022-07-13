kdtree
======

.. image:: http://nuclear.mutantstargoat.com/sw/kdtree/img/kdtree_logo.png

Overview
--------
kdtree is a simple, easy to use C library for working with kd-trees.

Kd-trees are an extension of binary search trees to k-dimensional data. They
facilitate very fast searching, and nearest-neighbor queries.

This particular implementation is designed to be efficient and very easy to
use. It is completely written in ANSI/ISO C, and thus completely
cross-platform. 

See under the ``doc/`` and ``examples/`` directories to find out how to use the
kdtree library.

Fork
----

This fork of the library has been modified to:
- Choose between using floats & doubles for coordinates with a compiler flag
- Run faster, as a result of several changes to...
- Use less memory
- Support application-provided memory buffers for most operations [0]
- Add an optional filter to skip nodes for with non-zero coordinates if the
  corresponding input coordinates have a value of zero.

[0] The only operations that still require internal allocations are the
nearest node queries that return a result set.  If you only care about
searching for the single nearest node, then you can use the library without
any internal allocations.

License
-------
Author: John Tsiombikas <nuclear@member.fsf.org>

kdtree is free software. You may use, modify, and redistribute it under the
terms of the 3-clause BSD license.

Download
--------
Latest release (0.5.6): http://nuclear.mutantstargoat.com/sw/kdtree/files/kdtree-0.5.6.tar.gz

You can find previous releases here:
http://nuclear.mutantstargoat.com/sw/kdtree/files/

You can also grab a copy of the source from github: https://github.com/jtsiomb/kdtree
