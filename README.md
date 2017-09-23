GGH15-Based Branching-Program Obfuscation
=========================================

An implementation of the [GGH15][1] Graded encoding scheme, and its
use for obfuscating read-once branching programs (BPs). A description
of many aspect of this implementation can be found [here][2].

To be able to handle anything more than just toy problems, we
developed a host of algorithmic and code-level optimizations. These include
new variants of discrete Gaussian sampler and lattice trapdoor sampler,
efficient matrix-manipulation routines, and many tradeoffs. Some of these
optimizations can likely find other uses in lattice-based cryptography.

This implementation is written in C++ and uses the [NTL mathematical library][3]
(version 10.4.0 or higher). It is distributed under the terms of the
[Apache License v2.0][4].

Content
-------
The top-level directory constains most of the implementations code,
and a Makefile to compile it. After cloning the repository, just run
`make depend` and then `make` in the top-level directory.

The `utils` subdirectory contains some utility functions (many of
which were lifted from [HElib][5] or [NTL][3]). The `programs`
subdirectory contains the obfuscation programs and some unit-test
programs. To build these programs, run `make <progName>_x` in the
top-level directory (e.g., `make initialize_x` to make the
initialization program).

Contributors
------------
This implementation was written by [Shai Halevi][6] from IBM,
[Tzipora Halevi][7] from Brooklyn College, and [Victor Shoup][8] and
[Noah Stephens-Davidowitz][9] from NYU. It was Supported by the
Defense Advanced Research Projects Agency (DARPA) and Army Research
Office(ARO) under Contract No. W911NF-15-C-0236.
 
  [1]: https://eprint.iacr.org/2014/645      "GGH15"
  [2]: https://eprint.iacr.org/2017/104      "HHSS17"
  [3]: http://www.shoup.net/ntl/             "NTL"
  [4]: http://www.apache.org/licenses/LICENSE-2.0  "Apache-v2.0"
  [5]: https://github.com/shaih/HElib        "HElib"
  [6]: https://shaih.github.io/              "Shai"
  [7]: http://www.brooklyn.cuny.edu/web/academics/faculty/faculty_profile.jsp?faculty=1333 "Tzipora"
  [8]: http://www.shoup.net/                 "Victor"
  [9]: http://www.noahsd.com/                "Noah"
