# tclfastme

Tcl wrapper around [FastME](http://www.atgc-montpellier.fr/fastme/) distance-based phylogeny
(balanced/OLS minimum evolution, NNI, BIONJ, NJ; Desper & Gascuel 2002).
Package `tclfastme` 0.1. It builds against Tcl 8.6 and Tcl 9 from the same sources.

```tcl
package require tclfastme
fastme namesList distList ?-init_tree bal|ols|bionj|nj? ?-branch_lengths bal|ols|none? ?-swap bal|ols|none?
```

`distList` is the full n×n distance matrix, flattened row by row, in the order of `namesList`.
Names must be unique. The command returns the tree as a Newick string.

## Build

TEA (autoconf) build, out of the source tree:

```sh
mkdir build-tcl8 && cd build-tcl8
../configure --with-tcl=/path/to/tcl8.6/lib && make     # -> libtclfastme0.1.so
mkdir build-tcl9 && cd build-tcl9
../configure --with-tcl=/path/to/tcl9/lib && make       # -> libtcl9tclfastme0.1.so
make install
```

On Windows, build with MinGW/MSYS2 the same way.
`configure` is already generated; run `autoconf` only after editing `configure.ac`.

## Tests

`generic/test.tcl` is a tcltest suite (51 tests). Run it from a build directory:

```sh
cd build-tcl9 && tclsh9.0 ../generic/test.tcl
```

## License

GPL-3.0-or-later, Copyright (c) Luc Moulinier. See [LICENSE](LICENSE).
The phylogeny code comes from FastME (Desper & Gascuel).
