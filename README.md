# HSSP

This program implements a branch and bound algorithm for the Hypervolume Subset Selection Problem (HSSP) for an arbitrary number of objectives. This implementation is a support material for the article: Implicit enumeration strategies for the hypervolume subset selection problem
R. Gomes, A. Guerreiro, T. Kuhn, L. Paquete Computers & Operations Research, 2018 (in press).

Building
--

```
cd src/
make
```

Usage
--

The program reads a set of points provided by a filename in the command line:
```
./hssp data
```

In the input file, each point is given in a separate line, and each coordinate within a line is separated by a whitespace.

The reference point can be given by the option `-r`.
```
./hssp -r "1 1 1" data
```
For the remainder options available, check the output of `hssp --help`.
