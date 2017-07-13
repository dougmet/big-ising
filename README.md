# big-ising

> Make a huge, really huge, Ising model at Tc

It is mostly a standard Wolff algorithm, but it does two clever things:

1. Stores spins bit-wise. So 8 spins per byte of memory.
2. The stack wraps around. So when your cluster gets huge (so long as you're clearing it quick enough) the stack starts overwriting itself. This allows the stack size to scale like $\sqrt{N}$ rather than $N$.

## Building

There are no big dependencies for this code. It's been tested on Mac with clang, Linux with gcc and Windows using MingGW-W64. You can get the Windows tool chain from [https://mingw-w64.org/](https://mingw-w64.org/). If you have the [Rtools](https://cran.r-project.org/bin/windows/Rtools/) tools set then it comes bundled with that. I used rtools34.exe.

I can't get the really really big models to work on Windows, only had success on Linux.

## Running

There are no input files. All variables are compiled in as `#define`s. If the `lattice.pos` file is present in the working directory then this will be loaded in and the simulation resumed from this point.

## Outputs

- `data` contains the magnetisation and energy for each MC sweep
- `lattice.pos` is a binary dump of the current spin configuration. This can be used for reloading or for making images.
- `swolf*.png` these files show the evolution of the configuration.

## Making Pngs

A big enough Ising model will be too big to plot. Therefore you need to block spins together to draw pixels. The `isingpng` programme is designed to help make pngs by choosing location and zoom level. The png library is from [Lode Vandevenne](http://lodev.org/lodepng/), please see its [license](https://github.com/lvandeve/lodepng/blob/master/lodepng.h) before distributing.

## License

MIT © 2016 Douglas Ashton, 2005-2016 Lode Vandevenne (lodepng)
