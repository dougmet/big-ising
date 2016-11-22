# big-ising

> Make a huge, really huge, Ising model at Tc

It is mostly a standard Wolff algorithm, but it does two clever things:

1. Stores spins bit-wise. So 8 spins per byte of memory.
2. The stack wraps around. So when your cluster gets huge (so long as you're clearing it quick enough) the stack starts overwriting itself.

## Building

There are no big dependencies for this code. It's been tested on Mac with clang, Linux with gcc and Windows using MingGW-W64. You can get the Windows tool chain from [https://mingw-w64.org/](https://mingw-w64.org/). If you have the [Rtools](https://cran.r-project.org/bin/windows/Rtools/) tools set then it comes bundled with that. I used rtools34.exe.

## License

MIT © Douglas Ashton, Lode Vandevenne
