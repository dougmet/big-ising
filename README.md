# big-ising

Some old code I wrote that did a huge, really huge, Ising model

Come back when I've got it working again.

It is mostly a standard Wolff algorithm, but it does two clever things:

1. Stores spins bit-wise. So 8 spins per byte of memory.
2. The stack wraps around. So when your cluster gets huge (so long as you're clearing it quick enough) the stack starts overwriting itself.

