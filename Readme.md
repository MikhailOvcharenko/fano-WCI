# The classification of smooth well-formed Fano weighted complete intersections

Author: Mikhail Ovcharenko.

We construct the classification of smooth well-formed Fano weighted complete intersections of any given variance, following ``The classification of smooth well-formed Fano weighted complete intersections'' (see [arXiv:2006.05666](https://arxiv.org/abs/2006.05666)).

For any details regarding this code, contact me at <ovcharenko@mi-ras.ru>.

The code runs using Sage 9.5 or more recent versions.

## How to run

The following command prints candidate generating families of smooth well-formed Fano weighted complete intersections of given variance r.

```
sage : load('fano-WCI.sage')
sage : bruteforce_WCI(r)
```

## Example for r = 2

```
[(1, 1, 1, 1, 1), (4,)]
[(1, 1, 1, 1, 3), (6,)]
[(1, 1, 1, 1, 1, 1, 1), (3, 3)]
[(1, 1, 1, 1, 1, 2, 2), (4, 4)]
[(1, 1, 1, 1, 1, 1, 2), (3, 4)]
[(1, 1, 1, 1, 2, 2, 3), (4, 6)]
[(1, 1, 1, 2, 2, 3, 3), (6, 6)]
```