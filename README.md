# qr-solve[WIP]

This module solves sparse linear systems, by finding the
QR-decomposition, and then using it to solve the system. It is
basically a Javascript port of
[Quern](https://www.cs.ubc.ca/~rbridson/quern/).

## Install

    npm install qr-solve

## Example

```javascript
var qrSolve = require('qr-solve')

/*
Below we specify the sparse matrix:

1.7  0    0    0    0    0    0    0    0.13 0
0    1.0  0    0    0.02 0    0    0    0    0.01
0    0    1.5  0    0    0    0    0    0    0
0    0    0    1.1  0    0    0    0    0    0
0    0.02 0    0    2.6  0    0.16 0.09 0.52 0.53
0    0    0    0    0    1.2  0    0    0    0
0    0    0    0    0.16 0    1.3  0    0    0.56
0    0    0    0    0.09 0    0    1.6  0.11 0
0.13 0    0    0    0.52 0    0    0.11 1.4  0
0    0.01 0    0    0.53 0    0.56 0    0    3.1
  */
var A = [
  [0, 0, +1.70],
  [0, 8, +0.13],

  [1, 1, +1.00],
  [1, 4, +0.02],
  [1, 9, +0.01],

  [2, 2, +1.50],

  [3, 3, +1.10],

  [4, 1, +0.02],
  [4, 4, +2.60],
  [4, 6, +0.16],
  [4, 7, +0.09],
  [4, 8, +0.52],
  [4, 9, +0.53],

  [5, 5, +1.20],

  [6, 4, +0.16],
  [6, 6, +1.30],
  [6, 9, +0.56],

  [7, 4, +0.09],
  [7, 7, +1.60],
  [7, 8, +0.11],

  [8, 0, +0.13],
  [8, 4, +0.52],
  [8, 7, +0.11],
  [8, 8, +1.40],

  [9, 1, +0.01],
  [9, 4, +0.53],
  [9, 6, +0.56],
  [9, 9, +3.10],
]

var b =  [0.287, 0.22, 0.45, 0.44, 2.486, 0.72, 1.55, 1.424, 1.621, 3.759]

var m = 10
var n = 10

var solve = qrSolve.prepare(A, m, n) // first decompose.

var solution = new Float64Array(n) // in here we put the solution
solve(b, solution)
console.log(solution) // then solve!
```

## API

### `require("qr-solve").prepare(A, m, n)`

Decomposes `A` into its QR-decomposition. A function is returned that
can be used to solve the equation `Ax = b`, for some given value of
`b`.

* `A` a list of the matrix coefficients of the sparse matrix
`A`. These MUST be given in the order of increasing rows and
columns numbers.
* `m` the number of rows in the matrix `A`
* `n` the number of columns in the matrix `A`

**Returns** A function that takes a first argument `b`, and a second
  argument `x`. The function
  put the solution to the equation `Ax = b` in the array `x`.

**NOTE** the module does no sanity checking on the input arguments. It is assumed that the user knows what he/she is doing!
