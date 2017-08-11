var qrSolve = require('./')
var test        = require('tape')
var almostEqual = require("almost-equal")

// deterministic RNG for generating test data.
var rng = new require('xorshift').constructor([1, 0, 2, 0]);

var eps = 1e-5

function qrSolveHelper(As, m, n, b) {
  var solve = qrSolve.prepare(As, m, n)
  return solve(b)
}

test('solve10x10matrix', function(t) {

  var As = [
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

  var foundSolution = qrSolveHelper(As, m, n, b)

  var expectedSolution = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

  for(var i=0; i< n; ++i) {
    t.assert(almostEqual(expectedSolution[i], foundSolution[i], eps, eps), "solution element " + i + ": "+ expectedSolution[i] + " = " + foundSolution[i])
  }

  t.end();
})

test('solve100x100matrix', function(t) {
  var L = []
  var n = 1000

  // first generate a lower-triangular matrix L
  for(var i = 0; i < n; ++i) {
    L[i] = []

    for(var j = 0; j <n; ++j) {
      L[i][j] = 0
    }

    for(var j = 0; j <= i; ++j) {
      if(rng.random() > 0.99 || i === j) {
        L[i][j] = Math.floor(rng.random() * 10)+1
      } else {
        L[i][j] = 0
      }
    }
  }

  // next, we simply multiply L with its transpose, and put the result in A.
  // the resulting matrix is symmetric, and positive definite,
  // thus it must have a cholesky decomposition.
  var A = []
  for(var i = 0; i < n; ++i) {
    A[i] = []
  }
  for(var i = 0; i < n; ++i) {
    for(var j = 0; j < n; ++j) {
      var s = 0.0
      for(var k = 0; k < n; ++k) {
        s += L[i][k] * L[j][k]
      }
      A[i][j] = s
    }
  }

  // now store A as a sparse matrix M.
  var M = []
  for(var row = 0; row < n; ++row) {
    for(var col = 0; col < n; ++col) {
      if(A[row][col] > 0.0001) {
        M.push([row, col, A[row][col]])
      }
    }
  }

  // In our test, we shall solve the equation
  // Mx = b
  // so randomly generate x.
  var x = []
  var b = []
  for(var i = 0; i < n; ++i) {
    x[i] = Math.floor(rng.random() * 9)
  }

  // Now compute b as b = Mx
  for(var i = 0; i < n; ++i) {
    var s = 0.0
    for(var k = 0; k < n; ++k) {
      s += A[i][k] * x[k]
    }
    b[i] = s
  }

  var foundSolution = qrSolveHelper(M, n, n, b)

  // check that the residual vector is 0.
  for(var i = 0; i < n; ++i) {
    var s = 0.0
    for(var k = 0; k < n; ++k) {
      s += A[i][k] * foundSolution[k]
    }
    var res = b[i] - s
    t.assert(almostEqual(0.0, res, eps, eps), "residual #" + i + ":" +  "0.0 = " + res)
  }

  t.end();
})
