var qrSolve = require('./')
var test        = require('tape')
var almostEqual = require("almost-equal")

// deterministic RNG for generating test data.
var rng = new require('xorshift').constructor([1, 0, 2, 0]);

var eps = 1e-5

function qrSolveHelper(M, b, m, n) {
  var solve = qrSolve.prepare(M, n, P)
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

  var solve = qrSolve.prepare(As, m, n)
  var foundSolution = solve(b)

  var expectedSolution = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

  for(var i=0; i< n; ++i) {
    t.assert(almostEqual(expectedSolution[i], foundSolution[i], eps, eps), "solution element " + i + ": "+ expectedSolution[i] + " = " + foundSolution[i])
  }

  t.end();
})
