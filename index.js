const QUERN_INPUT_ERROR = 1

function SparseEntry(/*int*/ index_ = 0, /* double*/ value_ = 0.0) {
  this.index = index_
  this.value = value_
}

function ListIterator(list_) {
  this.i = 0
  this.list = list_
}

ListIterator.prototype.inc = function() {
  ++this.i
}

ListIterator.prototype.arrow = function() {
  return this.list.data[this.i]
}

ListIterator.prototype.still_going = function() {
  return this.i != this.list.data.length;
}

function SparseVector() {
  this.data = []
}

SparseVector.prototype.begin = function(){
  return new ListIterator(this);
}

SparseVector.prototype.empty = function(){
  return this.data.length == 0
}

SparseVector.prototype.front = function(){
  return this.data[0]
}

SparseVector.prototype.swap = function(other){
  var temp =  this.data;
  this.data = other.data;
  other.data = temp
}

SparseVector.prototype.size = function()
{
  return this.data.length;
}

SparseVector.prototype.insert = function(/* ListIterator*/p, e)
{
  this.data.splice(p.i, 0, e)
}

SparseVector.prototype.erase = function(/*ListIterator<T>&*/ p) {
  this.data.splice(p.i, 1)
}

SparseVector.prototype.pop_front = function() {
  this.data.splice(0, 1)
}

SparseVector.prototype.push_front = function(e){
  this.data.splice(0, 0, e)
}

function encode(c, s)
{
  return [c, s]
}

function decode(tau) {
  return tau
}

//var proto = SparseVector.prototype

function copy_row(nnz,
                  offset,
                  A,
                  /*int* index, // const int*
                    double* value, // const double*
                  */

                  x) // x = sparseVector
{
  for(var i=nnz-1; i>=0; --i) if(A.value[i + offset]){
    x.push_front(new SparseEntry(A.column_index[i + offset], A.value[i + offset]));
  }
  return true;
}

function givens( a,
                 b)
{
  var c
  var s

  if(b==0){
    c=1.0;
    s=0.0;
  }else{
    if(Math.abs(b)>Math.abs(a)){
      var tau=-a/b;
      s=1.0/Math.sqrt(1+tau*tau);
      c=s*tau;
      if(c<0.0){ c=-c; s=-s; }
    }else{
      var tau=-b/a;
      c=1.0/Math.sqrt(1.0+tau*tau);
      s=c*tau;
    }
  }

  return [c, s]
}



function apply_givens(/*SparseVector&*/ x,
  /*SparseVector&*/ y,
  /*int*/ diagonal) // c, s
{
  // find the rotation we need
  var a=x.front().value, b=y.front().value;

  var ret = givens(a, b);
  var c= ret[0];
  var s= ret[1];

  // rotate the start of each list
  x.front().value=c*a-s*b;
  y.pop_front();
  // then update the rest (x_new = c*x-s*y and y_new = s*x+c*y)
  var p=x.begin(), q=y.begin();
  p.inc(); // skip the first value we already took care of
  while(p.still_going() && q.still_going()){
    if(p.arrow().index==q.arrow().index){ // two here.
      var xnew=c*p.arrow().value-s*q.arrow().value,
          ynew=s*p.arrow().value+c*q.arrow().value;
      if(xnew){
        p.arrow().value=xnew;
        p.inc();
      }else x.erase(p);
      if(ynew){
        q.arrow().value=ynew;
        q.inc();
      }else y.erase(q);
    }else if(p.arrow().index<q.arrow().index){
      var k=p.arrow().index;
      var xnew=c*p.arrow().value,
          ynew=s*p.arrow().value;
      p.arrow().value=xnew;
      p.inc();
      y.insert(q, new SparseEntry(k, ynew));
      q.inc();
    }else{
      var k=q.arrow().index;
      var xnew=-s*q.arrow().value,
          ynew=c*q.arrow().value;
      x.insert(p, new SparseEntry(k, xnew));
      p.inc();
      q.arrow().value=ynew;
      q.inc();
    }
  }
  if(p.still_going()){
    do{
      var k=p.arrow().index;
      var xnew=c*p.arrow().value,
          ynew=s*p.arrow().value;
      p.arrow().value=xnew;
      p.inc();
      y.insert(q, new SparseEntry(k, ynew));
      q.inc();
    }while(p.still_going());
  }else if(q.still_going()){
    do{
      var k=q.arrow().index;
      var xnew=-s*q.arrow().value,
          ynew=c*q.arrow().value;
      x.insert(p, new SparseEntry(k, xnew));
      p.inc();
      q.arrow().value=ynew;
      q.inc();
    }while(q.still_going());
  }
  return ret;

}



function QUERN_compute_qr(m,
                          n,

                          A,
                          out_Q,
                          out_R)
{
  if(m<=0 || n<=0 || m<n)
    // TODO: other checks, on A Q and R
    return QUERN_INPUT_ERROR;

  // set up lists for dynamically building Q and R

  //  SparseVector* Q=(SparseVector*)std::malloc(m*sizeof(SparseVector));
  var Q = []
  var R = []
  for(var i = 0; i < m; ++i) {
    Q[i] = new SparseVector()
    R[i] = new SparseVector()
  }

  // do the Givens QR
  var row = new SparseVector();
  var c, s;



  for(var a=0; a<m; ++a) {
    var i= a

    copy_row(
      A.row_start[i+1]-A.row_start[i],
      A.row_start[i],
      A,
      row)

    var q=Q[a].begin();
    while(!row.empty() && row.front().index<a && row.front().index<n) {
      var j=row.front().index;
      if(R[j].empty() || R[j].front().index>j){ // swap?
        R[j].swap(row);
        Q[a].insert(q, new SparseEntry(j, 1));
        q.inc();
      }else{ // use Givens

        var ret = apply_givens(R[j], row, j);

        Q[a].insert(q, new SparseEntry(j, encode(ret[0], ret[1])));
        q.inc();

      }
    }

    if(a<n){
      R[a].swap(row);
    }
  }

  //    var R_row_start=(int*)std::malloc((n+1)*sizeof(int));
  var R_row_start = []
  for(var i = 0; i < n+1; ++i) {
    R_row_start[i] = 0
  }

  R_row_start[0]=0;

  for(var i=0; i<n; ++i)
    R_row_start[i+1]=R_row_start[i]+R[i].size();
  var  Rnnz=R_row_start[n];

  ///    int* R_column_index=(int*)std::malloc(Rnnz*sizeof(int));
  var R_column_index = []
  for(var i = 0; i < Rnnz; ++i) {
    R_column_index[i] = 0
  }

  //    double* R_value=(double*)std::malloc(Rnnz*sizeof(double));
  var R_value = []
  for(var i = 0; i < Rnnz; ++i) {
    R_value[i] = 0
  }

  j=0;
  for(var i=0; i<n; ++i){
    var p;
    for(p=R[i].begin(); p.still_going(); p.inc()){
      R_column_index[j]=p.arrow().index;
      R_value[j]=p.arrow().value;
      ++j;
    }
  }


  var Q_row_start= []

  Q_row_start[0]=0;
  for(var i=0; i<m; ++i)
    Q_row_start[i+1]=Q_row_start[i]+Q[i].size();
  var Qnnz=Q_row_start[m];
  var Q_column_index= []

  var Q_value= []

  var j=0;
  for(var i=0; i<m; ++i){
    for(var q=Q[i].begin(); q.still_going(); q.inc()){
      Q_column_index[j]=q.arrow().index;
      Q_value[j]=q.arrow().value;
      ++j;
    }
  }
  // transfer R's lists to static CSR form
  var R_row_start= []

  R_row_start[0]=0;
  for(var i=0; i<n; ++i)
    R_row_start[i+1]=R_row_start[i]+R[i].size();
  var Rnnz=R_row_start[n];
  var R_column_index= []

  var R_value= []

  j=0;
  for(var i=0; i<n; ++i){
    for(var p=R[i].begin(); p.still_going(); p.inc()){
      R_column_index[j]=p.arrow().index;
      R_value[j]=p.arrow().value;
      ++j;
    }
  }


  out_Q.column_index = Q_column_index
  out_Q.row_start = Q_row_start
  out_Q.value = Q_value

  out_R.column_index = R_column_index
  out_R.row_start = R_row_start
  out_R.value = R_value
}

function QUERN_multiply_with_q_transpose(m,
                                         Q,
                                         x)
{
  if(m<=0)
    return QUERN_INPUT_ERROR;
  var i, j, k;
  var c, s;
  for(i=0; i<m; ++i){
    for(j=Q.row_start[i]; j<Q.row_start[i+1]; ++j){
      k=Q.column_index[j];
      if(Q.value[j]==1){ // swap?
        //                std::swap(x[i], x[k]);
        var temp = x[i]
        x[i] = x[k]
        x[k] = temp
      }else{
        var c = Q.value[j][0]
        var s = Q.value[j][1]

        var newxk=c*x[k]-s*x[i];
        x[i]=s*x[k]+c*x[i];
        x[k]=newxk;
      }
    }
  }
}

function QUERN_solve_with_r(intn,
                            R,
                            rhs,
                            result)
{
  // check input
  if(n<=0)
    return QUERN_INPUT_ERROR;
  // do the solve
  var x, rii;
  var i, j;
  for(i=n-1; i>=0; --i) {
    x=rhs[i];
    rii=0;
    j=R.row_start[i];
    if(j<R.row_start[i+1] && R.column_index[j]==i){
      rii=R.value[j];
      ++j;
    }
    if(rii){
      for(; j<R.row_start[i+1]; ++j){
        x-=R.value[j]*result[R.column_index[j]];
      }

      result[i]=x/rii;
    }else
      result[i]=0;
  }
}


function solve(m, n, A, b) {

  var Q = {
  }

  var R = {
  }

  QUERN_compute_qr(m, n,
                   A,
                   Q,
                   R
                  );

  QUERN_multiply_with_q_transpose(m,
                                  Q,
                                  b
                                 );
  //console.log(b)


  var solution = []

  QUERN_solve_with_r(n,
                     R,
                     b,
                     solution);


  return solution
  //console.log(A)
}

var A = {

  row_start: [
    0,
    3,
    6,
    9,
    11
  ],

  column_index: [
    0, 1, 2,
    0, 1, 2,
    0, 1, 2,

    0, 1
  ],

  value: [
    1, -1,  4,
    1,  4, -2,
    1,  4,  2,
    1,  -1
  ]
}
var b =  [11, 3, 15, -1]

var m = 4
var n = 3

console.log(solve(m, n, A ,b))