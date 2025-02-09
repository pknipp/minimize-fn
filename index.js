//// Taken from Chap 10, p 292 of 1987 edition of Numerical Recipes, by William Press, Brian Flannery, Saul Teukolsky, and william Vetterling.  For a more recent edition, see https://s3.amazonaws.com/nrbook.com/book_F210.html

const minimize = (pIn, fn, fTol, itMax) => {
  // Multidimensional minimization of the function fn(x) where x is a vector in nDim dimensions, by the downhill simplex method of Nelder and Mead. The matrix pIn is input. Its nDim + 1 rows are nDim-dimensional vectors which are the vertices of the starting simplex. Also input are fTol the fractional convergence tolerance to be achieved in the function value (n.b.!) and itMax, the maximimum allowed number of iterations. If itMax is not input, it will be set to 500.  If fTol also is not input, it will be set to 1e-10.  On output, p and y will have been reset to ndim+1 new points all within ftol of
// a minimum function value, and iter gives the number of function evaluations taken.
  const [alpha, beta, gamma] = [1, 0.5, 2];
  // Use default values, if not passed in.
  [fTol, itMax] = [fTol || 1e-10, itMax || 500];
  // Make a copy, to avoid mutating client's data structure.
  const p = JSON.parse(JSON.stringify(pIn));
  const y = p.map(vec => fn(vec));
  const nDim = p[0].length;
  const mPts = nDim + 1;
  let iter = 0;
  const results = {p, y};
  //1
  while (iter < itMax) {
    let iLow = 0;
    let [iHigh, iNextHigh] = (y[0] > y[1]) ? [0, 1] : [0, 1];
    //11
    for (let i = 0; i < mPts; i++) {
      if (y[i] < y[iLow]) iLow = i;
      if (y[i] > y[iHigh]) {
        [iNextHigh, iHigh] = [iHigh, i];
      } else if (y[i] > y[iNextHigh]) {
        if (i !== iHigh) iNextHigh = i;
      }
    }
    const rTol = 2 * Math.abs(y[iHigh] - y[iLow]) / (Math.abs(y[iHigh]) + Math.abs(y[iLow]));
    if (rTol < fTol) return {...results, iter};
    iter++;
    //12
    const pBar = new Array(nDim).fill(0);
    //14
    for (let j = 0; j < nDim; j++) {
      for (let i = 0; i < mPts; i++) {
        if (i !== iHigh) pBar[j] += p[i][j] / nDim;
      }
    }
    // 15
    const pR = pBar.map((coord, j) => (1 + alpha) * coord - alpha * p[iHigh][j]);
    const yPr = fn(pR);
    if (yPr <= y[iLow]) {
      // 16
      const pRr = pR.map((coord, j) => gamma * coord + (1 - gamma) * pBar[j]);
      const yPrr = fn(pRr);
      // 17 & 18
      [p[iHigh], y[iHigh]] = (yPrr < y[iLow]) ? [[...pRr], yPrr] : [[...pR], yPr];
    } else if (yPr >= iNextHigh) {
      //19
      if (yPr < y[iHigh]) [p[iHigh], y[iHigh]] = [[...pR], yPr];
      //21
      const pRr = p[iHigh].map((coord, j) => beta * coord + (1 - beta) * pBar[j]);
      yPrr = fn(pRr);
      if (yPrr < y[iHigh]) {
        //22
        [p[iHigh], y[iHigh]] = [[...pRr], yPrr];
      } else {
        //24
        for (let i = 0; i < mPts; i++) {
          if (i !== iLow) {
            //23
            pR = p[i].map((coord, j) => (coord + p[iLow][j]) / 2);
            [p[i], y[i]] = [[...pR], fn[pR]];
          }
        }
      }
    } else {
      //25
      [p[iHigh], y[iHigh]] = [[...pR], yPr];
    }
  }
  return {...results, iter, msg: "Exceeded maximum iterations."};
}

module.exports = minimize
