// Taken from Chap 10, p 292 of 1987 edition of Numerical Recipes, by William Press, Brian Flannery, Saul Teukolsky, and william Vetterling.  For a more recent edition, see https://s3.amazonaws.com/nrbook.com/book_F210.html

class Minimizer {
  // Multidimensional minimization of the function fn(x) where x is a vector in nDim dimensions, by the downhill simplex method of Nelder and Mead. The matrix pIn is input. Its nDim + 1 rows are nDim-dimensional vectors which are the vertices of the starting simplex. Also input are fTol the fractional convergence tolerance to be achieved in the function value (n.b.!) and itMax, the maximimum allowed number of iterations. If itMax is not input, it will be set to 500.  If fTol also is not input, it will be set to 1e-10.  After invoking the run method, the p and y fields will have been reset to nDim+1 new points all within ftol of a minimum function value, and iter gives the number of function evaluations taken.
  constructor (fn, nDimOrP, maxIter) {
    this.fn = fn;
    this.maxIter = maxIter || 500;
    // The following four fields will probably get mutated by an invocation of the run method.
    const isArray = Array.isArray(nDimOrP) && nDimOrP.length && Array.isArray(nDimOrP[0]) && nDimOrP[0].length;
    const isNumber = typeof nDimOrP === "number" && Number.isInteger(nDimOrP) && nDimOrP > 0;

    this.p = null;
    if (isArray) {
      this.p = JSON.parse(JSON.stringify(nDimOrP));
    }
    if (isNumber) {
      this.p = ((new Array(nDimOrP + 1).fill(null))).map((arr, i) => {
        return (new Array(nDimOrP).fill(0)).map((elem, j) => (i === j && i !== nDimOrP + 1) ? 1 : 0);
      });
    }

    this.y = null;
    this.iter = 0;
    const arg2 = `second argument (${JSON.stringify(nDimOrP)})`;

    this.error = null;
    if (typeof fn !== "function") {
      this.error = `The first argument (${fn}) of Minimizer must be a function not a ${typeof fn}.`;
    }
    if (!((typeof nDimOrP === "number" && Number.isInteger(nDimOrP) && nDimOrP > 0) || (isArray && nDimOrP.length > 1))) {
      this.error = `The ${arg2} is neither a positive integer nor an array whose length exceeds 1.`;
    }
    if (isArray && !(nDimOrP.every(arr => (Array.isArray(arr) && arr.length === nDimOrP.length - 1)))) {
      this.error = `Each element of the ${arg2} needs itself to be a ${nDimOrP.length - 1}-element array.`;
    }
    if (isArray && !(nDimOrP.every(arr => arr.every(element => typeof element === "number" && Number.isFinite(element) && !Number.isNaN(element))))) {
      this.error = `Each element of each element of the ${arg2} must be a number.`;
    }
    this.warnings = [];
  }

  run(fTol) {
    if (this.error) return this;
    fTol = fTol || 1e-10;
    const [alpha, beta, gamma] = [1, 0.5, 2];
    this.y = [];
    for (const vertex of this.p) {
      const result = this.fn(vertex);
      this.warnings = result.warnings;
      this.error = result.error;
      if (this.error) return this;
      this.y.push(result.value);
    }
    if (!this.y.every(element => typeof element === "number" && Number.isFinite(element) && !Number.isNaN(element))) {
      this.error = `The function is Infinite or NaN for at least one of the vertices ${JSON.stringify(this.y)} of the simplex.`;
      return this;
    }
    const mPts = this.p.length
    const nDim = mPts - 1;
    while (this.iter < this.maxIter) {
      //First determine which point is the highest (worst), next-highest, and lowest (best), by looping over the points in the simplex.
      let iLow = 0;
      let [iHigh, iNextHigh] = (this.y[0] > this.y[1]) ? [0, 1] : [0, 1];
      for (let i = 0; i < mPts; i++) {
        if (this.y[i] < this.y[iLow]) iLow = i;
        if (this.y[i] > this.y[iHigh]) {
          [iNextHigh, iHigh] = [iHigh, i];
        } else if (this.y[i] > this.y[iNextHigh]) {
          if (i !== iHigh) iNextHigh = i;
        }
      }
      // Compute the fractional range from highest to lowest and return if satisfactory.
      const rTol = 2 * Math.abs(this.y[iHigh] - this.y[iLow]) / (Math.abs(this.y[iHigh]) + Math.abs(this.y[iLow]));
      if (rTol < fTol) return this;
      this.iter++;
      // Begin a new iteration.  Computer the vector average of all points except for the highest, ie the center of the "face" of the simplex across from the high point.  We will subsequently explore along the ray from the high point through that center.
      const pBar = new Array(nDim).fill(0);
      for (let j = 0; j < nDim; j++) {
        for (let i = 0; i < mPts; i++) {
          if (i !== iHigh) pBar[j] += this.p[i][j] / nDim;
        }
      }
      // Extrapolate by a factor alpha through the face, ie reflect the simplex from the high point.
      let pR = pBar.map((coord, j) => (1 + alpha) * coord - alpha * this.p[iHigh][j]);
      // Evaluate the function at the reflected point.
      let result = this.fn(pR);
      this.warnings = result.warnings;
      this.error = result.error;
      if (this.error) return this;
      const yPr = result.value;
      if (yPr <= this.y[iLow]) {
        // This gives a better result than the high point, so try an additional extrapolation by a factor gamma.
        const pRr = pR.map((coord, j) => gamma * coord + (1 - gamma) * pBar[j]);
        // Check out the function after gamma-extrapolation.
        const result = this.fn(pRr);
        this.warnings = result.warnings;
        this.error = result.error;
        if (this.error) return this;
        const yPrr = result.value;
        // This ternary handles cases when the gamma-extrapolation is better or worse.
        [this.p[iHigh], this.y[iHigh]] = (yPrr < this.y[iLow]) ? [[...pRr], yPrr] : [[...pR], yPr];
      } else if (yPr >= iNextHigh) {
        // The reflected point is worse than the second highest.
        // The following line replaces the highest, if it's better.
        if (yPr < this.y[iHigh]) [this.p[iHigh], this.y[iHigh]] = [[...pR], yPr];
        // Look for an intermediate lower point.  In other words, perform a contraction of the simplex along one dimension ...
        const pRr = this.p[iHigh].map((coord, j) => beta * coord + (1 - beta) * pBar[j]);
        // ... and then evaluate the function.
        const result = this.fn(pRr);
        this.warnings = result.warnings;
        this.error = result.error;
        if (this.error) return this;
        const yPrr = result.value;
        if (yPrr < this.y[iHigh]) {
          // Contraction gives an improvement, so accept it.
          [this.p[iHigh], this.y[iHigh]] = [[...pRr], yPrr];
        } else {
          // We cannot seem to get rid of that high point, so we better contract along the lowest (best) point.
          for (let i = 0; i < mPts; i++) {
            if (i !== iLow) {
              pR = this.p[i].map((coord, j) => (coord + this.p[iLow][j]) / 2);
              const result = this.fn(pR);
              this.warnings = result.warnings;
              this.error = result.error;
              if (this.error) return this;
              [this.p[i], this.y[i]] = [[...pR], result.value];
            }
          }
        }
      } else {
        // We arrive here if the original reflection gives a middling point.  Replace the old high point and continue.
        [this.p[iHigh], this.y[iHigh]] = [[...pR], yPr];
      }
    }
    this.error = `Exceeded maximum iterations, which you set as ${this.maxIter}.`;
    return this;
  }
}

module.exports = Minimizer;
