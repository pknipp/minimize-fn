//amoeba.js

const minimize = (p, fn, fTol, itMax) => {
  const [alpha, beta, gamma] = [1, 0.5, 2];
  fTol = fTol || 1e-10;
  itMax = itMax || 500;
  const y = p.map(vec => fn(vec));
  const nDim = p[0].length;
  const mPts = nDim + 1;
  let iter = 0;
  const returnMe = {p, y};
  //1
  while (true) {
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
    if (rTol < fTol) return {...returnMe, iter};
    if (iter === itMax) return {...returnMe, iter, msg: "Exceeded maximum iterations."};
    iter++;
    //12
    const pBar = new Array(nDim).fill(0);
    //14
    for (let i = 0; i < mPts; i++) {
      if (i !== iHigh) {
        for (let j = 0; j < nDim; j++) {
          pBar[j] += p[i][j] / nDim;
        }
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
      if (yPr < y[iHigh]) {
        //19
        [p[iHigh], y[iHigh]] = [[...pR], yPr];
      }
      //21
      const pRr = p[iHigh].map((coord, j) => beta * coord + (1 - beta) * pBar[j]);
      yPrr = fn(pRr);
      if (yPrr < y[iHigh]) {
        //22
        p[iHigh] = [...pRr];
        y[iHigh] = yPrr;
      } else {
        //24
        for (let i = 0; i < mPts; i++) {
          if (i !== iLow) {
            //23
            for (let j = 0; j < nDim; j++) {
              p[iHigh][j] = pR[j];
            }
          }
        }
      }
    } else {
      //25
      p[iHigh] = [...pR];
      y[iHigh] = yPr;
    }
  }
}

module.exports = minimize
