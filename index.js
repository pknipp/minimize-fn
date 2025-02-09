//amoeba.js

const minimize = (p, fTol, fn) => {
  const [alpha, beta, gamma] = [1, 0.5, 2];
  const y = p.map(vec => fn(vec));
  const nDim = p[0].length;
  const mPts = nDim + 1;
  let iter = 0;
  let msg = "";
  const returnMe = {y, iter, msg};
  while (true) {
    const iLow = 0;
    let [iHigh, iNextHigh] = (y[0] > y[1]) ? [0, 1] : [0, 1];
    for (let i = 0; i < mPts; i++) {
      if (y[i] < y[iLow]) iLow = i;
      if (y[i] > y[iHigh]) {
        [iNextHigh, iHigh] = [iHigh, i];
      } else if (y[i] > y[iNextHigh]) {
        if (i !== iHigh) iNextHigh = i;
      }
      const rTol = 2 * Math.abs(y[iHigh] - y[iLow]) / (Math.abs(y[iHigh]) + Math.abs(y[iLow]));
      if (rTol < fTol) return returnMe;
      if (iter === itMax) return {...returnMe, msg: "Exceeded maximum iterations."};
      iter++;
      const pBar = new Array(nDim).fill(0);
      for (let i = 0; i < mPts; i++) {
        if (i !== iHigh) {
          for (let j = 0; j < nDim; j++) {
            pBar[j] += p[i][j];
          }
        }
      }
      const pR = [];
      // 12
      for (let i = 0; i < mPts; i++) {
        if (i !== iHigh) {
          for (let j = 0; j < nDim; j++) {
            pBar[j] += p[i][j];
          }
        }
      }
      //14
      for (let j = 0; j < nDim; j++) {
        pBar[j] /= nDim;
        pR.push((1 + alpha) * pBar[j] - alpha * p[iHigh][j]);
      }
      //15
      const yPr = fn(pR);
      if (yPr <= y[iLow]) {
        const pRr = [];
        for (let j = 0; j < nDim; j++) {
          pRr.push(gamma * pR[j] + (1 - gamma) * pBar[j]);
        }
        //16
        const yPrr = fn(pRr);
        if (yPrr < y[iLow]) {
          for (let j = 0; j < nDim; j++) {
            p[iHigh][j] = pRr[j];
          }
          //17
          y[iHigh] = yPrr;
        } else {
          for (let j = 0; j < nDim; j++) {
            p[iHigh][j] = pR[j];
          }
          //18
          y[iHigh] = yPr
        }
      } else if (yPr >= iNHigh) {}
    }
  }
}

module.exports = minimize
