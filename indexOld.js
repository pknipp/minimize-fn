// Taken from Chap 10, p 404 of https://s3.amazonaws.com/nrbook.com/book_F210.html
const amoeba = (p,fTol,func) => {
const nMax = 20;
const itMax = 5000;
const tiny = 1.e-10;
const y = p.map(vec => func(vec));
const nDim = p[0].length;

let iter = 0;
// Enter here when starting or have just overall contracted.
while (true) {
1 do 12 n=1,ndim
// Recompute psum.
const sum = 0;
for (let m = 0; m <= ndim; m++) {
    sum += p[m][n];
}
psum[n] = sum;
enddo 12
// Enter here when have just changed a single point.
2 const ilo = 1;
// Determine which point is the highest (worst), next-highest, and lowest (best) ...
if (y[0] > y[1]) {
    ihi = 1;
    inhi = 2;
} else {
    ihi = 2;
    inhi = 1;
}
// ... by looping over the points in the simplex.
for (let i = 0; i <= ndim; i++) {
    if(y[i] <= y[ilo]) ilo = i;
    if (y[i] > y[ihi]) {
        inhi = ihi;
        ihi = i;
    } else if (y[i] > y[inhi]) {
        if (i !== ihi) inhi = i;
    }
}

rtol = 2 * Math.abs(y[ihi] - y[ilo]) / (Math.abs(y[ihi]) + Math.abs(y[ilo]) + TINY);
// Compute the fractional range from highest to lowest and return if satisfactory. If returning, put best point and value in slot 1.
if (rtol < ftol) {
    swap = y[0];
    y[1] = y[ilo];
    y[ilo] = swap;
    for (let i = 0; i < ndim; i++) {
        swap = p[0][n];
        p[0, n] = p[ilo, n];
        p[ilo][n] = swap;
    }
    return;
}
if (iter >= ITMAX] pause ’ITMAX exceeded in amoeba’;
iter += 2;
Begin a new iteration. First extrapolate by a factor −1 through the face of the simplex across
from the high point, i.e., reflect the simplex from the high point.
ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0)
if (ytry.le.y(ilo)) then
Gives a result better than the best point, so try an additional extrapolation by a factor 2.
ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0)
else if (ytry.ge.y(inhi)) then
The reflected point is worse than the second-highest, so look for an intermediate lower point,
i.e., do a one-dimensional contraction.
ysave=y(ihi)
ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5)
// Can’t seem to get rid of that high point. Better contract around the lowest (best) point.
if (ytry.ge.ysave) then
do 16 i=1,ndim+1
if(i.ne.ilo)then
do 15 j=1,ndim
psum(j)=0.5*(p(i,j)+p(ilo,j))
p(i,j)=psum(j)
enddo 15
y(i)=funk(psum)
endif
enddo 16
// Keep track of function evaluations.
iter += ndim;
// Go back for the test of doneness and the next iteration.
goto 1
endif
else
// Correct the evaluation count.
iter++;
endif
goto 2
END
}

module.exports = minimize
