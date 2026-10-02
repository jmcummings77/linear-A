/** Real Hessenberg reduction and double-shift QR, adapted from public-domain
 * NIST/MathWorks JAMA 1.0.3 orthes/hqr2 (EISPACK/Martin-Wilkinson).
 * https://math.nist.gov/javanumerics/jama/
 * Workspace is private; eigenvectors are right columns, possibly dependent.
 */
function scalePower(value: number, exponent: number): number {
  while (exponent > 512) { value *= 2**512; exponent -= 512; }
  while (exponent < -512) { value *= 2**-512; exponent += 512; }
  return value*2**exponent;
}

function balance(work: number[][]): number[] {
  const n = work.length, powers = Array<number>(n).fill(0);
  const logNorm = (values: number[]) => {
    const largest = Math.max(...values.map(Math.abs),0);
    return largest === 0 ? null : Math.log2(largest)+Math.log2(values.reduce((sum,x)=>sum+Math.abs(x)/largest,0));
  };
  const logSum = (a: number,b: number) => Math.max(a,b)+Math.log2(1+2**(-Math.abs(a-b)));
  for (let sweep=0;sweep<64;sweep++) {
    let changed = false;
    for (let i=0;i<n;i++) {
      const row=work[i].filter((_,j)=>i!==j), col=work.filter((_,j)=>i!==j).map(row=>row[i]);
      const r=logNorm(row),c=logNorm(col);
      if (r===null || c===null) continue;
      const shift=Math.max(-512,Math.min(512,Math.round((r-c)/2)));
      if (!shift || logSum(r-shift,c+shift)>=logSum(r,c)+Math.log2(.95)) continue;
      const safe=(values:number[],power:number)=>values.every(x=>{
        const y=scalePower(x,power);
        return Number.isFinite(y) && (x===0 || y!==0) && scalePower(y,-power)===x;
      });
      if (!safe(row,-shift) || !safe(col,shift)) continue;
      for (let j=0;j<n;j++) if (i!==j) {
        work[i][j]=scalePower(work[i][j],-shift); work[j][i]=scalePower(work[j][i],shift);
      }
      powers[i]+=shift; changed=true;
    }
    if (!changed) break;
  }
  return powers;
}

export function generalEigen(input: Float64Array, size: number, maxIterations: number) {
  if (!Number.isSafeInteger(maxIterations) || maxIterations < 1 || maxIterations > 100000)
    throw new RangeError("maxIterations must be an integer between 1 and 100000");
  if (input.some(value => !Number.isFinite(value))) throw new RangeError("eigendecomposition requires finite inputs");
  let scale = 0, diagonal = true;
  for (let i = 0; i < size; i++) for (let j = 0; j < size; j++) {
    scale = Math.max(scale, Math.abs(input[i * size + j]));
    if (i !== j && input[i * size + j] !== 0) diagonal = false;
  }
  const n = size, d = Array<number>(n).fill(0), e = Array<number>(n).fill(0), ort = Array<number>(n).fill(0);
  const V = Array.from({length:n}, (_, i) => Array.from({length:n}, (_, j) => +(i === j)));
  let upper=true,lower=true;
  for (let i=0;i<n;i++) for (let j=0;j<n;j++) {
    if (i>j && input[i*n+j]!==0) upper=false;
    if (i<j && input[i*n+j]!==0) lower=false;
  }
  const reversedBasis=lower&&!upper;
  const H = Array.from({length:n}, (_,i) => Array.from({length:n}, (_,j) =>
    reversedBasis?input[(n-1-i)*n+n-1-j]:input[i*n+j]));
  const originalDiagonal=H.map((row,i)=>row[i]);
  const powers=upper||lower?Array<number>(n).fill(0):balance(H);
  scale=0; for (const row of H) for (const x of row) scale=Math.max(scale,Math.abs(x));
  scale=diagonal||scale===0?1:scale;
  for (const row of H) for(let j=0;j<n;j++) row[j]/=scale;
   function orthes() {
      let low = 0;
      let high = size-1;
      for (let m = low+1; m <= high-1; m++) {
         let scale = 0.0;
         for (let i = m; i <= high; i++) {
            scale = scale + Math.abs(H[i][m-1]);
         }
         if (scale != 0.0) {
            let h = 0.0;
            for (let i = high; i >= m; i--) {
               ort[i] = H[i][m-1]/scale;
               h += ort[i] * ort[i];
            }
            let g = Math.sqrt(h);
            if (ort[m] > 0) {
               g = -g;
            }
            h = h - ort[m] * g;
            ort[m] = ort[m] - g;
            for (let j = m; j < n; j++) {
               let f = 0.0;
               for (let i = high; i >= m; i--) {
                  f += ort[i]*H[i][j];
               }
               f = f/h;
               for (let i = m; i <= high; i++) {
                  H[i][j] -= f*ort[i];
               }
           }
           for (let i = 0; i <= high; i++) {
               let f = 0.0;
               for (let j = high; j >= m; j--) {
                  f += ort[j]*H[i][j];
               }
               f = f/h;
               for (let j = m; j <= high; j++) {
                  H[i][j] -= f*ort[j];
               }
            }
            ort[m] = scale*ort[m];
            H[m][m-1] = scale*g;
         }
      }
      for (let i = 0; i < n; i++) {
         for (let j = 0; j < n; j++) {
            V[i][j] = (i == j ? 1.0 : 0.0);
         }
      }
      for (let m = high-1; m >= low+1; m--) {
         if (H[m][m-1] != 0.0) {
            for (let i = m+1; i <= high; i++) {
               ort[i] = H[i][m-1];
            }
            for (let j = m; j <= high; j++) {
               let g = 0.0;
               for (let i = m; i <= high; i++) {
                  g += ort[i] * V[i][j];
               }
               g = (g / ort[m]) / H[m][m-1];
               for (let i = m; i <= high; i++) {
                  V[i][j] += g * ort[i];
               }
            }
         }
      }
   }
   let cdivr = 0, cdivi = 0;
   function cdiv(xr: number, xi: number, yr: number, yi: number) {
      let r = 0, d = 0;
      if (Math.abs(yr) > Math.abs(yi)) {
         r = yi/yr;
         d = yr + r*yi;
         cdivr = (xr + r*xi)/d;
         cdivi = (xi - r*xr)/d;
      } else {
         r = yr/yi;
         d = yi + r*yr;
         cdivr = (r*xr + xi)/d;
         cdivi = (r*xi - xr)/d;
      }
   }
   function hqr2() {
      let nn = size;
      let n = nn-1;
      let low = 0;
      let high = nn-1;
      let eps = Number.EPSILON;
      let exshift = 0.0;
      let p=0, q=0, r=0, s=0, z=0, t = 0, w = 0, x = 0, y = 0;
      let norm = 0.0;
      for (let i = 0; i < nn; i++) {
         if (i < low || i > high) {
            d[i] = H[i][i];
            e[i] = 0.0;
         }
         for (let j = Math.max(i-1,0); j < nn; j++) {
            norm = norm + Math.abs(H[i][j]);
         }
      }
      let iter = 0;
      while (n >= low) {
         let l = n;
         while (l > low) {
            s = Math.abs(H[l-1][l-1]) + Math.abs(H[l][l]);
            if (s == 0.0) {
               s = norm;
            }
            if (Math.abs(H[l][l-1]) <= eps * s) {
               break;
            }
            l--;
         }
         if (l == n) {
            H[n][n] = H[n][n] + exshift;
            d[n] = H[n][n];
            e[n] = 0.0;
            n--;
            iter = 0;
         } else if (l == n-1) {
            w = H[n][n-1] * H[n-1][n];
            p = (H[n-1][n-1] - H[n][n]) / 2.0;
            q = p * p + w;
            z = Math.sqrt(Math.abs(q));
            H[n][n] = H[n][n] + exshift;
            H[n-1][n-1] = H[n-1][n-1] + exshift;
            x = H[n][n];
            if (q >= 0) {
               if (p >= 0) {
                  z = p + z;
               } else {
                  z = p - z;
               }
               d[n-1] = x + z;
               d[n] = d[n-1];
               if (z != 0.0) {
                  d[n] = x - w / z;
               }
               e[n-1] = 0.0;
               e[n] = 0.0;
               x = H[n][n-1];
               s = Math.abs(x) + Math.abs(z);
               p = x / s;
               q = z / s;
               r = Math.sqrt(p * p+q * q);
               p = p / r;
               q = q / r;
               for (let j = n-1; j < nn; j++) {
                  z = H[n-1][j];
                  H[n-1][j] = q * z + p * H[n][j];
                  H[n][j] = q * H[n][j] - p * z;
               }
               for (let i = 0; i <= n; i++) {
                  z = H[i][n-1];
                  H[i][n-1] = q * z + p * H[i][n];
                  H[i][n] = q * H[i][n] - p * z;
               }
               for (let i = low; i <= high; i++) {
                  z = V[i][n-1];
                  V[i][n-1] = q * z + p * V[i][n];
                  V[i][n] = q * V[i][n] - p * z;
               }
            } else {
               d[n-1] = x + p;
               d[n] = x + p;
               e[n-1] = z;
               e[n] = -z;
            }
            n = n - 2;
            iter = 0;
         } else {
            x = H[n][n];
            y = 0.0;
            w = 0.0;
            if (l < n) {
               y = H[n-1][n-1];
               w = H[n][n-1] * H[n-1][n];
            }
            if (iter == 10) {
               exshift += x;
               for (let i = low; i <= n; i++) {
                  H[i][i] -= x;
               }
               s = Math.abs(H[n][n-1]) + Math.abs(H[n-1][n-2]);
               x = y = 0.75 * s;
               w = -0.4375 * s * s;
            }
            if (iter == 30) {
                s = (y - x) / 2.0;
                s = s * s + w;
                if (s > 0) {
                    s = Math.sqrt(s);
                    if (y < x) {
                       s = -s;
                    }
                    s = x - w / ((y - x) / 2.0 + s);
                    for (let i = low; i <= n; i++) {
                       H[i][i] -= s;
                    }
                    exshift += s;
                    x = y = w = 0.964;
                }
            }
            iter = iter + 1;
            if (iter > maxIterations) throw new RangeError("general eigendecomposition did not converge within maximum iterations");
            let m = n-2;
            while (m >= l) {
               z = H[m][m];
               r = x - z;
               s = y - z;
               p = (r * s - w) / H[m+1][m] + H[m][m+1];
               q = H[m+1][m+1] - z - r - s;
               r = H[m+2][m+1];
               s = Math.abs(p) + Math.abs(q) + Math.abs(r);
               p = p / s;
               q = q / s;
               r = r / s;
               if (m == l) {
                  break;
               }
               if (Math.abs(H[m][m-1]) * (Math.abs(q) + Math.abs(r)) <
                  eps * (Math.abs(p) * (Math.abs(H[m-1][m-1]) + Math.abs(z) +
                  Math.abs(H[m+1][m+1])))) {
                     break;
               }
               m--;
            }
            for (let i = m+2; i <= n; i++) {
               H[i][i-2] = 0.0;
               if (i > m+2) {
                  H[i][i-3] = 0.0;
               }
            }
            for (let k = m; k <= n-1; k++) {
               let notlast = (k != n-1);
               if (k != m) {
                  p = H[k][k-1];
                  q = H[k+1][k-1];
                  r = (notlast ? H[k+2][k-1] : 0.0);
                  x = Math.abs(p) + Math.abs(q) + Math.abs(r);
                  if (x == 0.0) {
                      continue;
                  }
                  p = p / x;
                  q = q / x;
                  r = r / x;
               }
               s = Math.sqrt(p * p + q * q + r * r);
               if (p < 0) {
                  s = -s;
               }
               if (s != 0) {
                  if (k != m) {
                     H[k][k-1] = -s * x;
                  } else if (l != m) {
                     H[k][k-1] = -H[k][k-1];
                  }
                  p = p + s;
                  x = p / s;
                  y = q / s;
                  z = r / s;
                  q = q / p;
                  r = r / p;
                  for (let j = k; j < nn; j++) {
                     p = H[k][j] + q * H[k+1][j];
                     if (notlast) {
                        p = p + r * H[k+2][j];
                        H[k+2][j] = H[k+2][j] - p * z;
                     }
                     H[k][j] = H[k][j] - p * x;
                     H[k+1][j] = H[k+1][j] - p * y;
                  }
                  for (let i = 0; i <= Math.min(n,k+3); i++) {
                     p = x * H[i][k] + y * H[i][k+1];
                     if (notlast) {
                        p = p + z * H[i][k+2];
                        H[i][k+2] = H[i][k+2] - p * r;
                     }
                     H[i][k] = H[i][k] - p;
                     H[i][k+1] = H[i][k+1] - p * q;
                  }
                  for (let i = low; i <= high; i++) {
                     p = x * V[i][k] + y * V[i][k+1];
                     if (notlast) {
                        p = p + z * V[i][k+2];
                        V[i][k+2] = V[i][k+2] - p * r;
                     }
                     V[i][k] = V[i][k] - p;
                     V[i][k+1] = V[i][k+1] - p * q;
                  }
               }
            }
         }
      }
      if (norm == 0.0) {
         return;
      }
      for (n = nn-1; n >= 0; n--) {
         p = d[n];
         q = e[n];
         if (q == 0) {
            let l = n;
            H[n][n] = 1.0;
            for (let i = n-1; i >= 0; i--) {
               w = H[i][i] - p;
               r = 0.0;
               for (let j = l; j <= n; j++) {
                  r = r + H[i][j] * H[j][n];
               }
               if (e[i] < 0.0) {
                  z = w;
                  s = r;
               } else {
                  l = i;
                  if (e[i] == 0.0) {
                     if (w != 0.0) {
                        H[i][n] = -r / w;
                     } else {
                        H[i][n] = -r / (eps * norm);
                     }
                  } else {
                     x = H[i][i+1];
                     y = H[i+1][i];
                     q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
                     t = (x * s - z * r) / q;
                     H[i][n] = t;
                     if (Math.abs(x) > Math.abs(z)) {
                        H[i+1][n] = (-r - w * t) / x;
                     } else {
                        H[i+1][n] = (-s - y * t) / z;
                     }
                  }
                  t = Math.abs(H[i][n]);
                  if ((eps * t) * t > 1) {
                     for (let j = i; j <= n; j++) {
                        H[j][n] = H[j][n] / t;
                     }
                  }
               }
            }
         } else if (q < 0) {
            let l = n-1;
            if (Math.abs(H[n][n-1]) > Math.abs(H[n-1][n])) {
               H[n-1][n-1] = q / H[n][n-1];
               H[n-1][n] = -(H[n][n] - p) / H[n][n-1];
            } else {
               cdiv(0.0,-H[n-1][n],H[n-1][n-1]-p,q);
               H[n-1][n-1] = cdivr;
               H[n-1][n] = cdivi;
            }
            H[n][n-1] = 0.0;
            H[n][n] = 1.0;
            for (let i = n-2; i >= 0; i--) {
               let ra = 0, sa = 0, vr = 0, vi = 0;
               ra = 0.0;
               sa = 0.0;
               for (let j = l; j <= n; j++) {
                  ra = ra + H[i][j] * H[j][n-1];
                  sa = sa + H[i][j] * H[j][n];
               }
               w = H[i][i] - p;
               if (e[i] < 0.0) {
                  z = w;
                  r = ra;
                  s = sa;
               } else {
                  l = i;
                  if (e[i] == 0) {
                     cdiv(-ra,-sa,w,q);
                     H[i][n-1] = cdivr;
                     H[i][n] = cdivi;
                  } else {
                     x = H[i][i+1];
                     y = H[i+1][i];
                     vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
                     vi = (d[i] - p) * 2.0 * q;
                     if (vr == 0.0 && vi == 0.0) {
                        vr = eps * norm * (Math.abs(w) + Math.abs(q) +
                        Math.abs(x) + Math.abs(y) + Math.abs(z));
                     }
                     cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi);
                     H[i][n-1] = cdivr;
                     H[i][n] = cdivi;
                     if (Math.abs(x) > (Math.abs(z) + Math.abs(q))) {
                        H[i+1][n-1] = (-ra - w * H[i][n-1] + q * H[i][n]) / x;
                        H[i+1][n] = (-sa - w * H[i][n] - q * H[i][n-1]) / x;
                     } else {
                        cdiv(-r-y*H[i][n-1],-s-y*H[i][n],z,q);
                        H[i+1][n-1] = cdivr;
                        H[i+1][n] = cdivi;
                     }
                  }
                  t = Math.max(Math.abs(H[i][n-1]),Math.abs(H[i][n]));
                  if ((eps * t) * t > 1) {
                     for (let j = i; j <= n; j++) {
                        H[j][n-1] = H[j][n-1] / t;
                        H[j][n] = H[j][n] / t;
                     }
                  }
               }
            }
         }
      }
      for (let i = 0; i < nn; i++) {
         if (i < low || i > high) {
            for (let j = i; j < nn; j++) {
               V[i][j] = H[i][j];
            }
         }
      }
      for (let j = nn-1; j >= low; j--) {
         for (let i = low; i <= high; i++) {
            z = 0.0;
            for (let k = low; k <= Math.min(j,high); k++) {
               z = z + V[i][k] * H[k][j];
            }
            V[i][j] = z;
         }
      }
   }

  if (diagonal) for (let i = 0; i < n; i++) d[i] = input[i*n+i];
  else { orthes(); hqr2(); }
  if (upper||lower) {
    for(let i=0;i<n;i++) { d[i]=originalDiagonal[i];e[i]=0; }
    scale=1;
  }
  if (reversedBasis) V.reverse();
  const order = Array.from({length:n}, (_,i) => i).sort((a,b) => d[a]-d[b] || e[a]-e[b]);
  const valuesReal = new Float64Array(n), valuesImag = new Float64Array(n);
  const vectorsReal = new Float64Array(n*n), vectorsImag = new Float64Array(n*n);
  for (let col = 0; col < n; col++) {
    const source = order[col], realCol = e[source] < 0 ? source-1 : source;
    valuesReal[col] = d[source]*scale; valuesImag[col] = e[source]*scale;
    const re = V.map(row => row[realCol]);
    const im = V.map(row => e[source] === 0 ? 0 : row[realCol+1]*(e[source] < 0 ? -1 : 1));
    if (powers.some(x=>x!==0)) {
      let exponent=-Infinity;
      for(let row=0;row<n;row++) if(re[row]!==0||im[row]!==0)
        exponent=Math.max(exponent,Math.floor(Math.log2(Math.max(Math.abs(re[row]),Math.abs(im[row]))))+powers[row]);
      for(let row=0;row<n;row++) {
        re[row]=scalePower(re[row],powers[row]-exponent); im[row]=scalePower(im[row],powers[row]-exponent);
      }
    }
    let magnitude = 0, pivot = 0;
    for (let row = 0; row < n; row++) {
      magnitude = Math.max(magnitude, Math.abs(re[row]), Math.abs(im[row]));
      if (Math.hypot(re[row], im[row]) > Math.hypot(re[pivot], im[pivot])) pivot = row;
    }
    if (!Number.isFinite(magnitude) || magnitude === 0) throw new RangeError("nonfinite or zero eigenvector");
    let norm = 0;
    for (let row = 0; row < n; row++) { re[row] /= magnitude; im[row] /= magnitude; norm = Math.hypot(norm,re[row],im[row]); }
    const phase = Math.hypot(re[pivot],im[pivot]), pr = re[pivot]/phase, pi = im[pivot]/phase;
    for (let row = 0; row < n; row++) {
      vectorsReal[row*n+col] = (re[row]*pr+im[row]*pi)/norm;
      vectorsImag[row*n+col] = (im[row]*pr-re[row]*pi)/norm;
    }
    if (!Number.isFinite(valuesReal[col]) || !Number.isFinite(valuesImag[col])) throw new RangeError("nonfinite eigenvalue");
  }
  return {valuesReal,valuesImag,vectorsReal,vectorsImag};
}
