#ifndef LINEAR_A_GENERAL_EIGEN_H
#define LINEAR_A_GENERAL_EIGEN_H

/* Hessenberg reduction and implicit real double-shift QR, adapted from JAMA's
 * public-domain EigenvalueDecomposition (orthes/hqr2), NIST/MathWorks:
 * https://math.nist.gov/javanumerics/jama/
 * The original algorithms are from Martin/Wilkinson and EISPACK.
 * Shared arithmetic kernel; the C and C++ APIs own and validate their storage.
 * Changes include a bounded QR iteration count, power-of-two scaling, explicit
 * complex columns, finite-result checks, normalization, and eigenpair sorting. */
#include <float.h>
#include <math.h>
#include <stddef.h>

#define H(i,j) hessenberg[(size_t)(i) * (size_t)width + (size_t)(j)]
#define V(i,j) v[(size_t)(i) * (size_t)width + (size_t)(j)]

static inline void la_eigen_cdiv(double xr, double xi, double yr, double yi,
                                double *real, double *imag) {
    double r, denominator;
    if (fabs(yr) > fabs(yi)) {
        r = yi / yr; denominator = yr + r * yi;
        *real = (xr + r * xi) / denominator;
        *imag = (xi - r * xr) / denominator;
    } else {
        r = yr / yi; denominator = yi + r * yr;
        *real = (r * xr + xi) / denominator;
        *imag = (r * xi - xr) / denominator;
    }
}

static inline void la_eigen_orthes(int width, double *hessenberg, double *v, double *ort) {
    const int n = width;


      //  This is derived from the Algol procedures orthes and ortran,
      //  by Martin and Wilkinson, Handbook for Auto. Comp.,
      //  Vol.ii-Linear Algebra, and the corresponding
      //  Fortran subroutines in EISPACK.

      int low = 0;
      int high = n-1;

      for (int m = low+1; m <= high-1; m++) {

         // Scale column.

         double scale = 0.0;
         for (int i = m; i <= high; i++) {
            scale = scale + fabs(H(i,m-1));
         }
         if (scale != 0.0) {

            // Compute Householder transformation.

            double h = 0.0;
            for (int i = high; i >= m; i--) {
               ort[i] = H(i,m-1)/scale;
               h += ort[i] * ort[i];
            }
            double g = sqrt(h);
            if (ort[m] > 0) {
               g = -g;
            }
            h = h - ort[m] * g;
            ort[m] = ort[m] - g;

            // Apply Householder similarity transformation
            // H = (I-u*u'/h)*H*(I-u*u')/h)

            for (int j = m; j < n; j++) {
               double f = 0.0;
               for (int i = high; i >= m; i--) {
                  f += ort[i]*H(i,j);
               }
               f = f/h;
               for (int i = m; i <= high; i++) {
                  H(i,j) -= f*ort[i];
               }
           }

           for (int i = 0; i <= high; i++) {
               double f = 0.0;
               for (int j = high; j >= m; j--) {
                  f += ort[j]*H(i,j);
               }
               f = f/h;
               for (int j = m; j <= high; j++) {
                  H(i,j) -= f*ort[j];
               }
            }
            ort[m] = scale*ort[m];
            H(m,m-1) = scale*g;
         }
      }

      // Accumulate transformations (Algol's ortran).

      for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
            V(i,j) = (i == j ? 1.0 : 0.0);
         }
      }

      for (int m = high-1; m >= low+1; m--) {
         if (H(m,m-1) != 0.0) {
            for (int i = m+1; i <= high; i++) {
               ort[i] = H(i,m-1);
            }
            for (int j = m; j <= high; j++) {
               double g = 0.0;
               for (int i = m; i <= high; i++) {
                  g += ort[i] * V(i,j);
               }
               // Double division avoids possible underflow
               g = (g / ort[m]) / H(m,m-1);
               for (int i = m; i <= high; i++) {
                  V(i,j) += g * ort[i];
               }
            }
         }
      }

}

/* Return 1 on iteration exhaustion; zero means QR/backsubstitution completed. */
static inline int la_eigen_hqr2(int width, double *hessenberg, double *v, double *d, double *e,
                               size_t max_iterations) {


      //  This is derived from the Algol procedure hqr2,
      //  by Martin and Wilkinson, Handbook for Auto. Comp.,
      //  Vol.ii-Linear Algebra, and the corresponding
      //  Fortran subroutine in EISPACK.

      // Initialize

      const int nn = width;
      double cdivr = 0, cdivi = 0;
      int n = nn-1;
      int low = 0;
      int high = nn-1;
      double eps = DBL_EPSILON;
      double exshift = 0.0;
      double p=0,q=0,r=0,s=0,z=0,t,w,x,y;

      // Store roots isolated by balanc and compute matrix norm

      double norm = 0.0;
      for (int i = 0; i < nn; i++) {
         if (i < low || i > high) {
            d[i] = H(i,i);
            e[i] = 0.0;
         }
         for (int j = fmax(i-1,0); j < nn; j++) {
            norm = norm + fabs(H(i,j));
         }
      }

      // Outer loop over eigenvalue index

      int iter = 0;
      while (n >= low) {

         // Look for single small sub-diagonal element

         int l = n;
         while (l > low) {
            s = fabs(H(l-1,l-1)) + fabs(H(l,l));
            if (s == 0.0) {
               s = norm;
            }
            if (fabs(H(l,l-1)) <= eps * s) {
               break;
            }
            l--;
         }

         // Check for convergence
         // One root found

         if (l == n) {
            H(n,n) = H(n,n) + exshift;
            d[n] = H(n,n);
            e[n] = 0.0;
            n--;
            iter = 0;

         // Two roots found

         } else if (l == n-1) {
            w = H(n,n-1) * H(n-1,n);
            p = (H(n-1,n-1) - H(n,n)) / 2.0;
            q = p * p + w;
            z = sqrt(fabs(q));
            H(n,n) = H(n,n) + exshift;
            H(n-1,n-1) = H(n-1,n-1) + exshift;
            x = H(n,n);

            // Real pair

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
               x = H(n,n-1);
               s = fabs(x) + fabs(z);
               p = x / s;
               q = z / s;
               r = sqrt(p * p+q * q);
               p = p / r;
               q = q / r;

               // Row modification

               for (int j = n-1; j < nn; j++) {
                  z = H(n-1,j);
                  H(n-1,j) = q * z + p * H(n,j);
                  H(n,j) = q * H(n,j) - p * z;
               }

               // Column modification

               for (int i = 0; i <= n; i++) {
                  z = H(i,n-1);
                  H(i,n-1) = q * z + p * H(i,n);
                  H(i,n) = q * H(i,n) - p * z;
               }

               // Accumulate transformations

               for (int i = low; i <= high; i++) {
                  z = V(i,n-1);
                  V(i,n-1) = q * z + p * V(i,n);
                  V(i,n) = q * V(i,n) - p * z;
               }

            // Complex pair

            } else {
               d[n-1] = x + p;
               d[n] = x + p;
               e[n-1] = z;
               e[n] = -z;
            }
            n = n - 2;
            iter = 0;

         // No convergence yet

         } else {

            // Form shift

            x = H(n,n);
            y = 0.0;
            w = 0.0;
            if (l < n) {
               y = H(n-1,n-1);
               w = H(n,n-1) * H(n-1,n);
            }

            // Wilkinson's original ad hoc shift

            if (iter == 10) {
               exshift += x;
               for (int i = low; i <= n; i++) {
                  H(i,i) -= x;
               }
               s = fabs(H(n,n-1)) + fabs(H(n-1,n-2));
               x = y = 0.75 * s;
               w = -0.4375 * s * s;
            }

            // MATLAB's new ad hoc shift

            if (iter == 30) {
                s = (y - x) / 2.0;
                s = s * s + w;
                if (s > 0) {
                    s = sqrt(s);
                    if (y < x) {
                       s = -s;
                    }
                    s = x - w / ((y - x) / 2.0 + s);
                    for (int i = low; i <= n; i++) {
                       H(i,i) -= s;
                    }
                    exshift += s;
                    x = y = w = 0.964;
                }
            }

            if ((size_t)iter >= max_iterations) return 1;
            iter = iter + 1;

            // Look for two consecutive small sub-diagonal elements

            int m = n-2;
            while (m >= l) {
               z = H(m,m);
               r = x - z;
               s = y - z;
               p = (r * s - w) / H(m+1,m) + H(m,m+1);
               q = H(m+1,m+1) - z - r - s;
               r = H(m+2,m+1);
               s = fabs(p) + fabs(q) + fabs(r);
               p = p / s;
               q = q / s;
               r = r / s;
               if (m == l) {
                  break;
               }
               if (fabs(H(m,m-1)) * (fabs(q) + fabs(r)) <
                  eps * (fabs(p) * (fabs(H(m-1,m-1)) + fabs(z) +
                  fabs(H(m+1,m+1))))) {
                     break;
               }
               m--;
            }

            for (int i = m+2; i <= n; i++) {
               H(i,i-2) = 0.0;
               if (i > m+2) {
                  H(i,i-3) = 0.0;
               }
            }

            // Double QR step involving rows l:n and columns m:n


            for (int k = m; k <= n-1; k++) {
               int notlast = (k != n-1);
               if (k != m) {
                  p = H(k,k-1);
                  q = H(k+1,k-1);
                  r = (notlast ? H(k+2,k-1) : 0.0);
                  x = fabs(p) + fabs(q) + fabs(r);
                  if (x == 0.0) {
                      continue;
                  }
                  p = p / x;
                  q = q / x;
                  r = r / x;
               }

               s = sqrt(p * p + q * q + r * r);
               if (p < 0) {
                  s = -s;
               }
               if (s != 0) {
                  if (k != m) {
                     H(k,k-1) = -s * x;
                  } else if (l != m) {
                     H(k,k-1) = -H(k,k-1);
                  }
                  p = p + s;
                  x = p / s;
                  y = q / s;
                  z = r / s;
                  q = q / p;
                  r = r / p;

                  // Row modification

                  for (int j = k; j < nn; j++) {
                     p = H(k,j) + q * H(k+1,j);
                     if (notlast) {
                        p = p + r * H(k+2,j);
                        H(k+2,j) = H(k+2,j) - p * z;
                     }
                     H(k,j) = H(k,j) - p * x;
                     H(k+1,j) = H(k+1,j) - p * y;
                  }

                  // Column modification

                  for (int i = 0; i <= fmin(n,k+3); i++) {
                     p = x * H(i,k) + y * H(i,k+1);
                     if (notlast) {
                        p = p + z * H(i,k+2);
                        H(i,k+2) = H(i,k+2) - p * r;
                     }
                     H(i,k) = H(i,k) - p;
                     H(i,k+1) = H(i,k+1) - p * q;
                  }

                  // Accumulate transformations

                  for (int i = low; i <= high; i++) {
                     p = x * V(i,k) + y * V(i,k+1);
                     if (notlast) {
                        p = p + z * V(i,k+2);
                        V(i,k+2) = V(i,k+2) - p * r;
                     }
                     V(i,k) = V(i,k) - p;
                     V(i,k+1) = V(i,k+1) - p * q;
                  }
               }  // (s != 0)
            }  // k loop
         }  // check convergence
      }  // while (n >= low)

      // Backsubstitute to find vectors of upper triangular form

      if (norm == 0.0) {
         return 0;
      }

      for (n = nn-1; n >= 0; n--) {
         p = d[n];
         q = e[n];

         // Real vector

         if (q == 0) {
            int l = n;
            H(n,n) = 1.0;
            for (int i = n-1; i >= 0; i--) {
               w = H(i,i) - p;
               r = 0.0;
               for (int j = l; j <= n; j++) {
                  r = r + H(i,j) * H(j,n);
               }
               if (e[i] < 0.0) {
                  z = w;
                  s = r;
               } else {
                  l = i;
                  if (e[i] == 0.0) {
                     if (w != 0.0) {
                        H(i,n) = -r / w;
                     } else {
                        H(i,n) = -r / (eps * norm);
                     }

                  // Solve real equations

                  } else {
                     x = H(i,i+1);
                     y = H(i+1,i);
                     q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
                     t = (x * s - z * r) / q;
                     H(i,n) = t;
                     if (fabs(x) > fabs(z)) {
                        H(i+1,n) = (-r - w * t) / x;
                     } else {
                        H(i+1,n) = (-s - y * t) / z;
                     }
                  }

                  // Overflow control

                  t = fabs(H(i,n));
                  if ((eps * t) * t > 1) {
                     for (int j = i; j <= n; j++) {
                        H(j,n) = H(j,n) / t;
                     }
                  }
               }
            }

         // Complex vector

         } else if (q < 0) {
            int l = n-1;

            // Last vector component imaginary so matrix is triangular

            if (fabs(H(n,n-1)) > fabs(H(n-1,n))) {
               H(n-1,n-1) = q / H(n,n-1);
               H(n-1,n) = -(H(n,n) - p) / H(n,n-1);
            } else {
               la_eigen_cdiv(0.0,-H(n-1,n),H(n-1,n-1)-p,q, &cdivr, &cdivi);
               H(n-1,n-1) = cdivr;
               H(n-1,n) = cdivi;
            }
            H(n,n-1) = 0.0;
            H(n,n) = 1.0;
            for (int i = n-2; i >= 0; i--) {
               double ra,sa,vr,vi;
               ra = 0.0;
               sa = 0.0;
               for (int j = l; j <= n; j++) {
                  ra = ra + H(i,j) * H(j,n-1);
                  sa = sa + H(i,j) * H(j,n);
               }
               w = H(i,i) - p;

               if (e[i] < 0.0) {
                  z = w;
                  r = ra;
                  s = sa;
               } else {
                  l = i;
                  if (e[i] == 0) {
                     la_eigen_cdiv(-ra,-sa,w,q, &cdivr, &cdivi);
                     H(i,n-1) = cdivr;
                     H(i,n) = cdivi;
                  } else {

                     // Solve complex equations

                     x = H(i,i+1);
                     y = H(i+1,i);
                     vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
                     vi = (d[i] - p) * 2.0 * q;
                     if (vr == 0.0 && vi == 0.0) {
                        vr = eps * norm * (fabs(w) + fabs(q) +
                        fabs(x) + fabs(y) + fabs(z));
                     }
                     la_eigen_cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi, &cdivr, &cdivi);
                     H(i,n-1) = cdivr;
                     H(i,n) = cdivi;
                     if (fabs(x) > (fabs(z) + fabs(q))) {
                        H(i+1,n-1) = (-ra - w * H(i,n-1) + q * H(i,n)) / x;
                        H(i+1,n) = (-sa - w * H(i,n) - q * H(i,n-1)) / x;
                     } else {
                        la_eigen_cdiv(-r-y*H(i,n-1),-s-y*H(i,n),z,q, &cdivr, &cdivi);
                        H(i+1,n-1) = cdivr;
                        H(i+1,n) = cdivi;
                     }
                  }

                  // Overflow control

                  t = fmax(fabs(H(i,n-1)),fabs(H(i,n)));
                  if ((eps * t) * t > 1) {
                     for (int j = i; j <= n; j++) {
                        H(j,n-1) = H(j,n-1) / t;
                        H(j,n) = H(j,n) / t;
                     }
                  }
               }
            }
         }
      }

      // Vectors of isolated roots

      for (int i = 0; i < nn; i++) {
         if (i < low || i > high) {
            for (int j = i; j < nn; j++) {
               V(i,j) = H(i,j);
            }
         }
      }

      // Back transformation to get eigenvectors of original matrix

      for (int j = nn-1; j >= low; j--) {
         for (int i = low; i <= high; i++) {
            z = 0.0;
            for (int k = low; k <= fmin(j,high); k++) {
               z = z + V(i,k) * H(k,j);
            }
            V(i,j) = z;
         }
      }

    return 0;
}

#undef H
#undef V

/* Power-of-two diagonal similarity balancing, B = D^-1 A D. Norms are
 * accumulated relative to a maximum, then compared in log space, so opposing
 * finite entries near the limits of double range are not multiplied together.
 * Reject any proposed step that would overflow or erase a nonzero entry. */
static inline void la_eigen_balance(size_t n, double *h, double *powers) {
    for (size_t pass = 0; pass < 64; pass++) {
        int changed = 0;
        for (size_t i = 0; i < n; i++) {
            double row_max = 0, col_max = 0;
            for (size_t j = 0; j < n; j++) if (j != i) {
                row_max = fmax(row_max, fabs(h[i*n+j]));
                col_max = fmax(col_max, fabs(h[j*n+i]));
            }
            if (row_max == 0 || col_max == 0) continue;
            double row_sum = 0, col_sum = 0;
            for (size_t j = 0; j < n; j++) if (j != i) {
                row_sum += fabs(h[i*n+j]) / row_max;
                col_sum += fabs(h[j*n+i]) / col_max;
            }
            const double row_log = log2(row_max) + log2(row_sum);
            const double col_log = log2(col_max) + log2(col_sum);
            int shift = (int)floor((row_log-col_log) / 2.0 + 0.5);
            if (shift > 512) shift = 512;
            if (shift < -512) shift = -512;
            if (shift == 0) continue;
            const double common = fmax(row_log, col_log);
            const double before = exp2(row_log-common) + exp2(col_log-common);
            const double after = exp2(row_log-shift-common) + exp2(col_log+shift-common);
            if (!(after < 0.95 * before)) continue;
            int safe = 1;
            for (size_t j = 0; j < n; j++) if (j != i) {
                const double row = scalbn(h[i*n+j], -shift), col = scalbn(h[j*n+i], shift);
                if (!isfinite(row) || !isfinite(col) ||
                    (row == 0 && h[i*n+j] != 0) || (col == 0 && h[j*n+i] != 0) ||
                    scalbn(row, shift) != h[i*n+j] || scalbn(col, -shift) != h[j*n+i]) safe = 0;
            }
            if (!safe) continue;
            for (size_t j = 0; j < n; j++) if (j != i) {
                h[i*n+j] = scalbn(h[i*n+j], -shift);
                h[j*n+i] = scalbn(h[j*n+i], shift);
            }
            powers[i] += shift;
            changed = 1;
        }
        if (!changed) break;
    }
}

/* Buffers: d/e/ort/powers have width entries; all others width*width entries.
 * Return 0 for success, 1 for nonconvergence, 2 for nonfinite arithmetic.
 * The caller validates finite input, dimensions, iteration budget, and memory. */
static inline int la_general_eigen(int width, const double *source, size_t max_iterations,
                                  double *d, double *e, double *real, double *imag,
                                  double *h, double *v, double *ort, double *powers) {
    const size_t n = (size_t)width;
    double maximum = 0;
    int diagonal = 1, upper = 1, lower = 1;
    for (size_t i = 0; i < n * n; i++) {
        maximum = fmax(maximum, fabs(source[i]));
        if (i / n != i % n && source[i] != 0) diagonal = 0;
        if (i / n > i % n && source[i] != 0) upper = 0;
        if (i / n < i % n && source[i] != 0) lower = 0;
    }
    if (diagonal) {
        for (size_t i = 0; i < n; i++) { d[i] = source[i * n + i]; real[i * n + i] = 1; }
    } else {
        const int reverse = lower && !upper;
        for (size_t row = 0; row < n; row++) for (size_t col = 0; col < n; col++)
            h[row*n+col] = reverse ? source[(n-1-row)*n+n-1-col] : source[row*n+col];
        if (!upper && !lower) la_eigen_balance(n, h, powers);
        maximum = 0;
        for (size_t i = 0; i < n * n; i++) maximum = fmax(maximum, fabs(h[i]));
        int exponent = 0;
        (void)frexp(maximum, &exponent);
        for (size_t i = 0; i < n * n; i++) h[i] = scalbn(h[i], -exponent);
        la_eigen_orthes(width, h, v, ort);
        if (la_eigen_hqr2(width, h, v, d, e, max_iterations)) return 1;
        for (size_t col = 0; col < n; col++) {
            const int positive = e[col] > 0, negative = e[col] < 0;
            if ((positive && col + 1 >= n) || (negative && col == 0)) return 2;
            for (size_t row = 0; row < n; row++) {
                const size_t input_row = reverse ? n - 1 - row : row;
                real[row * n + col] = v[input_row * n + (negative ? col - 1 : col)];
                imag[row * n + col] = positive ? v[input_row * n + col + 1] : negative ? -v[input_row * n + col] : 0;
            }
        }
        for (size_t col = 0; col < n; col++) {
            /* Triangular eigenvalues are known exactly before work scaling. */
            const size_t original = reverse ? n - 1 - col : col;
            d[col] = upper || lower ? source[original*n+original] : scalbn(d[col], exponent);
            e[col] = scalbn(e[col], exponent);
            if (!isfinite(d[col]) || !isfinite(e[col])) return 2;
            /* Undo D without overflowing: a common exponent per column cancels
             * during normalization, retaining representable tiny components. */
            int common_exponent = -2147483647;
            for (size_t row = 0; row < n; row++) {
                const size_t index = row * n + col;
                const double component = fmax(fabs(real[index]), fabs(imag[index]));
                if (!isfinite(component)) return 2;
                if (component != 0) {
                    int part_exponent = 0; (void)frexp(component, &part_exponent);
                    const int combined = part_exponent + (int)powers[row];
                    if (combined > common_exponent) common_exponent = combined;
                }
            }
            if (common_exponent == -2147483647) return 2;
            for (size_t row = 0; row < n; row++) {
                const size_t index = row * n + col;
                const int shift = (int)powers[row] - common_exponent;
                real[index] = scalbn(real[index], shift); imag[index] = scalbn(imag[index], shift);
            }
            double scale = 0, norm = 0;
            for (size_t row = 0; row < n; row++) {
                const size_t index = row * n + col;
                if (!isfinite(real[index]) || !isfinite(imag[index])) return 2;
                scale = fmax(scale, fmax(fabs(real[index]), fabs(imag[index])));
            }
            if (scale == 0) return 2;
            for (size_t row = 0; row < n; row++) {
                const size_t index = row * n + col;
                real[index] /= scale; imag[index] /= scale;
                norm = hypot(norm, hypot(real[index], imag[index]));
            }
            size_t largest = 0;
            for (size_t row = 0; row < n; row++) {
                const size_t index = row * n + col;
                real[index] /= norm; imag[index] /= norm;
                if (hypot(real[index], imag[index]) > hypot(real[largest * n + col], imag[largest * n + col])) largest = row;
            }
            /* Fix the otherwise arbitrary complex phase at the largest component. */
            const double magnitude = hypot(real[largest * n + col], imag[largest * n + col]);
            const double phase_real = real[largest * n + col] / magnitude;
            const double phase_imag = imag[largest * n + col] / magnitude;
            for (size_t row = 0; row < n; row++) {
                const size_t index = row * n + col;
                const double r = real[index], i = imag[index];
                real[index] = r * phase_real + i * phase_imag;
                imag[index] = i * phase_real - r * phase_imag;
                if (!isfinite(real[index]) || !isfinite(imag[index])) return 2;
            }
        }
    }
    for (size_t col = 0; col < n; col++) {
        size_t first = col;
        for (size_t j = col + 1; j < n; j++)
            if (d[j] < d[first] || (d[j] == d[first] && e[j] < e[first])) first = j;
        if (first != col) {
            double swap = d[col]; d[col] = d[first]; d[first] = swap;
            swap = e[col]; e[col] = e[first]; e[first] = swap;
            for (size_t row = 0; row < n; row++) {
                swap = real[row * n + col]; real[row * n + col] = real[row * n + first]; real[row * n + first] = swap;
                swap = imag[row * n + col]; imag[row * n + col] = imag[row * n + first]; imag[row * n + first] = swap;
            }
        }
    }
    return 0;
}
#endif
