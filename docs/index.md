- Semi-analytical alternative to finite difference methods which yield high degree of accuracy
- Canonically, restricted to periodic data but transform methods based on Chebyshev polynomial expansions can be used for non-periodic data

# Continuous Fourier series

A continuous periodic function $f(x)$ has a Fourier series representation:

$$f(x) = \sum^\infty_{k=-\infty} \hat{f}_k e^{ikx}$$
where $\hat{f}_k$ are the Fourier coefficients for wavenumber $k$, which are integer for a period of $2\pi$. Physically, $\hat{f}_k$ contains information on how much each harmonic contributes to $f(x)$. The Fourier series to the derivative, $f^\prime(x)$ results from straightforward differentiation:
$$f^\prime(x) = \sum_{k=-\infty}^{\infty}ik \hat{f}_k e^{ikx}$$
such that the Fourier coefficients of $f^\prime$ are simply $ik \hat{f}_k$.

# Discrete Fourier series

Suppose we have the periodic function $f_j \equiv f(x_j)$ for the $j = 0,\dots,N-1$ discrete grid points $x_0,\dots,x_{N-1}$. The DFT (also known as the forward transformation) of $f_0,\dots f_{N-1}$ is defined by:
$$f_j = \sum_{k=-N/2}^{N/2-1}\hat{f}_ke^{ikx_j}, \qquad j = 0,\dots,N-1$$
with discrete Fourier coefficients:
$$\hat{f}_{-N/2}, \hat{f}_{-N/2+1}, \dots, 0, \dots \hat{f}_{N/2-1}$$

Important notes:
- Assuming $N$ to be an **even** number here
- Assuming $f$ is periodic in $2\pi$, resulting in integer wavenumbers. Results in $f_0 = f_N$ such that there is no redundancy with sequence $f_0, \dots f_{N-1}$
- Grid spacing $h = 2\pi/N$ with $x_j = jh$.
- If more generally, the periodic was length $L$, then replace $k$ with $(2 \pi/L)k$ in above expressions, and $h = L/N$. Importantly, even with these substitutions, the DFT definition stays the same (although it does change for the derivative expression)

We have $N$ algebraic equations for the $N$ unknown Fourier coefficients, $\hat{f}_k$. Solve these by exploiting discrete orthogonality rather than, say, a usual linear-system solver method.

Consider the summation, $I$:

$$I = \sum_{j=0}^{N-1}e^{ikx_j}e^{-ik^\prime x_j} = \sum_{j=0}^{N-1}e^{ih(k-k^\prime)j}$$

If $h(k-k^\prime)$ is not a multiple of $2\pi$ then $I$ is like a geometric series. A geometric series is given by $I = \sum_{j=0}^{N-1}\alpha^j = (1-\alpha^N)/(1-\alpha)$  for $\alpha \neq 1$. So here, $\alpha \equiv e^{ih(k-k^\prime)}$ and thus the $h(k-k^\prime) \neq 1$ condition. The summation becomes:
$$I = \frac{1 - e^{ih(k - k^\prime)N}}{1 - e^{ih(k - k^\prime)}} = \frac{1 - e^{2\pi i (k - k^\prime)}}{1 - e^{2\pi i (k-k^\prime)/N}}$$
Where we used $h = 2\pi/N$. Recalling that the wavenumbers are any integers and $e^{2\pi in} = 1$ for another integer $n$, then the numerator becomes zero. Thus we have the following discrete orthogonality condition:

$$\sum_{j=0}^{N-1}e^{ikx_j}e^{-ik^\prime x_j} = \sum_{j=0}^{N-1}e^{ix_j(k - k^\prime)} =\begin{cases}
N, \qquad k = k^\prime + mN, m = 0,\pm1,\pm2,\dots \\
0, \qquad \text{otherwise}
\end{cases}$$
Now use this result to obtain an explicit expression for $\hat{f}_k$. Multiply the DFT definition by $e^{-ik^\prime x_j}$ and sum over $j=0,\dots,N-1$:
$$\sum_{j=0}^{N-1}f_j e^{-ik^\prime x_j} = \sum_{k=-N/2}^{N/2-1}\sum_{j=0}^{N-1}\hat{f}_k e^{ix_j(k-k^\prime)}$$
And invoking the orthogonality property:
$$\hat{f}_k = \frac{1}{N}\sum_{j=0}^{N-1}f_j e^{-ikx_j}, \qquad k = -\frac{N}{2},\dots,\frac{N}{2}-1$$
This defines the inverse DFT. Together with the DFT, they form the transform pair.

# Useful properties

## For real functions, the storage requirement in wavenumber space is reduced

If the function $f$ is purely real, the original $N$ $f_j$ datapoints can be represented with only $N/2$ complex Fourier coefficients. Simply set $k = -k$ in the IDFT definition and take its conjugate:
$$\hat{f}^*_{-k} = \left( \frac{1}{N}\sum_{j=0}^{N-1}f_j e^{ikx_j}\right)^* = \frac{1}{N}\sum_{j=0}^{N-1}f_j e^{-ikx_j}$$
When comparing with the original IDFT definition, we see that we simply have:
$$\hat{f}_{-k} = \hat{f}_k^*.$$

In 2D, we would have
$$\hat{f}^*_{-k_1,-k_2} = \hat{f}_{k_1,k_2}$$
Thus, **one half (not one quarter)** of the $(k_1,k_2)$ space are sufficient to determine the Fourier coefficients. 
## Odd and even functions

## Convolution

# Important concepts
## Aliasing

Consider the product $H(x) = f(x)g(x)$. Suppose we wish to represent $H$ with Fourier transforms of $f$ and $g$. Fourier coefficients of $H$ are given by:
$$\hat{H}_m = (\widehat{fg})_m = \frac{1}{N}\sum_{j=0}^{N-1}f_j g_j e^{imx_j}$$
Now sub in DFT representations of $f_j$ $g_j$:

$$\hat{H}_m = (\widehat{fg})_m = \frac{1}{N}\sum_{j=0}^{N-1}\sum_k \sum_{k^\prime} \hat{f}_k \hat{g}_{k^\prime} e^{ikx_j}e^{ik^\prime x_j}e^{imx_j}$$
The sum over $j$ is non-zero only when $k + k^\prime  = m$ or $m \pm N$ (recalling that $x_j = 2\pi Nj$). The part of the summation corresponding to $k + k^\prime = m \pm N$ is known as the **aliasing error**, and these terms should be discarded since the Fourier exponentials for these wavenumbers cannot be resolved on the grid size $N$.

Now according to the IDFT definition, this simplifies to

$$\hat{H}_m = \sum_{k = -N/2}^{N/2-1}\hat{f}_k \hat{g}_{m-k}$$
If we simply multiplied $f$ and $g$ at each grid point the transformed the result, this would not be equal to the result above due to contamination by aliasing errors.
We may circumvent the aliasing errors by computing the sum above directly, which is expensive $O(N^2)$ operations.
### Aliasing example


# Higher dimensions

Suppose we have $f(x,y)$ with $N_1$, $N_2$ grid points and it is doubly-periodic. Define $f_{m,l} \equiv f(x_m,y_l)$. The transform pairs are:
$$f(x_m,y_l) = \sum_{k_1 = -N_1/2}^{N_1/2-1}\sum_{k_2 = -N_2/2}^{N_2/2-1}\hat{f}_{k_1,k_2}e^{ik_1x_m}e^{ik_2 y_l}, \qquad m = 0,1,\dots,N_1-1, \quad l = 0,1,2,\dots,N_2-1$$
$$\hat{f}_{k_1,k_2} = \frac{1}{N_1}\frac{1}{N_2}\sum_{m=0}^{N_1}\sum_{l=0}^{N_2}f_{m,l}e^{-ik_1x_m}e^{-ik_2 y_l}, \qquad k_1 = -\frac{N_1}{2},\dots,\frac{N_1}{2}-1, \quad k_2 = -\frac{N_2}{2},\dots,\frac{N_2}{2}-1$$
# Discrete sine and cosine transforms

Not restricted to periodic signals.

Even function representations, $f(x) = f(-x)$ efficiently represented with cosines.

Odd functions, $f(x) = -f(-x)$, we can use sines.

# Using FFTW

Here, using general C/C++ syntax, and summarising from the FFTW manual. Link library with `-lfftw3` `-lm` compiler flags

Some summary notes:
- General workflow: initialise pointers to arrays, make plan, execute plan as many times as desired, free memory/destroy plan
- FFTW computes an unnormalised DFT. So FFT -> IFFT will produce the original input multiplied by the size of the transform, N
## Complex to complex example

Skeleton code for using FFTW looks like this 
```C
#include <fftw3.h>
...
{
fftw_complex *in, *out;
fftw_plan p;
...
in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
...
fftw_execute(p); /* repeat as needed */
...
fftw_destroy_plan(p);
fftw_free(in); fftw_free(out);
}
```

Steps are:
1) Allocate memory for input and output to FFTs. `fftw_malloc()` is recommended for this, can use the wrapper `fftw_alloc_complex(N)` for brevity as well
	- `fftw_complex` is a `double[2]` composed of real `in[i][0]` and imaginary `in[i][1]` components
2) Create a plan, which contains all the data needed to compute the FFT
	-  `fftw_plan fftw_plan_dft_1d(int n, fftw_complex *in, fftw_complex *out, int sign, unsigned flags);`
		- `n>0` is the size of the transform, `in` and `out` are pointers to the input/output arrays of the transformation. They can be equal, which would be an in-place transformation
		- `sign` can be either `FFT_FORWARD (-1)` or `FFT_BACKWARD (+1)`, which indicates the direction of the transform desired (sign of the exponent in the transform essentially)
		- `flags` usually set as `FFTW_MEASURE` or `FFT_ESTIMATE`. `FFTW_MEASURE` spends a few seconds calculating the most optimal way of computing the FFT for the given problem 
		- "You must create the plan before initializing the input, because `FFTW_MEASURE` overwrites the `in/out` arrays"
3) Once the plan is created, it can be used as many times as desired by `void fftw_execute(const fftw_plan plan);`
4) `void fftw_destroy_plan(fftw_plan plan);` once done with the plan.
5) "If you allocate an array with `fftw_malloc()` you must deallocate it with `fftw_free()`. Do not use `free()` or, heaven forbid, `delete`."

If we want to transform a *different* array of the same size, make another plan `fftw_plan_dft_1d` and FFTW automatically reuses information from the previous plan  if possible.

For the "guru" interface, we can assign particular plans.

## Transforms of real data

Hermitian redundancy $\hat{f}_i = \hat{f}_{N-i}$. Can take advantage of this for a $\approx$ factor of 2 improvement in memory usage and speed.

Speed/space advantage sacrificed for simplicity.
- Input/output arrays are of different *types* and *sizes*.
	- Input is of size `N` real numbers and transform result is size `N/2+1` complex numbers (division is rounded down)
- For "in-place" transforms (i.e. output memory = input memory) this requires a slight padding of the input array
- The inverse *complex to real* transform overwrites its input array by default

Snippet here for example
```C++
{
...
	int N = 8;
    int Nc = N/2 + 1;
    double* in;
    fftw_complex *out;
    fftw_plan p;
    fftw_plan invp;
    ...
    in = (double*) fftw_malloc(sizeof(double) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);
	...
    p = fftw_plan_dft_r2c_1d(N, in, out,FFTW_MEASURE);
    invp = fftw_plan_dft_c2r_1d(N, out, in,FFTW_MEASURE);

    fftw_execute(p); // Execute FFTW
    fftw_execute(invp); // Execute IFFT
	...
    fftw_destroy_plan(p);
    fftw_destroy_plan(invp);
    fftw_free(in); fftw_free(out);
    ...
}
```

- Still using `fftw_malloc()` for the `double` input
- There is no `sign` argument unlike with complex-to-complex
	- `r2c` is always `FFT_FORWARD` and `c2r` is always `FFT_BACKWARD`
# References
- [@moin2010]
- https://www.fftw.org/links.html
