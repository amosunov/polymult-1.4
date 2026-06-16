/*=============================================================================

    This file is part of POLYMULT.

    POLYMULT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    POLYMULT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with POLYMULT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2026 Anton Mosunov
 
******************************************************************************/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>

#include <flint/nmod_poly.h>

#include "init.h"





// returns ceil(sqrt(a/b))
// assumes that b is not zero
static ulong ceilsqrt(ulong a, ulong b)
{
    if (a == 0)
        return 0;

    ulong lo = 0;
    ulong hi = ((ulong) (-1));

    while (lo < hi) {
        ulong mid = lo + (hi - lo) / 2;

        __uint128_t lhs = (__uint128_t)mid * mid * b;

        if (lhs >= a)
            hi = mid;
        else
            lo = mid + 1;
    }

    return lo;
}





// returns the smallest x ≡ a (mod m) such that x >= min 
// assumes that m is positive
static ulong smallest_a_mod_m(const long a, const ulong m, const ulong min)
{
    /* normalize a into the range [0, m-1] */
    long am = a % (long)m;
    ulong r = (am < 0) ? (ulong)(am + (long)m) : (ulong)am;

    if (min <= r)
        return r;

    /* k = ceil((bound - r) / m) */
    ulong d = min - r;
    ulong k = (d + m - 1) / m;

    return r + k * m;
}





/*
* Initializes a block of coefficients of c*(q^r)*nabla_{a, m}(q^s) in the range [min, min + size)
* IMPORTANT: if m = 1 and a = 1, initializes (1/2)*c*(q^r)*nabla_{1, 1}(q^s), as opposed to c*(q^r)*nabla_{1, 1}(q^s)
* Must have gcd(a,m)=1, 0<=a<=m, and s*k*(m*k+a)/2 a non-zero integer for every integer k in order to work properly
*/
void init_block_nabla(int * block, const ulong size, const ulong min, const ulong c, const ulong r, const ulong s, const ulong a, const ulong m)
{
	ulong i, d;

	if (r >= min)
	{
		i = 0;
	}
	else
	{
		d = ceilsqrt(a*a*s + 8*m*(min-r), s);
		
        i = (d - a) / (2*m);

		if ((d - a) % (2*m) != 0)
		{
			i++;
		}
	}

	ulong t = r + ((s * i * (m * i + a)) >> 1) - min;

	for (; t < size; i++)
	{
		block[t] += c;
		t += ((s*(2*i*m + m + a)) >> 1);
	}

	if ((a != 1) || (m != 1))
	{
		if (r >= min)
		{
			i = 1;
		}
		else
		{
			i = (d + a) / (2 * m);

			if ((d + a) % (2*m) != 0)
			{
				i++;
			}
		}

		t = r + ((s * i * (m * i - a)) >> 1) - min;

		for (; t < size; i++)
		{
			block[t] += c;
			t += ((s*(2*i*m + (m - a))) >> 1);
		}
	}
}


/*
* Initializes a block of coefficients of c1*(q^r1)*nabla_{a1, m1}(q^s1) x c2*(q^r2)*nabla_{a2, m2}(q^s2) in the range [min, min + size)
* IMPORTANT: if m = 1 and a = 1, initializes (1/2)*c*(q^r)*nabla_{1, 1}(q^s), as opposed to c*(q^r)*nabla_{1, 1}(q^s)
* Must have gcd(a1,m1)=1, gcd(a2,m2)=1, 0<=a1<=m1, 0<=2<=m2, and s*k*(m*k+a)/2 a non-negative integer for every integer k in order to work properly
*/
void init_block_nabla_product(int * block, const ulong size, const ulong min, const ulong c1, const ulong r1, const ulong s1, const ulong a1, const ulong m1, const ulong c2, const ulong r2, const ulong s2, const ulong a2, const ulong m2) {
    const int c = (int) (c1 * c2);

    const ulong M = min;
    const ulong N = min + size - 1;

    ulong A = 8*m1*m2*M + s1*a1*a1*m2 + s2*a2*a2*m1;
    ulong B = 8*m1*m2*N + s1*a1*a1*m2 + s2*a2*a2*m1;
    const ulong C = 8*m1*m2*(r1 + r2);

    A = (A >= C) ? (A - C) : 0;
    B = (B >= C) ? (B - C) : 0;

    #ifdef DEBUG
    printf("c1=%lu, r1=%lu, s1=%lu, a1=%lu, m1=%lu\n", c1, r1, s1, a1, m1);
    printf("c2=%lu, r2=%lu, s2=%lu, a2=%lu, m2=%lu\n", c2, r2, s2, a2, m2);
    printf("A=%lu, B=%lu\n", A, B);
    #endif

    long b1 = (-1) * ((long) a1);
    long b2 = (-1) * ((long) a2);

    ulong square1 = smallest_a_mod_m(a1, 2*m1, 0);
    ulong k1 = (square1 - a1)/(2*m1);
    square1 *= square1;

    ulong k2, square2, index;

    while (square1 <= B/(s1*m2)) {
        square2 = smallest_a_mod_m(a2, 2*m2, (A >= s1*m2*square1) ? ceilsqrt(A - s1*m2*square1, s2*m1) : 0);
        k2 = (square2 - a2)/(2*m2);
        square2 *= square2;

        #ifdef DEBUG
        printf("CYCLE 1.1\n");
        printf("initial_square1=%lu, initial_square2=%lu\n", square1, square2);
        #endif

        index = (r1 + r2 + (s1*(square1 - a1*a1))/(8*m1) + (s2*(square2 - a2*a2))/(8*m2)) - M;
        while (index < size) {
            #ifdef DEBUG
            printf("tr1=%lu, tr2=%lu\n", (s1*(square1 - a1*a1))/(8*m1),  (s2*(square2 - a2*a2))/(8*m2)) ;
            printf("\tsquare1=%lu, square2=%lu, index=%lu\n", square1, square2, index);
            square2 += 4*m2*((2*k2+1)*m2 + a2);
            #endif
            block[index] += c;
            index += (s2*((2*k2+1)*m2+a2))/2;
            k2++;
        }

        if (a2 != m2) { // Notice that a2=m2 iff a2=m2=1. When a2=m2=1, we are generating sum q^(k*(k+1)/2) and not 2*sum q^(k*(k+1)/1)
            square2 = smallest_a_mod_m(b2, 2*m2, (A >= s1*m2*square1) ? ceilsqrt(A - s1*m2*square1, s2*m1) : 0);
            if (square2 == 0) square2 += 2; // Notice that square2=0 iff a2=0 and m2=1. This corresponds to k2=0, which should be counted only once
            k2 = (square2 + (-b2))/(2*m2);
            square2 *= square2;

            #ifdef DEBUG
            printf("CYCLE 1.2\n");
            printf("initial_square1=%lu, initial_square2=%lu\n", square1, square2);
            #endif

            index = (r1 + r2 + (s1*(square1 - a1*a1))/(8*m1) + (s2*(square2 - (-b2)*(-b2)))/(8*m2)) - M;
            while (index < size) {
                #ifdef DEBUG
                printf("tr1=%lu, tr2=%lu\n", (s1*(square1 - a1*a1))/(8*m1),  (s2*(square2 - b2*b2))/(8*m2)) ;
                printf("\tsquare1=%lu, square2=%lu, index=%lu\n", square1, square2, index);
                square2 += 4*m2*((2*k2+1)*m2 - (-b2));
                #endif
                block[index] += c;
                index += (s2*((2*k2+1)*m2-(-b2)))/2;
                k2++;
            }
        }

        square1 += 4*m1*((2*k1+1)*m1 + a1);
        k1++;
    }

    if (a1 == m1) return; // Notice that a1=m1 iff a1=m1=1. When a2=m2=1, we are generating sum q^(k*(k+1)/2) and not 2*sum q^(k*(k+1)/2)

    square1 = smallest_a_mod_m(b1, 2*m1, 0);
    if (square1 == 0) square1 += 2; // Notice that square1=0 iff a1=0 and m1=1. This corresponds to k1=0, which should be counted only once
    k1 = (square1 + (-b1))/(2*m1);
    square1 *= square1;

    while (square1 <= B/(s1*m2)) {
        square2 = smallest_a_mod_m(a2, 2*m2, (A >= s1*m2*square1) ? ceilsqrt(A - s1*m2*square1, s2*m1) : 0);
        k2 = (square2 - a2)/(2*m2);
        square2 *= square2;

        #ifdef DEBUG
        printf("CYCLE 2.1\n");
        printf("initial_square1=%lu, initial_square2=%lu\n", square1, square2);
        #endif

        index = (r1 + r2 + (s1*(square1 - (-b1)*(-b1)))/(8*m1) + (s2*(square2 - a2*a2))/(8*m2)) - M;
        while (index < size) {
            #ifdef DEBUG
            printf("tr1=%lu, tr2=%lu\n", (s1*(square1 - b1*b1))/(8*m1),  (s2*(square2 - a2*a2))/(8*m2)) ;
            printf("\tsquare1=%lu, square2=%lu, index=%lu\n", square1, square2, index);
            square2 += 4*m2*((2*k2+1)*m2 + a2);
            #endif
            block[index] += c;
            index += (s2*((2*k2+1)*m2+a2))/2;
            k2++;
        }

        if (a2 != m2) { // Notice that a2=m2 iff a2=m2=1. When a2=m2=1, we are generating sum q^(k*(k+1)/2) and not 2*sum q^(k*(k+1)/2)
            square2 = smallest_a_mod_m(b2, 2*m2, (A >= s1*m2*square1) ? ceilsqrt(A - s1*m2*square1, s2*m1) : 0);
            if (square2 == 0) square2 += 2; // Notice that square2=0 iff a2=0 and m2=1. This corresponds to k2=0, which should be counted only once
            k2 = (square2 + (-b2))/(2*m2);
            square2 *= square2;

            #ifdef DEBUG
            printf("CYCLE 2.2\n");
            printf("initial_square1=%lu, initial_square2=%lu\n", square1, square2);
            #endif

            index = (r1 + r2 + (s1*(square1 - (-b1)*(-b1)))/(8*m1) + (s2*(square2 - (-b2)*(-b2)))/(8*m2))  - M;
            while (index < size) {
                #ifdef DEBUG
                printf("tr1=%lu, tr2=%lu\n", (s1*(square1 - b1*b1))/(8*m1),  (s2*(square2 - b2*b2))/(8*m2)) ;
                printf("\tsquare1=%lu, square2=%lu, index=%lu\n", square1, square2, index);
                square2 += 4*m2*((2*k2+1)*m2 - (-b2));
                #endif
                block[index] += c;
                index += (s2*((2*k2+1)*m2-(-b2)))/2;
                k2++;
            }
        }

        square1 += 4*m1*((2*k1+1)*m1 - (-b1));
        k1++;
    }
}
