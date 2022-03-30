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

    Copyright (C) 2021 Anton Mosunov
 
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

// Block-by-block initialization


// returns floor(sqrt(x))
// the code is taken from https://dox.ipxe.org/isqrt_8c.html
ulong isqrt(ulong value) {
         ulong result = 0;
         ulong bit = ( 1UL << ( ( 8 * sizeof ( bit ) ) - 2 ) );
 
         while ( bit > value )
                 bit >>= 2;
         while ( bit ) {
                 if ( value >= ( result + bit ) ) {
                         value -= ( result + bit );
                         result = ( ( result >> 1 ) + bit );
                 } else {
                         result >>= 1;
                 }
                 bit >>= 2;
         }
         return result;
 }

// returns ceil(sqrt(a/b))
long ceilsqrt(ulong a, ulong b)
{
	ulong q = (a / b);
	ulong r = a - (q * b);
	ulong s = isqrt(q);
	if (r == 0 && (s * s) == q)
	{
		return(s);
	}
	return(s + 1);
}

/*
* Initializes a block of coefficients of c*(q^r)*nabla_{a, m}(q^s) in the range [min, min + size).
* IMPORTANT: if m = 1 and a = 1, we are initializing c*(1/2)*(q^r)*nabla_{1, 1}(q^s), as opposed to c*(q^r)*nabla_{1, 1}(q^s).
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
		printf("d=%lu\n", d);
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

void init_block_nabla_product(int * block, const ulong size, const ulong min, const ulong c1, const ulong r1, const ulong s1, const ulong a1, const ulong m1, const ulong c2, const ulong r2, const ulong s2, const ulong a2, const ulong m2)
{
	ulong i = 0, x = r1, j, y, t, d, c = (c1*c2);

	for (; x < min + size; i++)
	{
		if (x >= min)
		{
			j = 0;
		}
		else
		{
			d = ceilsqrt(a2*a2*s2 + 8*m2*(min-x), s2);
			j = (d - a2) / (2*m2);

			if ((d - a2) % (2*m2) != 0)
			{
				j++;
			}
		}
		
		y = r2 + ((s2*j*(m2*j+a2))>>1);
		t = x + y - min;

		for (; t < size; j++)
		{	
			block[t] += c;
			y += ((s2*(2*j*m2 + m2 + a2)) >> 1);
			t += ((s2*(2*j*m2 + m2 + a2)) >> 1);
		}

		if ((a2 != 1) || (m2 != 1))
		{
			if (x >= min)
			{
				j = 1;
			}
			else
			{
				d = ceilsqrt(a2*a2*s2 + 8*m2*(min-x), s2);
				j = (d + a2) / (2*m2);

				if ((d + a2) % (2*m2) != 0)
				{
					j++;
				}
			}

			y = r2 + ((s2*j*(m2*j-a2))>>1);
			t = x + y - min;

			for (; t < size; j++)
			{
				block[t] += c;
				y += ((s2*(2*j*m2 + (m2 - a2))) >> 1);
				t += ((s2*(2*j*m2 + (m2 - a2))) >> 1);
			}
		}

		x += ((s1*(2*i*m1 + m1 + a1)) >> 1);		
	}

	if ((a1 != 1) || (m1 != 1))
	{
		i = 1;
		x = r1 + ((s1*(m1-a1))>>1);

		for (; x < min + size; i++)
		{
			if (x >= min)
			{
				j = 0;
			}
			else
			{
				d = ceilsqrt(a2*a2*s2 + 8*m2*(min-x), s2);
				j = (d - a2) / (2*m2);

				if ((d - a2) % (2*m2) != 0)
				{
					j++;
				}
			}

			y = r2 + ((s2*j*(m2*j+a2))>>1);
			t = x + y - min;

			for (; t < size; j++)
			{
				block[t] += c;
				y += ((s2*(2*j*m2 + m2 + a2)) >> 1);
				t += ((s2*(2*j*m2 + m2 + a2)) >> 1);
			}

			if ((a2 != 1) || (m2 != 1))
			{
				if (x >= min)
				{       
					j = 1;
				}
				else
				{       
					d = ceilsqrt(a2*a2*s2 + 8*m2*(min-x), s2);
					j = (d - a2) / (2*m2);
    
					if ((d - a2) % (2*m2) != 0)
					{       
						j++;
					}
				}

				y = r2 + ((s2*j*(m2*j-a2))>>1);
				t = x + y - min;

				for (; t < size; j++)
				{
					block[t] += c;
					y += ((s2*(2*j*m2 + (m2 - a2))) >> 1);
					t += ((s2*(2*j*m2 + (m2 - a2))) >> 1);
				}
			}

			x += ((s1*(2*i*m1 + (m1 - a1))) >> 1);		
		}
	}
}
