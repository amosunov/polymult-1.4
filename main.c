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

#include "mult.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#ifndef uint
#define uint unsigned int
#endif

#ifndef ulong
#define ulong unsigned long
#endif

int main(int argc, char * argv[])
{

	//omp_set_num_threads(4);

	if ((argc < 22) || (argc % 5 != 2))
	{
		printf("Format: ./polymult [limit] [files] [bundle] [bound] [resultname] [folder] [c0] [r0] [s0] [a0] [m0] [c1] [r1] [s1] [a1] [m1] ...\n");
		exit(1);
	}

	const int series_total = (argc - 7)/5;
	ulong * series_params = malloc(sizeof(ulong) * series_total * 5);	

	printf("series_total = %lu\n", series_total);

	for (int i = 0; i < series_total; i++)
	{
		series_params[5*i] = atol(argv[7 + 5*i]);	
		series_params[5*i + 1] = atol(argv[7 + 5*i + 1]);
		series_params[5*i + 2] = atol(argv[7 + 5*i + 2]);
		series_params[5*i + 3] = atol(argv[7 + 5*i + 3]);
		series_params[5*i + 4] = atol(argv[7 + 5*i + 4]);
		printf("%d: (%lu, %lu, %lu, %lu, %lu)\n", i, series_params[5*i], series_params[5*i + 1], series_params[5*i + 2], series_params[5*i + 3], series_params[5*i + 4]);
	}

	const ulong limit = atol(argv[1]);
	const uint files = atoi(argv[2]);
	const uint bundle = atoi(argv[3]);
	const ulong bound = atol(argv[4]);
	const char * resultname = argv[5];
	const char * folder = argv[6];

	if ((limit / (files * bundle)) % (getpagesize() / (sizeof(ulong))) != 0)
	{
		printf("In order for MMAP to work, the quantity [limit]/([files] * [bundle]) = %lu should be divisible by the page size, which is %lu.\n", limit / (files * bundle), getpagesize() / sizeof(ulong));
		exit(1);
	}

	multiply(limit, files, bundle, bound, resultname, folder, series_total, series_params);

	free(series_params);
}
