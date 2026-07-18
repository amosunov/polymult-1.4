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

#ifndef INIT_H_
#define INIT_H_

#ifndef uint
#define uint unsigned int
#endif

#ifndef ulong
#define ulong unsigned long
#endif

#include <flint/nmod_poly.h>


// Block-by-block initialization

void init_block_nabla(int * block, const ulong size, const ulong min, const ulong c, const ulong r, const ulong s, const ulong a, const ulong m);

void init_block_nabla_product(int * block, const ulong size, const ulong min, const ulong c1, const ulong r1, const ulong s1, const ulong a1, const ulong m1, const ulong c2, const ulong r2, const ulong s2, const ulong a2, const ulong m2);



// FLINT polynomials

void nmod_poly_theta3(nmod_poly_t poly, const ulong size);

void nmod_poly_to_files(const nmod_poly_t poly, const ulong size, const uint files, const char * resultname);

#endif /* INIT_H_ */
