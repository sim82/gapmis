/**
    libgapmis: a library for pairwise sequence aligment with a single gap.
    Copyright (C) 2012 Nikos Alachiotis, Simon Berger, Tomas Flouri, and
    Solon P. Pissis. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#ifndef         ERRORS_H
#define         ERRORS_H

#define         LENGTH       		0x00000001
#define         MATRIX			0x00000002
#define         MAXGAP			0x00000004
#define         MALLOC			0x00000008
#define		NOGPU			0x00000016
#define         BADCHAR			0x00000032
#define         IO			0x00000064
#define		GPUERROR		0x00000128
#define		KERNEL			0x00000256
#define		GPUMALLOC		0x00000512
#define     	TEXTLEN   		0x00001024 /* texts of different length given to sse version */
#define     	THREADCOUNT   		0x00002048
#define     	UNHANDLED_INTERNAL 	0x00004096

#endif



