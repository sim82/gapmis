#ifndef         ERRORS_H
#define         ERRORS_H

/* REMARK: uuhm, using decimal values with 0x hex prefix seems a bit akward... are they supposed to be hex? */
#define         LENGTH       	0x00000001
#define         MATRIX		0x00000002
#define         MAXGAP		0x00000004
#define         MALLOC		0x00000008
#define		NOGPU		0x00000016
#define         BADCHAR		0x00000032
#define         IO		0x00000064
#define		GPUERROR	0x00000128
#define		KERNEL		0x00000256
#define		GPUMALLOC	0x00000512

#define     TEXTLEN   0x00001024 /* texts of different length given to sse version */
#define     THREADCOUNT   0x00002048

#endif



