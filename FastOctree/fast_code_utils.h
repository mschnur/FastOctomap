#ifndef __FAST_CODE_UTILS_H__
#define __FAST_CODE_UTILS_H__

#ifdef __cplusplus
extern "C" {
#endif

#ifndef USE_PMC
#define USE_PMC 0
#endif

#if USE_PMC
// source:
// https://software.intel.com/en-us/forums/software-tuning-performance-optimization-platform-monitoring/topic/595214
// rdpmc_actual_cycles uses a "fixed-function" performance counter to return the count of actual CPU core cycles
//       executed by the current core.  Core cycles are not accumulated while the processor is in the "HALT" state,
//       which is used when the operating system has no task(s) to run on a processor core.
__attribute__((always_inline)) inline unsigned long long rdtsc()
{
   static const unsigned c = (1<<30)+1;
   
   unsigned a, d;
   __asm__ volatile("rdpmc" : "=a" (a), "=d" (d) : "c" (c));
   return ((unsigned long long)a) | (((unsigned long long)d) << 32);;
}

#define NOMINAL_TO_BOOST_FREQUENCY_FACTOR (1.0)

#else
/*
Returns the processors clock value.
The clock ticks at nominal frequency,
which may differ from the clock frequency at which your programs execute.
*/
__attribute__((always_inline)) inline unsigned long long rdtsc() {
    unsigned a, d;
    __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));
    return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

#define NOMINAL_TO_BOOST_FREQUENCY_FACTOR (3.2 / 2.4)
#endif

#define INTEGER_ADD(dest, src) \
    __asm__ __volatile__( \
        "add %[rsrc], %[rdest]\n" \
        : [rdest] "+r"(dest) \
        : [rsrc] "r"(src) \
    );

#define INTEGER_MULTIPLY(dest, src) \
    __asm__ __volatile__( \
        "imul %[rsrc], %[rdest]\n" \
        : [rdest] "+r"(dest) \
        : [rsrc] "r"(src) \
    );

#define SIMD_ADD(dest, src) \
    __asm__ __volatile__( \
        "vaddpd %[rsrc], %[rdest], %[rdest]\n" \
        : [rdest] "+x"(dest) \
        : [rsrc] "x"(src) \
    );

#define SIMD_MULTIPLY(dest, src) \
    __asm__ __volatile__( \
        "vmulpd %[rsrc], %[rdest], %[rdest]\n" \
        : [rdest] "+x"(dest) \
        : [rsrc] "x"(src) \
    );
	
#define SIMD_FMA(dest, src1, src2) \
    __asm__ __volatile__( \
        "vfmadd231pd %[rsrc2], %[rsrc1], %[rdest] \n" \
        : [rdest] "+x"(dest) \
        : [rsrc1] "x"(src1), [rsrc2] "x"(src2) \
    );
	
#include <immintrin.h>
typedef union {
	__m256d simdBlock;
	double doubleVals[4];
} FourPackedDouble;

#ifndef DEBUG_PRINTS
#define DEBUG_PRINTS 1
#endif

#if DEBUG_PRINTS
#define dprintf(...) printf(__VA_ARGS__)
#else
#define dprintf(...) 
#endif

#ifndef STAT_PRINTS
#define STAT_PRINTS 0
#endif

#if STAT_PRINTS
#define statprintf(...) printf(__VA_ARGS__)
#else
#define statprintf(...) 
#endif

#define DO_TWO_TIMES(op) \
    op; \
	op;
	
#define DO_THREE_TIMES(op) \
    op; \
	op; \
	op; 

#define DO_FOUR_TIMES(op) \
    op; \
	op; \
	op; \
	op;
	
#define DO_FIVE_TIMES(op) \
    op; \
	op; \
	op; \
	op; \
	op;

#define DO_TEN_TIMES(op) \
    DO_FIVE_TIMES(op) \
	DO_FIVE_TIMES(op)

#define DO_HUNDRED_TIMES(op) DO_TEN_TIMES(DO_TEN_TIMES(op))

#define DO_FIVE_HUNDRED_TIMES(op) DO_FIVE_TIMES(DO_HUNDRED_TIMES(op))

#define DO_THOUSAND_TIMES(op) DO_HUNDRED_TIMES(DO_TEN_TIMES(op))

#define DO_TWO_THOUSAND_TIMES(op) \
	DO_THOUSAND_TIMES(op) \
	DO_THOUSAND_TIMES(op)

#define DO_FIVE_THOUSAND_TIMES(op) DO_FIVE_TIMES(DO_THOUSAND_TIMES(op))

#define DO_TEN_THOUSAND_TIMES(op) DO_TEN_TIMES(DO_THOUSAND_TIMES(op))

#define DO_HUNDRED_THOUSAND_TIMES(op) DO_TEN_TIMES(DO_TEN_THOUSAND_TIMES(op))
	
#ifdef __cplusplus
}
#endif	

#endif // #ifndef __FAST_CODE_UTILS_H__