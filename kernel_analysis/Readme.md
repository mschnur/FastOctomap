# Kernel Analysis
## Notes

goal: update leaf nodes?

Unroll comparisons for one ray - one comparison chain is an independent operation
each level needs 8 bits to represent nodes passed through

The kernel is:
- computing the nodes passed through for each ray.
- identify the node "state" (contains/before/after endpoint) and update accordingly

pseudo_code
    for each ray
        compute a
        compute txyz0, txyz1 (this can be vectorized)

    proc_subtree
        compute node state
    {
        compute txyzm
        compute next nodes (this is a dependent conditional chain)
        process all next next nodes
        
        update node
        if leaf node 
        if parent -> max of leaf nodes (maybe we don't worry about this - see goal)
    }

Parallelization:

- nodes at same level can be processed in parallel
- work queue with proc_subtree calls?
- synchronization issue with updating nodes - too many mutexes to track
- partition octree to different threads?
- not all threads might get used (partition never gets ray?)
- recursive partitioning?
- can i have a queue per thread that other threads add to?
- I can parallelize the node computation at each recursion level (OpenMP parallel for)

Terminal Logic:

- we can just not add a node to the list of nodes if it is "past" the ray endpoint
- 


## Plan of tasking

- extract node "computation" before recursive calls for one ray
- convert comparisons to arithmetic to take advantage of pipelining
- vectorize this extracted code

## Theoretical Peak thoughts

- Assume branch predicition always correct
- Compute worst case number of comparisons
- assume zero cost (or fixed cost) of function call 
- assume max number of nodes passed through (3?)

Instructions:
- FP MUL - 5 cycles, 1 instr/cycle p0
- FMA - 5 cycles, 2 instr/cycle, p01
- FP ADD/SUB - 3cycle lat, 1 instr/cycle p1

- OR/AND
    - 1 cyle lat, 4inst/cycle - p0156

- SAR
    - 1 cycle lat, 2inst/cycle - p06

FP ops will be bottleneck

Based on class: MULs/FMAs (if we use SIMD) will be bottleneck.
Need 5 rays at once to fill FP MUL/FMA pipelines (assuming the second FMA pipeline needs to be used by the other instructions);


## ray_parameter kernel
I have a bunch of the latency/throughputs/ports of these written down on paper, but I haven't gotten a chance to copy them to this document yet
and I haven't looked at bubbles or register allocation yet.

#### Micro-kernel 1 - Calculate direction
// I assume that the endpoint and origin will be packed into an aligned
// array of doubles, with the format [x0, y0, z0, 0, x1, y1, z1, 0, ...]
__m256d end = _mm256_load_pd(...)  // vmovapd
__m256d origin = _mm256_load_pd(...)  // vmovapd
sub end - origin = diff
mul diff * diff = d2
shuffle d2, d2, 5 -> d2_shuf  // swap elements 0 and 1 and elements 2 and 3
add d2 + d2_shuf = part_sum
vperm2f128 part_sum, part_sum, 1 -> perm_part_sum  // swap upper and lower 128 bits
add perm_part_sum + part_sum = magnitude
sqrt magnitude = magnitude
div diff / magnitude = direction

##### Issues with Micro-kernel 1 a written above
- sqrt (vsqrtpd) has a latency of 28-29 and a reciprocal throughput of 16-28
- div (vdivpd) has a latency of 19-35 and a reciprocal throughput of 16-28
- Solution 1: Use Fast Inverse Square Root (https://en.wikipedia.org/wiki/Fast_inverse_square_root) with double precision. The magic number for double precision is on page 33 of (https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf)
- Solution 2: Use single precision instead of double, and then we can use rsqrtps (single precision reciprocal square root) which has a latency of 5 and a reciprocal throughput of 1. Going to single precision also has the added benefit that we can use rcpps (single precision reciprocal - latency 5, reciprocal throughput 1) instead of vdivpd in Micro-kernel 3 when computing the reciprocal of the direction.

#### Micro-kernel 1 v2 - Calculate direction (single precision)
// I assume that the endpoint and origin will be packed into an aligned
// array of floats, with the format [x0, y0, z0, 0, x1, y1, z1, 0, ...]
__m256 end = _mm256_load_ps(...)  // vmovaps
__m256 origin = _mm256_load_ps(...)  // vmovaps
_mm256_sub_ps(end, origin) = diff  // vsubps
_mm256_mul_ps(diff, diff) = d2  // vmulps
// swap elements 0 and 1 and elements 2 and 3 of both the upper and lower halves
_mm256_shuffle_ps(d2, d2, 0xB1) -> d2_shuf  // vshufps
_mm256_add_ps(d2, d2_shuf) = part_sum  // vaddps
_mm256_shuffle_ps(part_sum, part_sum, 0x4E)  -> perm_part_sum  // (vshufps) swap upper and lower 64 bits of each 128-bit half
_mm256_add_ps(perm_part_sum, part_sum) = magnitude  // vaddps
add perm_part_sum + part_sum = magnitude
_mmm256_rsqrt_ps(magnitude) = recip_magnitude // vrsqrtps
_mm256_mul_ps(diff, recip_magnitude) = direction // vmulps

#### Micro-kernel 2 - Reflect negative directions and calculate "a"
union UShortsInInt64
{
    __Int64 int;
    uint16_t uShorts[4];
} intToShorts;
static const __m256 FourTwoOne = _mm256_set_ps(4.f, 2.f, 1.f, 0.f, 4.f, 2.f, 1.f, 0.f)  // What instruction?
static const __m256i signBits = _mm256_set1_epi32(0x80000000); // May be vpbroadcastd
// Cast signBits to __m256 (free - no instructions generated)
__m256 floatSignBits = _mm256_castsi256_ps(signBits)  // free (no instrucitons generated) 
// Handle reflections when direction is negative
and intDirection & floatSignBits = sign_bit_if_neg // _mm256_and_ps (vandps)
// toggle sign bit of origin
xor origin, sign_bit_if_neg = origin // _mm256_xor_ps (vxorps)
// clear sign bit from direction if direction was negative
not sign_bit_if_neg = all_but_sign_bit_if_neg
and direction & !sign_bit_if_neg = direction // _mm256_andnot_ps (vandnps)
// calculate "a"
and !sign_bit_if_neg & FourTwoOne = a_parts // _mm256_andnot_ps (vandnps)
// Convert a_parts into packed 32-bit signed integers
vcvttps2dq (_mm256_cvttps_epi32) a_part -> a_parts_int

// swap elements 0 and 1 of both the upper and lower halves
_mm256_shuffle_epi32(a_parts_int, a_parts_int, 0xE1) -> a_parts_int_shuff1  // vpshufd
_mm256_or_si256(a_parts_int, a_parts_int_shuff1) = partial_a_int  // vpor
// Put element 2 into slots 0 and 1 of both the upper and lower halves
_mm256_shuffle_epi32(a_parts_int, a_parts_int, 0xEA) -> a_parts_int_shuff2  // vpshufd
_mm256_or_si256(partial_a_int, a_parts_int_shuff2) = a_int  // vpor
// Now elements 0, 1, and 2 of both the upper and lower halves of a_int contain the "a" value
unsigned char a_low = (unsigned char) _mm256_cvtsi256_si32(a_int) // movd
_mm256_unpackhi_epi32(a_int, a_int) = a_int_hi // vpunpckhdq
unsigned char a_high = (unsigned char) _mm256_cvtss_f32(a_int_hi) // movd

#### Micro-kernel 3 - Compute t0 and t1
// All axes have the same minimum, so we don't actually need minX, minY, minZ
__m256 treeMins = _mm256_set1_ps((float) tree->min);
// Same with the maximums
__m256 treeMaxes = _mm256_set1_ps((float) tree->max);
_mm256_rcp_ps(direction) = dirInverse  // vrcpps
sub treeMins - origin = t0  // _mm256_sub_ps (vsubps)
mul t0 * dirInverse = t0  // _mm256_mul_ps (vmulps)
sub treeMaxes - origin = t1  // _mm256_sub_ps (vsubps)
mul t1 * dirInverse = t1  // _mm256_mul_ps (vmulps)
// Note: we need to convert the floats back to doubles when we store
// them in the output arrays (unless we want to change everything to
// use floats)


#### Micro-kernel 4 - Compare max of t0 and min of t1 
// swap elements 0 and 1 of both the upper and lower halves
_mm256_shuffle_ps(t0, t0, 0xE1) -> t0_shuff1  // vshufps
_mm256_max_ps(t0, t0_shuff1) = part_t0_max  // vmaxps
// Put element 2 into slots 0 and 1 of both the upper and lower halves
_mm256_shuffle_ps(t0, t0, 0xEA) -> t0_shuff2  // vshufps
_mm256_max_ps(part_t0_max, t0_shuff2) = t0_max  // vmaxps
// Now elements 0, 1, and 2 of t0_max contain the max t0 value

// swap elements 0 and 1 of both the upper and lower halves
_mm256_shuffle_ps(t1, t1, 0xE1) -> t1_shuff1  // vshufps
_mm256_min_ps(t1, t1_shuff1) = part_t1_min  // vminps
// Put element 2 into slots 0 and 1 of both the upper and lower halves
_mm256_shuffle_ps(t1, t1, 0xEA) -> t1_shuff2  // vshufps
_mm256_min_ps(part_t1_min, t1_shuff2) = t1_min  // vminps
// Now elements 0, 1, and 2 of t1_min contain the min t1 value

// each 32-bit element in cmp_result will contain 0xFFFFFFFF if max 
// is less than min for that index, otherwise it will contain 0 
_mm256_cmp_ps(max, min, 0x11) = cmp_result  // vcmpps (0x11 == _CMP_LT_OQ == less-than, ordered, non-signaling)
int low_result = (int) _mm256_cvtss_f32(cmp_result) // movss
_mm256_unpackhi_ps(cmp_result, cmp_result) = cmp_hi // vunpckhps
int high_result = (int) _mm256_cvtss_f32(cmp_hi) // movss


## Files

conversions.c
- contains the node computation functions with and without comparisons

extracted_node_computations.c
- computes the nodes passed through as a 8bit bit-vector.

remove_node_switch.c
- same as above, but replace the loop+switch statement with if statements
- in progress to replace with only bit manipulation

with_terminal_logic.c
- as one previous, but re-adding the terminal logic (with some adjustments)

with_terminal_logic_converted.c
- as one previous, but converting the conditional logic (also removed if statements around next_node)
