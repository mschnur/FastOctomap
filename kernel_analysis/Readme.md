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
- FP ADD/SUB - 3cycle lat, 1 instr/cycle p1

- OR/AND
    - 1 cyle lat, 4inst/cycle - p0156

- SAR
    - 1 cycle lat, 2inst/cycle - p06

FP ops will be bottleneck

Based on class: MULs/FMAs (if we use SIMD) will be bottleneck.
Need 5 rays at once to fill MUL/FMA pipelines.



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
