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


## Plan of tasking

- extract node "computation" before recursive calls for one ray
- convert comparisons to arithmetic to take advantage of pipelining
- vectorize this extracted code

## Theoretical Peak thoughts

- Assume branch predicition always correct
- Compute worst case number of comparisons
- assume zero cost (or fixed cost) of function call 
- assume max number of nodes passed through (3?)

## Files

conversions.c
- contains the node computation functions with and without comparisons

extracted_node_computations.c
- computes the nodes passed through as a 8bit bit-vector.

remove_node_switch.c
- same as above, but replace the loop+switch statement with if statements
