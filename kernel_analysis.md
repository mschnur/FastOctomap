goal: update leaf nodes?

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

observations and thoughts:

nodes at same level can be processed in parallel


work queue with proc_subtree calls? 
synchronization issue with updating nodes - too many mutexes to track
partition octree to different threads?
not all threads might get used (partition never gets ray?)

recursive partitioning?
can i have a queue per thread that other threads add to?
