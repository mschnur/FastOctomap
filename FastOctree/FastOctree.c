//
// Created by mschnur on 11/12/19.
//

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FastOctree.h"
#include "Stack.h"

#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))

#define TRUE 1
#define FALSE 0

//#define DEBUG_TRACE_FUNCTION_CALLS 1
//#define DEBUG_RAY_PARAMETER 1
//#define DEBUG_PROC_SUBTREE 1

const Node* const NO_CHILD = NULL;

const Node* const NO_CHILDREN[8] = {
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL
};

double CLAMPING_THRES_MIN;
double CLAMPING_THRES_MAX;
double PROB_HIT_LOG;
double PROB_MISS_LOG;
double OCC_PROB_THRES_LOG;

const double SIZE_LOOKUP_TABLE[MAX_DEPTH + 1] = {
        RESOLUTION * (double) (1u << (MAX_DEPTH - 0u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 1u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 2u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 3u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 4u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 5u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 6u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 7u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 8u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 9u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 10u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 11u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 12u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 13u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 14u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 15u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 16u))
};

static const size_t FastOctree_to_octomap_index_map[8] = {
        0, 4, 2, 6, 1, 5, 3, 7
};

//###########################################################################################
// Utility functions
//###########################################################################################

static inline int nodeHasAnyChildren(const Node* n)
{
    return (memcmp(n->children, NO_CHILDREN, 8 * sizeof(Node*)) != 0);
}

static inline int nodeIsPrunable(const Node* n)
{
    // If all of this node's children have the same log-likelihood and have no children themselves, then we can
    // prune these children, meaning we delete them.
    for (unsigned int i = 0; i < 8; ++i)
    {
        Node* child = n->children[i];
        if (child == NO_CHILD || child->logOdds != n->logOdds)
        {
            return FALSE;
        }
    }

    return TRUE;
}

static inline void deleteChild(Node* n, unsigned int childIndex)
{
    free(n->children[childIndex]);
    n->children[childIndex] = NULL;
}


static inline void deleteAllChildren(Node* n)
{
    for (unsigned int i = 0; i < 8; ++i)
    {
        deleteChild(n, i);
    }
}

static inline void createChild(Node* n, unsigned int childIndex)
{
    // In order for us to be able to free the nodes individually (which may happen during pruning) we need to
    // allocate each node separately, instead of the whole array at once. Using calloc instead of malloc
    // initializes the memory to zero, which means that the this new node's `children` array will be filled with
    // zeros (which is equivalent to filling it will NULL pointers).
    n->children[childIndex] = (Node*) calloc(1, sizeof(Node));
}

static inline int createChildIfItDoesntExist(Node* n, unsigned int childIndex)
{
    if (n->children[childIndex] == NO_CHILD)
    {
        createChild(n, childIndex);
        return TRUE;
    }

    return FALSE;
}

static inline double maxChildLogLikelihood(const Node* n)
{
    double maxLogLikelihood = -INFINITY;
    for (unsigned int i = 0; i < 8; ++i)
    {
        Node* child = n->children[i];
        if (child != NO_CHILD)
        {
            if (child->logOdds > maxLogLikelihood)
            {
                maxLogLikelihood = child->logOdds;
            }
#ifdef DEBUG_PROC_SUBTREE
            printf("    maxChildLogLikelihood: Child %u has log-odds of %lf\n", i, child->logOdds);
#endif
        } else {
#ifdef DEBUG_PROC_SUBTREE
            printf("    maxChildLogLikelihood: No child %u\n", i);
#endif
        }
    }

    return maxLogLikelihood;
}

static inline void expandPrunedNode(Node* n)
{
    assert(!nodeHasAnyChildren(n));

    for (unsigned int i = 0; i < 8; ++i)
    {
        createChild(n, i);

        // The children should all start with the same log odds as their parent.
        n->children[i]->logOdds = n->logOdds;
    }
}

static inline int is_less(double pointX, double pointY, double pointZ,
                   double minX, double minY, double minZ)
{
    return (minX > pointX && minY > pointY && minZ > pointZ);
}

static inline int is_greater(double pointX, double pointY, double pointZ,
                      double maxX, double maxY, double maxZ)
{
    return (pointX > maxX && pointY > maxY && pointZ > maxZ);
}

static inline int any_is_less(double pointX, double pointY, double pointZ,
                          double minX, double minY, double minZ)
{
    return (minX > pointX || minY > pointY || minZ > pointZ);
}

static inline int any_is_greater(double pointX, double pointY, double pointZ,
                             double maxX, double maxY, double maxZ)
{
    return (pointX > maxX || pointY > maxY || pointZ > maxZ);
}

static inline int is_between(double minX, double minY, double minZ,
                      double pointX, double pointY, double pointZ,
                      double maxX, double maxY, double maxZ)
{
    return (minX <= pointX && pointX <= maxX &&
            minY <= pointY && pointY <= maxY &&
            minZ <= pointZ && pointZ <= maxZ);
}

static inline int is_between_vector3d(const Vector3d* min,
                               const Vector3d* point,
                               const Vector3d* max)
{
    return (min->x <= point->x && point->x <= max->x &&
            min->y <= point->y && point->y <= max->y &&
            min->z <= point->z && point->z <= max->z);
}

//###########################################################################################
// Core algorithm functions
//###########################################################################################

int new_node(double txm, int x, double tym, int y, double tzm, int z)
{
    unsigned long long currentNode = 0;

    // txm = 0.45;
    // tym = 0.55;
    // tzm = 0.175;

    double tmp1 = txm - tym;
    double tmp2 = txm - tzm;
    double tmp3 = tym - tzm;

    //X
    currentNode |= (unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & 0x8000000000000000)>>(63-x);
    
    //Y
    currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>(63-y);
    
    //Z
    currentNode |= (unsigned long long)( ~*(unsigned long long*)(&tmp3) & ~*(unsigned long long*)(&tmp2) & 0x8000000000000000 )>>(63-z);

    return (int)currentNode;
}


void proc_subtree(double tx0, double ty0, double tz0,
                  double tx1, double ty1, double tz1,
                  unsigned int depth, int justCreated,
                  Node* n, Node *parent, unsigned char childIndex, unsigned char a, Ray* r)
{
#ifdef DEBUG_TRACE_FUNCTION_CALLS
    //printf("proc_subtree, depth=%u\n", depth);
#endif

#ifdef DEBUG_PROC_SUBTREE
    printf("proc_subtree called with depth = %u, a = %u, justCreated = %s, r->t_end = %lf\n"
           "                         tx0 = %lf, ty0 = %lf, tz0 = %lf,\n"
           "                         tx1 = %lf, ty1 = %lf, tz1 = %lf\n",
           depth, a, (justCreated ? "TRUE" : "FALSE"), r->t_end, tx0, ty0, tz0, tx1, ty1, tz1);
#endif

    double txm, tym, tzm;
    int currentNode = 0;

    if (tx1 < 0.0 || ty1 < 0.0 || tz1 < 0.0)
    {
#ifdef DEBUG_PROC_SUBTREE
        printf("tx1, ty1, or tz1 is negative; returning\n");
#endif
        return;
    }
    else if (n == NULL && parent != NULL)
    {
        if (any_is_less(r->t_end, r->t_end, r->t_end, tx0, ty0, tz0))
        {
#ifdef DEBUG_PROC_SUBTREE
            printf("r->t_end is less than at least one of tx0, ty0, and tz0; returning\n");
#endif
            // The ray endpoint happened before at least one of the minimum t values of this subtree, meaning the ray
            // ends before this subtree. Therefore, we don't want to do anything to this subtree.
            return;
        }

        justCreated = createChildIfItDoesntExist(parent, childIndex);
        n = parent->children[childIndex];
#ifdef DEBUG_PROC_SUBTREE
        printf("depth=%u: This node is the child with the index = %u, just created this node = %s\n",
                depth, childIndex, (justCreated ? "TRUE" : "FALSE"));
#endif
    }

    if (!nodeHasAnyChildren(n))
    {
        if (any_is_less(r->t_end, r->t_end, r->t_end, tx0, ty0, tz0))
        {
#ifdef DEBUG_PROC_SUBTREE
            printf("r->t_end is less than at least one of tx0, ty0, and tz0; returning\n");
#endif
            // The ray endpoint happened before at least one of the minimum t values of this subtree, meaning the ray
            // ends before this subtree. Therefore, we don't want to do anything to this subtree.
            return;
        }
        else if (depth == MAX_DEPTH)
        {
#ifdef DEBUG_PROC_SUBTREE
            printf("depth is equal to MAX_DEPTH, so updating the log likelihood of this node\n");
#endif

            double logLikelihoodUpdate = 0.0;

            if (any_is_greater(r->t_end, r->t_end, r->t_end, tx1, ty1, tz1))
            {
#ifdef DEBUG_PROC_SUBTREE
                printf("r->t_end is greater than any of tx1, ty1, and tz1; voting a MISS\n");
#endif
                // The ray endpoint does not occur in this subtree, but the ray passes through this subtree on its
                // way to the endpoint, and we're at our maximum depth. Therefore we need to give this node a vote
                // that it is free space.
                logLikelihoodUpdate = PROB_MISS_LOG;
            }
            else
            {
#ifdef DEBUG_PROC_SUBTREE
                printf("r->t_end is NOT greater than any of tx1, ty1, and tz1; voting a HIT\n");
#endif
                // The ray endpoint occurs within this subtree, and we're at our maximum depth. Therefore we need to
                // give this node a vote that it is free space.
                logLikelihoodUpdate = PROB_HIT_LOG;
            }

#ifdef DEBUG_PROC_SUBTREE
            double previousLogOdds = n->logOdds;
#endif

            // Do the update
            n->logOdds += logLikelihoodUpdate;

#ifdef DEBUG_PROC_SUBTREE
            printf("Log-odds update changes log-odds from %lf to %lf\n",
                    previousLogOdds, n->logOdds);
            double preClampLogOdds = n->logOdds;
#endif

            // Clamp the logOdds between the min/max
            n->logOdds = fmax(CLAMPING_THRES_MIN, fmin(n->logOdds, CLAMPING_THRES_MAX));

#ifdef DEBUG_PROC_SUBTREE
            printf("Clamping between bounds changes log-odds from %lf to %lf\n",
                   preClampLogOdds, n->logOdds);
#endif

#ifdef DEBUG_PROC_SUBTREE
            printf("Returning because we've updated the log odds of a node at the max depth\n");
#endif
            return;
        }
        else if (!justCreated)
        {
#ifdef DEBUG_PROC_SUBTREE
            printf("Not at our max depth, and this node was not just created, so expanding pruned node\n");
#endif
            // Either the ray endpoint occurs within this subtree or the ray passes through this subtree on its way to
            // the endpoint, but we're not at our maximum depth yet and this is a previously pruned node, so we want to
            // expand this node and then proceed as normal.
            expandPrunedNode(n);
        }
    }

    txm = 0.5 * (tx0 + tx1);
    tym = 0.5 * (ty0 + ty1);
    tzm = 0.5 * (tz0 + tz1);

#ifdef DEBUG_PROC_SUBTREE
    printf("txm = %lf, tym = %lf, tzm = %lf\n", txm, tym, tzm);
#endif

    double tmp1 = ty0 - tx0;
    double tmp2 = tz0 - tx0;
    double tmp3 = tym - tx0;
    double tmp4 = tzm - tx0;

    double tmp5 = tz0 - ty0;
    double tmp6 = txm - ty0;
    double tmp7 = tzm - ty0;

    double tmp8 = txm - tz0;
    double tmp9 = tym - tz0;


    currentNode |= (int)((unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>62);
    currentNode |= (int)((unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & *((unsigned long long*)(&tmp4)) & 0x8000000000000000)>>63);

    currentNode |= (int)((unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp6)) & 0x8000000000000000)>>61);
    currentNode |= (int)((unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp7)) & 0x8000000000000000)>>63);

    currentNode |= (int)((unsigned long long)(~*((unsigned long long*)(&tmp2)) & ~*((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp8)) & 0x8000000000000000)>>61);
    currentNode |= (int)((unsigned long long)(~*((unsigned long long*)(&tmp2)) & ~*((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp9)) & 0x8000000000000000)>>62);


#ifdef DEBUG_PROC_SUBTREE
    printf("depth=%u: first_node: %d\n", depth, currentNode);
#endif

    unsigned char nodes = 0;
    currentNode = 1u<<currentNode;


    int eq = ~(((1<<0) - currentNode)>>31);
    nodes |= currentNode & eq;
    currentNode = (new_node(txm, 4, tym, 2, tzm, 1) & eq) | (currentNode & ~eq);

    eq = ~(((1<<1) - currentNode)>>31);
    nodes |= currentNode & eq;
    currentNode = (new_node(txm, 5, tym, 3, tz1, 8) & eq)| (currentNode & ~eq);

    eq = ~(((1<<2) - currentNode)>>31);
    nodes |= currentNode & eq;
    currentNode = (new_node(txm, 6, ty1, 8, tzm, 3) & eq) | (currentNode & ~eq);

    eq = ~(((1<<3) - currentNode)>>31);
    nodes |= currentNode & eq;
    currentNode = (new_node(txm, 7, ty1, 8, tz1, 8) & eq) | (currentNode & ~eq);

    eq = ~(((1<<4) - currentNode)>>31);
    nodes |= currentNode & eq;
    currentNode = (new_node(tx1, 8, tym, 6, tzm, 5) & eq) | (currentNode & ~eq);

    eq = ~(((1<<5) - currentNode)>>31);
    nodes |= currentNode & eq;
    currentNode = (new_node(tx1, 8, tym, 7, tz1, 8) & eq) | (currentNode & ~eq);

    eq = ~(((1<<6) - currentNode)>>31);
    nodes |= currentNode & eq;
    currentNode = (new_node(tx1, 8, ty1, 8, tzm, 7) & eq) | (currentNode & ~eq);

    eq = ~(((1<<7) - currentNode)>>31);
    nodes |= currentNode & eq;

    for(int node = 0; node < 8; node++) {
        int createdChild = FALSE;
        if(nodes & (1u<<node)){
            switch (node) {
                case 0:
                    //createdChild = createChildIfItDoesntExist(n, a);
                    #ifdef DEBUG_PROC_SUBTREE
                                    //printf("depth=%u: current node = %d, child index = %u, created child = %s\n", depth, currentNode, a, (createdChild ? "TRUE" : "FALSE"));
                    #endif
                    proc_subtree(tx0, ty0, tz0, txm, tym, tzm, depth + 1, createdChild, n->children[a], n, a, a, r);
                    currentNode = new_node(txm, 4, tym, 2, tzm, 1);
                    #ifdef DEBUG_PROC_SUBTREE
                                    printf("depth=%u: previous node = 0, new_node = %d\n", depth, currentNode);
                    #endif
                    break;

                case 1:
                    //createdChild = createChildIfItDoesntExist(n, 1u^a);
                    #ifdef DEBUG_PROC_SUBTREE
                                    //printf("depth=%u: current node = %d, child index = %u, created child = %s\n", depth, currentNode, 1u^a, (createdChild ? "TRUE" : "FALSE"));
                    #endif
                    proc_subtree(tx0, ty0, tzm, txm, tym, tz1, depth + 1, createdChild, n->children[1u^a], n, 1u^a, a, r);
                    currentNode = new_node(txm, 5, tym, 3, tz1, 8);
                    #ifdef DEBUG_PROC_SUBTREE
                                    printf("depth=%u: previous node = 1, new_node = %d\n", depth, currentNode);
                    #endif
                    break;

                case 2:
                    //createdChild = createChildIfItDoesntExist(n, 2u^a);
                    #ifdef DEBUG_PROC_SUBTREE
                                    //printf("depth=%u: current node = %d, child index = %u, created child = %s\n", depth, currentNode, 2u^a, (createdChild ? "TRUE" : "FALSE"));
                    #endif
                    proc_subtree(tx0, tym, tz0, txm, ty1, tzm, depth + 1, createdChild, n->children[2u^a], n, 2u^a, a, r);
                    currentNode = new_node(txm, 6, ty1, 8, tzm, 3);
                    #ifdef DEBUG_PROC_SUBTREE
                                    printf("depth=%u: previous node = 2, new_node = %d\n", depth, currentNode);
                    #endif
                    break;

                case 3:
                    //createdChild = createChildIfItDoesntExist(n, 3u^a);
                    #ifdef DEBUG_PROC_SUBTREE
                                    //printf("depth=%u: current node = %d, child index = %u, created child = %s\n", depth, currentNode, 3u^a, (createdChild ? "TRUE" : "FALSE"));
                    #endif
                    proc_subtree(tx0, tym, tzm, txm, ty1, tz1, depth + 1, createdChild, n->children[3u^a], n, 3u^a, a, r);
                    currentNode = new_node(txm, 7, ty1, 8, tz1, 8);
                    #ifdef DEBUG_PROC_SUBTREE
                                    printf("depth=%u: previous node = 3, new_node = %d\n", depth, currentNode);
                    #endif
                    break;

                case 4:
                    //createdChild = createChildIfItDoesntExist(n, 4u^a);
                    #ifdef DEBUG_PROC_SUBTREE
                                    //printf("depth=%u: current node = %d, child index = %u, created child = %s\n", depth, currentNode, 4u^a, (createdChild ? "TRUE" : "FALSE"));
                    #endif
                    proc_subtree(txm, ty0, tz0, tx1, tym, tzm, depth + 1, createdChild, n->children[4u^a], n, 4u^a, a, r);
                    currentNode = new_node(tx1, 8, tym, 6, tzm, 5);
                    #ifdef DEBUG_PROC_SUBTREE
                                    printf("depth=%u: previous node = 4, new_node = %d\n", depth, currentNode);
                    #endif
                    break;

                case 5:
                    //createdChild = createChildIfItDoesntExist(n, 5u^a);
                    #ifdef DEBUG_PROC_SUBTREE
                                    //printf("depth=%u: current node = %d, child index = %u, created child = %s\n", depth, currentNode, 5u^a, (createdChild ? "TRUE" : "FALSE"));
                    #endif
                    proc_subtree(txm, ty0, tzm, tx1, tym, tz1, depth + 1, createdChild, n->children[5u^a], n, 5u^a, a, r);
                    currentNode = new_node(tx1, 8, tym, 7, tz1, 8);
                    #ifdef DEBUG_PROC_SUBTREE
                                    printf("depth=%u: previous node = 5, new_node = %d\n", depth, currentNode);
                    #endif
                    break;

                case 6:
                    //createdChild = createChildIfItDoesntExist(n, 6u^a);
                    #ifdef DEBUG_PROC_SUBTREE
                                    //printf("depth=%u: current node = %d, child index = %u, created child = %s\n", depth, currentNode, 6u^a, (createdChild ? "TRUE" : "FALSE"));
                    #endif
                    proc_subtree(txm, tym, tz0, tx1, ty1, tzm, depth + 1, createdChild, n->children[6u^a], n, 6u^a, a, r);
                    currentNode = new_node(tx1, 8, ty1, 8, tzm, 7);
                    #ifdef DEBUG_PROC_SUBTREE
                                    printf("depth=%u: previous node = 6, new_node = %d\n", depth, currentNode);
                    #endif
                    break;

                case 7:
                    //createdChild = createChildIfItDoesntExist(n, 7u^a);
                    #ifdef DEBUG_PROC_SUBTREE
                                    //printf("depth=%u: current node = %d, child index = %u, created child = %s\n", depth, currentNode, 7u^a, (createdChild ? "TRUE" : "FALSE"));
                    #endif
                    proc_subtree(txm, tym, tzm, tx1, ty1, tz1, depth + 1, createdChild, n->children[7u^a], n, 7u^a, a, r);
                    currentNode = 8;
                    #ifdef DEBUG_PROC_SUBTREE
                                    printf("depth=%u: previous node = 7, new_node = %d\n", depth, currentNode);
                    #endif
                    break;

                default:
                    assert(0);
            }
        }
    }

#ifdef DEBUG_PROC_SUBTREE
    double previousOdds = n->logOdds;
#endif
    n->logOdds = maxChildLogLikelihood(n);

#ifdef DEBUG_PROC_SUBTREE
    printf("Making the log-odds of this node equal to the max child log-odds changed the log-odds from %lf to %lf\n",
            previousOdds, n->logOdds);
#endif
}


void ray_parameter(Octree* tree, Ray* r) {
    int createdRoot = FALSE;
    if (tree->root == NULL) {
        // Using calloc instead of malloc initializes the memory to zero, which means that the the root's `children`
        // array will be filled with zeros (which is equivalent to filling it will NULL pointers). It also sets the
        // root's log odds to 0, which we want as our initial value.
        tree->root = (Node *) calloc(1, sizeof(Node));
        createdRoot = TRUE;
    }

#ifdef DEBUG_TRACE_FUNCTION_CALLS
    printf("ray_parameter\n");
#endif

    unsigned char a = 0;

#ifdef DEBUG_RAY_PARAMETER
    printf("Original ray origin: (%lf, %lf, %lf); direction: (%lf, %lf, %lf)\n",
            r->origin.x, r->origin.y, r->origin.z,
            r->direction.x, r->direction.y, r->direction.z);
#endif

    // Since the tree is centered at (0, 0, 0), reflecting the origins is as simple as negating them.
    if (r->direction.x < 0.0f) {
        r->origin.x = -r->origin.x;
        r->direction.x = -r->direction.x;
        a |= 4u;

#ifdef DEBUG_RAY_PARAMETER
        printf("X component of ray is negative, reflecting.\n");
#endif
    }

    if (r->direction.y < 0.0f) {
        r->origin.y = -r->origin.y;
        r->direction.y = -r->direction.y;
        a |= 2u;
#ifdef DEBUG_RAY_PARAMETER
        printf("Y component of ray is negative, reflecting.\n");
#endif
    }

    if (r->direction.z < 0.0f) {
        r->origin.z = -r->origin.z;
        r->direction.z = -r->direction.z;
        a |= 1u;
#ifdef DEBUG_RAY_PARAMETER
        printf("Z component of ray is negative, reflecting.\n");
#endif
    }

#ifdef DEBUG_RAY_PARAMETER
    printf("After potential reflection, ray origin: (%lf, %lf, %lf); direction: (%lf, %lf, %lf), a: %u\n",
           r->origin.x, r->origin.y, r->origin.z,
           r->direction.x, r->direction.y, r->direction.z, a);
#endif

    // Improve IEEE double stability
    double rdxInverse = 1.0 / r->direction.x;
    double rdyInverse = 1.0 / r->direction.y;
    double rdzInverse = 1.0 / r->direction.z;

    double tx0 = (tree->min.x - r->origin.x) * rdxInverse;
    double tx1 = (tree->max.x - r->origin.x) * rdxInverse;
    double ty0 = (tree->min.y - r->origin.y) * rdyInverse;
    double ty1 = (tree->max.y - r->origin.y) * rdyInverse;
    double tz0 = (tree->min.z - r->origin.z) * rdzInverse;
    double tz1 = (tree->max.z - r->origin.z) * rdzInverse;

#ifdef DEBUG_RAY_PARAMETER
    printf("txyz0: %lf %lf %lf\n", tx0, ty0, tz0);
    printf("txyz1: %lf %lf %lf\n", tx1, ty1, tz1);
#endif

    if (MAX(MAX(tx0, ty0), tz0) < MIN(MIN(tx1, ty1), tz1))
    {
        proc_subtree(tx0, ty0, tz0, tx1, ty1, tz1, 0, createdRoot, tree->root, NULL, 0, a, r);
    }
    else
    {
        printf("Ray outside of tree bounds\n");
    }
}

static inline void pruneNode(Node* node)
{
    if (nodeIsPrunable(node))
    {
        deleteAllChildren(node);
    }
}

void printNode(Node* node, size_t depth, size_t *number)
{
    char indentBuffer[50] = { 0 };
    for (size_t i = 0; i < 2 * depth; i += 2)
    {
        indentBuffer[i] = '|';
        indentBuffer[i + 1] = ' ';
    }

    indentBuffer[2 * depth] = ' ';

    char *allChildNodesNull = "";
    int noChildren = TRUE;
    if (node != NULL) {
        for (unsigned int i = 0; i < 8; ++i) {
            if (node->children[i] != NO_CHILD)
            {
                noChildren = FALSE;
                break;
            }
        }

        if (noChildren) {
            allChildNodesNull = " ; No children";
        }
    }

    *number = *number + 1;
    double val = (node == NULL) ? NAN : node->logOdds;
    if (node != NULL) {
        printf("%sd: %zu ; val: %lf ; addr: %p ; #: %zu%s\n", indentBuffer, depth, val, node, *number,
               allChildNodesNull);
    }

    if (!noChildren) {
        for (unsigned int i = 0; i < 8; ++i) {
            printNode(node->children[FastOctree_to_octomap_index_map[i]], depth + 1, number);
        }
    }
}

void printTree(Octree* tree)
{
    if (tree->root == NULL)
    {
        return;
    }

    size_t number = 0;
    printNode(tree->root, 0, &number);
}

size_t countNode(Node* node)
{
    if (node == NULL)
    {
        return 0;
    }

    size_t thisNodeTotal = 1;
    for (unsigned int i = 0; i < 8; ++i) {
        if (node->children[i] != NO_CHILD)
        {
            thisNodeTotal += countNode(node->children[i]);
        }
    }

    return thisNodeTotal;
}

size_t treeSize(Octree* tree)
{
    if (tree->root == NULL)
    {
        return 0;
    }

    return countNode(tree->root);
}

void pruneTree(Octree* tree)
{
    if (tree->root == NULL)
    {
        return;
    }

    static Stack* postOrderStack = NULL;

    if (postOrderStack == NULL)
    {
        postOrderStack = Stack__init(sizeof(Node*), 50000);
    }

    if (Stack__push(postOrderStack, &(tree->root)))
    {
        printf("Failed to push root element onto postOrderStack\n");
    }

    Node* lastNode = NULL;

    while (!Stack__isEmpty(postOrderStack)) {
        //printf("Stack size: %d\n", Stack__size(postOrderStack));

        Node *peeked = NULL;

        if (Stack__peek(postOrderStack, &peeked)) {
            printf("Failed to peek element from postOrderStack\n");
        }

        if (peeked == NULL)
        {
            printf("PEEKED IS NULL AFTER PEEK\n");
        }

        // check if lastNode is a child of peeked and if peeked has any children
        int lastNodeIsChild = FALSE;
        int anyChildren = FALSE;
        for (unsigned int i = 0; i < 8; ++i) {
            if (lastNode != NULL && peeked->children[i] == lastNode) {
                lastNodeIsChild = TRUE;
            }

            if (peeked->children[i] != NO_CHILD) {
                anyChildren = TRUE;
            }
        }

        // If current node has no children, or one child has been visited, then process and pop it
        if (!anyChildren || lastNodeIsChild) {
            pruneNode(peeked);

            if (Stack__pop(postOrderStack)) {
                printf("Failed to pop element from postOrderStack\n");
            }

            lastNode = peeked;
        }
        else // Otherwise, push the children onto the stack
        {
            for (int i = 7; i >= 0; --i) {
                Node *child = peeked->children[i];
                if (child != NO_CHILD) {
                    if (Stack__push(postOrderStack, &(peeked->children[i]))) {
                        printf("Failed to push element onto postOrderStack\n");
                    }
                }
            }
        }
    }
}

typedef struct NodeDetails
{
    Node* node;
    unsigned int depth;
    double centerX, centerY, centerZ, size, value;
} NodeDetails;

void createNodeCsv(Octree* tree, const char *filename)
{
    if (tree->root == NULL)
    {
        printf("Cannot create node CSV named %s because tree root is NULL\n", filename);
        return;
    }

    static Stack* treeDetailStack = NULL;

    if (treeDetailStack == NULL)
    {
        treeDetailStack = Stack__init(sizeof(NodeDetails*), 10000);
    }

    FILE *outFile = fopen(filename, "w");
    if (outFile == NULL)
    {
        printf("Failed to open file named %s\n", filename);
        return;
    }

    fprintf(outFile, "depth,centerX,centerY,centerZ,size,value\n");

    NodeDetails *rootNodeDetails = calloc(1, sizeof(NodeDetails));
    rootNodeDetails->node = tree->root;
    rootNodeDetails->depth = 0;
    rootNodeDetails->centerX = 0.0;
    rootNodeDetails->centerY = 0.0;
    rootNodeDetails->centerZ = 0.0;
    rootNodeDetails->size = SIZE_LOOKUP_TABLE[0];
    rootNodeDetails->value = tree->root->logOdds;

    if (Stack__push(treeDetailStack, &rootNodeDetails))
    {
        printf("Failed to push root element onto treeDetailStack\n");
    }

    while (!Stack__isEmpty(treeDetailStack))
    {
        NodeDetails *nodeDetails = NULL;

        if (Stack__peekAndPop(treeDetailStack, &nodeDetails))
        {
            printf("Failed to pop element from treeDetailStack\n");
        }

        fprintf(outFile, "%u,%lg,%lg,%lg,%lg,%lg\n",
                nodeDetails->depth, nodeDetails->centerX, nodeDetails->centerY, nodeDetails->centerZ,
                nodeDetails->size, nodeDetails->value);

        for (ssize_t childIndex = 7; childIndex >= 0; --childIndex) {
            size_t octomapEquivalentIndex = FastOctree_to_octomap_index_map[childIndex];
            if (nodeDetails->node->children[octomapEquivalentIndex] != NO_CHILD) {
                NodeDetails *newNodeDetails = calloc(1, sizeof(NodeDetails));
                newNodeDetails->node = nodeDetails->node->children[octomapEquivalentIndex];
                newNodeDetails->depth = nodeDetails->depth + 1;
                newNodeDetails->size = SIZE_LOOKUP_TABLE[newNodeDetails->depth];
                newNodeDetails->value = newNodeDetails->node->logOdds;
                double newNodeHalfSize = newNodeDetails->size / 2.0;

                switch (octomapEquivalentIndex) {
                    case 0:
                        newNodeDetails->centerX = nodeDetails->centerX - newNodeHalfSize;
                        newNodeDetails->centerY = nodeDetails->centerY - newNodeHalfSize;
                        newNodeDetails->centerZ = nodeDetails->centerZ - newNodeHalfSize;
                        break;

                    case 1:
                        newNodeDetails->centerX = nodeDetails->centerX - newNodeHalfSize;
                        newNodeDetails->centerY = nodeDetails->centerY - newNodeHalfSize;
                        newNodeDetails->centerZ = nodeDetails->centerZ + newNodeHalfSize;
                        break;

                    case 2:
                        newNodeDetails->centerX = nodeDetails->centerX - newNodeHalfSize;
                        newNodeDetails->centerY = nodeDetails->centerY + newNodeHalfSize;
                        newNodeDetails->centerZ = nodeDetails->centerZ - newNodeHalfSize;
                        break;

                    case 3:
                        newNodeDetails->centerX = nodeDetails->centerX - newNodeHalfSize;
                        newNodeDetails->centerY = nodeDetails->centerY + newNodeHalfSize;
                        newNodeDetails->centerZ = nodeDetails->centerZ + newNodeHalfSize;
                        break;

                    case 4:
                        newNodeDetails->centerX = nodeDetails->centerX + newNodeHalfSize;
                        newNodeDetails->centerY = nodeDetails->centerY - newNodeHalfSize;
                        newNodeDetails->centerZ = nodeDetails->centerZ - newNodeHalfSize;
                        break;

                    case 5:
                        newNodeDetails->centerX = nodeDetails->centerX + newNodeHalfSize;
                        newNodeDetails->centerY = nodeDetails->centerY - newNodeHalfSize;
                        newNodeDetails->centerZ = nodeDetails->centerZ + newNodeHalfSize;
                        break;

                    case 6:
                        newNodeDetails->centerX = nodeDetails->centerX + newNodeHalfSize;
                        newNodeDetails->centerY = nodeDetails->centerY + newNodeHalfSize;
                        newNodeDetails->centerZ = nodeDetails->centerZ - newNodeHalfSize;
                        break;

                    case 7:
                        newNodeDetails->centerX = nodeDetails->centerX + newNodeHalfSize;
                        newNodeDetails->centerY = nodeDetails->centerY + newNodeHalfSize;
                        newNodeDetails->centerZ = nodeDetails->centerZ + newNodeHalfSize;
                        break;

                    default:
                        printf("Unexpected childIndex: %zd\n", octomapEquivalentIndex);
                }

                if (Stack__push(treeDetailStack, &newNodeDetails))
                {
                    printf("Failed to push element onto treeDetailStack\n");
                }
            }
        }

        free(nodeDetails);
    }

    fclose(outFile);
}

void insertPointCloud(Octree* tree, Vector3d* points, size_t numPoints, Vector3d* sensorOrigin)
{
    static int notInitialized = TRUE;
    if (notInitialized) {
        CLAMPING_THRES_MIN = logodds(0.1192);
        CLAMPING_THRES_MAX = logodds(0.971);
        PROB_HIT_LOG = logodds(0.7);
        PROB_MISS_LOG = logodds(0.4);
        OCC_PROB_THRES_LOG = logodds(0.5);
        notInitialized = FALSE;
    }

    // for each point, create ray and call ray_parameter. After all points done, prune the tree
    for (size_t i = 0; i < numPoints; ++i)
    {
        Ray currentRay;
        initRay(&currentRay,
                sensorOrigin->x, sensorOrigin->y, sensorOrigin->z,
                points[i].x, points[i].y, points[i].z);

        ray_parameter(tree, &currentRay);

        //printTree(tree);
        //printf("Finished ray %zu. Press Any Key to Continue\n", i);
        //getchar();

        //pruneTree(tree);

        //printTree(tree);
        //printf("Pruned tree. Press Any Key to Continue\n", i);
        //getchar();

        //printf("Finished ray %zu\n", i);
    }

    printf("Done inserting all rays, now pruning\n");

    printf("tree size: %zu\n", treeSize(tree));
    //printTree(tree);

    //printf("Press Any Key to Continue\n");
    //getchar();

    pruneTree(tree);

    //printTree(tree);
    //printf("tree size: %zu\n", treeSize(tree));
    //printf("Press Any Key to Continue\n");
    //getchar();
}