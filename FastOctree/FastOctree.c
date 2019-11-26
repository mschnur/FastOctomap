//
// Created by mschnur on 11/12/19.
//

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pthread.h>

#include "FastOctree.h"
#include "NodeStack.h"

#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))

#define TRUE 1
#define FALSE 0

#define RESOLUTION 0.05
#define MAX_DEPTH 16u

const Node* const NO_CHILD = NULL;

const Node* const NO_CHILDREN[8] = {
        NO_CHILD, NO_CHILD, NO_CHILD, NO_CHILD, NO_CHILD, NO_CHILD, NO_CHILD, NO_CHILD
};

const double CLAMPING_THRES_MIN = -2.0;
const double CLAMPING_THRES_MAX = 3.5;
const double PROB_HIT_LOG = 0.85;
const double PROB_MISS_LOG = -0.4;
const double OCC_PROB_THRES_LOG = 0.0;

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

//###########################################################################################
// Utility functions
//###########################################################################################

inline int nodeHasAnyChildren(const Node* n)
{
    return (memcmp(n->children, NO_CHILDREN, 8 * sizeof(Node*)) != 0);
}

inline int nodeIsPrunable(const Node* n)
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

inline void deleteChild(Node* n, unsigned int childIndex)
{
    free(n->children[childIndex]);
    n->children[childIndex] = NULL;
}


inline void deleteAllChildren(Node* n)
{
    for (unsigned int i = 0; i < 8; ++i)
    {
        deleteChild(n, i);
    }
}

inline void createChild(Node* n, unsigned int childIndex)
{
    // In order for us to be able to free the nodes individually (which may happen during pruning) we need to
    // allocate each node separately, instead of the whole array at once. Using calloc instead of malloc
    // initializes the memory to zero, which means that the this new node's `children` array will be filled with
    // zeros (which is equivalent to filling it will NULL pointers).
    n->children[childIndex] = (Node*) calloc(1, sizeof(Node));

    // The children should all start with the same log odds as their parent.
    n->children[childIndex]->logOdds = n->logOdds;
}

inline void createChildIfItDoesntExist(Node* n, unsigned int childIndex)
{
    if (n->children[childIndex] == NO_CHILD)
    {
        createChild(n, childIndex);
    }
}

inline double maxChildLogLikelihood(const Node* n)
{
    double maxLogLikelihood = -INFINITY;
    for (unsigned int i = 0; i < 8; ++i)
    {
        Node* child = n->children[i];
        if (child != NO_CHILD && child->logOdds > maxLogLikelihood)
        {
            maxLogLikelihood = child->logOdds;
        }
    }

    return maxLogLikelihood;
}

inline void expandNode(Node* n)
{
    assert(!nodeHasAnyChildren(n));

    for (unsigned int i = 0; i < 8; ++i)
    {
        createChild(n, i);
    }
}

inline int is_less(double pointX, double pointY, double pointZ,
                   double minX, double minY, double minZ)
{
    return (minX > pointX && minY > pointY && minZ > pointZ);
}

inline int is_greater(double pointX, double pointY, double pointZ,
                      double maxX, double maxY, double maxZ)
{
    return (pointX > maxX && pointY > maxY && pointZ > maxZ);
}

inline int is_between(double minX, double minY, double minZ,
                      double pointX, double pointY, double pointZ,
                      double maxX, double maxY, double maxZ)
{
    return (minX <= pointX && pointX <= maxX &&
            minY <= pointY && pointY <= maxY &&
            minZ <= pointZ && pointZ <= maxZ);
}

inline int is_between_vector3d(const Vector3d* min,
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

int first_node(double tx0, double ty0, double tz0, double txm, double tym, double tzm)
{
    unsigned char answer = 0;	// initialize to 00000000
    // select the entry plane and set bits
    if(tx0 > ty0)
    {
        if(tx0 > tz0)
        {   // PLANE YZ
            if(tym < tx0) answer|=2u;	// set bit at position 1
            if(tzm < tx0) answer|=1u;	// set bit at position 0
            return (int) answer;
        }
    }
    else
    {
        if(ty0 > tz0)
        { // PLANE XZ
            if(txm < ty0) answer|=4u;	// set bit at position 2
            if(tzm < ty0) answer|=1u;	// set bit at position 0
            return (int) answer;
        }
    }
    // PLANE XY
    if(txm < tz0) answer|=4u;	// set bit at position 2
    if(tym < tz0) answer|=2u;	// set bit at position 1
    return (int) answer;
}


int new_node(double txm, int x, double tym, int y, double tzm, int z)
{
    if(txm < tym)
    {
        if(txm < tzm){return x;}  // YZ plane
    }
    else
    {
        if(tym < tzm){return y;} // XZ plane
    }

    return z; // XY plane;
}


void proc_subtree(double tx0, double ty0, double tz0,
                  double tx1, double ty1, double tz1,
                  unsigned int depth,
                  Node* n, unsigned char a, Ray* r)
{
    printf("proc_subtree\n");
    double txm, tym, tzm;
    int currentNode;

    if (tx1 < 0.0 || ty1 < 0.0 || tz1 < 0.0)
    {
        return;
    }

    if (!nodeHasAnyChildren(n))
    {
        if (is_less(r->t_end, r->t_end, r->t_end, tx0, ty0, tz0))
        {
            // The ray endpoint happened before at least one of the minimum t values of this subtree, meaning the ray
            // ends before this subtree. Therefore, we don't want to do anything to this subtree.
            return;
        }
        else if (depth == MAX_DEPTH)
        {
            double logLikelihoodUpdate = 0.0;

            if (is_greater(r->t_end, r->t_end, r->t_end, tx1, ty1, tz1))
            {
                // The ray endpoint does not occur in this subtree, but the ray passes through this subtree on its
                // way to the endpoint, and we're at our maximum depth. Therefore we need to give this node a vote
                // that it is free space.
                logLikelihoodUpdate = PROB_MISS_LOG;
            }
            else
            {
                // The ray endpoint occurs within this subtree, and we're at our maximum depth. Therefore we need to
                // give this node a vote that it is free space.
                logLikelihoodUpdate = PROB_HIT_LOG;
            }

            // Do the update
            n->logOdds += logLikelihoodUpdate;

            // Clamp the logOdds between the min/max
            n->logOdds = fmax(CLAMPING_THRES_MIN, fmin(n->logOdds, CLAMPING_THRES_MAX));
        }
        else
        {
            // Either the ray endpoint occurs within this subtree or the ray passes through this subtree on its way to
            // the endpoint, but we're not at our maximum depth yet so we want to expand this node and then proceed
            // as normal.
            expandNode(n);
        }
    }

    txm = 0.5 * (tx0 + tx1);
    tym = 0.5 * (ty0 + ty1);
    tzm = 0.5 * (tz0 + tz1);

    currentNode = first_node(tx0, ty0, tz0, txm, tym, tzm);
    do
    {
        switch (currentNode)
        {
            case 0:
                createChildIfItDoesntExist(n, a);
                proc_subtree(tx0, ty0, tz0, txm, tym, tzm, depth + 1, n->children[a], a, r);
                currentNode = new_node(txm, 4, tym, 2, tzm, 1);
                break;

            case 1:
                createChildIfItDoesntExist(n, 1u^a);
                proc_subtree(tx0, ty0, tzm, txm, tym, tz1, depth + 1, n->children[1u^a], a, r);
                currentNode = new_node(txm, 5, tym, 3, tz1, 8);
                break;

            case 2:
                createChildIfItDoesntExist(n, 2u^a);
                proc_subtree(tx0, tym, tz0, txm, ty1, tzm, depth + 1, n->children[2u^a], a, r);
                currentNode = new_node(txm, 6, ty1, 8, tzm, 3);
                break;

            case 3:
                createChildIfItDoesntExist(n, 3u^a);
                proc_subtree(tx0, tym, tzm, txm, ty1, tz1, depth + 1, n->children[3u^a], a, r);
                currentNode = new_node(txm, 7, ty1, 8, tz1, 8);
                break;

            case 4:
                createChildIfItDoesntExist(n, 4u^a);
                proc_subtree(txm, ty0, tz0, tx1, tym, tzm, depth + 1, n->children[4u^a], a, r);
                currentNode = new_node(tx1, 8, tym, 6, tzm, 5);
                break;

            case 5:
                createChildIfItDoesntExist(n, 5u^a);
                proc_subtree(txm, ty0, tzm, tx1, tym, tz1, depth + 1, n->children[5u^a], a, r);
                currentNode = new_node(tx1, 8, tym, 7, tz1, 8);
                break;

            case 6:
                createChildIfItDoesntExist(n, 6u^a);
                proc_subtree(txm, tym, tz0, tx1, ty1, tzm, depth + 1, n->children[6u^a], a, r);
                currentNode = new_node(tx1, 8, ty1, 8, tzm, 7);
                break;

            case 7:
                createChildIfItDoesntExist(n, 7u^a);
                proc_subtree(txm, tym, tzm, tx1, ty1, tz1, depth + 1, n->children[7u^a], a, r);
                currentNode = 8;
                break;

            default:
                assert(0);
        }
    } while (currentNode < 8);

    n->logOdds = maxChildLogLikelihood(n);
}


void ray_parameter(Octree* tree, Ray* r) {
    if (tree->root == NULL){
        // Using calloc instead of malloc initializes the memory to zero, which means that the the root's `children`
        // array will be filled with zeros (which is equivalent to filling it will NULL pointers). It also sets the
        // root's log odds to 0, which we want as our initial value.
        tree->root = (Node*) calloc(1, sizeof(Node));
    }

    printf("ray_parameter\n");
    unsigned char a = 0;

    if (r->direction.x < 0.0f) {
        r->origin.x = tree->size.x - r->origin.x;
        r->direction.x = -r->direction.x;
        a |= 4u;
    }

    if (r->direction.y < 0.0f) {
        r->origin.y = tree->size.y - r->origin.y;
        r->direction.y = -r->direction.y;
        a |= 2u;
    }

    if (r->direction.z < 0.0f) {
        r->origin.z = tree->size.z - r->origin.z;
        r->direction.z = -r->direction.z;
        a |= 1u;
    }

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

    printf("txyz0: %lf %lf %lf\n", tx0, ty0, tz0);
    printf("txyz1: %lf %lf %lf\n", tx1, ty1, tz1);

    if (MAX(MAX(tx0, ty0), tz0) < MIN(MIN(tx1, ty1), tz1))
    {
        proc_subtree(tx0, ty0, tz0, tx1, ty1, tz1, 0, tree->root, a, r);
    }
}

inline void pruneNode(Node* node)
{
    if (nodeIsPrunable(node))
    {
        deleteAllChildren(node);
    }
}

void pruneTree(Octree* tree)
{
    if (tree->root == NULL)
    {
        return;
    }

    static NodeStack preOrderStack = { 0 };
    static NodeStack postOrderStack = { 0 };

    init(&preOrderStack);
    init(&postOrderStack);

    push(&preOrderStack, tree->root);

    while (!isEmpty(&preOrderStack))
    {
        Node* popped;

        if (!pop(&preOrderStack, &popped))
        {
            printf("Failed to pop element from preOrderStack\n");
        }

        if (!push(&postOrderStack, popped))
        {
            printf("Failed to push element onto postOrderStack\n");
        }

        for (unsigned int i = 0; i < 8; ++i)
        {
            Node* child = popped->children[i];
            if (child != NO_CHILD)
            {
                if (!push(&preOrderStack, popped->children[i]))
                {
                    printf("Failed to push element onto preOrderStack\n");
                }
            }
        }
    }

    while (!isEmpty(&postOrderStack))
    {
        Node* popped;

        if (!pop(&postOrderStack, &popped))
        {
            printf("Failed to pop element from postOrderStack\n");
        }

        pruneNode(popped);
    }
}

void createNodeCsv(Octree* tree)
{
    if (tree->root == NULL)
    {
        return;
    }

    static NodeStack preOrderStack = { 0 };
    static NodeStack postOrderStack = { 0 };

    init(&preOrderStack);
    init(&postOrderStack);

    push(&preOrderStack, tree->root);

    while (!isEmpty(&preOrderStack))
    {
        Node* popped;

        if (!pop(&preOrderStack, &popped))
        {
            printf("Failed to pop element from preOrderStack\n");
        }

        if (!push(&postOrderStack, popped))
        {
            printf("Failed to push element onto postOrderStack\n");
        }

        for (unsigned int i = 0; i < 8; ++i)
        {
            Node* child = popped->children[i];
            if (child != NO_CHILD)
            {
                if (!push(&preOrderStack, popped->children[i]))
                {
                    printf("Failed to push element onto preOrderStack\n");
                }
            }
        }
    }

    while (!isEmpty(&postOrderStack))
    {
        Node* popped;

        if (!pop(&postOrderStack, &popped))
        {
            printf("Failed to pop element from postOrderStack\n");
        }

        pruneNode(popped);
    }
}

void insertPointCloud(Octree* tree, Vector3d* points, size_t numPoints)
{
    // for each point, create ray and call ray_parameter. After all points done, prune the tree
}