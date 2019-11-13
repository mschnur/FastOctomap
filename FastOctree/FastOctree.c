//
// Created by mschnur on 11/12/19.
//

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "FastOctree.h"

#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#define TRUE 1
#define FALSE 0

const Node* const NO_CHILD = NULL;

const Node* const NO_CHILDREN[8] = {
        NO_CHILD, NO_CHILD, NO_CHILD, NO_CHILD, NO_CHILD, NO_CHILD, NO_CHILD, NO_CHILD
};

const double resolution = 0.05;
const double clamping_thres_min = -2.0;
const double clamping_thres_max = 3.5;
const double prob_hit_log = 0.85;
const double prob_miss_log = -0.4;
const double occ_prob_thres_log = 0.0;
const double max_depth = 16;

//###########################################################################################
// Utility functions
//###########################################################################################

inline int nodeHasChildren(const Node* n)
{
    return (memcmp(n->children, NO_CHILDREN, 8 * sizeof(Node*)) != 0);
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
                  Vector3d* t_endpoint, unsigned int depth,
                  Node* n, unsigned char a, Ray* r)
{
    printf("proc_subtree\n");
    double txm, tym, tzm;
    int currentNode;

    if (tx1 < 0.0 || ty1 < 0.0 || tz1 < 0.0)
    {
        return;
    }

    if (!nodeHasChildren(n))
    {
        if (is_less(t_endpoint->x, t_endpoint->y, t_endpoint->z, tx0, ty0, tz0))
        {
            // The ray endpoint happened before at least one of the minimum t values of this subtree, meaning the ray
            // ends before this subtree. Therefore, we don't want to do anything to this subtree.
            return;
        }
        else if (depth == max_depth)
        {
            double logLikelihoodUpdate = 0.0;

            if (is_greater(t_endpoint->x, t_endpoint->y, t_endpoint->z, tx1, ty1, tz1))
            {
                // The ray endpoint does not occur in this subtree, but the ray passes through this subtree on its
                // way to the endpoint, and we're at our maximum depth. Therefore we need to give this node a vote
                // that it is free space.
                logLikelihoodUpdate = prob_miss_log;
            }
            else
            {
                // The ray endpoint occurs within this subtree, and we're at our maximum depth. Therefore we need to
                // give this node a vote that it is free space.
                logLikelihoodUpdate = prob_hit_log;
            }

            // Do the update
            n->logOdds += logLikelihoodUpdate;

            // Clamp the logOdds between the min/max
            n->logOdds = fmax(clamping_thres_min, fmin(n->logOdds, clamping_thres_max));
        }
        else
        {
            // Either the ray endpoint occurs within this subtree or the ray passes through this subtree on its way to
            // the endpoint, but we're not at our maximum depth yet so we want to expand this node and then proceed
            // as normal.

            // TODO: expand the node
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
                proc_subtree(tx0, ty0, tz0, txm, tym, tzm, n->children[a], a, r);
                currentNode = new_node(txm, 4, tym, 2, tzm, 1);
                break;

            case 1:
                proc_subtree(tx0, ty0, tzm, txm, tym, tz1, n->children[1u^a], a, r);
                currentNode = new_node(txm, 5, tym, 3, tz1, 8);
                break;

            case 2:
                proc_subtree(tx0, tym, tz0, txm, ty1, tzm, n->children[2u^a], a, r);
                currentNode = new_node(txm, 6, ty1, 8, tzm, 3);
                break;

            case 3:
                proc_subtree(tx0, tym, tzm, txm, ty1, tz1, n->children[3u^a], a, r);
                currentNode = new_node(txm, 7, ty1, 8, tz1, 8);
                break;

            case 4:
                proc_subtree(txm, ty0, tz0, tx1, tym, tzm, n->children[4u^a], a, r);
                currentNode = new_node(tx1, 8, tym, 6, tzm, 5);
                break;

            case 5:
                proc_subtree(txm, ty0, tzm, tx1, tym, tz1, n->children[5u^a], a, r);
                currentNode = new_node(tx1, 8, tym, 7, tz1, 8);
                break;

            case 6:
                proc_subtree(txm, tym, tz0, tx1, ty1, tzm, n->children[6u^a], a, r);
                currentNode = new_node(tx1, 8, ty1, 8, tzm, 7);
                break;

            case 7:
                proc_subtree(txm, tym, tzm, tx1, ty1, tz1, n->children[7u^a], a, r);
                currentNode = 8;
                break;

            default:
                assert(0);
        }
    } while (currentNode < 8);

    // TODO: If all of this node's children have the same log-likelihood and have no children themselves, then we can
    // TODO: prune these children, meaning we delete them.

    // TODO: Set the occupancy of this node to the maximum of the occupancy of its children
}



void ray_parameter(Octree* tree, Ray* r) {
    if (tree->root == NULL){
        // TODO: construct the root
    }

    printf("ray_parameter\n");
    unsigned char a = 0;

    if (r->direction.x < 0.0f) {
        r->origin.x = tree->size.x - r->origin.x;
        r->direction.x = -r->direction.x;
        a |= 4;
    }

    if (r->direction.y < 0.0f) {
        r->origin.y = tree->size.y - r->origin.y;
        r->direction.y = -r->direction.y;
        a |= 2;
    }

    if (r->direction.z < 0.0f) {
        r->origin.z = tree->size.z - r->origin.z;
        r->direction.z = -r->direction.z;
        a |= 1;
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
        proc_subtree(tx0, ty0, tz0, tx1, ty1, tz1, tree->root, a, r);
    }
}