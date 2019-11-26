//
// Created by mschnur on 11/12/19.
//

#ifndef FASTOCTREE_FASTOCTREE_H
#define FASTOCTREE_FASTOCTREE_H

extern const Node* const NO_CHILD = NULL;

extern const Node* const NO_CHILDREN[8];

extern const double resolution;
extern const double clamping_thres_min;
extern const double clamping_thres_max;
extern const double prob_hit_log;
extern const double prob_miss_log;
extern const double occ_prob_thres_log;

struct Node;

typedef struct Node {
    double logOdds;
    struct Node *children[8];
} Node;

typedef struct Vector3d {
    double x, y, z;
} Vector3d;

typedef struct Octree {
    Vector3d min, max, size;
    Node* root;


} Octree;

typedef struct Ray {
    Vector3d origin, end, direction;
} Ray;

#endif //FASTOCTREE_FASTOCTREE_H
