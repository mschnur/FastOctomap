//
// Created by mschnur on 11/12/19.
//

#ifndef FASTOCTREE_FASTOCTREE_H
#define FASTOCTREE_FASTOCTREE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>

struct Node;
typedef struct Node Node;

extern const Node *const NO_CHILD;

extern const Node *const NO_CHILDREN[8];

#define RESOLUTION 0.05
#define MAX_DEPTH 16u

extern double CLAMPING_THRES_MIN;
extern double CLAMPING_THRES_MAX;
extern double PROB_HIT_LOG;
extern double PROB_MISS_LOG;
extern double OCC_PROB_THRES_LOG;

extern const double SIZE_LOOKUP_TABLE[MAX_DEPTH + 1];

struct Node {
    double logOdds;
    struct Node *children[8];
};

typedef struct Vector3d {
    double x, y, z;
} Vector3d;

typedef struct Octree {
    Vector3d min, max, size;
    Node *root;


} Octree;

typedef struct Ray {
    Vector3d origin, end, direction;
    double t_end;
} Ray;

void pruneTree(Octree* tree);

void createNodeCsv(Octree *tree, const char *filename);

void insertPointCloud(Octree *tree, Vector3d *points, size_t numPoints, Vector3d *sensorOrigin);

int check_ray_parameter_kernel_correctness(Octree* tree, Vector3d* points, 
											size_t numPoints, Vector3d* sensorOrigin);

void compare_ray_parameter_times(Octree* tree, Vector3d* points,
                                 size_t numPoints, Vector3d* sensorOrigin);

/// compute log-odds from probability:
static inline double logodds(double probability){
    return log(probability / (1.0 - probability));
}

/// compute probability from logodds:
inline double probability(double logodds){
    return 1.0 - ( 1.0 / (1.0 + exp(logodds)));

}

static inline void initVector3d(Vector3d *v, double x, double y, double z) {
    v->x = x;
    v->y = y;
    v->z = z;
}

static inline void vectorSubtract(const Vector3d *a, const Vector3d *b, Vector3d *result) {
    result->x = a->x - b->x;
    result->y = a->y - b->y;
    result->z = a->z - b->z;
}

static inline void vectorAdd(const Vector3d *a, const Vector3d *b, Vector3d *result) {
    result->x = a->x + b->x;
    result->y = a->y + b->y;
    result->z = a->z + b->z;
}

static inline void vectorNormalizeInPlace(Vector3d *a) {
    double magnitude = sqrt((a->x * a->x) + (a->y * a->y) + (a->z * a->z));
    a->x /= magnitude;
    a->y /= magnitude;
    a->z /= magnitude;
}

static inline void initRay(Ray *r, double ox, double oy, double oz, double ex, double ey, double ez) {
    initVector3d(&(r->origin), ox, oy, oz);
    initVector3d(&(r->end), ex, ey, ez);
    vectorSubtract(&(r->end), &(r->origin), &(r->direction));
    vectorNormalizeInPlace(&(r->direction));
    r->t_end = (ex - ox) / r->direction.x;
}

static inline void initOctree(Octree *tree) {
    double min = -(SIZE_LOOKUP_TABLE[0] / 2.0);
    double max = -min;
    initVector3d(&(tree->min), min, min, min);
    initVector3d(&(tree->max), max, max, max);
    initVector3d(&(tree->size), SIZE_LOOKUP_TABLE[0], SIZE_LOOKUP_TABLE[0], SIZE_LOOKUP_TABLE[0]);
    tree->root = NULL;
}

#ifdef __cplusplus
}
#endif

#endif //FASTOCTREE_FASTOCTREE_H
