//
// Created by mschnur on 11/12/19.
//

#ifndef FASTOCTREE_FASTOCTREE_H
#define FASTOCTREE_FASTOCTREE_H

#include <math.h>

struct Node;
typedef struct Node Node;

extern const Node* const NO_CHILD;

extern const Node* const NO_CHILDREN[8];

extern const double resolution;
extern const double CLAMPING_THRES_MIN;
extern const double CLAMPING_THRES_MAX;
extern const double prob_hit_log;
extern const double PROB_MISS_LOG;
extern const double OCC_PROB_THRES_LOG;

struct Node {
    double logOdds;
    struct Node *children[8];
};

typedef struct Vector3d {
    double x, y, z;
} Vector3d;

typedef struct Octree {
    Vector3d min, max, size;
    Node* root;


} Octree;

typedef struct Ray {
    Vector3d origin, end, direction;
    double t_end;
} Ray;

void createNodeCsv(Octree* tree);

void insertPointCloud(Octree* tree, Vector3d* points, size_t numPoints);

inline void initVector3d(Vector3d* v, double x, double y, double z)
{
    v->x = x;
    v->y = y;
    v->z = z;
}

inline void vectorSubtract(const Vector3d* a, const Vector3d* b, Vector3d* result)
{
    result->x = a->x - b->x;
    result->y = a->y - b->y;
    result->z = a->z - b->z;
}

inline void vectorAdd(const Vector3d* a, const Vector3d* b, Vector3d* result)
{
    result->x = a->x + b->x;
    result->y = a->y + b->y;
    result->z = a->z + b->z;
}

inline void vectorNormalizeInPlace(Vector3d* a)
{
    double magnitude = sqrt((a->x * a->x) + (a->y * a->y) + (a->z * a->z));
    a->x /= magnitude;
    a->y /= magnitude;
    a->z /= magnitude;
}

inline void initRay(Ray* r, double ox, double oy, double oz, double ex, double ey, double ez)
{
    initVector3d(&(r->origin), ox, oy, oz);
    initVector3d(&(r->end), ex, ey, ez);
    vectorSubtract(r->end, )
}



#endif //FASTOCTREE_FASTOCTREE_H