#include <stdio.h>
#include <stdlib.h>

typedef struct Node {
    double logOdds;
    struct Node *children[8];
} Node;

typedef struct Vector3d {
    double x, y, z;
} Vector3d;

typedef struct Ray {
    Vector3d origin, end, direction;
    double t_end;
} Ray;

#define MAX_DEPTH 2
#define N 5
#include "lists_of_params1.c"

#define num_rays 2
int main(){
    double tx0[num_rays] = {0.0,0.1};
    double ty0[num_rays] = {0.0,0.1};
    double tz0[num_rays] = {0.0,0.1};
    double tx1[num_rays] = {1.0,0.8};
    double ty1[num_rays] = {1.0,1.0};
    double tz1[num_rays] = {1.0,0.25};

    double endpoints[num_rays] = {0.9,0.9};

    proc_subtree(tx0,ty0,tz0,
                 tx1,ty1,tz1,
                 num_rays,0,
                 NULL,0,endpoints);
}