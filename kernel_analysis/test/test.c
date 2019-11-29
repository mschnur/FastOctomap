#include <stdio.h>
#include <stdlib.h>

#define MAX_DEPTH

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

// #include "original.c"
// #include "original_converted.c"
// #include "extracted_node_computation.c"
//#include "remove_node_switch.c"
#include "with_terminal_logic_converted.c"

int main(){
    double tx0 = 0.1;//0.0;
    double ty0 = 0.1;//0.0;
    double tz0 = 0.1;//0.0;
    double tx1 = 0.8;//1.0;
    double ty1 = 1.0;//1.0;
    double tz1 = 0.25;//1.0;


    Ray r = {
        {0,0,0},
        {0,0,0},
        {0,0,0},
        0.9
    };

    proc_subtree(tx0,ty0,tz0,tx1,ty1,tz1,0,NULL,0,&r);
}