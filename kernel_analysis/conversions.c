#include <stdio.h>

unsigned char first_node_one(double tx0, double ty0, double tz0, double txm, double tym, double tzm){
    unsigned char currentNode = 0;
    
    // select the entry plane and set bits
    if(tx0 > ty0) {
        if(tx0 > tz0) {   // PLANE YZ
            if(tym < tx0) currentNode|=2u;	// set bit at position 1
            if(tzm < tx0) currentNode|=1u;	// set bit at position 0
        }
    }
    else if(ty0 > tz0) { 
            // PLANE XZ
            if(txm < ty0) currentNode|=4u;	// set bit at position 2
            if(tzm < ty0) currentNode|=1u;	// set bit at position 0
    }
    else {
        // PLANE XY
        if(txm < tz0) currentNode|=4u;	// set bit at position 2
        if(tym < tz0) currentNode|=2u;	// set bit at position 1
    }
    return currentNode;
}

unsigned char first_node_two(double tx0, double ty0, double tz0, double txm, double tym, double tzm){
    unsigned long long currentNode = 0;
    
    // select the entry plane and set bits
    
    double tmp1 = ty0 - tx0;
    double tmp2 = tz0 - tx0;
    double tmp3 = tym - tx0;
    double tmp4 = tzm - tx0;

    double tmp5 = tz0 - ty0;
    double tmp6 = txm - ty0;
    double tmp7 = tzm - ty0;

    double tmp8 = txm - tz0;
    double tmp9 = tym - tz0;

    currentNode |= (unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>62;
    currentNode |= (unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & *((unsigned long long*)(&tmp4)) & 0x8000000000000000)>>63;

    currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp6)) & 0x8000000000000000)>>61;
    currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp7)) & 0x8000000000000000)>>63;

    currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & ~*((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp8)) & 0x8000000000000000)>>61;
    currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & ~*((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp9)) & 0x8000000000000000)>>62;

    return (unsigned char)currentNode;
}

unsigned int new_node_one(double txm, double tym, double tzm, unsigned char x, unsigned char y, unsigned char z){
    if(txm < tym)
    {
        if(txm < tzm){return (1u<<x);}  // YZ plane
    }
    else
    {
        if(tym < tzm){return (1u<<y);} // XZ plane
    }

    return (1u<<z); // XY plane;
}

unsigned int new_node_two(double txm, double tym, double tzm, unsigned char x, unsigned char y, unsigned char z) {
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

    return (unsigned int)currentNode;
}


int main(){
    unsigned char answer1 = 0, answer2 = 0;
    double tx0,ty0,tz0,txm,tym,tzm;

    //first_node tests
    printf("first_node tests\n");
    {
        {    
            tx0 = 2.0;
            ty0 = 1.0;
            tz0 = 0.0;

            txm = 0.5;
            tym = 0.5;
            tzm = 0.5;

            answer1 = first_node_one(tx0,ty0,tz0,txm,tym,tzm);
            answer2 = first_node_two(tx0,ty0,tz0,txm,tym,tzm);

            printf("%d:%d\n", (int)answer1,(int)answer2);
        }

        {
            tx0 = 2.0;
            ty0 = 1.0;
            tz0 = 0.0;

            txm = 0.5;
            tym = 0.5;
            tzm = 2.0;

            answer1 = first_node_one(tx0,ty0,tz0,txm,tym,tzm);
            answer2 = first_node_two(tx0,ty0,tz0,txm,tym,tzm);

            printf("%d:%d\n", (int)answer1,(int)answer2);
        }

        {
            tx0 = 2.0;
            ty0 = 2.0;
            tz0 = 0.0;

            txm = 0.5;
            tym = 0.5;
            tzm = 0.5;

            answer1 = first_node_one(tx0,ty0,tz0,txm,tym,tzm);
            answer2 = first_node_two(tx0,ty0,tz0,txm,tym,tzm);

            printf("%d:%d\n", (int)answer1,(int)answer2);
        }

        {
            tx0 = 2.0;
            ty0 = 2.0;
            tz0 = 0.0;

            txm = 0.5;
            tym = 0.5;
            tzm = 2.0;

            answer1 = first_node_one(tx0,ty0,tz0,txm,tym,tzm);
            answer2 = first_node_two(tx0,ty0,tz0,txm,tym,tzm);

            printf("%d:%d\n", (int)answer1,(int)answer2);
        }

        {
            tx0 = 2.0;
            ty0 = 2.0;
            tz0 = 2.0;

            txm = 0.5;
            tym = 0.5;
            tzm = 0.5;

            answer1 = first_node_one(tx0,ty0,tz0,txm,tym,tzm);
            answer2 = first_node_two(tx0,ty0,tz0,txm,tym,tzm);

            printf("%d:%d\n", (int)answer1,(int)answer2);
        }

        {
            tx0 = 2.0;
            ty0 = 2.0;
            tz0 = 2.0;

            txm = 2.0;
            tym = 0.5;
            tzm = 0.5;

            answer1 = first_node_one(tx0,ty0,tz0,txm,tym,tzm);
            answer2 = first_node_two(tx0,ty0,tz0,txm,tym,tzm);

            printf("%d:%d\n", (int)answer1,(int)answer2);
        }
    }
    printf("\n");

    //next_node tests
    printf("next_node tests\n");
    {
        {
            txm = 0.0;
            tym = 1.0;
            tzm = 1.0;

            answer1 = new_node_one(txm,tym,tzm,1,2,3);
            answer2 = new_node_two(txm,tym,tzm,1,2,3);

            printf("%d:%d\n", (int)answer1,(int)answer2);
        }

        {
            txm = 0.0;
            tym = 0.0;
            tzm = 1.0;

            answer1 = new_node_one(txm,tym,tzm,1,2,3);
            answer2 = new_node_two(txm,tym,tzm,1,2,3);

            printf("%d:%d\n", (int)answer1,(int)answer2);
        }

        {
            txm = 0.45;
            tym = 0.55;
            tzm = 0.175;

            answer1 = new_node_one(txm,tym,tzm,1,2,3);
            answer2 = new_node_two(txm,tym,tzm,1,2,3);
            
            printf("%d:%d\n", (int)answer1,(int)answer2);
        }

    }

}