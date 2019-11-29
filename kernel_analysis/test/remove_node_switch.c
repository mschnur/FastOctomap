#include <assert.h>

int new_node(double txm, unsigned char x, double tym, unsigned char y, double tzm, unsigned char z)
{
    unsigned long long currentNode = 0;

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


void proc_subtree(double tx0, double ty0, double tz0,
                  double tx1, double ty1, double tz1,
                  unsigned int depth,
                  Node* n, unsigned char a, Ray* r) {
    printf("proc_subtree\n");
    double txm, tym, tzm;
    int currentNode = 0;
    unsigned char nodes = 0; //each bit represents a node


    if (tx1 < 0.0 || ty1 < 0.0 || tz1 < 0.0) {
        return;
    }

    /****
     * LOGIC REMOVED HERE
     * //////

    /******
     * Compute all nodes the ray passes through
     ******/

    {
        txm = 0.5 * (tx0 + tx1);
        tym = 0.5 * (ty0 + ty1);
        tzm = 0.5 * (tz0 + tz1);


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
        
        currentNode |= 1<<currentNode;


        int eq = ~(((1<<0) - currentNode)>>31);
        nodes |= currentNode & eq;
        currentNode = (new_node(txm, 4, tym, 2, tzm, 1) & eq) | (currentNode & ~eq);

        eq = ~(((1<<1) - currentNode)>>31);
        nodes |= currentNode & eq;
        currentNode = (new_node(txm, 5, tym, 3, tz1, 8) & eq); // | (currentNode & ~eq);

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
    }

    /******
     * Process all computed nodes
     ******/
    printf("processing nodes\n");
    for (int node = 0; node < 8; node++){
        if(nodes & (1u<<node)){
            printf("%d\n", node);
            switch (node) {
                case 0:
                    // createChildIfItDoesntExist(n, a);
                    // proc_subtree(tx0, ty0, tz0, txm, tym, tzm, depth + 1, n->children[a], a, r);
                    break;

                case 1:
                    // createChildIfItDoesntExist(n, 1u^a);
                    // proc_subtree(tx0, ty0, tzm, txm, tym, tz1, depth + 1, n->children[1u^a], a, r);
                    break;

                case 2:
                    // createChildIfItDoesntExist(n, 2u^a);
                    // proc_subtree(tx0, tym, tz0, txm, ty1, tzm, depth + 1, n->children[2u^a], a, r);
                    break;

                case 3:
                    // createChildIfItDoesntExist(n, 3u^a);
                    // proc_subtree(tx0, tym, tzm, txm, ty1, tz1, depth + 1, n->children[3u^a], a, r);
                    break;

                case 4:
                    // createChildIfItDoesntExist(n, 4u^a);
                    // proc_subtree(txm, ty0, tz0, tx1, tym, tzm, depth + 1, n->children[4u^a], a, r);
                    break;

                case 5:
                    // createChildIfItDoesntExist(n, 5u^a);
                    // proc_subtree(txm, ty0, tzm, tx1, tym, tz1, depth + 1, n->children[5u^a], a, r);
                    break;

                case 6:
                    // createChildIfItDoesntExist(n, 6u^a);
                    // proc_subtree(txm, tym, tz0, tx1, ty1, tzm, depth + 1, n->children[6u^a], a, r);
                    break;

                case 7:
                    // createChildIfItDoesntExist(n, 7u^a);
                    // proc_subtree(txm, tym, tzm, tx1, ty1, tz1, depth + 1, n->children[7u^a], a, r);
                    break;

                default:
                    assert(0);
            }
        }
    }
}