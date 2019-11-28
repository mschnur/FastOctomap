

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
    currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & ~*((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>(63-z);

    return (unsigned int)currentNode;
}


void proc_subtree(double tx0, double ty0, double tz0,
                  double tx1, double ty1, double tz1,
                  unsigned int depth,
                  Node* n, unsigned char a, Ray* r) {
    printf("proc_subtree\n");
    double txm, tym, tzm;
    unsigned int currentNode = 0;
    unsigned char nodes = 0; //each bit represents a node


    if (tx1 < 0.0 || ty1 < 0.0 || tz1 < 0.0) {
        return;
    }

    //doesn't matter what depth we are at or if we have children
    if (is_less(r->t_end, r->t_end, r->t_end, tx0, ty0, tz0))
    {
        // The ray endpoint happened before at least one of the minimum t values of this subtree, meaning the ray
        // ends before this subtree. Therefore, we don't want to do anything to this subtree.
        return;
    }

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


        nodes |= 1u<<currentNode;
        currentNode = nodes;

        if(currentNode == (1u<<0)) {
            currentNode = new_node(txm, 4, tym, 2, tzm, 1);
            nodes |= currentNode;
        }

        if(currentNode == (1u<<1)) {
            currentNode = new_node(txm, 5, tym, 3, tz1, 8);
            nodes |= currentNode;
        }

        if(currentNode == (1u<<2)) {
            currentNode = new_node(txm, 6, ty1, 8, tzm, 3);
            nodes |= currentNode;
        }

        if(currentNode == (1u<<3)) {
            currentNode = new_node(txm, 7, ty1, 8, tz1, 8);
            nodes |= currentNode;
        }

        if(currentNode == (1u<<4)) {
            currentNode = new_node(tx1, 8, tym, 6, tzm, 5);
            nodes |= currentNode;
        }

        if(currentNode == (1u<<5)) {
            currentNode = new_node(tx1, 8, tym, 7, tz1, 8);
            nodes |= currentNode;
        }

        if(currentNode == (1u<<6)) {
            currentNode = new_node(tx1, 8, ty1, 8, tzm, 7);
            nodes |= currentNode;
        }

    }

    /******
     * Process all computed nodes
     ******/
    for (int node = 0; node < 8; node++){
        if(nodes && (1u<<node)){
            switch (node) {
                case 0:
                    createChildIfItDoesntExist(n, a);
                    proc_subtree(tx0, ty0, tz0, txm, tym, tzm, depth + 1, n->children[a], a, r);
                    break;

                case 1:
                    createChildIfItDoesntExist(n, 1u^a);
                    proc_subtree(tx0, ty0, tzm, txm, tym, tz1, depth + 1, n->children[1u^a], a, r);
                    break;

                case 2:
                    createChildIfItDoesntExist(n, 2u^a);
                    proc_subtree(tx0, tym, tz0, txm, ty1, tzm, depth + 1, n->children[2u^a], a, r);
                    break;

                case 3:
                    createChildIfItDoesntExist(n, 3u^a);
                    proc_subtree(tx0, tym, tzm, txm, ty1, tz1, depth + 1, n->children[3u^a], a, r);
                    break;

                case 4:
                    createChildIfItDoesntExist(n, 4u^a);
                    proc_subtree(txm, ty0, tz0, tx1, tym, tzm, depth + 1, n->children[4u^a], a, r);
                    break;

                case 5:
                    createChildIfItDoesntExist(n, 5u^a);
                    proc_subtree(txm, ty0, tzm, tx1, tym, tz1, depth + 1, n->children[5u^a], a, r);
                    break;

                case 6:
                    createChildIfItDoesntExist(n, 6u^a);
                    proc_subtree(txm, tym, tz0, tx1, ty1, tzm, depth + 1, n->children[6u^a], a, r);
                    break;

                case 7:
                    createChildIfItDoesntExist(n, 7u^a);
                    proc_subtree(txm, tym, tzm, tx1, ty1, tz1, depth + 1, n->children[7u^a], a, r);
                    break;

                default:
                    assert(0);
            }
        }
    }


    //doesn't matter if we have children if we are at max depth
    if (depth == MAX_DEPTH)
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
    //we should always have a child unless we are at max depth
    else {
        n->logOdds = maxChildLogLikelihood(n);
    }

}