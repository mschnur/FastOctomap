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

void proc_subtree(double* tx0, double* ty0, double* tz0,
                  double* tx1, double* ty1, double* tz1,
                  int num_rays, unsigned int depth,
                  Node* n, unsigned char a, double* endpoints) {
    printf("proc_subtree; depth: %d\n", depth);

    //assume we know how many are in this
    //consider this implementation psuedo-code so far
    int cur_index[8] = {0};
    double new_tx0[8][N] = {{0},{0},{0},{0},{0},{0},{0},{0}};
    double new_ty0[8][N] = {{0},{0},{0},{0},{0},{0},{0},{0}};   
    double new_tz0[8][N] = {{0},{0},{0},{0},{0},{0},{0},{0}};
    double new_tx1[8][N] = {{0},{0},{0},{0},{0},{0},{0},{0}};
    double new_ty1[8][N] = {{0},{0},{0},{0},{0},{0},{0},{0}};
    double new_tz1[8][N] = {{0},{0},{0},{0},{0},{0},{0},{0}};
    double new_endpoints[8][N] = {{0},{0},{0},{0},{0},{0},{0},{0}};
    //unsigned char* new_a[8][N] = {{0},{0},{0},{0},{0},{0},{0},{0}};


    // this should never happen?
    // if (tx1 < 0.0 || ty1 < 0.0 || tz1 < 0.0) {
    //     return;
    // }

    //if we are at max depth there wont be any children
    if (depth == MAX_DEPTH) {
        /*
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
            // give this node a vote that it is occupied.
            logLikelihoodUpdate = PROB_HIT_LOG;
        }

        // Do the update
        n->logOdds += logLikelihoodUpdate;

        // Clamp the logOdds between the min/max
        n->logOdds = fmax(CLAMPING_THRES_MIN, fmin(n->logOdds, CLAMPING_THRES_MAX));
        */
        return;
    }

    /******
     * Compute all nodes the ray passes through
     ******/

    for(int i = 0; i < num_rays; i++) {
        double txm, tym, tzm;
        int currentNode = 0;

        double t1,t2,t3;
        int eq;
        long long valid_node;

        printf("ray %d\n", i);
        txm = 0.5 * (tx0[i] + tx1[i]);
        tym = 0.5 * (ty0[i] + ty1[i]);
        tzm = 0.5 * (tz0[i] + tz1[i]);

        //first node should always be valid, assuming we are in a valid node
        // select the entry plane and set bits
        double tmp1 = ty0[i] - tx0[i];
        double tmp2 = tz0[i] - tx0[i];
        double tmp3 = tym - tx0[i];
        double tmp4 = tzm - tx0[i];

        double tmp5 = tz0[i] - ty0[i];
        double tmp6 = txm - ty0[i];
        double tmp7 = tzm - ty0[i];

        double tmp8 = txm - tz0[i];
        double tmp9 = tym - tz0[i];

        currentNode |= (unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>62;
        currentNode |= (unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & *((unsigned long long*)(&tmp4)) & 0x8000000000000000)>>63;

        currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp6)) & 0x8000000000000000)>>61;
        currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp7)) & 0x8000000000000000)>>63;

        currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & ~*((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp8)) & 0x8000000000000000)>>61;
        currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & ~*((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp9)) & 0x8000000000000000)>>62;


        currentNode |= 1u<<currentNode;

        eq = ~(((1<<0) - currentNode)>>31);
        //this should compute if the currentNode is valid.
        t1 = tx0[i] - endpoints[i];
        t2 = ty0[i] - endpoints[i];
        t3 = tz0[i] - endpoints[i];
        valid_node = (long long)(*((unsigned long long*)(&t1)) & *((unsigned long long*)(&t2)) & *((unsigned long long*)(&t3)) & 0x8000000000000000)>>63;
        new_tx0[0][cur_index[0]] = tx0[i];
        new_ty0[0][cur_index[0]] = ty0[i];
        new_tz0[0][cur_index[0]] = tz0[i];
        new_tx1[0][cur_index[0]] = txm;
        new_ty1[0][cur_index[0]] = tym;
        new_tz1[0][cur_index[0]] = tzm;
        new_endpoints[0][cur_index[0]] = endpoints[i];
        cur_index[0] = ((cur_index[0])+1 & valid_node & eq) | (cur_index[0] & valid_node & ~eq);
        currentNode = (new_node(txm, 4, tym, 2, tzm, 1) & eq) | (currentNode & ~eq);

        eq = ~(((1<<1) - currentNode)>>31);
        t1 = tx0[i] - endpoints[i];
        t2 = ty0[i] - endpoints[i];
        t3 = tzm - endpoints[i];
        valid_node = (long long)(*((unsigned long long*)(&t1)) & *((unsigned long long*)(&t2)) & *((unsigned long long*)(&t3)) & 0x8000000000000000)>>63;
        new_tx0[1][cur_index[1]] = tx0[i];
        new_ty0[1][cur_index[1]] = ty0[i];
        new_tz0[1][cur_index[1]] = tzm;
        new_tx1[1][cur_index[1]] = txm;
        new_ty1[1][cur_index[1]] = tym;
        new_tz1[1][cur_index[1]] = tz1[i];
        new_endpoints[1][cur_index[1]] = endpoints[i];
        cur_index[1] = ((cur_index[1])+1 & valid_node & eq) | (cur_index[1] & valid_node & ~eq);
        currentNode = (new_node(txm, 5, tym, 3, tz1[i], 8) & eq) | (currentNode & ~eq);

        eq = ~(((1<<2) - currentNode)>>31);
        t1 = tx0[i] - endpoints[i];
        t2 = tym - endpoints[i];
        t3 = tz0[i] - endpoints[i];
        valid_node = (long long)(*((unsigned long long*)(&t1)) & *((unsigned long long*)(&t2)) & *((unsigned long long*)(&t3)) & 0x8000000000000000)>>63;
        new_tx0[2][cur_index[2]] = tx0[i];
        new_ty0[2][cur_index[2]] = tym;
        new_tz0[2][cur_index[2]] = tz0[i];
        new_tx1[2][cur_index[2]] = txm;
        new_ty1[2][cur_index[2]] = ty1[i];
        new_tz1[2][cur_index[2]] = tzm;
        new_endpoints[2][cur_index[2]] = endpoints[i];
        cur_index[2] = ((cur_index[2])+1 & valid_node & eq) | (cur_index[2] & valid_node & ~eq);
        currentNode = (new_node(txm, 6, ty1[i], 8, tzm, 3) & eq) | (currentNode & ~eq);

        eq = ~(((1<<3) - currentNode)>>31);
        t1 = tx0[i] - endpoints[i];
        t2 = tym - endpoints[i];
        t3 = tzm - endpoints[i];
        valid_node = (long long)(*((unsigned long long*)(&t1)) & *((unsigned long long*)(&t2)) & *((unsigned long long*)(&t3)) & 0x8000000000000000)>>63;
        new_tx0[3][cur_index[3]] = tx0[i];
        new_ty0[3][cur_index[3]] = tym;
        new_tz0[3][cur_index[3]] = tzm;
        new_tx1[3][cur_index[3]] = txm;
        new_ty1[3][cur_index[3]] = ty1[i];
        new_tz1[3][cur_index[3]] = tz1[i];
        new_endpoints[3][cur_index[3]] = endpoints[i];
        cur_index[3] = ((cur_index[3])+1 & valid_node & eq) | (cur_index[3] & valid_node & ~eq);
        currentNode = (new_node(txm, 7, ty1[i], 8, tz1[i], 8) & eq) | (currentNode & ~eq);

        eq = ~(((1<<4) - currentNode)>>31);
        t1 = txm - endpoints[i];
        t2 = ty0[i] - endpoints[i];
        t3 = tz0[i] - endpoints[i];
        valid_node = (long long)(*((unsigned long long*)(&t1)) & *((unsigned long long*)(&t2)) & *((unsigned long long*)(&t3)) & 0x8000000000000000)>>63;
        new_tx0[4][cur_index[4]] = txm;
        new_ty0[4][cur_index[4]] = ty0[i];
        new_tz0[4][cur_index[4]] = tz0[i];
        new_tx1[4][cur_index[4]] = tx1[i];
        new_ty1[4][cur_index[4]] = tym;
        new_tz1[4][cur_index[4]] = tzm;
        new_endpoints[4][cur_index[4]] = endpoints[i];
        cur_index[4] = ((cur_index[4])+1 & valid_node & eq) | (cur_index[4] & valid_node & ~eq);
        currentNode = (new_node(tx1[i], 8, tym, 6, tzm, 5) & eq) | (currentNode & ~eq);

        eq = ~(((1<<5) - currentNode)>>31);
        t1 = txm - endpoints[i];
        t2 = ty0[i] - endpoints[i];
        t3 = tzm - endpoints[i];
        valid_node = (long long)(*((unsigned long long*)(&t1)) & *((unsigned long long*)(&t2)) & *((unsigned long long*)(&t3)) & 0x8000000000000000)>>63;
        new_tx0[5][cur_index[5]] = txm;
        new_ty0[5][cur_index[5]] = ty0[i];
        new_tz0[5][cur_index[5]] = tzm;
        new_tx1[5][cur_index[5]] = tx1[i];
        new_ty1[5][cur_index[5]] = tym;
        new_tz1[5][cur_index[5]] = tz1[i];
        new_endpoints[5][cur_index[5]] = endpoints[i];
        cur_index[5] = ((cur_index[5])+1 & valid_node & eq) | (cur_index[5] & valid_node & ~eq);
        currentNode = (new_node(tx1[i], 8, tym, 7, tz1[i], 8) & eq) | (currentNode & ~eq);

        eq = ~(((1<<6) - currentNode)>>31);
        t1 = txm - endpoints[i];
        t2 = tym - endpoints[i];
        t3 = tz0[i] - endpoints[i];
        valid_node = (long long)(*((unsigned long long*)(&t1)) & *((unsigned long long*)(&t2)) & *((unsigned long long*)(&t3)) & 0x8000000000000000)>>63;
        new_tx0[6][cur_index[6]] = txm;
        new_ty0[6][cur_index[6]] = tym;
        new_tz0[6][cur_index[6]] = tz0[i];
        new_tx1[6][cur_index[6]] = tx1[i];
        new_ty1[6][cur_index[6]] = ty1[i];
        new_tz1[6][cur_index[6]] = tzm;
        new_endpoints[6][cur_index[6]] = endpoints[i];
        cur_index[6] = ((cur_index[6])+1 & valid_node & eq) | (cur_index[6] & valid_node & ~eq);
        currentNode = (new_node(tx1[i], 8, ty1[i], 8, tzm, 7) & eq) | (currentNode & ~eq);

        eq = ~(((1<<7) - currentNode)>>31);
        t1 = txm - endpoints[i];
        t2 = tym - endpoints[i];
        t3 = tzm - endpoints[i];
        valid_node = (long long)(*((unsigned long long*)(&t1)) & *((unsigned long long*)(&t2)) & *((unsigned long long*)(&t3)) & 0x8000000000000000)>>63;
        new_tx0[7][cur_index[7]] = txm;
        new_ty0[7][cur_index[7]] = tym;
        new_tz0[7][cur_index[7]] = tzm;
        new_tx1[7][cur_index[7]] = tx1[i];
        new_ty1[7][cur_index[7]] = ty1[i];
        new_tz1[7][cur_index[7]] = tz1[i];
        new_endpoints[7][cur_index[7]] = endpoints[i];
        cur_index[7] = ((cur_index[7])+1 & valid_node & eq) | (cur_index[7] & valid_node & ~eq);
    }


    /******
     * Process all computed nodes
     ******/
    for (int node = 0; node < 8; node++){
        if(cur_index[node] > 0){
            printf("cur_index[%d] = %d\n", node, cur_index[node]);
            // createChildIfItDoesntExist(n, (unsigned)node^a);
            proc_subtree(new_tx0[node], new_ty0[node], new_tz0[node],
                         new_tx1[node], new_ty1[node], new_tz1[node], 
                         cur_index[node], depth + 1,
                         NULL, a, new_endpoints[node]);
        }
    }

    //this should only happen if call was recusive
    //n->logOdds = maxChildLogLikelihood(n);

}