//
// Created by mschnur on 11/16/19.
//

#ifndef FASTOCTREE_NODESTACK_H
#define FASTOCTREE_NODESTACK_H

#include "FastOctree.h"

#define NODE_STACK_SIZE 10000

typedef struct NodeStack
{
    int topIndex;
    Node* items[NODE_STACK_SIZE];
} NodeStack;

void init(NodeStack* stack);

int size(const NodeStack* stack);

int isEmpty(const NodeStack* stack);

int isFull(const NodeStack* stack);

int push(NodeStack* stack, Node* node);

int peek(const NodeStack* stack, Node** poppedElement);

int pop(NodeStack* stack, Node** poppedElement);

#endif //FASTOCTREE_NODESTACK_H
