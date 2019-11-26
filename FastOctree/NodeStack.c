//
// Created by mschnur on 11/16/19.
//

#include "NodeStack.h"

#include <stdio.h>
#include <stdlib.h>

void init(NodeStack* stack)
{
    stack->topIndex = -1;
}

int size(const NodeStack* stack)
{
    return stack->topIndex + 1;
}

int isEmpty(const NodeStack* stack)
{
    return stack->topIndex == -1;
}

int isFull(const NodeStack* stack)
{
    return stack->topIndex == NODE_STACK_SIZE - 1;
}

int push(NodeStack* stack, Node* node)
{
    if (isFull(stack))
    {
        return 1;
    }

    // Pre-increment the index, then add the new node pointer to the stack
    stack->items[++(stack->topIndex)] = node;
}

int peek(const NodeStack* stack, Node** poppedElement)
{
    if (isEmpty(stack))
    {
        return 1;
    }

    *poppedElement = stack->items[stack->topIndex];
}

int pop(NodeStack* stack, Node** poppedElement)
{
    if (isEmpty(stack))
    {
        return 1;
    }

    // Get the element at topIndex, then decrement it.
    *poppedElement = stack->items[stack->topIndex--];
}