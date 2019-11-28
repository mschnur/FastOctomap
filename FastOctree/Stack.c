//
// Created by mschnur on 11/16/19.
//

#include "Stack.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

Stack* Stack__init(int elementSize, int numElements)
{
    Stack *s = calloc(1, sizeof(Stack));
    s->topIndex = -1;
    s->elementSize = elementSize;
    s->capacity = numElements;
    s->data = calloc(numElements, elementSize);
    return s;
}

void Stack__destroy(Stack* stack)
{
    free(stack->data);
    free(stack);
}

int Stack__size(const Stack* stack)
{
    return stack->topIndex + 1;
}

int Stack__isEmpty(const Stack* stack)
{
    return stack->topIndex == -1;
}

int Stack__isFull(const Stack* stack)
{
    return stack->topIndex == (stack->capacity - 1);
}

int Stack__push(Stack* stack, void* element)
{
    if (Stack__isFull(stack))
    {
        return -1;
    }

    stack->topIndex++;
    void* target = ((uint8_t*) stack->data) + (stack->topIndex * stack->elementSize);
    memcpy(target, element, stack->elementSize);
    return 0;
}

int Stack__peek(const Stack* stack, void* peekedElement)
{
    if (Stack__isEmpty(stack))
    {
        return -1;
    }

    void* source = ((uint8_t*) stack->data) + (stack->topIndex * stack->elementSize);
    memcpy(peekedElement, source, stack->elementSize);
    return 0;
}

int Stack__peekAndPop(Stack* stack, void* poppedElement)
{
    if (Stack__isEmpty(stack))
    {
        return -1;
    }

    void* source = ((uint8_t*) stack->data) + (stack->topIndex * stack->elementSize);
    memcpy(poppedElement, source, stack->elementSize);
    stack->topIndex--;
    return 0;
}

int Stack__pop(Stack* stack)
{
    if (Stack__isEmpty(stack))
    {
        return -1;
    }

    stack->topIndex--;
    return 0;
}