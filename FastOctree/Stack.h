//
// Created by mschnur on 11/16/19.
//

#ifndef FASTOCTREE_STACK_H
#define FASTOCTREE_STACK_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

typedef struct Stack {
    int topIndex;
    size_t elementSize;
    size_t capacity;
    void *data;
} Stack;

Stack *Stack__init(size_t elementSize, size_t numElements);

void Stack__destroy(Stack *stack);

int Stack__size(const Stack *stack);

int Stack__isEmpty(const Stack *stack);

int Stack__isFull(const Stack *stack);

int Stack__push(Stack *stack, void *element);

int Stack__peek(const Stack *stack, void *peekedElement);

int Stack__pop(Stack *stack, void *poppedElement);

#ifdef __cplusplus
}
#endif

#endif //FASTOCTREE_STACK_H
