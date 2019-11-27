#include <stdio.h>

void main(){
    unsigned int one = (1<<8);
    unsigned char two = 0;

    two |= one;

    printf("%u\n", one);
    printf("%u\n", two);
}