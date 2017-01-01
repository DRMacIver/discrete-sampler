#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

void *sampler_family_new(size_t capacity, uint64_t seed);
void sampler_family_free(void *samplers);
size_t sampler_family_sample(void *samplers, size_t n_items, double *weights);



#define SIZE 12


int main(int argc, char **argv){
    void *family = sampler_family_new(1, 2);
    size_t accumulator = 0;

    double weights[SIZE];

     for(size_t c = 32767; c < (1L << (2 * SIZE)); c++){
        size_t d = c;
        for(int i = 0; i < SIZE; i++){
            weights[i] = d & 1;
            d >>= 1; 
            weights[i] += d & 1;
            d >>= 1; 
        }
        fflush(stdout);
        size_t i = sampler_family_sample(family, SIZE, weights);
        assert(i < SIZE);
        assert(weights[i] > 0.0);
        accumulator += i;
     }
    printf("Accumulator %d\n", (int)accumulator);
    sampler_family_free(family);
}
