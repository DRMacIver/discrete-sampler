#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

void *sampler_family_new(size_t capacity, uint64_t seed);
void sampler_family_free(void *samplers);
size_t sampler_family_sample(void *samplers, size_t n_items, double *weights);


size_t accumulator;

void sample_test(size_t n, size_t k, void *samplers){
    for(size_t i = 0; i < n; i++){
        double *weights = malloc(k * sizeof(double));
        for(size_t j = 0; j < k; j++) weights[j] = 1.0;
        accumulator += sampler_family_sample(samplers, k, weights);        
        free(weights);
    }
}


int main(int argc, char **argv){
    void *family = sampler_family_new(2048, 2);
    accumulator = 0;

     for(int c = 0; c < 100000; c++){
        sample_test(14, 256, family);
        sample_test(1, 20, family);
        sample_test(7, 1, family);
        sample_test(1, 2, family);
        sample_test(7, 1, family);
        sample_test(1, 8, family);
     }
    printf("Accumulator %d\n", (int)accumulator);
    sampler_family_free(family);
}
