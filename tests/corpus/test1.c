
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

void *sampler_family_new(size_t capacity, uint64_t seed);
void sampler_family_free(void *samplers);
size_t sampler_family_sample(void *samplers, size_t n_items, double *weights);

void do_sample(void *samplers, size_t n_items, double *weights){
    size_t result = sampler_family_sample(samplers, n_items, weights);
    if(weights[result] <= 0){
        for(size_t i = 0; i < n_items; i++) assert(weights[i] <= 0);
    }
    printf("%zu\n", result);
}

int main(int argc, char **argv){
    void *family = sampler_family_new(1, 0ul);
    double weights0[2] = {0.0, 1.0};
    for(int i = 0; i < 3; i++){
        do_sample(family, 2, weights0);    }
    sampler_family_free(family);
    return 0;
}
