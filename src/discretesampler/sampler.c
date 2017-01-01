#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#define mt_PARAM_N 312

struct mt {
    uint64_t mt[mt_PARAM_N]; /* the array for the state vector  */
    int16_t mti;
};


typedef struct {
  size_t capacity;
  size_t n_items;
  size_t n_viable_items;
  size_t item_mask;
  uint8_t n_bits_needed;
  double total;
  double max;
  double *weights;

  size_t table_size;
  size_t *alias_table;
} random_sampler;


static void mt_seed(struct mt *mt, uint64_t seed);
static uint64_t genrand64_int64(struct mt *r);


struct mt *mersenne_twister_new(uint64_t seed){
    struct mt *result = malloc(sizeof(struct mt));
    mt_seed(result, seed);
    return result;
}

void mersenne_twister_free(struct mt *mt){
    free(mt);
}


static uint8_t highest_set_bit(size_t i){
    // Note: Not very efficient, but not even close to being the hot path for
    // where we use it.
    uint8_t result = 0;
    while(i){
        result++;
        i >>= 1;
    }
    return result;
}

static random_sampler *random_sampler_create(size_t capacity){
  void *data = malloc(sizeof(random_sampler) + capacity * sizeof(size_t) * 2 + capacity * sizeof(double));
  random_sampler *result = (random_sampler*)data;
  result->capacity = capacity;
  result->n_items = 0;
  result->weights = (double*)(data + sizeof(random_sampler));
  result->alias_table =  (size_t*)((void*)result->weights + capacity * sizeof(double));
  return result;
}


static size_t right_fill(size_t mask){
  mask |= (mask >> 1);
  mask |= (mask >> 2);
  mask |= (mask >> 4);
  mask |= (mask >> 8);
  mask |= (mask >> 16);
  mask |= (mask >> 32);
  return mask;
}

static void init_mask(random_sampler *result){
  size_t mask = right_fill(result->table_size);
  result->item_mask = mask;
  result->n_bits_needed = highest_set_bit(mask);
}

static void init_full_table(random_sampler *result){
  assert(result->n_viable_items > 1);
  assert(result->n_viable_items == result->table_size);

  result->table_size = 0;

  for(size_t i = 0; i < result->n_items; i++){
    if(result->weights[i] / result->total <= 0) continue;
    double w = floor(result->n_viable_items * result->weights[i] / result->total) + 1;
    for(int j = 0; j < w; j++){
        result->alias_table[result->table_size++] = i;
        assert(result->table_size <= 2 * result->n_viable_items);
    }
  }
  init_mask(result);
}


static void random_sampler_initialize(random_sampler *result, size_t n_items, double *weights){
  assert(n_items <= result->capacity);

  memcpy(result->weights, weights, sizeof(double) * n_items);
   
  result->n_items = n_items;
  double min = INFINITY;
  double max = -INFINITY;
  double total = 0.0;
  result->table_size = 0;
  result->n_viable_items = 0;

  for(size_t i = 0; i < n_items; i++){
    double x = weights[i];
    if(x > 0){
      result->alias_table[result->table_size++] = i;
      result->n_viable_items++;
      if(x < min) min = x;
      if(x > max) max = x;
      total += x;
    }
  }

  assert(min < max);
  assert(max > 0);

  init_mask(result);
  result->total = total;
  init_full_table(result); 
}

random_sampler *random_sampler_new(size_t n_items, double *weights){
  random_sampler *result = random_sampler_create(n_items);
  random_sampler_initialize(result, n_items, weights);
  return result;
}

void random_sampler_free(random_sampler *sampler){
  free(sampler);
}

static double mt_random_double(struct mt *mt);

size_t random_sampler_sample(random_sampler *sampler, struct mt *mt){
  assert(sampler->table_size > 0);
  if(sampler->table_size == 1) return sampler->alias_table[0];
  size_t i = sampler->n_items;
  while(true){
    uint64_t probe = genrand64_int64(mt);
    for(uint8_t t = 0; t < 64; t += sampler->n_bits_needed){
        i = probe & sampler->item_mask;
        if(i < sampler->table_size){
          size_t result = sampler->alias_table[i];
          assert(result < sampler->n_items);
          if((i > 0) && (result == sampler->alias_table[i - 1])){ 
            return result;
          }
          double q = sampler->n_viable_items * sampler->weights[result] / sampler->total;
          double threshold = 1.0 - (q - floor(q));
          assert(!isnan(threshold));
          if((threshold < 1.0) && (mt_random_double(mt) > threshold)){
            return result;
          }
        }
        probe >>= sampler->n_bits_needed;
    }
  }
}


#define RECENCY 4


typedef struct {
    random_sampler *sampler;
    uint32_t hash;
    size_t access_date;
} sampler_entry;


typedef struct {
    size_t *buffer;
    size_t buffer_capacity;

    size_t capacity;
    size_t max_items;
    sampler_entry *entries;
    struct mt mersenne_twister;
    uint64_t generation;
    size_t last_index;
    size_t recent[RECENCY];
} sampler_family;


uint64_t hash64(uint64_t key){
    // Thomas Wang's 64 bit integer hash function
    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

uint64_t hash_double(uint64_t seed, double x){
    uint64_t key;
    memcpy(&key, &x, sizeof(double));
    return hash64(seed ^ key);
}

void MurmurHash3_x86_32 ( const void * key, int len,
                          uint32_t seed, void * out );

static uint32_t memhash(size_t n, char *bytes){
    uint32_t result;
    MurmurHash3_x86_32(bytes, n, 0, &result);
    return result;
}

static uint64_t hash_doubles(size_t n_items, double* weights){
    return memhash(n_items * sizeof(double), (char*)weights);
}

static void push_index(sampler_family *family, size_t index){
    assert(family->last_index < RECENCY);
    family->last_index = (family->last_index + 1) & (RECENCY - 1);
    family->recent[family->last_index] = index; 
    //printf("Writing %zu to %zu\n", index, family->last_index);
}

#define PROBE_MAX 8

static random_sampler *lookup_sampler(
    sampler_family *family, size_t n_items, double* weights
){
    if(n_items > family->max_items) family->max_items = n_items;

    for(size_t _i = 0; _i < RECENCY; _i++){
        size_t i = (_i + family->last_index) % RECENCY;

        sampler_entry *recent_entry = family->entries + family->recent[i];

        if(
            (recent_entry->sampler != NULL) &&
            (recent_entry->sampler->n_items == n_items) &&
            (memcmp(
                recent_entry->sampler->weights, weights, n_items * sizeof(double)) == 0)
        ){
            recent_entry->access_date = ++(family->generation);
            family->last_index = i;
            return recent_entry->sampler;
        }
    }

    uint64_t hash = family->capacity > PROBE_MAX ? hash_doubles(n_items, weights) : 0;

    size_t probe = (hash % (uint64_t)family->capacity);
    size_t target = probe;
    uint64_t target_date = family->entries[probe].access_date;
    for(size_t _i = 0; _i < PROBE_MAX; _i++){
        size_t i = (probe + _i) % family->capacity;
        sampler_entry *existing = family->entries + i;
        if(existing->sampler == NULL){
            target = i;
            break;
        }
        if((existing->hash == hash) && (existing->sampler->n_items == n_items)){
            assert(existing->sampler->capacity >= n_items);
            existing->access_date = ++(family->generation);
            if(memcmp(
                existing->sampler->weights, weights, n_items * sizeof(double)
            ) == 0){
                //printf("Cache hit after %d\n", (int)_i);
                push_index(family, i);
                return existing->sampler;
            }
        }
        if(existing->access_date < target_date){
            target = i;
            target_date = existing->access_date;
        }
    }
    // We didn't find an existing sampler, so it's time to create a new one.

    //printf("Cache miss\n");
    push_index(family, target);
    sampler_entry *result = family->entries + target;
   
    if((result->sampler == NULL) || (result->sampler->capacity < n_items)){
        random_sampler_free(result->sampler);
        result->sampler = random_sampler_create(family->max_items);
    }
    random_sampler_initialize(result->sampler, n_items, weights);
    result->hash = hash;
    result->access_date = ++(family->generation);
    return result->sampler;
}


size_t sampler_family_sample(sampler_family *samplers, size_t n_items, double *weights){
    if(n_items <= 1) return 0;

    if(samplers->buffer_capacity < n_items){
        free(samplers->buffer);
        samplers->buffer_capacity = right_fill(n_items) + 1;
        samplers->buffer = malloc(sizeof(size_t) * samplers->buffer_capacity);
    }

    double max = -INFINITY;
    size_t fill = 0;
    for(size_t i = 0; i < n_items; i++){
        double w = weights[i];
        if(w > 0) samplers->buffer[fill++] = i; 
        if(w > max) max = w;
    }

    struct mt *mt = &samplers->mersenne_twister;


    if(fill == 2){
        size_t i = samplers->buffer[0];
        size_t j = samplers->buffer[1];
        assert(i < j);
        if(mt_random_double(mt) <= weights[i] / (weights[i] + weights[j])){
            return i;
        } else { 
            return j;
        }
    };

    size_t mask = right_fill(fill);
    size_t n_bits_needed = highest_set_bit(mask);

    int rejections = 0;

    while(true){
      uint64_t probe = genrand64_int64(mt);
      for(uint8_t t = 0; t < 64; t += n_bits_needed){
        size_t i = probe & mask;
        if(max <= 0){
            if(i < n_items) return i;
        } else if(i < fill){
          size_t result = samplers->buffer[i];
          double w = weights[result];
          assert(w > 0);
          if(w < max){
              if(mt_random_double(mt) <= w / max){
                  return result;
              } else {
                  rejections++;
                  if(rejections > 1) goto fallback;
              }
          } else {
              return result;
          }
        }
        probe >>= n_bits_needed;
      }
    }

    fallback:;

    if(samplers->capacity == 0){
        random_sampler *sampler = random_sampler_new(n_items, weights);
        size_t result = random_sampler_sample(sampler, &samplers->mersenne_twister);
        random_sampler_free(sampler);
        return result;
    }
    random_sampler *sampler = lookup_sampler(samplers, n_items, weights);
    return random_sampler_sample(sampler, &(samplers->mersenne_twister));
}


sampler_family *sampler_family_new(size_t capacity, uint64_t seed){
    sampler_family *result = (sampler_family*)calloc(1, sizeof(sampler_family));
    result->last_index = 0;
    result->generation = 0;
    result->capacity = capacity;
    result->entries = calloc(capacity, sizeof(sampler_entry));
    mt_seed(&(result->mersenne_twister), seed);
    return result;
}

void sampler_family_free(sampler_family *family){
    for(size_t i = 0; i < family->capacity; i++){
        sampler_entry *entry = family->entries + i;
        random_sampler_free(entry->sampler);
    }
    free(family->entries);
    free(family);
}


/* 
   A C-program for MT19937-64 (2004/9/29 version).
   Coded by Takuji Nishimura and Makoto Matsumoto.

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_genrand64(seed)  
   or init_by_array64(init_key, key_length).

   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   References:
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
     ACM Transactions on Modeling and 
     Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
     ``Mersenne Twister: a 623-dimensionally equidistributed
       uniform pseudorandom number generator''
     ACM Transactions on Modeling and 
     Computer Simulation 8. (Jan. 1998) 3--30.

   Any feedback is very welcome.
   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
*/

/* The code has been modified to store internal state in heap/stack
 * allocated memory, rather than statically allocated memory, to allow
 * multiple instances running in the same address space. */

#define NN mt_PARAM_N
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */

static uint64_t genrand64_int64(struct mt *r);

/* initializes mt[NN] with a seed */
static void mt_seed(struct mt *mt, uint64_t seed)
{
    mt->mt[0] = seed;
    uint16_t mti = 0;
    for (mti=1; mti<NN; mti++) {
        mt->mt[mti] = (6364136223846793005ULL *
            (mt->mt[mti-1] ^ (mt->mt[mti-1] >> 62)) + mti);
    }
    mt->mti = mti;
}

/* Generate a random number on [0,1]-real-interval. */
static double mt_random_double(struct mt *mt)
{
    return (genrand64_int64(mt) >> 11) * (1.0/9007199254740991.0);
}

/* generates a random number on [0, 2^64-1]-interval */
static uint64_t genrand64_int64(struct mt *r)
{
    int i;
    uint64_t x;
    static uint64_t mag01[2]={0ULL, MATRIX_A};

    if (r->mti >= NN) { /* generate NN words at one time */

        /* if init has not been called, */
        /* a default initial seed is used */
        if (r->mti == NN+1)
            mt_seed(r, 5489ULL);

        for (i=0;i<NN-MM;i++) {
            x = (r->mt[i]&UM)|(r->mt[i+1]&LM);
            r->mt[i] = r->mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        for (;i<NN-1;i++) {
            x = (r->mt[i]&UM)|(r->mt[i+1]&LM);
            r->mt[i] = r->mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        x = (r->mt[NN-1]&UM)|(r->mt[0]&LM);
        r->mt[NN-1] = r->mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];

        r->mti = 0;
    }
  
    x = r->mt[r->mti++];

    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);

    return x;
}

// Following code is from MurmurHash3, which has been placed into the public
// domain.

uint32_t fmix32 ( uint32_t h )
{
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;

  return h;
}

uint32_t rotl32 ( uint32_t x, int8_t r )
{
  return (x << r) | (x >> (32 - r));
}

uint32_t getblock32 ( const uint32_t * p, int i )
{
  return p[i];
}


void MurmurHash3_x86_32 ( const void * key, int len,
                          uint32_t seed, void * out )
{
  const uint8_t * data = (const uint8_t*)key;
  const int nblocks = len / 4;

  uint32_t h1 = seed;

  const uint32_t c1 = 0xcc9e2d51;
  const uint32_t c2 = 0x1b873593;

  //----------
  // body

  const uint32_t * blocks = (const uint32_t *)(data + nblocks*4);

  for(int i = -nblocks; i; i++)
  {
    uint32_t k1 = getblock32(blocks,i);

    k1 *= c1;
    k1 = rotl32(k1,15);
    k1 *= c2;
    
    h1 ^= k1;
    h1 = rotl32(h1,13); 
    h1 = h1*5+0xe6546b64;
  }

  //----------
  // tail

  const uint8_t * tail = (const uint8_t*)(data + nblocks*4);

  uint32_t k1 = 0;

  switch(len & 3)
  {
  case 3: k1 ^= tail[2] << 16;
  case 2: k1 ^= tail[1] << 8;
  case 1: k1 ^= tail[0];
          k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1;
  };

  //----------
  // finalization

  h1 ^= len;

  h1 = fmix32(h1);

  *(uint32_t*)out = h1;
} 
