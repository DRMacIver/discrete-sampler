from hypothesis import given, strategies as st, settings

import subprocess
import os


PREAMBLE = """
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
    printf("%%zu\\n", result);
}

int main(int argc, char **argv){
    void *family = sampler_family_new(%d, %dul);
"""

FINISH = """
    sampler_family_free(family);
    return 0;
}
"""[1:]


simple_floats = st.integers(0, 1000).map(float)


def write_file(fn, capacity, seed, weights, draws):
    with open(fn, 'w') as o:
        o.write(PREAMBLE % (capacity, seed))
        used = {i for _, i in draws}

        for i, w in enumerate(weights):
            if i not in used:
                continue
            name = "weights%d" % (i,)
            o.write("    double %s[%d] = {%s};" % (
                name, len(w), ', '.join(map(str, w))
            ))
            o.write('\n')

        for n, i in draws:
            name = "weights%d" % (i,)
            w = weights[i]
            if n == 1:
                o.write("    do_sample(family, %d, %s);" % (
                    len(w), name))
            else:
                o.write("    for(int i = 0; i < %d; i++){\n" % (n,))
                o.write("        do_sample(family, %d, %s);" % (
                    len(w), name))
                o.write("    }")
            o.write('\n')
        o.write(FINISH)


def run_file(fn):
    subprocess.check_call((
        "gcc -Wall -Werror --std=c99 -O2 %s "
        "src/discretesampler/sampler.c -lm -otest") % (fn,),
        shell=True, env=os.environ)

    results = subprocess.check_output("./test", shell=True, env=os.environ)

    return [
        int(l.strip().decode('ascii')) for l in results.split(b'\n')
        if l.strip()
    ]


@settings(timeout=-1, max_examples=10**6, max_iterations=10**6)
@given(st.data())
def test_always_runs_validly(data):
    capacity = data.draw(st.integers(1, 128))
    seed = data.draw(st.integers(0, 2 ** 64 - 1))
    weights = data.draw(
        st.lists(
            st.lists(simple_floats, min_size=1),
            min_size=1, average_size=20))
    draws = data.draw(st.lists(
        st.tuples(st.integers(1, 1000), st.integers(0, len(weights) - 1)),
        average_size=20
    ))
    write_file('test.c', capacity, seed, weights, draws)
    write_file('test2.c', 0, seed, weights, draws)

    try:
        results = run_file('test.c')
        assert len(results) == sum(n for n, _ in draws)
        assert results == run_file('test2.c')

        counter = 0
        for n, i in draws:
            w = weights[i]
            trivial = not any(w)
            for _ in range(n):
                if not trivial:
                    assert w[results[counter]] != 0
                counter += 1
    except Exception:
        with open("test.c") as i:
            data = i.read()
        with open("failing-test.c", "w") as o:
            o.write(data)
        raise

if __name__ == '__main__':
    try:
        test_always_runs_validly()
    except Exception:
        i = 1
        while True:
            fname = "tests/corpus/test%d.c" % (i,)
            if not os.path.exists(fname):
                os.rename("test.c", fname)
                print("Wrote test file to %s" % (fname,))
                raise
            i += 1
