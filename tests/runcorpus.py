import glob
import subprocess
import os


if __name__ == '__main__':
    files = sorted(glob.glob('tests/corpus/test*.c'))
    for f in files:
        print(f)
        name = os.path.basename(f)[:-2]
        subprocess.check_call((
            "gcc -fprofile-arcs -ftest-coverage -Wall -Werror --std=c99 -O2 %s"
            " src/discretesampler/sampler.c -lm -o%s") % (f, name),
            shell=True)

        results = subprocess.check_output("./%s" % (name,), shell=True)

    subprocess.check_call([
        "geninfo", ".", "-o", "./corpus.info" % (name,)
    ])


