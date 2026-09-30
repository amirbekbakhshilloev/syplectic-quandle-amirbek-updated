import sys, json, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reproduce import axis0_sym
from symq.links import L1
out = {}
for p in [int(x) for x in sys.argv[1].split(',')]:
    t = time.time(); a0 = axis0_sym(L1, p)
    out[p] = dict(axis0=a0, p4=p**4, extra=a0 - p**4, H=2 * p**4 - p**2, sec=round(time.time() - t, 1))
    print(p, out[p], flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'primes_extra_L1.json'), 'w'), indent=1)
