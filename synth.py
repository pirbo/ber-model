import random
import bisect

default_values = {
    "base_state": "OK",
    "dg": None,
    "ape": None,
    "pol": None,
    "xrc": None,
    "dnm": None,
    "lig": None,
    "meth": "false"
}

def compute_cdf(items):
    s, fn, run_sum = sum(items.values()), [], 0
    for k, w in items.items():
        run_sum += w
        fn.append( (run_sum/s, k) ) 
    return ([x[0] for x in fn], [x[1] for x in fn])

def draw(cdf):
    x = random.random()
    idx = bisect.bisect(cdf[0], x)
    return cdf[1][idx]

def dual(x):
    return {"C": "G", "A": "T", "G": "C", "T": "A"}[x]

def node(base, base_type, up, prev, nxt):
    s = "DNA(" + ", ".join([ "{}~{}".format(k, v) if v else k for k, v in  default_values.items()])
    pn = ("e3" if up else "e5") + "~lig"
    nn = ("e5" if up else "e3") + "~lig"
    if prev > 0: pn = pn + "!{}".format(prev)
    if nxt: nn = nn + "!{}".format(nxt)
    s = s + ", {}, base~{}!{}, init~{}".format(pn, base_type, base, base_type) + \
        ", " + nn + ")"
    return s
    
def pair(base_type, prev, nxt = True):
    p1, p2 = prev
    base = max(p1, p2) + 1
    n1, n2 = base + 1, base + 2
    return (node(base, base_type, True, p1, n1 if nxt else None) + ", \\\n" + \
        node(base, dual(base_type), False, p2, n2 if nxt else False), n1, n2)

def strand(track):
    n1, n2, l = 0, 0, []
    for x in track[:-1]:
        s, n1, n2 = pair(x, (n1, n2))
        l.append(s)
    l.append(pair(track[-1], (n1, n2), nxt = False)[0])
    return "%init: 1 " + ", \\\n\t".join(l)

def gen_strand(nb, items):
    cdf = compute_cdf(items)
    return strand([draw(cdf) for i in range(nb)])
