# build_templates.py -- builds resources/basepair-idealized (RUBRIC's class templates).
#
# Saved 2026-11-11 exactly as run on 2026-11-03 to produce the installed template set
# (recovered from the session log: original script plus its one recorded patch that adds the
# majority hydrogen bonds and writes symmetric classes in both residue orders). Logic unchanged.
#
# How it was run:  cd PARSE-RNA-Analysis && python <this file>
# Inputs it expects (paths are hard-coded below and must be adjusted to rerun):
#   S/bgsu/<class>.json         BGSU/NAKB base-pair catalog listings,
#                               https://nakb.org/ndbmodule/bp-catalog/static/data/<class>.json
#   S/impact/res_old/           a copy of the engine resources/ directory (output is written to
#                               S/impact/res_full/basepair-idealized, then copied into resources/)
#   PARSE-RNA-Analysis data/     annotations of the six other programs, curated reference set,
#                               structures
#
# Reproducibility note (verified 2026-11-11): rerunning reproduces 146 of the 178 template files
# byte for byte. The other 32 are class-sequence combinations with more than 300 pairs, where a
# random sample of 300 is drawn from a list whose order depends on Python's set iteration order
# (hash seed), so the medoid can change. 31 of them move by <= 0.25 A RMSD; cWW/UU moves by 1.3 A.
# Setting PYTHONHASHSEED or sorting the pair keys before sampling would make the build deterministic.

"""Build the complete consensus template set (5 tiers) into a fresh resource copy."""
import sys, json, collections, shutil, urllib.request
from pathlib import Path
import numpy as np, gemmi
S = Path('/private/tmp/claude-531458859/-Users-kdewan2-Desktop-Projects-find-pair-score-rna/a6f8054d-8094-45be-b248-f3b5238b0f2e/scratchpad')
ROOT = Path('/Users/kdewan2/Desktop/Projects/find-pair-score-rna/PARSE-RNA-Analysis')
sys.path.insert(0, str(ROOT / 'scripts/validation'))
from figure2_comparison import Comparison
OTHERS = ['dssr', 'fr3d', 'rnaview', 'barnaba', 'rnapolis', 'bpnet']
CLASSES = ["cWW","tWW","cWH","tWH","cWS","tWS","cHH","tHH","cHS","tHS","cSS","tSS"]
ATOMS = {"A": ["C1'","N9","C8","N7","C5","C6","N6","N1","C2","N3","C4"],
         "G": ["C1'","N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"],
         "C": ["C1'","N1","C2","O2","N3","C4","N4","C5","C6"],
         "U": ["C1'","N1","C2","O2","N3","C4","O4","C5","C6"]}
canon = lambda lw: lw[0]+lw[2]+lw[1] if {'W':0,'H':1,'S':2}[lw[1]] > {'W':0,'H':1,'S':2}[lw[2]] else lw
sym = lambda lw: lw[1] == lw[2]
def combo(lw, seq): return (lw, min(seq, seq[::-1])) if sym(lw) else (lw, seq)
UNIVERSE = sorted({combo(lw, a+b) for lw in CLASSES for a in 'ACGU' for b in 'ACGU'})

comp = Comparison(ROOT / 'data'); tools = {t: comp.tool(t) for t in OTHERS}; names = comp.residue_names
curated = set(json.load(open(ROOT / 'data/reference/reference_sets/reference_set.json'))['pdbs'])
pdbs = sorted(set.intersection(*(set(v) for v in tools.values())) & curated)
pairs = collections.defaultdict(list)   # combo -> [(votes, pdb, resA, resB)] residues ordered to match combo seq
for pdb in pdbs:
    for key in set().union(*(tools[t][pdb].keys() for t in OTHERS)):
        votes = collections.Counter()
        for t in OTHERS:
            s = tools[t][pdb].get(key)
            if s and len(s) == 1:
                lw = next(iter(s))
                if len(lw) == 3 and lw[0] in 'ct' and lw[1] in 'WHS' and lw[2] in 'WHS': votes[lw] += 1
        if not votes: continue
        (lw, v), *rest = votes.most_common()
        if v < 4 or (rest and rest[0][1] == v): continue
        b = [names.get(pdb, {}).get(f"{c}|{n}|{i}") for c, n, i in key]
        if None in b or not set(b) <= set('ACGU'): continue
        a_, b_ = (key, b) if lw == canon(lw) else (key[::-1], b[::-1])
        cb = combo(canon(lw), ''.join(b_))
        if cb[1] != ''.join(b_): a_, b_ = a_[::-1], b_[::-1]
        pairs[cb].append((v, pdb, a_[0], a_[1]))

def superpose(p, q):
    pc, qc = p.mean(0), q.mean(0)
    u, s, vt = np.linalg.svd((q - qc).T @ (p - pc))
    rot = u @ np.diag([1, 1, np.sign(np.linalg.det(u @ vt))]) @ vt
    moved = (q - qc) @ rot + pc
    return moved, float(np.sqrt(((moved - p) ** 2).sum(1).mean()))
def swap(x, seq):  # residue swap for same-base pairs in symmetric classes
    n = len(ATOMS[seq[0]]); return np.vstack([x[n:], x[:n]])
def best_fit(ref, x, lw, seq):
    out = superpose(ref, x)
    if sym(lw) and seq[0] == seq[1]:
        alt = superpose(ref, swap(x, seq))
        if alt[1] < out[1]: out = alt
    return out
def medoid(shapes, lw, seq):
    mean = shapes[0]
    for _ in range(3): mean = np.mean([best_fit(mean, s, lw, seq)[0] for s in shapes], axis=0)
    r = [best_fit(mean, s, lw, seq)[1] for s in shapes]
    i = int(np.argmin(r)); return shapes[i], float(np.median(r))

rng = np.random.default_rng(0); plan = {}
for cb in UNIVERSE:
    pool = pairs.get(cb, [])
    for tier, need, label in ((1, 6, 'six of six'), (2, 5, 'five of six'), (3, 4, 'four of six')):
        sel = [p for p in pool if p[0] >= need]
        if len(sel) >= 5: plan[cb] = (tier, label, sel); break
CATALOG = {('cHH','GG'), ('tHH','GG'), ('cHS','GA'), ('cHS','UA'), ('cHS','UU'), ('cSS','CC'), ('cSS','CU'),
           ('cSS','UU'), ('tSS','AA'), ('tSS','AC'), ('tSS','AU'), ('tSS','GU')}
for cb in UNIVERSE:
    if cb in plan: continue
    if cb in CATALOG: plan[cb] = (4, 'catalog exemplar', None)
    elif pairs.get(cb): plan[cb] = (5, f'four or more of six, {len(pairs[cb])} pair(s)', pairs[cb])
need_coords = collections.defaultdict(list)
for cb, (tier, label, sel) in plan.items():
    if sel is None: continue
    if len(sel) > 300: sel = [sel[i] for i in rng.choice(len(sel), 300, replace=False)]
    plan[cb] = (tier, label, sel)
    for v, pdb, a, b in sel: need_coords[pdb].append((cb, a, b))
shapes = collections.defaultdict(list)
def residue(model, rid):
    ch = model.find_chain(rid[0])
    return next((r for r in ch or [] if r.seqid.num == rid[1] and r.seqid.icode.strip() == rid[2]), None)
def coords(res, base):
    pts = []
    for nm in ATOMS[base]:
        at = res.find_atom(nm, '*')
        if at is None: return None
        pts.append(at.pos.tolist())
    return np.array(pts)
for pdb, items in need_coords.items():
    model = gemmi.read_structure(str(ROOT / f'data/structures/{pdb}.cif'))[0]
    for cb, a, b in items:
        ra, rb = residue(model, a), residue(model, b)
        if ra is None or rb is None: continue
        xa, xb = coords(ra, cb[1][0]), coords(rb, cb[1][1])
        if xa is not None and xb is not None: shapes[cb].append(np.vstack([xa, xb]))
templates = {}
for cb, (tier, label, sel) in plan.items():
    if sel is None: continue
    if not shapes[cb]: continue
    xyz, spread = medoid(shapes[cb], *cb) if len(shapes[cb]) > 1 else (shapes[cb][0], 0.0)
    templates[cb] = (xyz, tier, label, len(shapes[cb]), spread)

DONORS = {"A": {"N6"}, "G": {"N1", "N2"}, "C": {"N4"}, "U": {"N3"}}
ACCEPTORS = {"A": {"N1", "N3", "N7"}, "G": {"O6", "N3", "N7"}, "C": {"O2", "N3"}, "U": {"O2", "O4"}}
def contacts(x, seq, cut=3.5):
    na = len(ATOMS[seq[0]]); res = (dict(zip(ATOMS[seq[0]], x[:na])), dict(zip(ATOMS[seq[1]], x[na:])))
    out = set()
    for i, j in ((0, 1), (1, 0)):
        for d in DONORS[seq[i]]:
            for a in ACCEPTORS[seq[j]]:
                if d in res[i] and a in res[j] and np.linalg.norm(res[i][d] - res[j][a]) <= cut:
                    out.add((i, d, a, float(np.linalg.norm(res[i][d] - res[j][a]))))
    return out
def majority_hbonds(shape_list, seq):
    cnt = collections.Counter(); dist = collections.defaultdict(list)
    for x in shape_list:
        for i, d, a, r in contacts(x, seq):
            cnt[(i, d, a)] += 1; dist[(i, d, a)].append(r)
    return [(k, float(np.median(dist[k]))) for k, c in cnt.items() if c >= 0.5 * len(shape_list)]
HBONDS = {cb: majority_hbonds(shapes[cb], cb[1]) for cb in templates}

# catalog exemplars
cat_dir = S / 'bgsu'
def catalog_xyz(lw, seq):
    items = json.load(open(cat_dir / f'{lw}.json'))['items']
    for k, v in items.items():
        if not v.get('coordinates_exist') or v['pdb'] == 'Modeled': continue
        if not ((k.upper() in (seq, seq[::-1])) if sym(lw) else k.upper() == seq): continue
        fam = v['family']; case = '0'
        if fam[1].lower() == fam[2].lower() and fam[1] != fam[2]: case = '1' if fam[1] > fam[2] else '2'
        url = f"https://nakb.org/ndbmodule/bp-catalog/static/data/{lw}/{k.upper()}-{case}.pdb"
        f = cat_dir / f"{lw}_{k.upper()}-{case}.pdb"
        if not f.exists(): urllib.request.urlretrieve(url, f)
        st = gemmi.read_structure(str(f)); res = [r for ch in st[0] for r in ch]
        if len(res) != 2: continue
        bases = ''.join(r.name.strip()[-1] for r in res)
        order = res if bases == k.upper() else res[::-1]
        if ''.join(r.name.strip()[-1] for r in order) != k.upper(): continue
        xa, xb = coords(order[0], k.upper()[0]), coords(order[1], k.upper()[1])
        if xa is None or xb is None: continue
        x = np.vstack([xa, xb])
        if k.upper() != seq: x = swap(x, k.upper())     # symmetric class, reversed listing
        return x, f"catalog exemplar {v['pdb']} {','.join(v['units'])}", url
    return None
# orientation check: catalog vs consensus on combos both cover (asymmetric classes)
check = []
for cb in [c for c in templates if not sym(c[0])][:200]:
    got = catalog_xyz(*cb) if (cb[0], cb[1]) else None
    if got is None: continue
    x = got[0]; t = templates[cb][0]
    as_is = superpose(t, x)[1]
    na = len(ATOMS[cb[1][0]]); flipped_seq = cb[1][::-1]
    check.append((cb, round(as_is, 2)))
    if len(check) >= 12: break
print('catalog vs consensus (asymmetric classes, as listed):', check)
for cb in CATALOG:
    if cb in templates: continue
    got = catalog_xyz(*cb)
    if got: templates[cb] = (got[0], 4, got[1], 1, 0.0)

out = S / 'impact' / 'res_full'
if out.exists(): shutil.rmtree(out)
shutil.copytree(S / 'impact' / 'res_old', out)
ideal = out / 'basepair-idealized'
for d in ideal.iterdir():
    if d.is_dir(): shutil.rmtree(d)
(ideal / 'hbond_patterns_cache.json').unlink(missing_ok=True)
summary = []
def write(lw, seq, xyz, tier, label, n, spread, hb):
    names_ = [(seq[0], 1, a) for a in ATOMS[seq[0]]] + [(seq[1], 2, a) for a in ATOMS[seq[1]]]
    lines = [f"REMARK   1 Template: {seq[0]}-{seq[1]}-{lw}", f"REMARK   2 Tier {tier}: {label}",
             f"REMARK   3 Pairs: {n}; median RMSD to mean shape {spread:.2f} A"]
    for (i, d, a), r in sorted(hb):   # donor on residue i+1, acceptor on the other residue
        lines.append(f"REMARK   7   {d}-{a}: {r:.2f} A (present in most consensus pairs)")
    for i, ((base, rn, atom), (x, y, z)) in enumerate(zip(names_, xyz), 1):
        nm = atom if len(atom) == 4 else ' ' + atom
        lines.append(f"ATOM  {i:5d} {nm:<4s}   {base} {'A' if rn == 1 else 'B'}{rn:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {atom[0]}")
    (ideal / lw).mkdir(exist_ok=True)
    (ideal / lw / f'{seq}.pdb').write_text('\n'.join(lines + ['END']) + '\n')
for (lw, seq), (xyz, tier, label, n, spread) in sorted(templates.items()):
    hb = HBONDS.get((lw, seq)) or [(k, r) for (*k, r) in [(i, d, a, r) for i, d, a, r in contacts(xyz, seq)]]
    write(lw, seq, xyz, tier, label, n, spread, hb)
    summary.append({'class': lw, 'sequence': seq, 'tier': tier, 'source': label, 'pairs': n, 'spread': round(spread, 2),
                    'hbonds': [f"{d}-{a}" for (i, d, a), r in sorted(hb)]})
    if sym(lw) and seq[0] != seq[1]:
        write(lw, seq[::-1], swap(xyz, seq), tier, label, n, spread, [((1 - i, d, a), r) for (i, d, a), r in hb])
json.dump(summary, open(ideal / 'summary.json', 'w'), indent=1)
tc = collections.Counter(s['tier'] for s in summary)
print(f"written {len(summary)} of {len(UNIVERSE)} combinations; by tier {dict(sorted(tc.items()))}; missing {sorted(set(UNIVERSE) - set(templates))}")
