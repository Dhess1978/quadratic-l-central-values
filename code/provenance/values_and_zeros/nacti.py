"""Read study V blocks, verify SHA-256 against manifest, save per-family arrays (D, L) as npz."""
import json, gzip, hashlib, os, sys, numpy as np
B = sys.argv[1]; OUT = sys.argv[2]
man = json.load(open(os.path.join(B, 'manifest.json')))
print('format', man['format'], 'xmax', man.get('xmax'), 'divz', man.get('divz'))
print('prostredi', man.get('prostredi'), man.get('prostredi_posledni'))
tot = 0; rep = {}
for F, r in man['rodiny'].items():
    Ds = []; Ls = []; ok = bad = 0; nL = 0
    for k in sorted(r['bloky'], key=lambda s: int(s)):
        info = r['bloky'][k]; raw = open(os.path.join(B, F, info['soubor']), 'rb').read()
        if hashlib.sha256(raw).hexdigest() != info['sha256']: bad += 1; continue
        ok += 1; cnt = 0
        for line in gzip.decompress(raw).decode().splitlines():
            if line.strip():
                o = json.loads(line); Ds.append(o['D']); Ls.append(o['L']); cnt += 1
        assert cnt == info['nL'], (F, k, cnt, info['nL'])
        nL += cnt
    D = np.array(Ds, dtype=np.int64); L = np.array(Ls)
    assert len(np.unique(D)) == len(D), F + ' duplicate D'
    np.savez_compressed(os.path.join(OUT, F + '.npz'), D=D, L=L)
    rep[F] = dict(bloky=r['bloku_celkem'], sha_ok=ok, sha_bad=bad, nL=nL, spec=r['spec'])
    tot += nL
    print(F, 'blocks', ok, '/', r['bloku_celkem'], 'bad', bad, 'L', nL, 'D range', D.min(), D.max(), 'pos', int((D > 0).sum()), 'min L', L.min(), 'neg', int((L <= 0).sum()), 'nan', int(np.isnan(L).sum()))
json.dump(rep, open(os.path.join(OUT, 'kontrola_sha.json'), 'w'), indent=1, ensure_ascii=False)
print('total L', tot)
