"""34.85a: kontrola přenesených bloků nul rodiny S7 proti manifestu z Macu (SHA-256 každého bloku, počty L a nul, otisk skriptu).
usage: python3 -I s7_sha.py SLOZKA_BEHU(manifest.json + S7/) VYSTUP.json"""
import sys, os, json, hashlib, gzip
B, OUT = sys.argv[1].rstrip('/'), sys.argv[2]; m = json.load(open(B + '/manifest.json')); r = m['rodiny']['S7']; bl = r['bloky']
ok = bad = nL = nn = 0; sek = 0.0; spatne = []
for k in sorted(bl, key=int):
    b = bl[k]; p = os.path.join(B, 'S7', b['soubor'])
    if not os.path.exists(p): bad += 1; spatne.append([b['soubor'], 'chybí']); continue
    h = hashlib.sha256(open(p, 'rb').read()).hexdigest()
    if h == b['sha256']: ok += 1
    else: bad += 1; spatne.append([b['soubor'], 'sha'])
    nL += b['nL']; nn += b['nnul']; sek += b.get('sekund', 0)
navic = sorted(set(f for f in os.listdir(B + '/S7') if f.endswith('.jsonl.gz')) - set(b['soubor'] for b in bl.values()))
out = dict(format=m['format'], xmax=m['xmax'], divz=m['divz'], prostredi=m['prostredi'], skript_sha256=m['prostredi']['skript_sha256'], skript_sha_ocekavany='1b53fe0e093e4c998e2eda9c63594e297ab08666da4cfec8644bc29497a103e8',
           bloku_celkem=r['bloku_celkem'], bloku_v_manifestu=len(bl), sha_ok=ok, sha_spatne=bad, spatne=spatne[:20], souboru_navic=navic[:20], nL=nL, nnul=nn, procesorovych_hodin=sek / 3600, spec=r['spec'])
out['v_poradku'] = (ok == r['bloku_celkem'] == len(bl) and bad == 0 and not navic and out['skript_sha256'] == out['skript_sha_ocekavany'])
json.dump(out, open(OUT, 'w'), indent=1, ensure_ascii=False); print({k: v for k, v in out.items() if k not in ('spatne', 'souboru_navic', 'spec', 'prostredi')})
