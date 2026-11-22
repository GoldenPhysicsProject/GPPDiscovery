"""Task 3, step 0: inventory of Platt's tables.  Per file: name, size, blocks, t-range, zero count, first and last
zero, contiguity with the next file; zero #1 against mpmath.zetazero(1); one mid-file zero per file against
mpmath.zetazero at 25 digits.  Usage: python inventory.py ZETA_FOLDER OUT_DIR"""
import sys, csv, math
from pathlib import Path
import mpmath
from mpmath import mp, mpf
import plattio

zeta = Path(sys.argv[1]); out = Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
data = zeta / 'data'
files = sorted(data.glob('zeros_*.dat'), key=lambda q: int(q.stem.split('_')[1]))
mp.prec = 400
EPS = mpf(2) ** -101


def zero_from(t0, ints):
    return mpf(t0) + mpf(ints) * EPS


def digits(x, n=30):
    mp.dps = 40
    s = mpmath.nstr(x, n)
    mp.prec = 400
    return s


rows = []
prev_N1 = prev_t1 = None
total = 0
for p in files:
    size = p.stat().st_size
    blks = list(plattio.blocks(p))
    count = sum(b[5] - b[4] for b in blks)
    expect = 8 + sum(32 + 13 * (b[5] - b[4]) for b in blks)
    first_b, last_b = blks[0], blks[-1]
    arr0 = plattio.read_block(p, first_b[1], 1)
    first = zero_from(first_b[2], plattio.entry_int(arr0, 0))
    arrL = plattio.read_block(p, last_b[1], last_b[5] - last_b[4])
    last = zero_from(last_b[2], plattio.delta_sum_int(arrL))
    # mid-file zero: index N (1-based), found by block
    N_target = first_b[4] + count // 2 + 1
    for b in blks:
        if b[4] < N_target <= b[5]:
            arr = plattio.read_block(p, b[1], b[5] - b[4])
            i = N_target - b[4] - 1
            mid = zero_from(b[2], plattio.delta_sum_int(arr[:i + 1]))
            break
    mp.dps = 40
    ref = mpmath.zetazero(N_target).imag
    mp.prec = 400
    diff = abs(mid - ref)
    agree = float(-mpmath.log10(diff / mid)) if diff > 0 else 99.0
    rows.append(dict(file=p.name, size_bytes=size, size_expected=expect, blocks=len(blks), t_lo=first_b[2],
                     t_hi=last_b[3], N_lo=first_b[4], N_hi=last_b[5], zeros=count, first_zero=digits(first),
                     last_zero=digits(last), contiguous_with_prev=(prev_N1 is None or prev_N1 == first_b[4]),
                     t_contiguous=(prev_t1 is None or prev_t1 == first_b[2]), check_N=N_target,
                     check_zero=digits(mid), check_digits_agree_with_mpmath=round(agree, 1)))
    prev_N1, prev_t1 = last_b[5], last_b[3]
    total += count
    print(p.name, count, flush=True)

mp.dps = 40
z1 = rows[0]['first_zero']
ref1 = mpmath.nstr(mpmath.zetazero(1).imag, 30)
ok1 = z1[:27] == ref1[:27]
with open(out / 'inventory.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
with open(out / 'inventory.md', 'w') as f:
    f.write("# Platt zero tables: inventory (NUMERICAL data audit)\n\n")
    f.write(f"Files: {len(rows)}. Total zeros: **{total:,}** (N from {rows[0]['N_lo']} to {rows[-1]['N_hi']}). "
            f"Height range: {rows[0]['t_lo']} to {rows[-1]['t_hi']}.\n\n")
    f.write(f"Zero #1 from the table: `{z1}`; mpmath.zetazero(1): `{ref1}`; first 27 digits agree: **{ok1}**.\n\n")
    f.write("| file | size (bytes) | blocks | t range | zeros | first zero | last zero | contiguous | mid-file zero N | digits agreeing with mpmath.zetazero |\n|---|---|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        f.write(f"| {r['file']} | {r['size_bytes']:,}{'' if r['size_bytes']==r['size_expected'] else ' (expected '+str(r['size_expected'])+')'} | {r['blocks']} | "
                f"[{r['t_lo']:g}, {r['t_hi']:g}) | {r['zeros']:,} | {r['first_zero'][:22]} | {r['last_zero'][:22]} | "
                f"{r['contiguous_with_prev'] and r['t_contiguous']} | {r['check_N']:,} | {r['check_digits_agree_with_mpmath']} |\n")
bad = [r['file'] for r in rows if not (r['contiguous_with_prev'] and r['t_contiguous']) or r['check_digits_agree_with_mpmath'] < 25]
print("zero #1 ok:", ok1, "problems:", bad)
sys.exit(0 if ok1 and not bad else 2)
