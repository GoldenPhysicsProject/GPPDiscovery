"""Build index.db for platt_zeros.py (the file that ships without it).
Schema is read off the SQL in platt_zeros.py: table zero_index with columns (t, N, filename, offset, block_number)
in that order (`select *` is unpacked as t0, N0, filename, offset, block_number; queries filter on t and N).
One row per block of every data file: t = block t0, N = Nt0 (zeros below t0), offset = byte offset of the block header,
block_number = index of the block in its file.
Usage: python build_index.py ZETA_FOLDER   (expects ZETA_FOLDER/data/zeros_*.dat; writes ZETA_FOLDER/index.db)"""
import sqlite3
import sys
from pathlib import Path
import plattio

zeta = Path(sys.argv[1])
data = zeta / 'data'
db_path = zeta / 'index.db'
if db_path.exists():
    db_path.unlink()
db = sqlite3.connect(db_path)
db.execute("CREATE TABLE zero_index (t REAL, N INTEGER, filename TEXT, offset INTEGER, block_number INTEGER)")
rows = 0
for p in sorted(data.glob('zeros_*.dat'), key=lambda q: int(q.stem.split('_')[1])):
    batch = []
    for b, off, t0, t1, N0, N1 in plattio.blocks(p):
        batch.append((t0, N0, p.name, off, b))
    db.executemany("INSERT INTO zero_index VALUES (?,?,?,?,?)", batch)
    rows += len(batch)
db.execute("CREATE INDEX zero_index_t ON zero_index (t)")
db.execute("CREATE INDEX zero_index_N ON zero_index (N)")
db.commit()
print(f"index.db: {rows} blocks from {len(list(data.glob('zeros_*.dat')))} files")
