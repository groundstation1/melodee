"""Print the named functions' disassembly from build/melodee.dis, with each one's size."""
import re, sys
want = set(sys.argv[2:]); cur = None; out = {}
for line in open(sys.argv[1]):
    m = re.match(r"^([0-9a-f]+) <([A-Za-z_][A-Za-z_0-9.]*)>:$", line.strip()) or re.match(r"^([A-Za-z_][A-Za-z_0-9.]*):$", line.strip())
    if m:
        name = m.group(m.lastindex)
        if not name.startswith("."): cur = name if name in want else None
        if cur: out.setdefault(cur, [])
        continue
    if cur: out[cur].append(line.rstrip())
for n, lines in out.items():
    addrs = [int(l.split(":")[0], 16) for l in lines if re.match(r"^\s*[0-9a-f]+:", l)]
    print(f"==== {n}: {len(addrs)} instructions, {addrs[-1] - addrs[0] + 4 if addrs else 0} bytes")
    print("\n".join(lines))
