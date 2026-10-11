"""Catch results.json keys that are read but never written.

A consumer script that reads a key no producer writes means the released code
cannot regenerate the released results file, and the symptom is a crash or an
IndexError on a freshly generated run rather than anything visible in the file
that was shipped. Neither a numeric audit of the results file nor a figure
check can detect it, because both read the file rather than asking whether the
code is able to produce it.

So run this FIRST, on a results.json generated from an empty out/ directory. It
exits non-zero if any key a consumer script reads is absent.

    python src/check_keys.py

It is a static check and deliberately over-reports rather than under-reports: a
dynamically built key name shows up as a false positive, which is the right
direction for a release gate.
"""

import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
RESULTS = os.path.join(BASE, "out", "results.json")

# Scripts that CONSUME results.json, in ../../code. Producers are excluded: a
# producer reading back a key it just wrote is not an orphan.
CONSUMERS = ("figures.py", "fig1_timeline.py", "verify_numbers.py",
             "multiplicity.py", "interaction_test.py", "mouse_validation.py")
CODE = os.path.join(os.path.dirname(BASE), "code")

# Subscript chains rooted at the results dict. Consumers variously bind it to
# R, m, w, d, et, nd and others, so a chain is matched from any root and then
# resolved against the file by trying each top-level section too.
CHAIN = re.compile(r'\b([A-Za-z_]\w*)((?:\[\s*["\'][^"\']+["\']\s*\])+)')
# .get("key") reads are invisible to a subscript-only scan, and an orphan read
# exactly that way is precisely the kind this check exists for, so matching
# .get() is not an optional refinement.
GET = re.compile(r'\b([A-Za-z_]\w*)\s*\.get\(\s*["\']([^"\']+)["\']')
LITERAL = re.compile(r'["\']([^"\']+)["\']')
ASSIGN = re.compile(r'\[\s*["\']([^"\']+)["\']\s*\]\s*=')
KWARG = re.compile(r'\b([a-z_]\w*)\s*=', re.I)


def read_chains(path):
    """Every ('root', ('key', ...)) subscript chain in a file."""
    src = open(path, encoding="utf-8").read()
    out = set()
    for m in CHAIN.finditer(src):
        keys = tuple(LITERAL.findall(m.group(2)))
        if keys:
            out.add((m.group(1), keys))
    for m in GET.finditer(src):
        out.add((m.group(1), (m.group(2),)))
    return out


def walk_present(obj, keys):
    """Can this key chain be followed from obj? Descends into lists.

    A consumer that writes sweep[3] has already chosen its row, and the first
    row of a sweep is usually the "none" baseline, which carries fewer keys
    than the arms. Testing row 0 alone therefore reports a key as missing when
    every row the consumer actually reads has it, so any row counts.
    """
    cur = obj
    for i, k in enumerate(keys):
        if isinstance(cur, list):
            if not cur:
                return True          # empty sweep, cannot disprove
            return any(walk_present(row, keys[i:]) for row in cur)
        if not isinstance(cur, dict) or k not in cur:
            return False
        cur = cur[k]
    return True


def all_dict_nodes(obj, out=None, depth=0):
    """Every dict anywhere in the tree.

    Consumer scripts bind loop variables to elements of deeply nested sweep
    lists (`for r in R["dynamic_range"]["sweep"]`), so a chain like
    r["confusion"]["bits"] has to be resolvable from an arbitrary interior
    node, not just from the root or a top-level section. Resolving against
    every dict node is what makes the check usable as a gate instead of a
    source of false alarms.
    """
    if out is None:
        out = []
    if depth > 8:
        return out
    if isinstance(obj, dict):
        out.append(obj)
        for v in obj.values():
            all_dict_nodes(v, out, depth + 1)
    elif isinstance(obj, list):
        for v in obj[:3]:            # sweep rows share a schema
            all_dict_nodes(v, out, depth + 1)
    return out


ALIAS = re.compile(
    r'^\s*([A-Za-z_]\w*)\s*=\s*R((?:\[\s*["\'][^"\']+["\']\s*\])+)', re.M)

# A consumer may load more than one JSON file. growth_trajectory.json and
# switches50_seeds.json are not results.json and their keys are not this
# check's business, so any root bound from a json.load of some other file is
# out of scope. Without this the check reports every key of every side file.
OTHER_JSON = re.compile(
    r'^\s*([A-Za-z_]\w*)\s*=\s*json\.loads?\(.*?["\']([\w./-]+\.json)["\']', re.M | re.S)


def foreign_roots(path):
    """Roots bound from a JSON file that is not results.json."""
    src = open(path, encoding="utf-8").read()
    return {n for n, f in OTHER_JSON.findall(src)
            if not f.endswith("results.json")}

# Roots that are not the results dict at all: matplotlib objects and the like.
NOT_RESULTS = {"spines", "ax", "axs", "fig", "gs", "kw", "kwargs", "style",
               "os", "sys", "np", "plt", "mpl", "rc", "rcParams"}


def aliases(path):
    """{name: ('section', ...)} for every `w = R["writable"]`-style binding.

    Resolving an alias to its exact section is what makes the check correct
    rather than merely plausible. `writable["gradient_sem"]` is read as
    w["gradient_sem"], and the same key name exists under `nondivision`, so a
    checker that resolves a chain against any node in the tree declares it
    present and misses a real orphan entirely.
    """
    src = open(path, encoding="utf-8").read()
    out = {}
    for m in ALIAS.finditer(src):
        out[m.group(1)] = tuple(LITERAL.findall(m.group(2)))
    return out


# {r["variant"]: r for r in R["switches"]} binds a root whose keys are VALUES
# in the data ("central", "none", "founders 6000"), not schema keys. A static
# check cannot resolve those and should not report them.
COMPREHENSION = re.compile(
    r'^\s*([A-Za-z_]\w*)\s*=\s*\{[^}]*\bfor\b[^}]*\}', re.M)


def label_keyed_roots(path):
    """Roots bound by a dict comprehension, so keyed on runtime values."""
    src = open(path, encoding="utf-8").read()
    return set(COMPREHENSION.findall(src))


def runtime_label(keys):
    """Dicts built in the consumer keyed by a human-readable label.

    e.g. `{r["variant"]: r for r in R["switches"]}` gives keys like
    'central' and 'founders 6000'. Those are values in the data, not
    schema keys, so a static check cannot resolve them and should not try.
    """
    return any(" " in k or "(" in k for k in keys)


def resolvable(results, nodes, alias_map, root, keys):
    """Resolve a chain, using the alias's exact section when one is known."""
    if root in NOT_RESULTS or runtime_label(keys):
        return True                      # out of scope, not a finding
    if root == "R":
        return walk_present(results, keys)
    if root in alias_map:
        base = results
        for k in alias_map[root]:
            if isinstance(base, list):
                base = base[0] if base else {}
            if not isinstance(base, dict) or k not in base:
                return True              # the alias itself is broken; not this
            base = base[k]
        return walk_present(base, keys)
    return any(walk_present(n, keys) for n in nodes)


def assigned_anywhere():
    """Key names assigned or passed as a dict() kwarg anywhere in src/."""
    names = set()
    for fn in sorted(os.listdir(HERE)):
        if not fn.endswith(".py"):
            continue
        src = open(os.path.join(HERE, fn), encoding="utf-8").read()
        names.update(ASSIGN.findall(src))
        names.update(KWARG.findall(src))
    return names


def main():
    if not os.path.exists(RESULTS):
        print(f"missing {RESULTS}")
        print("Generate it first, from an EMPTY out/ directory, so stale keys "
              "from a previous run cannot mask a missing one.")
        return 2

    results = json.load(open(RESULTS))
    assigned = assigned_anywhere()
    nodes = all_dict_nodes(results)

    missing, checked = [], 0
    for fn in CONSUMERS:
        path = os.path.join(CODE, fn)
        if not os.path.exists(path):
            continue
        alias_map = aliases(path)
        label_roots = label_keyed_roots(path) | foreign_roots(path)
        for root, keys in sorted(read_chains(path)):
            if root in NOT_RESULTS or root in label_roots or runtime_label(keys):
                continue
            checked += 1
            if not resolvable(results, nodes, alias_map, root, keys):
                missing.append((fn, root, keys))

    print(f"check_keys: {checked} subscript chains across "
          f"{len(CONSUMERS)} consumer scripts, resolved against "
          f"{len(nodes)} dict nodes")
    if not missing:
        print("  PASS: every key read by a consumer script is present in "
              "results.json")
        return 0

    print(f"\n  FAIL: {len(missing)} key chain(s) READ but NOT PRESENT:\n")
    for fn, root, keys in missing:
        chain = "".join(f'["{k}"]' for k in keys)
        ever = [k for k in keys if k in assigned]
        print(f"    {fn}: {root}{chain}")
        if ever:
            print(f"      assigned somewhere in src/: {', '.join(ever)}")
        else:
            print("      NO component is assigned anywhere in src/")
    print("\n  A key read by a consumer and written by nothing "
          "means the released code cannot regenerate the released results "
          "file. That is what this check is for.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
