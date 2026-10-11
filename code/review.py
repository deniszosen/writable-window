# -*- coding: utf-8 -*-
"""Every check this manuscript has to pass, in one place.

Run it, fix what it reports, run it again. Revision introduces defects that
reading does not catch: a duplicated paragraph, a stale panel reference, a
number that no longer matches the data. These checks do.
"""
import json, re, subprocess, sys, itertools, difflib, collections
from pathlib import Path

HERE = Path(__file__).resolve().parent
MD = HERE.parent / "markdown/manuscript.md"
SUP = HERE.parent / "markdown/supplementary.md"

fails = []
def chk(name, ok, note=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}{'  ' + note if note else ''}")
    if not ok:
        fails.append(name)


def main():
    t = MD.read_text()
    sup = SUP.read_text()
    body = t[:t.index("## References")]
    paras = [p.strip() for p in body.split("\n\n") if p.strip()]

    print("STRUCTURE")
    heads = re.findall(r"^#{2,3} (.+)$", body, re.M)
    secs = {h.split()[0] for h in heads if re.match(r"^\d\.\d ", h)}
    chk("7 Results subsections", len([s for s in secs if s.startswith("2.")]) == 7)
    chk("8 Methods subsections", len([s for s in secs if s.startswith("4.")]) == 8)
    refs = set(re.findall(r"[Ss]ection (\d\.\d)", body))
    chk("no dangling section reference", refs <= secs,
        f"({sorted(refs - secs) or 'none'})")
    caps = re.findall(r"\*\*Figure (\d)\.", body)
    chk("figures captioned in order", caps == [str(i) for i in range(1, len(caps) + 1)],
        f"({caps})")
    inline = set(re.findall(r"figure (\d)", body))
    chk("no dangling figure reference", inline <= set(caps), f"({sorted(inline)})")
    links = re.findall(r"!\[Figure (\d)\]", body)
    chk("one image per caption", links == caps)
    for cn in caps:
        cap = next(p for p in paras if p.startswith(f"**Figure {cn}."))
        # the caption writes panels as (**a**), which r"\((\w)\)" cannot
        # match, so this check silently skipped every figure until now
        panels = set(re.findall(r"\(\*\*(\w)\*\*\)", cap))
        # "figures 1b and 1c" is as valid as two separate references, and a
        # pattern anchored on the singular "figure 1b" sees neither of them
        used = set()
        for m in re.finditer(
                rf"figures?\s+{cn}([a-e])\b((?:\s*(?:,|and)\s*(?:{cn})?[a-e]\b)*)",
                body):
            used.add(m.group(1))
            used.update(re.findall(rf"(?:,|and)\s*(?:{cn})?([a-e])\b", m.group(2)))
        if panels:
            chk(f"figure {cn}: every panel referenced", used == panels,
                f"(caption {sorted(panels)}, text {sorted(used)})")
    tabs = set(re.findall(r"Table (S\d)", t))
    have = set(re.findall(r"## (Table S\d)", sup)) | {"S" + str(i) for i in (1, 2, 3)
                                                      if f"Table S{i}" in sup}
    chk("supplementary tables exist", all(x in sup for x in tabs), f"({sorted(tabs)})")

    print("\nDUPLICATION")
    sents = [s.strip() for s in re.split(r"(?<=[.!?])\s+", re.sub(r"\n+", " ", body))
             if len(s.split()) >= 10]
    dupes = []
    for a, b in itertools.combinations(range(len(sents)), 2):
        if abs(len(sents[a]) - len(sents[b])) > 60:
            continue
        r = difflib.SequenceMatcher(None, sents[a], sents[b]).ratio()
        if r > 0.72:
            dupes.append((round(r, 2), sents[a][:85], sents[b][:85]))
    chk("no near-duplicate sentences", not dupes, f"({len(dupes)})")
    for r, a, b in dupes[:6]:
        print(f"        [{r}] A: {a}\n              B: {b}")
    pc = collections.Counter(p for p in paras if len(p.split()) > 12)
    chk("no repeated paragraph", max(pc.values(), default=0) == 1)

    print("\nREFERENCES")
    nrefs = len(re.findall(r"^\d+\. [A-Z]", t[t.index("## References"):], re.M))
    cited = set(int(x) for g in re.findall(r"\[([\d,\s]+)\]", t + sup)
                for x in g.split(",") if x.strip().isdigit())
    chk(f"all {nrefs} references cited", all(i in cited for i in range(1, nrefs + 1)),
        f"({sorted(i for i in range(1, nrefs + 1) if i not in cited) or 'none missing'})")
    chk("no citation beyond the list", max(cited) <= nrefs, f"(max [{max(cited)}])")
    titled = [l for l in re.findall(r"^\d+\. (.+)$", t[t.index("## References"):], re.M)
              if re.search(r"[.?]\s\*", l)]
    chk("every reference carries a title", len(titled) == nrefs, f"({len(titled)}/{nrefs})")

    print("\nTERMINOLOGY AND STYLE")
    stray = [m.group(0) for m in
             re.finditer(r"(?<!ticker )(?<!ticker-)\btapes?\b", body)]
    chk("'array' not 'tape'", not stray,
        "(ticker tape is a device name and is allowed)")
    chk("no bare 'editor'", not re.search(r"(?<!prime )(?<!Prime )\beditor\b", body))
    chk("arms vs mechanisms consistent", "candidate mechanisms" not in body)
    for ab, full in (("dNTP", "deoxynucleoside"), ("FEN1", "flap endonuclease"),
                     ("AAV", "adeno-associated"), ("LIG1", "DNA ligase")):
        i = body.find(ab)
        chk(f"{ab} expanded at first use", i < 0 or full in body[max(0, i - 160):i + 60])
    chk("no markdown horizontal rules", "\n---\n" not in body)

    print("\nLENGTH")
    cap_w = sum(len(p.split()) for p in paras if p.startswith("**Figure "))
    head_w = sum(len(p.split()) for p in paras if p.startswith("#"))
    ab = t[t.index("## Abstract") + 12:t.index("**Keywords")]
    bw = len(body.split()) - cap_w - head_w - len(ab.split()) - 40
    # Physical Biology: "Research Papers ... not normally more than 10,000
    # words". IOP does not say whether Methods counts, and the two conventions
    # give very different answers for this paper, so both are reported and the
    # cover letter quotes both rather than picking the flattering one. The
    # check passes if the paper is inside the guideline on either convention.
    mstart = body.index("## 4. Methods")
    meth = body[mstart:]
    mw = (len(meth.split())
          - sum(len(p.split()) for p in paras if p.startswith("#") and p in meth))
    chk("body text <= 10,000 words on at least one counting convention",
        min(bw, bw - mw) <= 10000,
        f"(with Methods {bw:,}; without Methods {bw - mw:,})")
    chk("abstract <= 300 words", len(ab.split()) <= 300, f"({len(ab.split())})")
    L = [len(s.split()) for s in re.split(r"(?<=[.!?])\s+", re.sub(r"\n+", " ", body))
         if s.strip()]
    chk("no sentence over 75 words", max(L) <= 75, f"(longest {max(L)})")

    # The supplementary was never checked, which is how two bare "editor"s and
    # an uncited Table S4 survived. Same rules, applied to it.
    print("\nSUPPLEMENTARY")
    import re as _re
    sup_tables = _re.findall(r"^## (Table S\d)", sup, _re.M)
    cited = set(_re.findall(r"Table (S\d)", t)) | set(_re.findall(r"Table (S\d)", sup))
    uncited = [x for x in sup_tables if x.split()[1] not in cited]
    chk("every supplementary table is cited", not uncited, f"({uncited})" if uncited else "")
    sup_figs = _re.findall(r"\*\*Figure (S\d)", sup)
    nofig = [f for f in sup_figs if f"figure {f}".lower() not in t.lower()]
    chk("every supplementary figure is cited in the text", not nofig,
        f"({nofig})" if nofig else "")
    bare = _re.findall(r"(?<!prime )(?<!Prime )\beditor\b", sup)
    chk("no bare 'editor' in the supplementary", not bare, f"({len(bare)})" if bare else "")
    us = _re.findall(r"\b\w*(?:ization|ized|izing|yze|behavior|color|modeling)\w*\b", sup)
    chk("supplementary uses British spelling", not us, f"({set(us)})" if us else "")
    chk("supplementary title <= 30 characters",
        len(_re.search(r"\*\*Title:\*\* (.+)", sup).group(1).strip()) <= 30)
    chk("supplementary description <= 30 words",
        len(_re.search(r"\*\*Description:\*\* (.+)", sup).group(1).split()) <= 30)

    print("\nNUMBERS")
    r = subprocess.run([sys.executable, str(HERE / "verify_numbers.py")],
                       capture_output=True, text=True)
    chk("every mapped claim matches the data",
        "every mapped claim matches" in r.stdout, "")
    if "MISMATCH" in r.stdout:
        print("        " + r.stdout.split("MISMATCH")[1][:400])

    print()
    print("RESULT:", "clean" if not fails else f"{len(fails)} to fix: {fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
