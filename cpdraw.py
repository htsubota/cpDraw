#!/usr/bin/env python3
"""
cpDraw v1.0.0

A stand-alone Python program for generating editable circular
plastid genome maps from annotated GenBank files.

Input:
    GenBank flat file

Output:
    SVG genome map
    Optional label-placement diagnostic CSV

Requirement:
    Biopython

Copyright (c) 2026 Hiromi Tsubota
Licensed under the MIT License.
See the LICENSE file for details.
"""

from __future__ import annotations
import argparse, json, math, re, sys
from collections import Counter
from pathlib import Path
from xml.sax.saxutils import escape

VERSION = "1.0.0"
PROGRAM_NAME = "cpDraw"
LEGACY_PROGRAM = "plastome_draw"
LEGACY_VERSION = "1.1.8"

PALETTE = {
    "photosystem I": "#00733E", "photosystem II": "#2E8B2D",
    "cytochrome b/f complex": "#7C9E16", "ATP synthase": "#9FC400",
    "NADH dehydrogenase": "#FFE600", "RubisCO large subunit": "#00A63C",
    "photosystem assembly/stability factors": "#FFFDE0", "RNA polymerase": "#C81020",
    "ribosomal proteins (SSU)": "#E3B77D", "ribosomal proteins (LSU)": "#9A7040",
    "transfer RNAs": "#173A91", "ribosomal RNAs": "#ED001C",
    "clpP, matK": "#F56600", "other genes": "#AF2694",
    "hypothetical chloroplast reading frames (ycf)": "#D9D9D9", "ORFs": "#55B7A8",
    "origin of replication": "#F47B86", "polycistronic transcripts": "#E000E8",
}
LEGEND_ORDER = list(PALETTE)

def category(name, feature_type):
    n = name.lower()
    if feature_type == "tRNA" or n.startswith("trn"): return "transfer RNAs"
    if feature_type == "rRNA" or n.startswith("rrn"): return "ribosomal RNAs"
    if n.startswith("psa"): return "photosystem I"
    if n.startswith("psb"): return "photosystem II"
    if n.startswith("pet"): return "cytochrome b/f complex"
    if n.startswith("atp"): return "ATP synthase"
    if n.startswith("ndh"): return "NADH dehydrogenase"
    if n == "rbcl": return "RubisCO large subunit"
    if n in {"pbf1","psai","psaj","ycf3","ycf4"}: return "photosystem assembly/stability factors"
    if n.startswith("rpo"): return "RNA polymerase"
    if n.startswith("rps"): return "ribosomal proteins (SSU)"
    if n.startswith("rpl"): return "ribosomal proteins (LSU)"
    if n in {"clpp","matk"}: return "clpP, matK"
    if n.startswith("ycf"): return "hypothetical chloroplast reading frames (ycf)"
    if n.startswith("orf"): return "ORFs"
    return "other genes"

# ---------- GenBank parsing (Biopython imported lazily so the drawing code can be unit-tested) ----------
def qualifier_first(feature, *keys):
    for key in keys:
        values = feature.qualifiers.get(key)
        if values:
            return str(values[0])
    return ""

def feature_parts(feature, with_strand=False):
    """Return location parts; optionally retain the strand of each part.

    Biopython CompoundLocation parts can carry their own strand.  Keeping that
    information is essential for trans-spliced genes such as rps12.
    """
    from Bio.SeqFeature import CompoundLocation
    loc = feature.location
    src = loc.parts if isinstance(loc, CompoundLocation) else [loc]
    if with_strand:
        fallback = loc.strand or 1
        return [(int(p.start), int(p.end), p.strand if p.strand is not None else fallback)
                for p in src]
    return [(int(p.start), int(p.end)) for p in src]

def collect_features(record, trans_gap=10000, show_pseudogene=False):
    """Collect non-redundant drawable gene features.

    Typed records (tRNA, rRNA, CDS) are preferred to generic gene records when
    they describe the same biological feature. Auxiliary single-part gene
    records that overlap multipart components or exactly fill an intron interval
    are suppressed. Separate genomic copies remain distinct.
    """
    accepted = {"gene", "CDS", "tRNA", "rRNA"}
    priority = {"tRNA": 0, "rRNA": 0, "CDS": 1, "gene": 2}
    entries = []
    for order, f in enumerate(record.features):
        if f.type not in accepted:
            continue
        name = qualifier_first(f, "gene", "locus_tag", "product")
        if not name:
            continue
        raw_parts = sorted(feature_parts(f, with_strand=True),
                           key=lambda q: (q[0], q[1]))
        parts = [(a, b) for a, b, _ in raw_parts]
        part_strands = [st for _, _, st in raw_parts]
        strand = f.location.strand or (part_strands[0] if part_strands else 1) or 1
        entries.append({"feature": f, "name": name, "norm": name.lower(),
                        "order": order, "priority": priority[f.type],
                        "raw_parts": raw_parts, "parts": parts,
                        "part_strands": part_strands, "strand": strand})

    # Exact structural duplicates: keep the biologically typed feature.
    exact = {}
    for e in entries:
        key = (e["norm"], tuple(e["raw_parts"]), e["strand"])
        old = exact.get(key)
        if old is None or (e["priority"], e["order"]) < (old["priority"], old["order"]):
            exact[key] = e
    entries = sorted(exact.values(), key=lambda e: e["order"])

    by_name = {}
    for e in entries:
        by_name.setdefault(e["norm"], []).append(e)

    def strand_compatible(a, b):
        if a["strand"] == b["strand"]:
            return True
        return any(st == a["strand"] for st in b["part_strands"])

    def typed_supersedes_gene(gene, typed):
        if gene["feature"].type != "gene" or typed["feature"].type == "gene":
            return False
        if not strand_compatible(gene, typed):
            return False
        # Every typed component must be contained in a component envelope of the
        # generic gene. This supports unequal exon granularity (e.g. rps12).
        return all(any(ga <= a and b <= gb for ga, gb in gene["parts"])
                   for a, b in typed["parts"])

    def auxiliary_fragment(gene, multipart):
        if gene["feature"].type != "gene" or len(gene["parts"]) != 1:
            return False
        if len(multipart["parts"]) < 2 or not strand_compatible(gene, multipart):
            return False
        a, b = gene["parts"][0]
        # A component-overlapping generic feature is redundant.
        if any(max(a, x) < min(b, y) for x, y in multipart["parts"]):
            return True
        # A feature exactly matching a gap between components is an intron
        # annotation, not a second copy of the gene.
        for (_, left_end), (right_start, _) in zip(multipart["parts"],
                                                    multipart["parts"][1:]):
            if a == left_end and b == right_start:
                return True
        return False

    filtered = []
    for e in entries:
        peers = by_name[e["norm"]]
        suppress = any(typed_supersedes_gene(e, p) for p in peers if p is not e)
        if not suppress:
            suppress = any(auxiliary_fragment(e, p) for p in peers if p is not e)
        if not suppress:
            filtered.append(e)

    out = []
    for e in filtered:
        f = e["feature"]
        out.append(build_feature(e["name"], f.type, e["strand"], e["parts"],
                                 "pseudo" in f.qualifiers, trans_gap,
                                 part_strands=e["part_strands"],
                                 show_pseudogene=show_pseudogene))
    return out

def build_feature(name, ftype, strand, ps, pseudo=False, trans_gap=10000,
                  part_strands=None, show_pseudogene=False):
    """Turn sorted 0-based parts into a drawable feature dict.

    ``part_strands`` parallels ``ps`` and is used only when separate parts are
    drawn. Ordinary cis-spliced features retain the v0.8.5 behaviour.
    """
    if part_strands is None:
        part_strands = [strand] * len(ps)
    paired = sorted(zip(ps, part_strands), key=lambda x: (x[0][0], x[0][1]))
    ps = [p for p, _ in paired]
    part_strands = [s if s is not None else strand for _, s in paired]
    gaps = [ps[i+1][0]-ps[i][1] for i in range(len(ps)-1)]
    is_trans = bool(gaps and max(gaps) >= trans_gap)
    is_spliced = len(ps) > 1
    draw_parts = ps if is_trans else [(ps[0][0], ps[-1][1])]
    draw_part_strands = part_strands if is_trans else [strand]
    label = (("psi " + name) if (pseudo and show_pseudogene) else name) + ("*" if is_spliced else "")
    return {"name": name, "label": label, "type": ftype, "strand": strand,
            "parts": draw_parts, "part_strands": draw_part_strands,
            "source_parts": ps, "source_part_strands": part_strands, "is_trans": is_trans,
            "pseudo": pseudo, "category": category(name, ftype)}

def collect_regions(record):
    regions = []
    for f in record.features:
        if f.type not in {"misc_feature", "repeat_region"}:
            continue
        note = " ".join(f.qualifiers.get("note", [])).lower()
        lab = None
        if "large single copy" in note or re.search(r"\blsc\b", note): lab = "LSC"
        elif "small single copy" in note or re.search(r"\bssc\b", note): lab = "SSC"
        elif "inverted repeat" in note or re.search(r"\bir[ab]?\b", note):
            lab = "IRA" if re.search(r"\bira\b", note) else ("IRB" if re.search(r"\birb\b", note) else "IR")
        if lab:
            p = feature_parts(f)[0]
            regions.append((lab, p[0], p[1]))
    return regions

# ---------- geometry ----------
def polar(cx, cy, r, deg):
    a = math.radians(deg - 90)
    return cx + r*math.cos(a), cy + r*math.sin(a)

def ang_sep(a, b):
    d = abs(a-b) % 360
    return min(d, 360-d)

def annulus(cx, cy, r0, r1, a0, a1):
    p1 = polar(cx,cy,r1,a0); p2 = polar(cx,cy,r1,a1)
    p3 = polar(cx,cy,r0,a1); p4 = polar(cx,cy,r0,a0)
    large = 1 if (a1-a0) % 360 > 180 else 0
    return (f"M{p1[0]:.2f},{p1[1]:.2f} A{r1:.2f},{r1:.2f} 0 {large} 1 {p2[0]:.2f},{p2[1]:.2f} "
            f"L{p3[0]:.2f},{p3[1]:.2f} A{r0:.2f},{r0:.2f} 0 {large} 0 {p4[0]:.2f},{p4[1]:.2f} Z")

def orientation_offset(regions, L):
    for lab, s, e in regions:
        if lab == "LSC":
            return 360.0 * ((s + e) / 2.0) / L
    return 0.0

def mapped_angle(pos, L, offset):
    return (offset - 360.0 * pos / L) % 360.0

def gc_sliding(seq, window, step):
    """Circular sliding-window GC. Returns [(center_bp, gc)], genome_mean.
    center_bp is NOT reduced mod n so angles stay monotonic; wrap-around is used
    so terminal windows keep the full width. No clipping anywhere."""
    s = str(seq).upper(); n = len(s)
    if n == 0: return [], 0.0
    ext = s + s[:window]
    res = []; i = 0
    while i < n:
        chunk = ext[i:i+window]
        valid = chunk.count("A")+chunk.count("C")+chunk.count("G")+chunk.count("T")
        gc = (chunk.count("G")+chunk.count("C"))/valid if valid else 0.0
        res.append((i + window/2.0, gc)); i += step
    den = s.count("A")+s.count("C")+s.count("G")+s.count("T")
    mean = (s.count("G")+s.count("C"))/den if den else 0.0
    return res, mean

# ---------- label placement (EPS-derived) ----------
def estimated_text_width(label, font_size):
    """Conservative Arial-width estimate used only for collision detection."""
    narrow = set("ilI1.,:;|'")
    wide = set("mwMWQO@%")
    units = 0.0
    for ch in label:
        if ch in narrow: units += 0.30
        elif ch in wide: units += 0.86
        elif ch.isupper(): units += 0.66
        else: units += 0.56
    return max(font_size * units, font_size * 1.2)


def label_box(z, cx, cy, font_size):
    """Return the four corners of the rendered radial label.

    Unlike the old axis-aligned estimate, this model respects SVG start/end
    anchoring: outer labels extend away from the circle and inner labels extend
    toward the centre.  This fixes missed overlaps between neighbouring labels.
    """
    ang = z["angle"] % 360.0
    x, y = polar(cx, cy, z["radius"], ang)
    w = estimated_text_width(z["f"]["label"], font_size)
    h = font_size * 1.18
    a = math.radians(ang - 90.0)
    ux, uy = math.cos(a), math.sin(a)       # outward radial unit vector
    tx, ty = -uy, ux                       # tangential unit vector
    # SVG anchor lies at the block-facing end of the text.
    if z["strand"] > 0:                   # outer text extends outward
        q0, q1 = 0.0, w
    else:                                  # inner text extends inward
        q0, q1 = -w, 0.0
    return [
        (x + q0*ux - h/2*tx, y + q0*uy - h/2*ty),
        (x + q1*ux - h/2*tx, y + q1*uy - h/2*ty),
        (x + q1*ux + h/2*tx, y + q1*uy + h/2*ty),
        (x + q0*ux + h/2*tx, y + q0*uy + h/2*ty),
    ]


def boxes_overlap(a, b, padding=0.8):
    """Separating-axis test for two rotated label rectangles."""
    def axes(poly):
        for i in range(4):
            x1,y1=poly[i]; x2,y2=poly[(i+1)%4]
            dx,dy=x2-x1,y2-y1; n=math.hypot(dx,dy)
            if n: yield (-dy/n, dx/n)
    for ax,ay in list(axes(a))+list(axes(b)):
        pa=[x*ax+y*ay for x,y in a]; pb=[x*ax+y*ay for x,y in b]
        if max(pa)+padding <= min(pb) or max(pb)+padding <= min(pa):
            return False
    return True


def place_labels(group, natural_radius, displaced_radius, font_size, cx, cy,
                 collision_angle_step=0.18, max_angle_shift=7.0,
                 collision_padding=1.2):
    """Two-ring packing with final rotated-box verification (v0.8.5).

    Natural labels are retained whenever possible. A conflicting label moves to
    the OGDRAW displaced ring. Displaced labels are packed in genomic order,
    which prevents line crossings, and every candidate is tested against all
    already placed labels. A final audit repairs any residual overlap.
    """
    if not group:
        return
    group.sort(key=lambda z: z["natural"] % 360.0)
    accepted=[]; displaced=[]
    for z in group:
        z.update(angle=z["natural"]%360.0, radius=natural_radius, lane=0,
                 leader=False, leader_reason="natural", angular_shift=0.0)
        poly=label_box(z,cx,cy,font_size)
        if any(boxes_overlap(poly,p,collision_padding) for _,p in accepted):
            z.update(radius=displaced_radius, lane=1, leader=True,
                     leader_reason="collision-displaced")
            displaced.append(z)
        else:
            accepted.append((z,poly))

    # Pack displaced labels. Prefer tiny shifts, but preserve circular order.
    placed=[]; previous=None
    step=max(0.05,collision_angle_step)
    steps=max(1,int(max_angle_shift/step))
    for z in displaced:
        base=z["natural"]%360.0
        candidates=[0.0]
        for n in range(1,steps+1): candidates.extend((n*step,-n*step))
        best=None
        for shift in candidates:
            cand=(base+shift)%360.0
            # Unwrap relative to previous displaced label and preserve order.
            unwrap=cand
            if previous is not None:
                while unwrap <= previous: unwrap += 360.0
                if unwrap-previous > 180.0: unwrap -= 360.0
                if unwrap <= previous: continue
            z["angle"]=cand
            poly=label_box(z,cx,cy,font_size)
            if not any(boxes_overlap(poly,p,collision_padding) for _,p in accepted+placed):
                best=(shift,unwrap,poly); break
        if best is None:
            # Deterministic outward sweep beyond the nominal limit.
            for n in range(steps+1,steps+61):
                shift=n*step; cand=(base+shift)%360.0; z["angle"]=cand
                poly=label_box(z,cx,cy,font_size)
                unwrap=cand
                if previous is not None:
                    while unwrap <= previous: unwrap += 360.0
                if not any(boxes_overlap(poly,p,collision_padding) for _,p in accepted+placed):
                    best=(shift,unwrap,poly); break
        if best is None:
            z["angle"]=base; best=(0.0,base,label_box(z,cx,cy,font_size))
            z["leader_reason"]="collision-displaced+unresolved"
        else:
            z["leader_reason"]=("collision-displaced" if abs(best[0])<1e-9
                                else "collision-displaced+angle")
        z["angular_shift"]=abs(best[0]); previous=best[1]
        placed.append((z,best[2]))

    # Final audit: if an overlap survived because of the circular seam, advance
    # only the later displaced label until the polygons separate.
    all_items=accepted+placed
    for i in range(len(all_items)):
        zi,pi=all_items[i]
        for j in range(i):
            if not boxes_overlap(pi,all_items[j][1],collision_padding): continue
            if not zi.get("leader"):
                zi.update(radius=displaced_radius,lane=1,leader=True,
                          leader_reason="final-audit-displaced")
            for n in range(1,81):
                zi["angle"]=(zi["angle"]+step)%360.0
                pi=label_box(zi,cx,cy,font_size)
                if not any(boxes_overlap(pi,all_items[q][1],collision_padding)
                           for q in range(i)):
                    zi["angular_shift"]=ang_sep(zi["angle"],zi["natural"])
                    break
            all_items[i]=(zi,pi)

def _orient(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def _segments_intersect(a,b,c,d):
    o1,o2,o3,o4=_orient(a,b,c),_orient(a,b,d),_orient(c,d,a),_orient(c,d,b)
    return o1*o2 < -1e-9 and o3*o4 < -1e-9


def _point_segment_distance(p,a,b):
    vx,vy=b[0]-a[0],b[1]-a[1]; wx,wy=p[0]-a[0],p[1]-a[1]
    vv=vx*vx+vy*vy
    if vv == 0: return math.hypot(wx,wy)
    t=max(0.0,min(1.0,(wx*vx+wy*vy)/vv))
    q=(a[0]+t*vx,a[1]+t*vy)
    return math.hypot(p[0]-q[0],p[1]-q[1])


def segment_polygon_clearance(a,b,poly):
    """Minimum distance from a segment to a convex label polygon (0 if crossing)."""
    best=min(_point_segment_distance(p,a,b) for p in poly)
    for i in range(len(poly)):
        c,d=poly[i],poly[(i+1)%len(poly)]
        if _segments_intersect(a,b,c,d): return 0.0
        best=min(best,_point_segment_distance(a,c,d),_point_segment_distance(b,c,d))
    return best


def segment_segment_clearance(a,b,c,d):
    if _segments_intersect(a,b,c,d): return 0.0
    return min(_point_segment_distance(a,c,d),_point_segment_distance(b,c,d),
               _point_segment_distance(c,a,b),_point_segment_distance(d,a,b))


def route_leader(z, all_labels, prior_segments, prior_elbows, cx, cy, k,
                 label_clearance=4.2, line_clearance=2.0, elbow_clearance=5.0):
    """Choose a compact two-segment leader that avoids every other label.

    Several OGDRAW-like elbow radii and angular fractions are tested. The
    shortest clean candidate is selected; prior leader segments are also kept
    apart. This addresses lines that pass through or graze neighbouring text.
    """
    if z["strand"] > 0:
        anchor_r, end_r = 413*k, 451*k
        elbow_rs=[438*k,440*k,442*k,436*k,444*k]
    else:
        anchor_r, end_r = 362*k, 324*k
        elbow_rs=[345*k,347*k,349*k,343*k,351*k]
    p0=polar(cx,cy,anchor_r,z["natural"])
    p2=polar(cx,cy,end_r,z["angle"])
    delta=((z["angle"]-z["natural"]+180.0)%360.0)-180.0
    fractions=[0.58,0.42,0.72,0.25,0.85,0.0,1.0,-0.22,1.22,-0.45,1.45]
    own=id(z); candidates=[]

    # v1.0.4: Prefer a genuinely straight leader whenever the direct segment is
    # clean.  The direct route must satisfy the same label and prior-line
    # clearances as elbow routes.  If it fails any test, retain the established
    # two-segment OGDRAW-like routing below.
    straight_clean=True; straight_min_lab=1e9; straight_min_line=1e9
    for other,poly in all_labels:
        if id(other)==own:
            continue
        cl=segment_polygon_clearance(p0,p2,poly)
        straight_min_lab=min(straight_min_lab,cl)
        if cl < label_clearance:
            straight_clean=False
            break
    if straight_clean:
        for c,d in prior_segments:
            cl=segment_segment_clearance(p0,p2,c,d)
            straight_min_line=min(straight_min_line,cl)
            if cl < line_clearance:
                straight_clean=False
                break
    if straight_clean:
        # None marks a one-segment route; callers must not create a phantom elbow.
        return (False, math.dist(p0,p2), p0, None, p2,
                straight_min_lab, straight_min_line)

    # trnfM: prefer an exactly radial first leg from the block, followed by one
    # oblique segment to the label. No circumferential dogleg is introduced.
    trnfm_radial = (z["f"].get("name","").lower() == "trnfm")
    if trnfm_radial:
        fractions = [0.0]
    for er in elbow_rs:
        for frac in fractions:
            ea=(z["natural"]+frac*delta)%360.0
            p1=polar(cx,cy,er,ea)
            segs=[(p0,p1),(p1,p2)]
            clean=True; min_lab=1e9; min_line=1e9
            for other,poly in all_labels:
                if id(other)==own: continue
                for a,b in segs:
                    cl=segment_polygon_clearance(a,b,poly); min_lab=min(min_lab,cl)
                    if cl < label_clearance: clean=False; break
                if not clean: break
            if clean and any(math.dist(p1, ep) < elbow_clearance for ep in prior_elbows):
                clean = False
            if clean:
                for a,b in segs:
                    for c,d in prior_segments:
                        cl=segment_segment_clearance(a,b,c,d); min_line=min(min_line,cl)
                        if cl < line_clearance: clean=False; break
                    if not clean: break
            length=math.dist(p0,p1)+math.dist(p1,p2)
            # Prefer visible clearance over a marginally shorter line. This makes a
            # diagonal detour win when the radial route grazes a neighbouring name.
            clearance_penalty = 18.0 / max(min(min_lab, 20.0), 0.25)
            line_penalty = 3.0 / max(min(min_line, 10.0), 0.25) if prior_segments else 0.0
            elbow_penalty = sum(max(0.0, elbow_clearance*1.8-math.dist(p1,ep))*2.0
                                for ep in prior_elbows)
            score = length + clearance_penalty + line_penalty + elbow_penalty + 0.08*abs(frac-.58)
            candidates.append((not clean,score,p0,p1,p2,min_lab,min_line))
    candidates.sort(key=lambda q:(q[0],q[1]))
    return candidates[0]


def promote_labels_obscured_by_leaders(labels, cx, cy, k, font_size,
                                      visual_clearance=6.0):
    """Promote natural labels that sit in a displaced label's leader corridor.

    This handles cases such as trnY -> trnE: although the text boxes may not
    mathematically overlap, the line behind the neighbouring label makes the
    assignment visually ambiguous. The obstructing natural label is moved to
    the displaced ring and receives its own leader.
    """
    promoted=[]
    polygons=[(z,label_box(z,cx,cy,font_size)) for z in labels]
    leaders=[z for z in labels if z.get("leader")]
    naturals=[z for z in labels if not z.get("leader")]
    for lead in leaders:
        if lead["strand"] > 0:
            anchor_r, end_r, displaced_r = 413*k, 451*k, 453*k
        else:
            anchor_r, end_r, displaced_r = 362*k, 324*k, 322*k
        p0=polar(cx,cy,anchor_r,lead["natural"])
        p2=polar(cx,cy,end_r,lead["angle"])
        delta=((lead["angle"]-lead["natural"]+180.0)%360.0)-180.0
        p1=polar(cx,cy,(440 if lead["strand"]>0 else 347)*k,
                 (lead["natural"]+0.58*delta)%360.0)
        for z,poly in polygons:
            if z not in naturals or z["strand"] != lead["strand"]: continue
            clearance=min(segment_polygon_clearance(p0,p1,poly),
                          segment_polygon_clearance(p1,p2,poly))
            if clearance < visual_clearance*k:
                z.update(radius=displaced_r,lane=1,leader=True,
                         leader_reason="leader-corridor-promoted")
                promoted.append(z); naturals.remove(z)
    return promoted


def enforce_non_crossing_order(group, min_gap=0.06):
    """Preserve circular gene order in final label angles.

    Leader lines between two concentric rings cannot cross when the angular
    order of their block endpoints and label endpoints is the same.  This
    function keeps all existing collision decisions, radii, and shifts, and
    only repairs angular inversions with the smallest possible adjustment.
    The circular seam is placed at the largest natural-angle gap, where no
    dense label cluster is split.
    """
    if len(group) < 2:
        return 0
    natural_sorted=sorted(group,key=lambda z:z["natural"]%360.0)
    gaps=[]
    for i,z in enumerate(natural_sorted):
        a=z["natural"]%360.0
        b=(natural_sorted[(i+1)%len(natural_sorted)]["natural"]%360.0)+(360.0 if i==len(natural_sorted)-1 else 0.0)
        gaps.append(b-a)
    cut=(max(range(len(gaps)),key=lambda i:gaps[i])+1)%len(natural_sorted)
    ordered=natural_sorted[cut:]+natural_sorted[:cut]
    natural=[]
    for i,z in enumerate(ordered):
        a=z["natural"]%360.0
        if i:
            while a<=natural[-1]: a+=360.0
        natural.append(a)
    final=[]
    for a,z in zip(natural,ordered):
        f=z["angle"]%360.0
        while f-a>180.0: f-=360.0
        while f-a<=-180.0: f+=360.0
        final.append(f)
    repaired=0
    # Forward pass is sufficient to remove adjacent inversions; a short
    # backward balancing pass shares the correction and limits drift.
    for i in range(1,len(final)):
        required=final[i-1]+min_gap
        if final[i]<required:
            overlap=required-final[i]
            final[i-1]-=overlap/2.0
            final[i]+=overlap/2.0
            repaired+=1
    for i in range(len(final)-2,-1,-1):
        maximum=final[i+1]-min_gap
        if final[i]>maximum:
            overlap=final[i]-maximum
            final[i]-=overlap/2.0
            final[i+1]+=overlap/2.0
            repaired+=1
    for a,z,f in zip(natural,ordered,final):
        old=z["angle"]%360.0
        z["angle"]=f%360.0
        extra=abs(((z["angle"]-old+180.0)%360.0)-180.0)
        if extra>1e-7:
            z["leader"]=True
            previous=z.get("leader_reason","natural")
            z["leader_reason"]=(previous+"+uncross") if previous!="natural" else "uncross"
        z["angular_shift"]=abs(((z["angle"]-z["natural"]+180.0)%360.0)-180.0)
    return repaired


def repair_local_order_inversions(group, cluster_gap=8.0, max_extra_shift=7.0):
    """Repair only local final-angle inversions within one strand/lane.

    Labels are clustered by adjacent natural-angle gaps <= cluster_gap. Existing
    final-angle slots are reassigned only inside a cluster whose final order is
    actually inverted. A repair is rejected if any label would move more than
    max_extra_shift degrees from its v1.0.0 final angle.
    """
    if len(group) < 2:
        return 0
    ordered=sorted(group,key=lambda z:z["natural"]%360.0)
    # rotate at largest natural gap to avoid splitting a seam-spanning cluster
    gaps=[]
    for i,z in enumerate(ordered):
        a=z["natural"]%360.0
        b=ordered[(i+1)%len(ordered)]["natural"]%360.0
        if i==len(ordered)-1: b+=360.0
        gaps.append(b-a)
    cut=(max(range(len(gaps)),key=lambda i:gaps[i])+1)%len(ordered)
    ordered=ordered[cut:]+ordered[:cut]
    nat=[]
    for i,z in enumerate(ordered):
        a=z["natural"]%360.0
        if i:
            while a<=nat[-1]: a+=360.0
        nat.append(a)
    clusters=[]; start=0
    for i in range(1,len(ordered)):
        if nat[i]-nat[i-1] > cluster_gap:
            clusters.append(ordered[start:i]); start=i
    clusters.append(ordered[start:])
    repaired=0
    for cluster in clusters:
        if len(cluster)<2: continue
        centre=sum(z["natural"]%360.0 for z in cluster)/len(cluster)
        def unwrap_final(z):
            f=z["angle"]%360.0
            while f-centre>180.0: f-=360.0
            while f-centre<=-180.0: f+=360.0
            return f
        natural_order=sorted(cluster,key=lambda z:z["natural"]%360.0)
        final_order=sorted(cluster,key=unwrap_final)
        if [id(z) for z in natural_order] == [id(z) for z in final_order]:
            continue
        slots=sorted(unwrap_final(z) for z in cluster)
        candidate=list(zip(natural_order,slots))
        if any(ang_sep(slot%360.0,z["angle"]%360.0)>max_extra_shift for z,slot in candidate):
            continue
        for z,slot in candidate:
            new=slot%360.0
            if ang_sep(new,z["angle"]%360.0)>1e-7:
                z["angle"]=new
                z["angular_shift"]=ang_sep(new,z["natural"])
                z["leader_reason"]=z.get("leader_reason","natural")+"+local-order-preserved"
                z["order_adjusted"]=True; repaired+=1
            z["order_preserved"]=True
    return repaired

# ---------- SVG ----------
def make_svg(features, regions, seq, rec_id, organism, args):
    S = args.size; k = S/1200.0; cx = cy = S/2.0; L = len(seq)
    orient = orientation_offset(regions, L)
    backbone = 387*k; block = 24*k
    out_text_r = 414*k; in_text_r = 360*k
    # OGDRAW text-anchor radii measured from the original vector output.
    # With start/end anchoring these correspond to observed text-centre radii
    # of about 426/349 (natural) and 466/309 (displaced).
    out_displaced_r = 453*k; in_displaced_r = 322*k
    region_r = 276*k
    # OGDRAW EPS measurement: GC bars terminate at r=250 on a 1200-unit canvas.
    gc_outer = 260*k; gc_inner = 200*k
    leader_w = args.leader_width*k
    out = []; A = out.append
    A(f'<svg xmlns="http://www.w3.org/2000/svg" width="{S}" height="{S}" viewBox="0 0 {S} {S}">')
    # v1.0.9: embed compact provenance metadata without changing the drawing.
    svg_metadata = {
        "generator": PROGRAM_NAME,
        "version": VERSION,
        "legacy_program": LEGACY_PROGRAM,
        "legacy_version": LEGACY_VERSION,
        "record_id": str(rec_id),
        "organism": str(organism),
        "genome_length_bp": L,
        "input_file": Path(getattr(args, "genbank", "sequence.gb")).name,
    }
    A(f'<metadata id="cpDraw_metadata">{escape(json.dumps(svg_metadata, ensure_ascii=False, separators=(",", ":")))}</metadata>')
    A('<defs><marker id="greyArrow" markerWidth="5.2" markerHeight="4.2" refX="3.8" refY="2.1" '
      'orient="auto" markerUnits="strokeWidth"><path d="M0,0 L5.2,2.1 L0,4.2 Z" '
      'fill="#B7B7B7" stroke="none"/></marker></defs>')
    A('<rect width="100%" height="100%" fill="white"/>')
    A(f'<style>text{{font-family:Arial,Helvetica,sans-serif;fill:#111}} '
      f'.gene{{font-size:{args.font_size}px;paint-order:stroke fill;stroke:#fff;stroke-width:{1.8*k:.2f}px;stroke-linejoin:round}} .legend{{font-size:{args.font_size+1}px}}</style>')

    # ---- GC content ring (v0.9.8: 0-100% radial GC bars + reference circles) ----
    # Radial domain is widened to r=260..200 at 1200 px: 0% GC maps to r=260,
    # 100% GC to r=200. 50% and genome-wide mean GC are drawn behind the bars.
    seq_u = str(seq).upper()
    valid = sum(1 for b in seq_u if b in "ACGT")
    gc_total = sum(1 for b in seq_u if b in "GC")
    mean = gc_total / valid if valid else 0.0
    n_bars = args.gc_bars
    gc_span = gc_outer - gc_inner
    A('<g id="gc_content">')
    A(f'<circle cx="{cx}" cy="{cy}" r="{(gc_outer+gc_inner)/2:.2f}" fill="none" '
      f'stroke="#E4E4E4" stroke-width="{gc_outer-gc_inner:.2f}"/>')
    # Reference radii are calculated here but drawn after the bars so that
    # the black 50% and mean-GC guides remain fully visible in the foreground.
    r_gc50 = gc_outer - 0.5 * gc_span
    r_gcmean = gc_outer - mean * gc_span
    A('<g id="gc_bars">')
    for i in range(n_bars):
        st = (i * L) // n_bars
        en = ((i + 1) * L) // n_bars
        fragment = seq_u[st:en]
        gc_count = fragment.count("G") + fragment.count("C")
        valid_count = sum(fragment.count(b) for b in "ACGT")
        gc_fraction = (gc_count / valid_count) if valid_count else 0.0
        pos = (st + en) / 2.0 if en > st else float(st)
        ang = mapped_angle(pos, L, orient)
        # 0% GC = r=260; 100% GC = r=200 at the default 1200-px canvas.
        r0 = gc_outer - gc_fraction * gc_span
        p0 = polar(cx, cy, r0, ang)
        p1 = polar(cx, cy, gc_outer, ang)
        A(f'<path d="M{p0[0]:.2f},{p0[1]:.2f} L{p1[0]:.2f},{p1[1]:.2f}" '
          f'fill="none" stroke="#B4B4B4" stroke-width="{args.gc_bar_width*k:g}" '
          f'stroke-linecap="butt"/>')
    A('</g>')
    A(f'<circle cx="{cx}" cy="{cy}" r="{gc_outer:.2f}" fill="none" '
      f'stroke="#B4B4B4" stroke-width="{0.6*k:.2f}"/>')
    # Foreground GC reference guides (v0.9.9): black 50% solid, black mean dashed.
    A(f'<circle cx="{cx}" cy="{cy}" r="{r_gc50:.2f}" fill="none" '
      f'stroke="#000" stroke-width="{0.65*k:.2f}"/>')
    A(f'<circle cx="{cx}" cy="{cy}" r="{r_gcmean:.2f}" fill="none" '
      f'stroke="#000" stroke-width="{0.7*k:.2f}" stroke-dasharray="{3*k:.2f},{3*k:.2f}"/>')
    A('</g>')

    # ---- region circle + boundary ticks + tangential region names ----
    # v0.7.4: restore the successful v0.6.4-style region-label orientation.
    # Gene-label placement is unchanged from v0.7.3.
    A(f'<circle cx="{cx}" cy="{cy}" r="{region_r:.2f}" fill="none" stroke="#111" stroke-width="3"/>')
    # v0.9.1: overlay IRA/IRB segments at a configurable width.
    # Zero disables the overlay; the underlying region circle remains unchanged.
    if args.ir_region_width > 0:
        for lab, st, en in regions:
            if lab not in {"IRA", "IRB"}:
                continue
            a0 = mapped_angle(en, L, orient)
            a1 = mapped_angle(st, L, orient)
            p0 = polar(cx, cy, region_r, a0)
            p1 = polar(cx, cy, region_r, a1)
            large = 1 if (a1-a0) % 360 > 180 else 0
            A(f'<path d="M{p0[0]:.2f},{p0[1]:.2f} A{region_r:.2f},{region_r:.2f} 0 {large} 1 {p1[0]:.2f},{p1[1]:.2f}" '
              f'fill="none" stroke="#111" stroke-width="{args.ir_region_width:g}"/>')
    bounds = set()
    A('<g id="regions">')
    for lab, st, en in regions:
        bounds.update((st, en)); mid = mapped_angle((st+en)/2, L, orient)
        x, y = polar(cx, cy, region_r, mid)
        # Tangential orientation: the label baseline follows the region-circle arc.
        # Flip labels on the lower half so the lettering stays upright.
        rot = mid
        if 90 < (mid % 360) < 270:
            rot = mid + 180
        A(f'<rect x="{x-16:.2f}" y="{y-8:.2f}" width="32" height="16" fill="white" '
          f'transform="rotate({rot:.2f} {x:.2f} {y:.2f})"/>')
        A(f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="middle" dominant-baseline="middle" font-size="{10}" '
          f'transform="rotate({rot:.2f} {x:.2f} {y:.2f})">{lab}</text>')
    for bp in bounds:
        ang = mapped_angle(bp, L, orient)
        p0 = polar(cx, cy, region_r-10*k, ang); p1 = polar(cx, cy, region_r+10*k, ang)
        A(f'<path d="M{p0[0]:.2f},{p0[1]:.2f} L{p1[0]:.2f},{p1[1]:.2f}" stroke="#111" stroke-width="2"/>')
    A('</g>')

    # ---- optional genomic position scale (v1.1.2) ----
    # A separate overlay only: no existing map geometry is changed.
    if args.show_scale:
        tick_base_r = region_r + 2*k
        A('<g id="genome_scale" stroke="#555" fill="#000">')
        for bp in range(0, L, 1000):
            if bp == 0 or bp % 10000 == 0:
                tick_len, tick_width = 8*k, 1.0*k
            elif bp % 5000 == 0:
                tick_len, tick_width = 5*k, 0.8*k
            else:
                tick_len, tick_width = 3*k, 0.6*k
            ang = mapped_angle(bp, L, orient)
            p0 = polar(cx, cy, tick_base_r, ang)
            p1 = polar(cx, cy, tick_base_r + tick_len, ang)
            A(f'<path d="M{p0[0]:.2f},{p0[1]:.2f} L{p1[0]:.2f},{p1[1]:.2f}" '
              f'stroke="#555" stroke-width="{tick_width:.2f}"/>')
            if bp and bp % 10000 == 0:
                text_r = tick_base_r + 13*k
                x, y = polar(cx, cy, text_r, ang)
                rot = ang if not (90 < (ang % 360) < 270) else ang + 180
                A(f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="middle" '
                  f'dominant-baseline="middle" font-size="{8*k:.2f}" '
                  f'transform="rotate({rot:.2f} {x:.2f} {y:.2f})">{bp//1000} kb</text>')
        A('</g>')

    # ---- central text ----
    A(f'<text x="{cx}" y="{cy-44*k}" text-anchor="middle" font-size="{29}" font-style="italic">{escape(organism)}</text>')
    A(f'<text x="{cx}" y="{cy-6*k}" text-anchor="middle" font-size="{24}">chloroplast genome</text>')
    A(f'<text x="{cx}" y="{cy+36*k}" text-anchor="middle" font-size="{22}">{L:,} bp</text>')
    A(f'<text x="{cx}" y="{cy+88*k}" text-anchor="middle" font-size="{34}">{escape(args.center_text or rec_id)}</text>')

    # ---- gene blocks around the backbone ----
    counts = Counter(); labels = []
    # v0.8.7: a trans-spliced exon may be shared by two gene definitions
    # (notably the first exon of plastid rps12). Draw/label an identical physical
    # part only once when gene name, coordinates, and strand are all identical.
    # This is deliberately limited to rendered parts; distinct IR copies remain.
    drawn_gene_parts = set()
    A('<g id="genes">')
    for f in features:
        counts[f["category"]] += 1
        part_strands = f.get("part_strands", [f["strand"]] * len(f["parts"]))
        # Blocks remain part-specific, but labels are feature-specific. For each
        # strand represented by a feature, use only its longest rendered part as
        # the label anchor. This prevents one multipart feature from producing a
        # repeated label for every exon/component.
        label_anchor_by_strand = {}
        for (st, en), part_strand in zip(f["parts"], part_strands):
            draw_strand = part_strand if part_strand is not None else f["strand"]
            draw_sign = 1 if draw_strand >= 0 else -1
            part_key = (f["name"].lower(), st, en, draw_sign)
            if part_key in drawn_gene_parts:
                continue
            drawn_gene_parts.add(part_key)
            r0, r1 = ((backbone, backbone+block) if draw_sign > 0
                      else (backbone-block, backbone))
            a0 = mapped_angle(en, L, orient); a1 = mapped_angle(st, L, orient)
            A(f'<path d="{annulus(cx,cy,r0,r1,a0,a1)}" fill="{PALETTE[f["category"]]}" '
              f'stroke="#111" stroke-width="{args.block_outline_width:g}"><title>{escape(f["name"])} {st+1}..{en}</title></path>')
            span = en - st
            old_anchor = label_anchor_by_strand.get(draw_sign)
            if old_anchor is None or span > old_anchor[0]:
                label_anchor_by_strand[draw_sign] = (span, st, en)
        for draw_sign, (span, st, en) in label_anchor_by_strand.items():
            nat = mapped_angle((st+en)/2.0, L, orient)
            labels.append({"f": f, "natural": nat, "angle": nat,
                           "strand": draw_sign, "anchor_span": span})
    A('</g>')
    # v1.0.6: collapse near-identical labels generated by overlapping GenBank
    # feature definitions. Only same-name, same-strand labels within 3 degrees
    # are merged; distant IR copies remain independent. The longest-part anchor
    # is retained, with a typed feature preferred in a tie.
    deduped_labels = []
    label_priority = {"tRNA": 0, "rRNA": 0, "CDS": 1, "gene": 2}
    for z in labels:
        duplicate = None
        for old in deduped_labels:
            if (old["f"]["name"].lower() == z["f"]["name"].lower()
                    and old["strand"] == z["strand"]
                    and ang_sep(old["natural"], z["natural"]) <= 3.0):
                duplicate = old
                break
        if duplicate is None:
            deduped_labels.append(z)
            continue
        old_rank = (-duplicate.get("anchor_span", 0),
                    label_priority.get(duplicate["f"].get("type", "gene"), 9))
        new_rank = (-z.get("anchor_span", 0),
                    label_priority.get(z["f"].get("type", "gene"), 9))
        if new_rank < old_rank:
            deduped_labels[deduped_labels.index(duplicate)] = z
    labels = deduped_labels
    # v0.9.0: draw the backbone after gene blocks so the ring is in front
    # (higher SVG z-order), matching the OGDRAW layer order.
    A(f'<circle cx="{cx}" cy="{cy}" r="{backbone:.2f}" fill="none" stroke="#111" stroke-width="4"/>')
    # v0.9.2: optionally thicken only the IRA/IRB portions of the outer backbone.
    # The complete backbone remains width 4; this overlay changes only IR segments.
    if args.ir_backbone_width > 0:
        for lab, st, en in regions:
            if lab not in {"IRA", "IRB"}:
                continue
            a0 = mapped_angle(en, L, orient)
            a1 = mapped_angle(st, L, orient)
            p0 = polar(cx, cy, backbone, a0)
            p1 = polar(cx, cy, backbone, a1)
            large = 1 if (a1-a0) % 360 > 180 else 0
            A(f'<path d="M{p0[0]:.2f},{p0[1]:.2f} A{backbone:.2f},{backbone:.2f} 0 {large} 1 {p1[0]:.2f},{p1[1]:.2f}" '
              f'fill="none" stroke="#111" stroke-width="{args.ir_backbone_width:g}"/>')

    # ---- place forward labels outside, reverse labels inside (independent engines) ----
    outer = [z for z in labels if z["strand"] > 0]
    inner = [z for z in labels if z["strand"] < 0]
    place_labels(outer, out_text_r, out_displaced_r, args.font_size, cx, cy,
                 args.collision_angle_step, args.max_label_shift, args.collision_padding)
    place_labels(inner, in_text_r, in_displaced_r, args.font_size, cx, cy,
                 args.collision_angle_step, args.max_label_shift, args.collision_padding)
    # Keep gene order and label order identical on each strand. This prevents
    # leader-line crossings without changing natural labels or adding lanes.
    outer_uncrossed = 0
    inner_uncrossed = 0

    # If an existing leader visually passes behind a natural label, promote the
    # obstructing label as well. Then relax all displaced labels once more.
    promoted = promote_labels_obscured_by_leaders(
        labels, cx, cy, k, args.font_size, args.leader_promotion_clearance)
    if promoted:
        for group, dr in ((outer, out_displaced_r), (inner, in_displaced_r)):
            moved=[z for z in group if z.get("leader")]
            fixed=[(z,label_box(z,cx,cy,args.font_size)) for z in group if not z.get("leader")]
            moved.sort(key=lambda z:z["natural"]%360.0)
            accepted=[]
            for z in moved:
                base=z["natural"]%360.0; found=False
                for n in range(0,81):
                    for sign in ((1,) if n==0 else (1,-1)):
                        z["angle"]=(base+sign*n*args.collision_angle_step)%360.0
                        z["radius"]=dr
                        poly=label_box(z,cx,cy,args.font_size)
                        if not any(boxes_overlap(poly,p,args.collision_padding)
                                   for _,p in fixed+accepted):
                            accepted.append((z,poly)); found=True; break
                    if found: break

    # v1.0.1: repair proven order inversions only within local (<=8 deg)
    # clusters on the same strand and displaced lane. Normal v1.0.0 layouts are
    # untouched when their local natural/final orders already agree.
    for z in labels:
        z["order_preserved"]=True; z["order_adjusted"]=False
    outer_uncrossed = repair_local_order_inversions(
        [z for z in outer if z.get("leader") and z.get("lane")==1],
        cluster_gap=8.0, max_extra_shift=args.max_label_shift)
    inner_uncrossed = repair_local_order_inversions(
        [z for z in inner if z.get("leader") and z.get("lane")==1],
        cluster_gap=8.0, max_extra_shift=args.max_label_shift)

    # ---- selective leaders with global fan-out routing ----
    label_diagnostics = []
    label_polygons = [(z, label_box(z, cx, cy, args.font_size)) for z in labels]
    prior_leader_segments, prior_elbows, routes = [], [], {}
    route_order = []
    for strand in (-1, 1):
        side = [z for z in labels if z.get("leader") and z["strand"] == strand]
        side.sort(key=lambda z: z["angle"] % 360.0)
        route_order.extend(side)

    # v1.0.5: jointly route local leader clusters before greedy routing. Labels
    # within 8 degrees on the same strand are tested as a complete straight fan.
    # The fan is accepted only when every segment clears all other labels, all
    # segments clear one another, and the fan clears already accepted clusters.
    def local_route_clusters(side, gap=8.0):
        if not side:
            return []
        ordered=sorted(side,key=lambda z:z["angle"]%360.0)
        # Break the cycle at the largest final-angle gap.
        gaps=[]
        for i,z in enumerate(ordered):
            a=z["angle"]%360.0
            b=ordered[(i+1)%len(ordered)]["angle"]%360.0
            if i==len(ordered)-1: b+=360.0
            gaps.append(b-a)
        cut=(max(range(len(gaps)),key=lambda i:gaps[i])+1)%len(ordered)
        ordered=ordered[cut:]+ordered[:cut]
        vals=[]
        for i,z in enumerate(ordered):
            a=z["angle"]%360.0
            if i:
                while a<=vals[-1]: a+=360.0
            vals.append(a)
        clusters=[]; start=0
        for i in range(1,len(ordered)):
            if vals[i]-vals[i-1] > gap:
                clusters.append(ordered[start:i]); start=i
        clusters.append(ordered[start:])
        return [c for c in clusters if len(c)>=2]

    clustered=set()
    for strand in (-1,1):
        side=[z for z in route_order if z["strand"]==strand]
        for cluster in local_route_clusters(side):
            proposed=[]; ok=True
            for z in cluster:
                if z["strand"]>0: anchor_r,end_r=413*k,451*k
                else: anchor_r,end_r=362*k,324*k
                p0=polar(cx,cy,anchor_r,z["natural"])
                p2=polar(cx,cy,end_r,z["angle"])
                min_lab=1e9
                for other,poly in label_polygons:
                    if id(other)==id(z): continue
                    cl=segment_polygon_clearance(p0,p2,poly); min_lab=min(min_lab,cl)
                    if cl < args.leader_label_clearance*k:
                        ok=False; break
                if not ok: break
                min_line=1e9
                for c,d in prior_leader_segments:
                    cl=segment_segment_clearance(p0,p2,c,d); min_line=min(min_line,cl)
                    if cl < args.leader_line_clearance*k:
                        ok=False; break
                if not ok: break
                proposed.append((z,p0,p2,min_lab,min_line))
            if ok:
                for i in range(len(proposed)):
                    for j in range(i):
                        cl=segment_segment_clearance(proposed[i][1],proposed[i][2],
                                                     proposed[j][1],proposed[j][2])
                        if cl < args.leader_line_clearance*k:
                            ok=False; break
                    if not ok: break
            if not ok:
                continue
            for z,p0,p2,min_lab,min_line in proposed:
                routes[id(z)]=(p0,None,p2)
                z["leader_route_clean"]=True
                z["leader_label_clearance"]=min_lab
                z["leader_route_style"]="cluster-straight"
                prior_leader_segments.append((p0,p2))
                clustered.add(id(z))

    # Labels not accepted as a clean local fan retain v1.0.4 routing.
    for z in route_order:
        if id(z) in clustered:
            continue
        blocked, score, p0, p1, p2, min_lab, min_line = route_leader(
            z, label_polygons, prior_leader_segments, prior_elbows, cx, cy, k,
            args.leader_label_clearance*k, args.leader_line_clearance*k,
            args.leader_elbow_clearance*k)
        z["leader_route_clean"] = not blocked
        z["leader_label_clearance"] = min_lab
        routes[id(z)] = (p0,p1,p2)
        z["leader_route_style"] = "straight" if p1 is None else "elbow"
        if p1 is None:
            prior_leader_segments.append((p0,p2))
        else:
            prior_elbows.append(p1)
            prior_leader_segments.extend([(p0,p1),(p1,p2)])

    # ---- transcription-direction arrows (behind leaders and gene labels) ----
    # v1.0.0 changes only z-order: v0.9.9 arrow geometry/style is preserved.
    o_r = backbone+80*k; i_r = backbone-66*k
    o_ang = 8.0; i_ang = o_ang*o_r/i_r
    op0 = polar(cx,cy,o_r,0); op1 = polar(cx,cy,o_r,-o_ang)
    ip0 = polar(cx,cy,i_r,0); ip1 = polar(cx,cy,i_r,i_ang)
    A('<g id="transcription_arrows">')
    A(f'<path d="M{op0[0]:.2f},{op0[1]:.2f} A{o_r:.2f},{o_r:.2f} 0 0 0 {op1[0]:.2f},{op1[1]:.2f}" '
      f'fill="none" stroke="#B7B7B7" stroke-width="3.3" stroke-linecap="butt" stroke-linejoin="miter" marker-end="url(#greyArrow)"/>')
    A(f'<path d="M{ip0[0]:.2f},{ip0[1]:.2f} A{i_r:.2f},{i_r:.2f} 0 0 1 {ip1[0]:.2f},{ip1[1]:.2f}" '
      f'fill="none" stroke="#B7B7B7" stroke-width="3.3" stroke-linecap="butt" stroke-linejoin="miter" marker-end="url(#greyArrow)"/>')
    A('</g>')

    A(f'<g id="leaders" stroke="#303030" stroke-width="{leader_w:.3f}" fill="none" '
      'stroke-linecap="round" stroke-linejoin="round">')
    for z in labels:
        z.setdefault("leader_route_clean", True)
        z.setdefault("leader_label_clearance", None)
        if z.get("leader",False):
            p0,p1,p2 = routes[id(z)]
            if p1 is None:
                A(f'<path d="M{p0[0]:.2f},{p0[1]:.2f} L{p2[0]:.2f},{p2[1]:.2f}"/>')
            else:
                A(f'<path d="M{p0[0]:.2f},{p0[1]:.2f} L{p1[0]:.2f},{p1[1]:.2f} '
                  f'L{p2[0]:.2f},{p2[1]:.2f}"/>')
        label_diagnostics.append({
            "gene": z["f"]["name"], "strand": "outer" if z["strand"] > 0 else "inner",
            "natural_angle": z["natural"] % 360, "final_angle": z["angle"] % 360,
            "angular_shift": z.get("angular_shift",0.0), "lane": z.get("lane",0),
            "radius": z["radius"], "leader": z.get("leader",False),
            "reason": z.get("leader_reason","natural"),
            "leader_route_clean": z.get("leader_route_clean",True),
            "leader_label_clearance": z.get("leader_label_clearance"),
            "leader_route_style": z.get("leader_route_style", "none"),
            "order_preserved": z.get("order_preserved",True),
            "order_adjusted": z.get("order_adjusted",False),
        })
    A('</g>')

    # ---- label text (forward reads outward, reverse reads inward) ----
    A('<g id="labels">')
    for z in labels:
        ang = z["angle"] % 360; x, y = polar(cx, cy, z["radius"], ang); rot = ang-90
        if z["strand"] > 0:
            anchor = "start"
            if 180 < ang < 360: rot += 180; anchor = "end"
        else:
            anchor = "end"
            if 180 < ang < 360: rot += 180; anchor = "start"
        A(f'<text class="gene" x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" '
          f'dominant-baseline="middle" transform="rotate({rot:.2f} {x:.2f} {y:.2f})">'
          f'{escape(z["f"]["label"])}</text>')
    A('</g>')


    # ---- legend (lower-left, vertical, only used categories) ----
    A('<g id="legend">'); x = 42*k; y = 902*k
    for name in LEGEND_ORDER:
        A(f'<rect x="{x:.2f}" y="{y-9:.2f}" width="14" height="11" fill="{PALETTE[name]}" '
          f'stroke="#111" stroke-width=".5"/><text class="legend" x="{x+17:.2f}" y="{y:.2f}">{escape(name)}</text>')
        y += 16*k
    A('</g></svg>')
    return "\n".join(out), mean, label_diagnostics

# ---------- CLI ----------
def main():
    examples = f"""
Examples:
  python cpdraw.py
      Read sequence.gb and create sequence.svg plus sequence_label_report.csv.

  python cpdraw.py plastome.gb
      Draw plastome.gb using the default settings.

  python cpdraw.py plastome.gb --block-outline-width 0
      Draw gene blocks without outlines.

  python cpdraw.py plastome.gb --show-scale
      Add a genomic position ruler with 1-kb ticks and 10-kb labels.

  python cpdraw.py plastome.gb --label-report
      Also create plastome_label_report.csv.

  python cpdraw.py plastome.gb --label-report diagnostics.csv
      Write the diagnostic CSV as diagnostics.csv.

Standard output:
  SVG map (*.svg)

Optional output:
  Label-placement diagnostic report CSV
    --label-report          use <input>_label_report.csv
    --label-report FILE     use the specified file name

Requirement:
  Biopython

Use --help to list all public options.
"""
    ap = argparse.ArgumentParser(
        description=f"{PROGRAM_NAME} {VERSION}: OGDRAW-like circular plastome map generator",
        epilog=examples,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genbank", nargs="?", default="sequence.gb",
                    help="input GenBank file (default: sequence.gb)")
    ap.add_argument("-o", "--output")
    ap.add_argument("--size", type=int, default=1200)
    ap.add_argument("--font-size", type=float, default=10.0)
    ap.add_argument("--block-outline-width", type=float, default=0.25, help="gene-block outline width (default: 0.25; 0 disables outlines)")
    ap.add_argument("--ir-region-width", type=float, default=3.0,
                    help="IRA/IRB segment width on the inner region ring (default: 3.0; 0 disables IR overlay)")
    ap.add_argument("--ir-backbone-width", type=float, default=6.0,
                    help="IRA/IRB segment width on the outer gene-block backbone (default: 6.0; 0 disables IR thickening)")
    ap.add_argument("--show-pseudogene", action="store_true",
                    help="prefix pseudogene labels with 'psi' (default: hidden, OGDRAW-like)")
    ap.add_argument("--show-scale", action="store_true",
                    help="show genomic ruler: 1-kb ticks, 5-kb medium ticks, 10-kb labeled ticks")
    ap.add_argument("--radial-nudge", type=float, default=0.0, help=argparse.SUPPRESS)
    ap.add_argument("--leader-width", type=float, default=0.25, help="leader stroke width (OGDRAW default: 0.25 units @1200px)")
    ap.add_argument("--leader-label-clearance", type=float, default=4.2, help="minimum clearance between a leader and another label")
    ap.add_argument("--leader-line-clearance", type=float, default=2.0, help="minimum clearance between separate leader lines")
    ap.add_argument("--leader-promotion-clearance", type=float, default=6.0, help="promote a natural label when another leader passes within this visual corridor")
    ap.add_argument("--leader-elbow-clearance", type=float, default=5.0, help="minimum separation between leader elbows; prevents ambiguous shared junctions")
    ap.add_argument("--leader-angle-threshold", type=float, default=999.0, help=argparse.SUPPRESS)
    ap.add_argument(
        "--label-report", nargs="?", const="__AUTO__", default=None, metavar="FILE",
        help=("optionally write a diagnostic CSV; if FILE is omitted, use "
              "<input>_label_report.csv (default: no CSV)"))
    ap.add_argument("--collision-angle-step", type=float, default=0.18, help="tangential search step in degrees for colliding labels")
    ap.add_argument("--max-label-shift", type=float, default=7.0, help="maximum tangential shift before radial fallback")
    ap.add_argument("--collision-radial-offset", type=float, default=39.0, help=argparse.SUPPRESS)
    ap.add_argument("--collision-padding", type=float, default=1.2, help="padding used by rotated-label collision detection")
    ap.add_argument("--gc-bars", type=int, default=3600,
                    help="number of OGDRAW-like radial GC bars (default: 3600)")
    ap.add_argument("--gc-percent-span", type=float, default=60.0, help=argparse.SUPPRESS)
    ap.add_argument("--gc-bar-width", type=float, default=1.0,
                    help="GC bar stroke width at 1200px (default: 1.0)")
    ap.add_argument("--gc-window", type=int, default=500, help=argparse.SUPPRESS)
    ap.add_argument("--gc-step", type=int, default=50, help=argparse.SUPPRESS)
    ap.add_argument("--trans-gap", type=int, default=10000)
    ap.add_argument("--center-text")
    args = ap.parse_args()
    if args.block_outline_width < 0:
        ap.error("--block-outline-width must be 0 or greater")
    if args.ir_region_width < 0:
        ap.error("--ir-region-width must be 0 or greater")
    if args.ir_backbone_width < 0:
        ap.error("--ir-backbone-width must be 0 or greater")
    if args.gc_bars < 1:
        ap.error("--gc-bars must be 1 or greater")
    if args.gc_percent_span < 0 or args.gc_bar_width < 0:
        ap.error("--gc-percent-span and --gc-bar-width must be 0 or greater")

    input_path = Path(args.genbank)
    if not input_path.is_file():
        ap.error(f"GenBank input file not found: {args.genbank}")
    try:
        from Bio import SeqIO
    except ImportError:
        print(
            "ERROR: Biopython is required to read GenBank files.\n"
            f"Install it with one of the following commands and run {PROGRAM_NAME} again:\n"
            "  python -m pip install biopython\n"
            "  py -m pip install biopython    (Windows)",
            file=sys.stderr)
        raise SystemExit(1)
    rec = SeqIO.read(str(input_path), "genbank")
    feats = collect_features(rec, args.trans_gap, args.show_pseudogene)
    regs = collect_regions(rec)
    organism = (rec.annotations.get("organism")
                or (rec.description.split(" chloroplast")[0] if rec.description else rec.id))
    svg, mean, label_diagnostics = make_svg(feats, regs, str(rec.seq), rec.id, organism, args)
    base = Path(args.output) if args.output else Path(args.genbank).with_suffix("")
    if base.suffix.lower() in {".gb", ".gbk"}: base = base.with_suffix("")
    svg_path = base.with_suffix(".svg"); svg_path.write_text(svg, encoding="utf-8")
    created = [svg_path]
    # v1.1.4: diagnostic CSV is opt-in. When requested, every data row records
    # the generating program and exact version for provenance.
    if args.label_report:
        report_path = (base.with_name(base.name + "_label_report.csv")
                       if args.label_report == "__AUTO__"
                       else Path(args.label_report))
        import csv
        with report_path.open("w", newline="", encoding="utf-8-sig") as fh:
            fields = ["program", "version", "legacy_program", "legacy_version",
                      "gene", "strand", "natural_angle", "final_angle",
                      "angular_shift", "lane", "radius", "leader",
                      "reason", "leader_route_clean", "leader_label_clearance",
                      "leader_route_style", "order_preserved", "order_adjusted"]
            writer = csv.DictWriter(fh, fieldnames=fields)
            writer.writeheader()
            for row in label_diagnostics:
                cooked = dict(row)
                cooked["program"] = PROGRAM_NAME
                cooked["version"] = VERSION
                cooked["legacy_program"] = LEGACY_PROGRAM
                cooked["legacy_version"] = LEGACY_VERSION
                for key in ("natural_angle", "final_angle", "angular_shift", "radius"):
                    cooked[key] = f"{cooked[key]:.4f}"
                if cooked.get("leader_label_clearance") is not None:
                    cooked["leader_label_clearance"] = f"{cooked['leader_label_clearance']:.4f}"
                writer.writerow(cooked)
        created.append(report_path)
    natural_count = sum(1 for row in label_diagnostics if not row["leader"])
    leader_count = len(label_diagnostics) - natural_count
    order_adjusted_count=sum(1 for row in label_diagnostics if row.get("order_adjusted"))
    print(f"Label analysis: natural/no leader={natural_count}; collision-adjusted/leader={leader_count}; order-adjusted={order_adjusted_count}; threshold={args.leader_angle_threshold:.2f} deg")
    print(f"{PROGRAM_NAME} {VERSION}: {rec.id}, {len(rec.seq):,} bp, GC {mean*100:.2f}%")
    print("Label geometry @1200px: natural outer/inner anchors=414/360; displaced anchors=453/322; leaders=25+13 units")
    print(f"GC bars={args.gc_bars}; 0-100% domain=r260..r200; 50% solid; mean GC={mean*100:.2f}% dashed; bar width={args.gc_bar_width:g}")
    print(f"Features: {len(feats)}; regions: {', '.join(r[0] for r in regs) or 'none'}")
    print("Created: " + ", ".join(str(p) for p in created))

if __name__ == "__main__":
    main()
