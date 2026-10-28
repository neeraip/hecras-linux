"""RASMapper's Delaunay TIN, ported (RasMapperLib DelaunayTINFactory,
TINFactory, DevelopmentTin, DevelopmentTriangle, CircleM, PlaneM,
RasterizeTriangle.ComputeCellsZ, HEC-RAS 7.0.2).

Only the path RASMapper takes for terrain stitches is ported: a TIN built
from stitch segments alone (no perimeter, breaklines or mass points), whose
triangles bridge narrow gaps in the terrain and are painted onto a raster.
The algorithm is followed step for step, including the order points and
triangles are visited in, because ties (four points on one circle, which grid
points nearly always are) are decided by that order.
"""
from __future__ import annotations

import math
import struct

SUPER_POINT_COUNT = 8

# DevelopmentTriangle flags
E0_BREAK, E1_BREAK, E2_BREAK = 1, 2, 4
E0_PERIM, E1_PERIM, E2_PERIM = 8, 16, 32
DELETED = 64


def _div(a: float, b: float) -> float:
    """IEEE division as .NET does it: x/0 is ±inf, 0/0 is NaN."""
    try:
        return a / b
    except ZeroDivisionError:
        if a == 0 or math.isnan(a):
            return math.nan
        return math.copysign(math.inf, a) * math.copysign(1.0, b)


class Tri:
    __slots__ = ("p0", "p1", "p2", "at0", "at1", "at2", "flags")

    def __init__(self, p0, p1, p2, at0, at1, at2):
        self.p0, self.p1, self.p2 = p0, p1, p2
        self.at0, self.at1, self.at2 = at0, at1, at2
        self.flags = 0

    def deleted(self):
        return self.flags & DELETED != 0

    def mark_deleted(self):
        self.flags |= DELETED

    def set_edge(self, e, perimeter):
        self.flags |= ((E0_PERIM, E1_PERIM, E2_PERIM) if perimeter else (E0_BREAK, E1_BREAK, E2_BREAK))[min(e, 2)]

    def edge_bl_or_perim(self, e):
        return self.flags & ((E0_BREAK | E0_PERIM), (E1_BREAK | E1_PERIM), (E2_BREAK | E2_PERIM))[e] != 0

    def change_at(self, src, dst):
        if self.at0 == src:
            self.at0 = dst
        elif self.at1 == src:
            self.at1 = dst
        elif self.at2 == src:
            self.at2 = dst
        else:
            raise ValueError(f"triangle {src} is not adjacent")

    def cw(self, pt):
        if self.p0 == pt:
            return self.at0, self.edge_bl_or_perim(0)
        if self.p1 == pt:
            return self.at1, self.edge_bl_or_perim(1)
        if self.p2 == pt:
            return self.at2, self.edge_bl_or_perim(2)
        return -1, False

    def ccw(self, pt):
        if self.p0 == pt:
            return self.at2, self.edge_bl_or_perim(2)
        if self.p1 == pt:
            return self.at0, self.edge_bl_or_perim(0)
        if self.p2 == pt:
            return self.at1, self.edge_bl_or_perim(1)
        return -1, False


def _close_to(a, bx, by, d):
    """PointM.CloseTo."""
    dx = abs(a[0] - bx)
    if dx > d:
        return False
    dy = abs(a[1] - by)
    if dy > d:
        return False
    return not (math.sqrt(dx ** 2 + dy ** 2) > d)


def _dist(a, bx, by):
    return math.sqrt((a[0] - bx) ** 2 + (a[1] - by) ** 2)


def _cross(a, b, c):
    """TINFactory.CrossProduct."""
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _circumcircle(a, b, c):
    """CircleM.Compute."""
    n1, n2 = b[0] - a[0], b[1] - a[1]
    n3 = n1 * n1 + n2 * n2
    n4, n5 = c[0] - a[0], c[1] - a[1]
    n6 = n4 * n4 + n5 * n5
    n7 = 2.0 * (n1 * n5 - n2 * n4)
    if n7 == 0.0:
        raise ArithmeticError("Colinear points")
    x = a[0] + (n5 * n3 - n2 * n6) / n7
    y = a[1] + (n1 * n6 - n4 * n3) / n7
    return x, y, _dist(a, x, y)


def _seg_intersect(a, b, c, d, eps=1e-6):
    """SegmentM.Intersects of a-b with c-d: (hit, x, y, t1, t2)."""
    if not (min(c[0], d[0]) <= max(a[0], b[0]) + eps and max(c[0], d[0]) >= min(a[0], b[0]) - eps
            and min(c[1], d[1]) <= max(a[1], b[1]) + eps and max(c[1], d[1]) >= min(a[1], b[1]) - eps):
        return False, None, None, None, None     # outputs left as they were
    dx, dy = b[0] - a[0], b[1] - a[1]
    ex, ey = d[0] - c[0], d[1] - c[1]
    den = ey * dx - ex * dy
    if den == 0.0:
        return False, None, None, None, None
    n6, n7 = a[1] - c[1], a[0] - c[0]
    t1 = (ex * n6 - ey * n7) / den
    t2 = (dx * n6 - dy * n7) / den
    x, y = a[0] + t1 * dx, a[1] + t1 * dy
    if 0.0 <= t1 <= 1.0 and 0.0 <= t2 <= 1.0:
        return True, x, y, t1, t2
    if -eps < t1 < 1.0 + eps and -eps < t2 < 1.0 + eps:
        return True, x, y, min(max(t1, 0.0), 1.0), min(max(t2, 0.0), 1.0)
    return False, x, y, t1, t2


def _point_at(a, b, t):
    """SegmentM.PointAtFraction (x, y, z, m)."""
    return (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]),
            a[2] + (b[2] - a[2]) * t, a[3] + (b[3] - a[3]) * t)


def _dist_to_segment(a, b, p):
    """SegmentM.DistanceToPoint (distance only)."""
    length = math.sqrt((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2)
    if length == 0.0:
        return _dist(a, p[0], p[1])
    t = ((b[0] - a[0]) * (p[0] - a[0]) + (b[1] - a[1]) * (p[1] - a[1])) / length
    if t < 0.0:
        return _dist(a, p[0], p[1])
    if t > length:
        return _dist(b, p[0], p[1])
    d = _dist(a, p[0], p[1])
    r = d * d - t * t
    return math.sqrt(r) if r > 0.0 else 0.0


def segment_hash(a) -> int:
    """.NET's default GetHashCode of a SegmentM: the hash of its first field,
    point A, which (four doubles, no padding) XORs its eight 32-bit words."""
    h = 0
    for w in struct.unpack("<8i", struct.pack("<4d", *a)):
        h ^= w
    return h


class StitchTin:
    """DelaunayTINFactory.ComputeTIN(null, null, null, StitchCache) with
    RemoveDeletedTriangles: the triangles left, as (x, y, z) corner triples."""

    def __init__(self, stitches, tolerance=0.01):
        self.tol = tolerance
        self.stitches = [tuple(tuple(float(v) for v in p) for p in sg) for sg in stitches]  # [((x,y,z,m), (x,y,z,m))]
        self.points = []                          # [(x, y, z, m)]
        self.tris = []
        self.free = []                            # _deletedTrianglesStack
        self.last_t = 0                           # FindTriangleContainingPoint's static t

    # -- TINFactory ------------------------------------------------------------
    def add_triangle(self, p0, p1, p2, at0, at1, at2):
        t = Tri(p0, p1, p2, at0, at1, at2)
        if self.free:
            n = self.free.pop()
            self.tris[n] = t
        else:
            n = len(self.tris)
            self.tris.append(t)
        return n

    def find_triangle(self, pt):
        t = self.last_t
        if t < 0 or t > len(self.tris) - 1:
            t = 0
        prev, count = -1, 0
        P = self.points
        while True:
            tr = self.tris[t]
            c0 = _cross(P[tr.p0], P[tr.p1], pt)
            c1 = _cross(P[tr.p1], P[tr.p2], pt)
            c2 = _cross(P[tr.p2], P[tr.p0], pt)
            k = (0 if c0 < c2 else 2) if c0 < c1 else (1 if c1 < c2 else 2)   # Math.Min3Index
            if (c0, c1, c2)[k] >= 0.0:
                self.last_t = t
                return t
            nxt = (tr.at0, tr.at1, tr.at2)[k]
            if prev == nxt:
                self.last_t = t
                return t
            prev, t = t, nxt
            count += 1
            if t == -1 or count == len(self.tris):
                self.last_t = t
                return -1

    def _closest(self, pt, p0, p1, p2):
        P = self.points
        d0, d1, d2 = _dist(P[p0], pt[0], pt[1]), _dist(P[p1], pt[0], pt[1]), _dist(P[p2], pt[0], pt[1])
        if d0 < d1:
            return p0 if d0 < d2 else p2
        return p1 if d1 < d2 else p2

    def find_triangle_and_point(self, pt):
        t = self.find_triangle(pt)
        if t == -1:
            return t, -1
        tr = self.tris[t]
        P = self.points
        if any(_close_to(pt, P[q][0], P[q][1], self.tol) for q in (tr.p0, tr.p1, tr.p2)):
            return t, self._closest(pt, tr.p0, tr.p1, tr.p2)
        return t, -1

    def find_triangle_and_nearest(self, pt):
        t = self.find_triangle(pt)
        if t == -1:
            return t, -1
        tr = self.tris[t]
        return t, self._closest(pt, tr.p0, tr.p1, tr.p2)

    def update_zm(self, p, pt):
        x, y, z, m = self.points[p]
        if not math.isnan(pt[2]):
            z = pt[2]
        if not math.isnan(pt[3]):
            m = pt[3]
        self.points[p] = (x, y, z, m)

    def set_edge_type(self, t, a, b, perimeter):
        tr = self.tris[t]
        nxt = -1
        for e, (u, v) in enumerate(((tr.p0, tr.p1), (tr.p1, tr.p2), (tr.p2, tr.p0))):
            if (u == a and v == b) or (u == b and v == a):
                tr.set_edge(e, perimeter)
                nxt = (tr.at0, tr.at1, tr.at2)[e]
                break
        if nxt != -1:
            tr = self.tris[nxt]
            for e, (u, v) in enumerate(((tr.p0, tr.p1), (tr.p1, tr.p2), (tr.p2, tr.p0))):
                if (u == a and v == b) or (u == b and v == a):
                    tr.set_edge(e, perimeter)
                    break

    def edge(self, t, e):
        tr = self.tris[t]
        P = self.points
        return ((P[tr.p0], P[tr.p1]), (P[tr.p1], P[tr.p2]), (P[tr.p2], P[tr.p0]))[e]

    def edge_first(self, t, e):
        tr = self.tris[t]
        return (tr.p0, tr.p1, tr.p2)[e]

    def edge_last(self, t, e):
        tr = self.tris[t]
        return (tr.p1, tr.p2, tr.p0)[e]

    # -- DelaunayTINFactory ---------------------------------------------------------
    def compute(self):
        segs = self.stitches
        ext = None
        for a, b in segs:                         # ComputeMaximalExtent
            lo_x, hi_x = min(a[0], b[0]), max(a[0], b[0])
            lo_y, hi_y = min(a[1], b[1]), max(a[1], b[1])
            if ext is None:
                ext = [lo_x, hi_x, lo_y, hi_y]
            elif not (ext[0] <= a[0] <= ext[1] and ext[2] <= a[1] <= ext[3]
                      and ext[0] <= b[0] <= ext[1] and ext[2] <= b[1] <= ext[3]):
                ext = [min(ext[0], lo_x), max(ext[1], hi_x), min(ext[2], lo_y), max(ext[3], hi_y)]
        ext = _enforce_min_dims(ext, self.tol * 10.0, self.tol * 10.0)
        dx, dy = ext[1] - ext[0], ext[3] - ext[2]
        if dx > dy:
            if dy < dx * 0.125:
                ext = _enforce_min_dims(ext, 0.0, dx * 0.125)
        elif dx < dy * 0.125:
            ext = _enforce_min_dims(ext, dy * 0.125, 0.0)
        self.offset = (-ext[0], -ext[2])
        self._add_octagon(ext)
        for a, b in segs:                         # InsertStitchSegments
            self.insert_point(self._translate(a))
            self.insert_point(self._translate(b))
        for a, b in segs:
            self.enforce_edge(self._translate(a), self._translate(b))
        self._mark_bad()
        ox, oy = self.offset
        self.points = [(x - ox, y - oy, z, m) for x, y, z, m in self.points]
        out = []                                  # RemoveDeletedTriangles keeps order
        for tr in self.tris:
            if not tr.deleted():
                out.append(tuple((float(self.points[p][0]), float(self.points[p][1]), float(self.points[p][2]))
                                 for p in (tr.p0, tr.p1, tr.p2)))
        return out

    def _translate(self, p):
        return (p[0] + self.offset[0], p[1] + self.offset[1], p[2], p[3])

    def _add_octagon(self, ext):
        ox, oy = self.offset
        min_x, max_x, min_y, max_y = ext[0] + ox, ext[1] + ox, ext[2] + oy, ext[3] + oy
        n = max(max_x - min_x, max_y - min_y) * 0.05
        mid_x, mid_y = (min_x + max_x) * 0.5, (max_y + min_y) * 0.5
        nan = math.nan
        for x, y in ((mid_x, max_y + n * 2.0), (max_x + n, max_y + n), (max_x + n * 2.0, mid_y),
                     (max_x + n, min_y - n), (mid_x, min_y + (-n) * 2.0), (min_x - n, min_y - n),
                     (min_x + (-n) * 2.0, mid_y), (min_x - n, max_y + n)):
            self.points.append((x, y, nan, nan))
        if (ext[1] - ext[0]) > (ext[3] - ext[2]):
            spec = [(0, 7, 6, -1, -1, 4), (6, 5, 4, -1, -1, 4), (4, 3, 2, -1, -1, 5),
                    (2, 1, 0, -1, -1, 5), (0, 6, 4, 0, 1, 5), (0, 4, 2, 4, 2, 3)]
        else:
            spec = [(0, 7, 6, -1, -1, 4), (6, 5, 4, -1, -1, 5), (4, 3, 2, -1, -1, 5),
                    (2, 1, 0, -1, -1, 4), (0, 6, 2, 0, 5, 3), (6, 4, 2, 1, 2, 4)]
        for s in spec:
            self.add_triangle(*s)

    def insert_point(self, pt):
        """DelaunayInsertPoint (point already translated)."""
        t, p = self.find_triangle_and_point(pt)
        if p >= 0:
            self.update_zm(p, pt)
            return False
        p = len(self.points)
        self.points.append(pt)
        stack = [t]
        while stack:
            t = stack.pop()
            self.free.append(t)
            tr = self.tris[t]
            tr.mark_deleted()
            for at in (tr.at0, tr.at1, tr.at2):
                if self._inside_hole(at, pt):
                    stack.append(at)
        self._fill_hole(p)
        return True

    def _inside_hole(self, at, pt):
        if at == -1:
            return False
        tr = self.tris[at]
        if tr.deleted():
            return False
        P = self.points
        a, b, c = P[tr.p0], P[tr.p1], P[tr.p2]
        x, y, r = _circumcircle(a, b, c)
        if not _close_to(pt, x, y, r):
            return False
        for u, v in ((a, b), (a, c), (b, c)):
            if u[0] == v[0] == pt[0] and ((pt[1] > u[1] and pt[1] > v[1]) or (pt[1] < u[1] and pt[1] < v[1])):
                return False
        for u, v in ((a, b), (a, c), (b, c)):
            if u[1] == v[1] == pt[1] and ((pt[0] > u[0] and pt[0] > v[0]) or (pt[0] < u[0] and pt[0] < v[0])):
                return False
        return True

    def _fill_hole(self, p):
        """MakeInsertionHoleFromDeletedTrianglesAndFill."""
        edges = []
        for t in reversed(self.free):             # foreach over a Stack: newest first
            tr = self.tris[t]
            for (u, v, at) in ((tr.p0, tr.p1, tr.at0), (tr.p1, tr.p2, tr.at1), (tr.p2, tr.p0, tr.at2)):
                if at == -1 or not self.tris[at].deleted():
                    edges.append([u, v, at, -1])
        around = {}
        for edge in edges:
            u, v, outside, _ = edge
            n = self.add_triangle(u, v, p, outside, -1, -1)
            edge[3] = n
            around.setdefault(u, [0, 0])[1] = n   # [triangleLeft, triangleRight]
            around.setdefault(v, [0, 0])[0] = n
            if outside >= 0:
                o = self.tris[outside]
                bl = perim = False
                if o.p0 == v:
                    o.at0 = n
                    bl, perim = o.flags & E0_BREAK != 0, o.flags & E0_PERIM != 0
                elif o.p1 == v:
                    o.at1 = n
                    bl, perim = o.flags & E1_BREAK != 0, o.flags & E1_PERIM != 0
                elif o.p2 == v:
                    o.at2 = n
                    bl, perim = o.flags & E2_BREAK != 0, o.flags & E2_PERIM != 0
                if bl or perim:
                    self.tris[n].set_edge(0, perim)
        for u, v, _, n in edges:
            tr = self.tris[n]
            tr.at1 = around[v][1]
            tr.at2 = around[u][0]
        return edges

    def enforce_edge(self, a, b, perimeter=False):
        """EnforceEdgeSegment (segment already translated)."""
        P = self.points
        _, p = self.find_triangle_and_nearest(b)
        b = P[p]
        t3 = 0.0                                  # kept across iterations, as in RASMapper
        while True:
            t2, p2 = self.find_triangle_and_nearest(a)
            a = self.points[p2]
            while True:
                tr = self.tris[t2]
                if p2 == tr.p0:
                    e, n_cw, n_ccw = 1, tr.at0, tr.at2
                elif p2 == tr.p1:
                    e, n_cw, n_ccw = 2, tr.at1, tr.at0
                else:
                    e, n_cw, n_ccw = 0, tr.at2, tr.at1
                sa, sb = self.edge(t2, e)
                if _cross(self.points[p2], sa, b) < 0.0:
                    t2 = n_cw
                    continue
                if not (_cross(self.points[p2], sb, b) > 0.0):
                    break
                t2 = n_ccw
            if p in (tr.p0, tr.p1, tr.p2):
                self.set_edge_type(t2, p2, p, perimeter)
                return
            _, _, _, t, _ = _seg_intersect(a, b, sa, sb, 0.0)
            if t is not None:
                t3 = t
            pm = _point_at(a, b, t3)
            first, last = self.edge_first(t2, e), self.edge_last(t2, e)
            if _close_to(pm, sa[0], sa[1], self.tol):
                self.set_edge_type(t2, p2, first, perimeter)
                p2 = first
                self.update_zm(p2, pm)
                a = self.points[p2]
                continue
            if _close_to(pm, sb[0], sb[1], self.tol):
                self.set_edge_type(t2, p2, last, perimeter)
                p2 = last
                self.update_zm(p2, pm)
                a = self.points[p2]
                continue
            count = len(self.points)
            self.points.append(pm)
            adj = (tr.at0, tr.at1, tr.at2)[e]
            bl = tr.flags & (E0_BREAK, E1_BREAK, E2_BREAK)[e] != 0
            perim = tr.flags & (E0_PERIM, E1_PERIM, E2_PERIM)[e] != 0
            self.free.append(t2)
            self.free.append(adj)
            self.tris[t2].mark_deleted()
            self.tris[adj].mark_deleted()
            edges = self._fill_hole(count)
            if bl or perim:
                for edge in edges:
                    self.set_edge_type(edge[3], count, first, perim)
                    self.set_edge_type(edge[3], count, last, perim)
            self.set_edge_type(len(self.tris) - 1, p2, count, perimeter)
            a = pm

    def _is_super(self, t):
        tr = self.tris[t]
        return tr.p0 < SUPER_POINT_COUNT or tr.p1 < SUPER_POINT_COUNT or tr.p2 < SUPER_POINT_COUNT

    def _mark_bad(self):
        """MarkBadTrianglesAsDeleted_SinglePart, with stitches and no perimeter."""
        for tr in self.tris:
            if tr.p0 < SUPER_POINT_COUNT or tr.p1 < SUPER_POINT_COUNT or tr.p2 < SUPER_POINT_COUNT:
                tr.mark_deleted()
        m = 0.0
        for a, _ in self.stitches:
            if m < a[3]:
                m = a[3]
        limit = m * 2.0001
        P = self.points
        for tr in self.tris:
            if tr.deleted():
                continue
            a, b, c = P[tr.p0], P[tr.p1], P[tr.p2]
            if (abs(a[0] - b[0]) > limit or abs(a[1] - b[1]) > limit or abs(b[0] - c[0]) > limit
                    or abs(b[1] - c[1]) > limit or abs(c[0] - a[0]) > limit or abs(c[1] - a[1]) > limit):
                tr.mark_deleted()
        for a, b in self.stitches:
            self._delete_left_of(self._translate(a), self._translate(b))

    def _delete_left_of(self, a, b):
        """DeleteTrianglesLeftOfStitchEdge (segment already translated)."""
        t, p = self.find_triangle_and_point(a)
        if p == -1:
            p = self.tris[t].p1
        seen = set()
        while True:
            while True:
                tr = self.tris[t]
                if p == tr.p0:
                    e, n_cw, n_ccw = 1, tr.at0, tr.at2
                elif p == tr.p1:
                    e, n_cw, n_ccw = 2, tr.at1, tr.at0
                else:
                    e, n_cw, n_ccw = 0, tr.at2, tr.at1
                sa, sb = self.edge(t, e)
                if _cross(self.points[p], sa, b) < 0.0:
                    t = n_cw
                    continue
                if not (_cross(self.points[p], sb, b) > 0.0):
                    break
                t = n_ccw
            if t in seen:
                return
            seen.add(t)
            if _dist_to_segment(a, b, sa) < _dist_to_segment(a, b, sb):
                p = self.edge_first(t, e)
                e = (2, 0, 1)[e]
                n = (self.tris[t].at0, self.tris[t].at1, self.tris[t].at2)[e]
                if n == -1:
                    raise IndexError("stitch edge runs off the TIN")
                self.tris[n].mark_deleted()
                self._rotate_delete(n)
            else:
                p = self.edge_last(t, e)
                self.tris[t].mark_deleted()
                self._rotate_delete(t)
            if _close_to(self.points[p], b[0], b[1], self.tol):
                return

    def _rotate_delete(self, t0):
        tr = self.tris[t0]
        for pt in (tr.p0, tr.p1, tr.p2):
            cur = self.tris[t0]
            cur.mark_deleted()
            while True:
                nxt, stop = cur.cw(pt)
                if stop or nxt == -1 or nxt == t0 or self._is_super(nxt):
                    break
                cur = self.tris[nxt]
                cur.mark_deleted()
            while True:
                nxt, stop = cur.ccw(pt)
                if stop or nxt == -1 or nxt == t0 or self._is_super(nxt):
                    break
                cur = self.tris[nxt]
                cur.mark_deleted()


def _enforce_min_dims(ext, width, height):
    """Extent.EnforceMinDimensions on [min_x, max_x, min_y, max_y]."""
    min_x, max_x, min_y, max_y = ext
    if width > 0.0 and (max_x - min_x) < width:
        cx = (min_x + max_x) * 0.5
        min_x, max_x = cx - width / 2.0, cx + width / 2.0
    if height > 0.0 and (max_y - min_y) < height:
        cy = (min_y + max_y) * 0.5
        min_y, max_y = cy - height / 2.0, cy + height / 2.0
    return [min_x, max_x, min_y, max_y]


def paint(tris, values, rows, cols, min_x, max_y, cs, f32):
    """DevelopmentTin.Rasterize → RasterizeTriangle.ComputeCellsZ: each
    triangle's plane onto the raster cells whose centre it covers, in order."""
    if not tris:
        return
    min_x, max_y, cs = float(min_x), float(max_y), float(cs)
    max_x, min_y = min_x + cols * cs, max_y - rows * cs
    xs = [p[0] for t in tris for p in t]
    ys = [p[1] for t in tris for p in t]
    if min(xs) > max_x or max(xs) < min_x or min(ys) > max_y or max(ys) < min_y:
        return                                    # raster.Overlaps(_extent)
    def _index(v):                                # RasterM.Row / Col
        if v >= 2147483647.0:
            return 2147483647
        if v <= -2147483648.0:
            return -2147483648
        return int(v)
    row = lambda y: _index(math.floor((max_y - y) / cs))
    col = lambda x: _index(math.floor((x - min_x) / cs))
    cy = lambda r: max_y - cs * (r + 0.5)
    cx = lambda c: min_x + cs * (c + 0.5)
    for p0, p1, p2 in tris:
        if ((p0[0] < min_x and p1[0] < min_x and p2[0] < min_x) or (p0[0] > max_x and p1[0] > max_x and p2[0] > max_x)
                or (p0[1] < min_y and p1[1] < min_y and p2[1] < min_y) or (p0[1] > max_y and p1[1] > max_y and p2[1] > max_y)):
            continue
        # PlaneM.SetCoeffientsPtPtPt
        d1 = (p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2])
        d2 = (p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2])
        A = d1[1] * d2[2] - d1[2] * d2[1]
        B = d1[2] * d2[0] - d1[0] * d2[2]
        C = d1[0] * d2[1] - d1[1] * d2[0]
        D = (0.0 - A) * p0[0] - B * p0[1] - C * p0[2]
        step = _div(0.0 - A, C) * cs
        a, b, c = _sort3_by_y(p0, p1, p2)
        if a[1] == c[1]:
            continue
        if c[1] > max_y:
            r0 = 0
        else:
            r0 = row(c[1])
            if c[1] < cy(r0):
                r0 += 1
        if a[1] < min_y:
            r1 = rows - 1
        else:
            r1 = row(a[1])
            if a[1] > cy(r1):
                r1 -= 1
        if r0 > r1:
            continue
        s4 = _div(c[0] - a[0], c[1] - a[1])
        s5 = _div(c[0] - b[0], c[1] - b[1])
        s6 = _div(b[0] - a[0], b[1] - a[1])
        y = cy(r0)
        x4 = c[0] + s4 * (y - c[1])
        x5 = c[0] + s5 * (y - c[1])
        x6 = b[0] + s6 * (y - b[1])
        for r in range(r0, r1 + 1):
            n1 = x4
            n2 = x6 if (not (y > b[1]) and not math.isnan(x6)) else x5
            if n1 > n2:
                n1, n2 = n2, n1
            c0 = col(n1)
            if cx(c0) < n1:
                c0 += 1
            c1 = col(n2)
            if cx(c1) > n2:
                c1 -= 1
            if c1 >= 0 and c0 < cols:
                c0 = max(c0, 0)
                c1 = min(c1, cols - 1)
                z = (0.0 if C == 0.0 else (0.0 - D - A * cx(c0) - B * y) / C)
                for cc in range(c0, c1 + 1):
                    values[r * cols + cc] = f32(z)
                    z += step
            y -= cs
            x4 -= s4 * cs
            x5 -= s5 * cs
            x6 -= s6 * cs


def _sort3_by_y(a, b, c):
    """PointM.Sort3PointsByY."""
    if not (a[1] <= b[1]) or not (b[1] <= c[1]):
        if a[1] <= c[1] <= b[1]:
            b, c = c, b
        elif b[1] <= a[1] <= c[1]:
            a, b = b, a
        elif b[1] <= c[1] <= a[1]:
            a, b = b, a
            b, c = c, b
        elif c[1] <= a[1] <= b[1]:
            a, b = b, a
            a, c = c, a
        else:
            a, c = c, a
    return a, b, c


class PolygonTin:
    """Geospatial.Vectors DelaunayTinFactory.ComputeDevelopmentTin of one
    polygon (perimeter only, deleted triangles kept) and Tin.ContainsPoint,
    which ParameterResampler uses for a calibration region. A point on the
    region's edge is inside when the walk reaches an inside triangle first,
    and the walk starts where the last query ended, so queries must come in
    RASMapper's order."""

    def __init__(self, ring):
        pts = [(float(x), float(y)) for x, y in ring]       # open ring
        area = 0.0                                           # Polygon.ComputeArea
        for i in range(len(pts)):
            (ax, ay), (bx, by) = pts[i], pts[(i + 1) % len(pts)]
            area += (ay + by) * (bx - ax)
        self.tol = math.sqrt(abs(area * 0.5)) / 1000000.0
        self.points = []
        self.tris = []
        self.free = []
        self.last_t = 0                                      # _previousSearchTriangle
        self._ring_count = len(pts)
        self._build(pts)
        xs = [q[0] for q in self.points]
        ys = [q[1] for q in self.points]
        self.extent = (min(xs), max(xs), min(ys), max(ys))   # Tin.Extent, cached

    # -- Tin --------------------------------------------------------------------
    def add_triangle(self, p0, p1, p2, at0, at1, at2):
        t = Tri(p0, p1, p2, at0, at1, at2)
        if self.free:
            n = self.free.pop()
            self.tris[n] = t
        else:
            n = len(self.tris)
            self.tris.append(t)
        return n

    def find_triangle(self, pt):
        """Tin.FindTriangleContainingPoint."""
        t = self.last_t
        if t < 0 or t > len(self.tris) - 1:
            t = 0
        count, prev = 0, -1
        P = self.points
        while t != -1 and count != len(self.tris):
            tr = self.tris[t]
            c0 = _cross(P[tr.p0], P[tr.p1], pt)
            c1 = _cross(P[tr.p1], P[tr.p2], pt)
            c2 = _cross(P[tr.p2], P[tr.p0], pt)
            k = (0 if c0 < c2 else 2) if c0 < c1 else (1 if c1 < c2 else 2)   # Math.Min3Index
            if (c0, c1, c2)[k] >= 0.0:
                self.last_t = t
                return t
            nxt = (tr.at0, tr.at1, tr.at2)[k]
            if prev == nxt:
                return min(nxt, t)
            prev, t = t, nxt
            count += 1
        return -1

    def _closest(self, pt, p0, p1, p2):
        P = self.points
        d0, d1, d2 = (_hypot(P[q], pt) for q in (p0, p1, p2))
        if d0 < d1:
            return p0 if d0 < d2 else p2
        return p1 if d1 < d2 else p2

    def _near(self, a, b):
        """Point.CloseToXY."""
        return abs(a[0] - b[0]) < self.tol and abs(a[1] - b[1]) < self.tol

    def contains(self, pt) -> bool:
        """Tin.ContainsPoint (the extent check, then the walk)."""
        min_x, max_x, min_y, max_y = self.extent
        if not (min_x <= pt[0] < max_x and min_y < pt[1] <= max_y):
            return False
        t = self.find_triangle(pt)
        return t != -1 and not self.tris[t].deleted()

    def _set_edge_type(self, t, a, b, perimeter):
        """Tin.SetEdgeFlags(GetEdgeDefinition(t, a, b)): GetEdge012 answers
        edge 2 for an edge the triangle doesn't have."""
        tr = self.tris[t]
        e = 0 if {a, b} == {tr.p0, tr.p1} else (1 if {a, b} == {tr.p1, tr.p2} else 2)
        tr.set_edge(e, perimeter)
        nb = (tr.at0, tr.at1, tr.at2)[e]
        if nb != -1:
            o = self.tris[nb]
            o.set_edge(0 if o.at0 == t else (1 if o.at1 == t else 2), perimeter)

    # -- DelaunayTinFactory ---------------------------------------------------------
    def _build(self, ring):
        min_x, max_x = min(p[0] for p in ring), max(p[0] for p in ring)
        min_y, max_y = min(p[1] for p in ring), max(p[1] for p in ring)
        ext = _enforce_min_dims([min_x, max_x, min_y, max_y], self.tol * 10.0, self.tol * 10.0)
        dx, dy = ext[1] - ext[0], ext[3] - ext[2]
        if dx > dy:
            if dy < dx * 0.125:
                ext = _enforce_min_dims(ext, 0.0, dx * 0.125)
        elif dx < dy * 0.125:
            ext = _enforce_min_dims(ext, dy * 0.125, 0.0)
        ox, oy = 0.0 - ext[0], 0.0 - ext[2]
        lo_x, hi_x, lo_y, hi_y = ext[0] + ox, ext[1] + ox, ext[2] + oy, ext[3] + oy
        n = max(hi_x - lo_x, hi_y - lo_y) * 0.05
        mid_x, mid_y = (hi_x + lo_x) * 0.5, (hi_y + lo_y) * 0.5
        for x, y in ((mid_x, hi_y + n * 2.0), (hi_x + n, hi_y + n), (hi_x + n * 2.0, mid_y),
                     (hi_x + n, lo_y - n), (mid_x, lo_y + (-n) * 2.0), (lo_x - n, lo_y - n),
                     (lo_x + (-n) * 2.0, mid_y), (lo_x - n, hi_y + n)):
            self.points.append((x, y))
        if (ext[1] - ext[0]) > (ext[3] - ext[2]):
            spec = [(0, 7, 6, -1, -1, 4), (6, 5, 4, -1, -1, 4), (4, 3, 2, -1, -1, 5),
                    (2, 1, 0, -1, -1, 5), (0, 6, 4, 0, 1, 5), (0, 4, 2, 4, 2, 3)]
        else:
            spec = [(0, 7, 6, -1, -1, 4), (6, 5, 4, -1, -1, 5), (4, 3, 2, -1, -1, 5),
                    (2, 1, 0, -1, -1, 4), (0, 6, 2, 0, 5, 3), (6, 4, 2, 1, 2, 4)]
        for s in spec:
            self.add_triangle(*s)
        moved = [(x + ox, y + oy) for x, y in ring]
        for pt in moved:                                     # InsertPerimeterPoints
            self._insert(pt)
        for i in range(len(moved)):                          # EnforcePerimeterEdges
            self._enforce(moved[i], moved[(i + 1) % len(moved)])
        self._mark_bad()
        self.points = [(x + -ox, y + -oy) for x, y in self.points]

    def _insert(self, pt):
        """DelaunayInsertPoint (FindTriangleCheckNearbyCorners, then the hole)."""
        t = self.find_triangle(pt)
        if t != -1:
            tr = self.tris[t]
            P = self.points
            if any(self._near(pt, P[q]) for q in (tr.p0, tr.p1, tr.p2)):
                return False
        p = len(self.points)
        self.points.append(pt)
        stack = [t]
        while stack:
            t = stack.pop()
            self.free.append(t)
            tr = self.tris[t]
            tr.mark_deleted()
            for at in (tr.at0, tr.at1, tr.at2):
                if self._inside_hole(at, pt):
                    stack.append(at)
        self._fill_hole(p)
        return True

    def _inside_hole(self, at, pt):
        if at == -1:
            return False
        tr = self.tris[at]
        if tr.deleted():
            return False
        P = self.points
        a, b, c = P[tr.p0], P[tr.p1], P[tr.p2]
        x, y, _ = _circumcircle(a, b, c)
        if not (_hypot(pt, (x, y)) <= _hypot(a, (x, y))):   # Circle.Contains
            return False
        for u, v in ((a, b), (a, c), (b, c)):
            if u[0] == v[0] == pt[0] and ((pt[1] > u[1] and pt[1] > v[1]) or (pt[1] < u[1] and pt[1] < v[1])):
                return False
        for u, v in ((a, b), (a, c), (b, c)):
            if u[1] == v[1] == pt[1] and ((pt[0] > u[0] and pt[0] > v[0]) or (pt[0] < u[0] and pt[0] < v[0])):
                return False
        return True

    _fill_hole = StitchTin._fill_hole

    def _enforce(self, a, b):
        """EnforceEdgeSegment for a perimeter segment (already translated)."""
        P = self.points
        t = self.find_triangle(b)
        tr = self.tris[t]
        p = self._closest(b, tr.p0, tr.p1, tr.p2)
        b = P[p]
        while True:
            t2 = self.find_triangle(a)
            tr = self.tris[t2]
            p2 = self._closest(a, tr.p0, tr.p1, tr.p2)
            a = self.points[p2]
            while True:                                      # SpinFindTriangleToBisect
                tr = self.tris[t2]
                if p2 == tr.p0:
                    e, n_cw, n_ccw = 1, tr.at0, tr.at2
                elif p2 == tr.p1:
                    e, n_cw, n_ccw = 2, tr.at1, tr.at0
                else:
                    e, n_cw, n_ccw = 0, tr.at2, tr.at1
                sa, sb = ((self.points[tr.p0], self.points[tr.p1]), (self.points[tr.p1], self.points[tr.p2]),
                          (self.points[tr.p2], self.points[tr.p0]))[e]
                if _cross(self.points[p2], sa, b) < 0.0:
                    t2 = n_cw
                    continue
                if _cross(self.points[p2], sb, b) > 0.0:
                    t2 = n_ccw
                    continue
                break
            if p in (tr.p0, tr.p1, tr.p2):
                break
            t = _geo_intersect_t(a, b, sa, sb)
            pm = (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))
            first, last = (tr.p0, tr.p1, tr.p2)[e], (tr.p1, tr.p2, tr.p0)[e]
            if self._near(pm, sa):
                self._set_edge_type(t2, p2, first, True)
                p2 = first
                a = self.points[p2]
                continue
            if self._near(pm, sb):
                self._set_edge_type(t2, p2, last, True)
                p2 = last
                a = self.points[p2]
                continue
            count = len(self.points)
            self.points.append(pm)
            adj = (tr.at0, tr.at1, tr.at2)[e]
            bl = tr.flags & (E0_BREAK, E1_BREAK, E2_BREAK)[e] != 0
            perim = tr.flags & (E0_PERIM, E1_PERIM, E2_PERIM)[e] != 0
            self.free.append(t2)
            self.free.append(adj)
            self.tris[t2].mark_deleted()
            self.tris[adj].mark_deleted()
            edges = self._fill_hole(count)
            if bl or perim:
                for edge in edges:
                    self._set_edge_type(edge[3], count, first, perim)
                    self._set_edge_type(edge[3], count, last, perim)
            self._set_edge_type(len(self.tris) - 1, p2, count, True)
            a = pm
        self._set_edge_type(t2, p2, p, True)

    def _mark_bad(self):
        """MarkBadTrianglesAsDeleted_SinglePart: flood from the super
        triangles, stopping at perimeter edges."""
        stack = []
        for i, tr in enumerate(self.tris):
            if tr.p0 < SUPER_POINT_COUNT or tr.p1 < SUPER_POINT_COUNT or tr.p2 < SUPER_POINT_COUNT:
                tr.mark_deleted()
                stack.append(i)
        if self._ring_count <= 3:
            return
        while stack:
            t = stack.pop()
            tr = self.tris[t]
            tr.mark_deleted()
            for at, flag in ((tr.at0, E0_PERIM), (tr.at1, E1_PERIM), (tr.at2, E2_PERIM)):
                if at != -1 and not self.tris[at].deleted() and not tr.flags & flag:
                    stack.append(at)


def _hypot(a, b):
    """Point.DistanceTo."""
    dx, dy = a[0] - b[0], a[1] - b[1]
    return math.sqrt(dx * dx + dy * dy)


def _geo_intersect_t(a, b, c, d):
    """Geospatial Segment.Intersects(seg, ..., epsilon 0): t along a-b, NaN
    when the segments' boxes don't meet or they are parallel."""
    if (max(a[0], b[0]) < min(c[0], d[0]) or min(a[0], b[0]) > max(c[0], d[0])
            or max(a[1], b[1]) < min(c[1], d[1]) or min(a[1], b[1]) > max(c[1], d[1])):
        return math.nan
    dx, dy = b[0] - a[0], b[1] - a[1]
    ex, ey = d[0] - c[0], d[1] - c[1]
    den = ey * dx - ex * dy
    if den == 0.0:
        return math.nan
    return (ex * (a[1] - c[1]) - ey * (a[0] - c[0])) / den


def paint_bary(tris, values, rows, cols, min_x, max_x, min_y, max_y, cs):
    """Tin.Rasterize → RasterizeTriangle.ComputeCellsZNew: each triangle's
    plane onto the cells whose centre it covers, edges included, stepping
    barycentric weights across the raster as RASMapper does."""
    def _col(x):
        return int(math.floor((x - min_x) / cs))

    def _row(y):
        return int(math.floor((max_y - y) / cs))

    for p0, p1, p2 in tris:
        if (p0[1] == p1[1] == p2[1]) or (p0[0] == p1[0] == p2[0]):
            continue
        lo_x, hi_x = min(p0[0], p1[0], p2[0]), max(p0[0], p1[0], p2[0])
        lo_y, hi_y = min(p0[1], p1[1], p2[1]), max(p0[1], p1[1], p2[1])
        if lo_x > max_x or hi_x < min_x or lo_y > max_y or hi_y < min_y:
            continue
        c0 = 0 if lo_x < min_x else _col(lo_x)
        r0 = rows - 1 if lo_y < min_y else _row(lo_y)
        c1 = cols - 1 if hi_x > max_x else _col(hi_x)
        r1 = 0 if hi_y > max_y else _row(hi_y)
        ox, oy = min_x + cs * (c0 + 0.5), max_y - cs * (r0 + 0.5)
        q = [(p[0] - ox, p[1] - oy) for p in (p0, p1, p2)]

        def bary(v0, v1):                                   # GetBary at the origin
            dx, dy = v0[1] - v1[1], v1[0] - v0[0]
            return [dx * 0.0 + dy * 0.0 + (v0[0] * v1[1] - v0[1] * v1[0]), dx, dy]

        b = [bary(q[1], q[2]), bary(q[2], q[0]), bary(q[0], q[1])]
        sx = [w[1] * cs for w in b]
        sy = [w[2] * cs for w in b]
        total = b[0][0] + b[1][0] + b[2][0]
        for r in range(r0, r1 - 1, -1):
            w0, w1, w2 = b[0][0], b[1][0], b[2][0]
            for c in range(c0, c1 + 1):
                if w0 >= 0.0 and w1 >= 0.0 and w2 >= 0.0:
                    values[r * cols + c] = (w0 * p0[2] + w1 * p1[2] + w2 * p2[2]) / total
                w0 += sx[0]
                w1 += sx[1]
                w2 += sx[2]
            for k in range(3):
                b[k][0] += sy[k]


class DelaunayTin(StitchTin):
    """DelaunayTINFactory.ComputeTIN with a perimeter and/or mass points (the
    path RASMapper takes for a 2D area's perimeter TIN and for the TIN of its
    cell centres), deleted triangles kept: a DevelopmentTin. Points include
    the eight super points; tris include deleted ones, with adjacency."""

    def __init__(self, tolerance, perimeter=None, points=None, mbr=None):
        super().__init__([], tolerance)
        # perimeter: a closed ring (first point repeated last), as RASMapper's Polygon
        self.perimeter = [(float(x), float(y), math.nan, math.nan) for x, y in perimeter] if perimeter is not None else None
        self.mass = [(float(x), float(y), math.nan, math.nan) for x, y in points] if points is not None else None
        self.mbr = mbr                             # (min_x, max_x, min_y, max_y) or None
        self.duplicates = []                       # mass point indexes already in the TIN

    def compute(self):
        xs, ys = [], []
        for pts in (self.perimeter, self.mass):
            if pts:
                xs += [p[0] for p in pts]
                ys += [p[1] for p in pts]
        ext = _enforce_min_dims([min(xs), max(xs), min(ys), max(ys)], self.tol * 10.0, self.tol * 10.0)
        dx, dy = ext[1] - ext[0], ext[3] - ext[2]
        if dx > dy:
            if dy < dx * 0.125:
                ext = _enforce_min_dims(ext, 0.0, dx * 0.125)
        elif dx < dy * 0.125:
            ext = _enforce_min_dims(ext, dy * 0.125, 0.0)
        if self.mbr is not None:                   # Extent.UnionXY
            ext = [min(ext[0], self.mbr[0]), max(ext[1], self.mbr[1]),
                   min(ext[2], self.mbr[2]), max(ext[3], self.mbr[3])]
        self.extent = ext
        self.offset = (-ext[0], -ext[2])
        self._add_octagon(ext)
        for p in self.perimeter or ():             # InsertPerimeterPoints
            self.insert_point(self._translate(p))
        self.mass_start = len(self.points)         # MassPointInsertionsStartHere
        for i, p in enumerate(self.mass or ()):    # InsertMassPoints
            if not self.insert_point(self._translate(p)):
                self.duplicates.append(i)
        if self.perimeter:                         # EnforcePerimeterEdges
            for a, b in zip(self.perimeter[:-1], self.perimeter[1:]):
                self.enforce_edge(self._translate(a), self._translate(b), perimeter=True)
        self._mark_bad()
        ox, oy = self.offset
        self.points = [(x + -ox, y + -oy, z, m) for x, y, z, m in self.points]
        self.last_t = 0                            # the DevelopmentTin's own search starts afresh
        return self

    def _mark_bad(self):
        """MarkBadTrianglesAsDeleted_SinglePart without stitches: the super
        triangles, then, with a perimeter, everything reached from them
        without crossing a perimeter edge."""
        stack = []
        for i, tr in enumerate(self.tris):
            if tr.p0 < SUPER_POINT_COUNT or tr.p1 < SUPER_POINT_COUNT or tr.p2 < SUPER_POINT_COUNT:
                tr.mark_deleted()
                stack.append(i)
        if not self.perimeter or len(self.perimeter) <= 3:
            return
        while stack:
            tr = self.tris[stack.pop()]
            tr.mark_deleted()
            for at, flag in ((tr.at0, E0_PERIM), (tr.at1, E1_PERIM), (tr.at2, E2_PERIM)):
                if at != -1 and not self.tris[at].deleted() and not tr.flags & flag:
                    stack.append(at)

    def contains(self, pt) -> bool:
        """DevelopmentTin.ContainsPoint: inside the TIN's extent, then in a
        triangle that isn't deleted."""
        x, y = float(pt[0]), float(pt[1])
        min_x, max_x, min_y, max_y = self.extent
        if x < min_x or x > max_x or y < min_y or y > max_y:
            return False
        t = self.find_triangle((x, y))
        return t != -1 and not self.tris[t].deleted()
