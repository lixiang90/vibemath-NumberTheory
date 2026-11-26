"""Enclose the contact and test smaller targets in Liu's Theorem 7.33.

Requires python-flint 0.9.0. The input is the published differential equation,
initial value and active branch. All ODE steps and comparisons use Arb balls.
"""

from flint import arb, arb_series, ctx


def curve(r):
    radical = (7*r**4*(r-arb(2)/3)).root(5)
    return arb(5)/24 + (19*r-5*radical)/16


def flow(y, r, radical):
    """y derivatives of the cutoff radius and its positive fifth root."""
    curve_value = arb(5)/24 + (19*r-5*radical)/16
    width = (curve_value-arb(1)/2-y)/2
    width_rate = (width+arb(1)/3-(y+r-curve_value)/(5*(y+r-1)))/y
    radical_r = radical*(4/r+1/(r-arb(2)/3))/5
    radius_rate = (1+2*width_rate)*16/(19-5*radical_r)
    return radius_rate, radical_r*radius_rate


def series(y, r, radical, order):
    rs, hs = [r], [radical]
    for k in range(order):
        yy = arb_series([y, 1], prec=k+1)
        dr, dh = flow(yy, arb_series(rs, prec=k+1), arb_series(hs, prec=k+1))
        rs.append(dr[k]/(k+1))
        hs.append(dh[k]/(k+1))
    return rs, hs


def step(y, r, radical, length, order):
    # The interval Picard bound encloses the trajectory before the Taylor
    # remainder is evaluated. length is signed and may itself be an interval.
    radius = abs(length).upper()
    time_box = y.union(y+length)
    r_box = r+arb(0, 6*radius)
    h_box = radical+arb(0, 16*radius)
    dr, dh = flow(time_box, r_box, h_box)
    if not (abs(dr)<6 and abs(dh)<16):
        raise ArithmeticError("Picard tube fails to enclose the flow")
    coefficients = series(y, r, radical, order-1)
    remainder = series(time_box, r_box, h_box, order)
    result = []
    for poly, bounds in zip(coefficients, remainder, strict=True):
        value = bounds[-1]
        for coeff in reversed(poly):
            value = coeff+length*value
        result.append(value)
    return tuple(result)


class Width:
    def __init__(self, division=4096, order=28):
        self.division = division
        self.order = order
        self.start = arb(7)/18
        self.nodes = [(arb(7)/9, arb(7)/9)]

    def __call__(self, y):
        if not arb(".369") < y < self.start:
            raise ArithmeticError("outside the published active-branch range")
        index = int(((self.start-y)*self.division).floor().unique_fmpz())
        length = -arb(1)/self.division
        while len(self.nodes)<=index:
            k = len(self.nodes)-1
            self.nodes.append(step(self.start+k*length, *self.nodes[k], length, self.order))
        base = self.start+index*length
        r, radical = step(base, *self.nodes[index], y-base, self.order)
        if not (radical**5).overlaps(7*r**4*(r-arb(2)/3)):
            raise ArithmeticError("inconsistent fifth-root enclosure")
        width = (arb(5)/24+(19*r-5*radical)/16-arb(1)/2-y)/2
        if not 2*y < r:
            raise ArithmeticError("wrong cutoff branch at the requested point")
        return width


def radius_at_target(target):
    lo, hi = arb(".738"), arb(7)/9
    if not curve(lo) < target < curve(hi):
        raise ArithmeticError("critical radius outside the curved range")
    for _ in range(168):
        mid = (lo+hi)/2
        val = curve(mid)
        if val < target:
            lo = mid
        elif val > target:
            hi = mid
        else:
            # The derivative of the curve exceeds 0.35 on this bracket.
            return mid+arb(0, abs(val-target)/arb(".35"))
    return lo.union(hi)


def contact_data(theta, width):
    length = arb(11)/3-4*theta
    d = (theta-arb(1)/2-length/2)/(arb(1)/2+length)
    r = radius_at_target(4/(9*d)-arb(1)/6)
    z = (arb(7)/18+d/6)/(arb(5)/6-d)-r
    f = width(z)
    required = 2/(9*d)-z/2-arb(1)/3
    deficit = required-f
    return z, f, deficit, -d*deficit/2


def main():
    ctx.prec = 272
    ctx.cap = 40
    width = Width()
    lo, hi = arb(".8733326374459049352"), arb(".8733326374459051924")
    if not (contact_data(lo, width)[-1]<0<contact_data(hi, width)[-1]):
        raise ArithmeticError("published contact bracket was not reproduced")
    for _ in range(46):
        mid = (lo+hi)/2
        margin = contact_data(mid, width)[-1]
        if margin<0:
            lo = mid
        elif margin>0:
            hi = mid
        else:
            break
    print("contact", lo.union(hi).str(42), "width of bracket", (hi-lo).str(6))
    reported_lo = arb(".87333263744590506611657458483")
    reported_hi = arb(".87333263744590506611657458485")
    if not reported_lo < lo < hi < reported_hi:
        raise ArithmeticError("the bisection is outside the reported enclosure")
    # A different mesh and truncation order certify the signs directly at
    # the two published rational endpoints, without using the bisection.
    second = Width(8192, 32)
    for engine in (width, second):
        signs = [contact_data(t, engine)[-1] for t in (reported_lo, reported_hi)]
        if not signs[0] < 0 < signs[1]:
            raise ArithmeticError("endpoint signs failed")
        print("endpoint margins", *(x.str(12) for x in signs))
    targets = ((".8733326", "1.1328391e-6", "1.1328393e-6"),
               (".873332", "1.9273564e-5", "1.9273566e-5"),
               (".87333", "7.9594571e-5", "7.9594573e-5"),
               (".8733", ".0009569133", ".0009569134"))
    for theta, missing_lo, missing_hi in targets:
        z, f, deficit, margin = contact_data(arb(theta), width)
        if not (margin<0 and arb(missing_lo)<deficit<arb(missing_hi)):
            raise ArithmeticError("failed to enclose the width deficit")
        if not second(z).overlaps(f):
            raise ArithmeticError("independent ODE meshes disagree at the target")
        print(theta, "short length", z.str(14), "ordered width", f.str(14),
              "missing width", deficit.str(14), "margin", margin.str(14))


if __name__ == "__main__":
    main()
