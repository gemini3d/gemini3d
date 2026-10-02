"""IGRF-14 geocentric field oracle for qualification, not a GEMINI grid generator.

North/east/down components are in nT; latitude is geocentric degrees; radius
is geocentric km. Declared domain: 2025..2030, -89.9..89.9 latitude, radius
6371.2..8371.2 km. Full degree-13 field; no centered-dipole substitution.
"""
from math import sin, cos, sqrt, factorial, radians, isfinite
from pathlib import Path


def coefficients(path: Path):
    rows = {}
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) < 5 or fields[0] not in {"g", "h"}:
            continue
        kind, n, m = fields[:3]
        rows[kind, int(n), int(m)] = float(fields[-2]), float(fields[-1])
    if len(rows) != 195:
        raise ValueError("Expected all 195 IGRF degree-13 coefficients")
    return rows


def field(rows, year, radius_km, latitude_deg, longitude_deg):
    if not all(map(isfinite, [year, radius_km, latitude_deg, longitude_deg])):
        raise ValueError("Nonfinite IGRF input")
    if not (2025 <= year <= 2030 and 6371.2 <= radius_km <= 8371.2
            and -89.9 <= latitude_deg <= 89.9):
        raise ValueError("Outside declared IGRF qualification domain")
    theta, lon = radians(90-latitude_deg), radians(longitude_deg % 360)
    st, ct = sin(theta), cos(theta)
    p, dp = {(0, 0): 1.0}, {(0, 0): 0.0}
    br = bt = bp = 0.0
    for n in range(1, 14):
        for m in range(n+1):
            if m == n:
                p[n,m] = (2*n-1)*st*p[n-1,m-1]
                dp[n,m] = (2*n-1)*(ct*p[n-1,m-1]+st*dp[n-1,m-1])
            else:
                p[n,m] = ((2*n-1)*ct*p[n-1,m]-(n+m-1)*p.get((n-2,m),0))/(n-m)
                dp[n,m] = ((2*n-1)*(-st*p[n-1,m]+ct*dp[n-1,m])
                           -(n+m-1)*dp.get((n-2,m),0))/(n-m)
            norm = sqrt((2 if m else 1)*factorial(n-m)/factorial(n+m))
            g0, gv = rows['g',n,m]
            h0, hv = rows.get(('h',n,m),(0,0))
            g, h = g0+(year-2025)*gv, h0+(year-2025)*hv
            trig = g*cos(m*lon)+h*sin(m*lon)
            radial = (6371.2/radius_km)**(n+2)
            br += (n+1)*radial*trig*p[n,m]*norm
            bt -= radial*trig*dp[n,m]*norm
            bp += radial*m*(g*sin(m*lon)-h*cos(m*lon))*p[n,m]*norm/st
    return -bt, bp, -br
