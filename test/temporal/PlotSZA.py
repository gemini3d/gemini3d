import astropy.coordinates as ac
import astropy
import astropy.units as u
from argparse import ArgumentParser


def solarzenithangle(t, glat, glon, alt_m):
    """

    Parameters
    ----------

    t : datetime
      time of observation
    glat : float
      latitude
    glon : float
      longitude
    alt_m : float
      observer altitude [meters]

    Returns
    -------
    sza : float
      solar zenith angle [degrees]
    """

    obs = ac.EarthLocation(lat=glat * u.deg, lon=glon * u.deg, height=alt_m * u.m)
    times = astropy.time.Time(t, scale="ut1")
    sun = ac.get_sun(times)
    sunobs = sun.transform_to(ac.AltAz(obstime=times, location=obs))

    return 90.0 - sunobs.alt.degree


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("time", help="UTC")
    p.add_argument("latlon", nargs=2, type=float)
    p.add_argument("altitude", help="meters", type=float)
    a = p.parse_args()

    sza = solarzenithangle(a.time, a.latlon[0], a.latlon[1], a.altitude)

    print(sza)
