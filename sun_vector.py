import math
import numpy as np
from datetime import datetime, timezone

# ---------- Utilities ----------
def _normalize(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("Zero-length vector.")
    return v / n

def _rodrigues(v, k, ang):
    k = _normalize(k); v = np.asarray(v, dtype=float)
    c, s = math.cos(ang), math.sin(ang)
    return v*c + np.cross(k, v)*s + k*np.dot(k, v)*(1-c)

def _rotmat_to_quat_wxyz(R):
    """Convert proper rotation matrix to quaternion (w,x,y,z)."""
    m = np.asarray(R, dtype=float)
    tr = m[0,0] + m[1,1] + m[2,2]
    if tr > 0:
        S = math.sqrt(tr + 1.0) * 2.0
        w = 0.25 * S
        x = (m[2,1] - m[1,2]) / S
        y = (m[0,2] - m[2,0]) / S
        z = (m[1,0] - m[0,1]) / S
    else:
        if m[0,0] > m[1,1] and m[0,0] > m[2,2]:
            S = math.sqrt(1.0 + m[0,0] - m[1,1] - m[2,2]) * 2.0
            w = (m[2,1] - m[1,2]) / S
            x = 0.25 * S
            y = (m[0,1] + m[1,0]) / S
            z = (m[0,2] + m[2,0]) / S
        elif m[1,1] > m[2,2]:
            S = math.sqrt(1.0 + m[1,1] - m[0,0] - m[2,2]) * 2.0
            w = (m[0,2] - m[2,0]) / S
            x = (m[0,1] + m[1,0]) / S
            y = 0.25 * S
            z = (m[1,2] + m[2,1]) / S
        else:
            S = math.sqrt(1.0 + m[2,2] - m[0,0] - m[1,1]) * 2.0
            w = (m[1,0] - m[0,1]) / S
            x = (m[0,2] + m[2,0]) / S
            y = (m[1,2] + m[2,1]) / S
            z = 0.25 * S
    q = np.array([w, x, y, z], dtype=float)
    return q / np.linalg.norm(q)

# ---------- Sun vector (ICRF/J2000) ----------
def sun_vector_icrf(utc_datetime):
    """
    Low-precision Sun direction in ICRF/J2000 (unit vector).
    Input: timezone-aware datetime in UTC.
    Source: standard solar position approximation (Meeus-style).
    """
    if utc_datetime.tzinfo is None:
        # assume UTC if naive
        utc_datetime = utc_datetime.replace(tzinfo=timezone.utc)
    else:
        utc_datetime = utc_datetime.astimezone(timezone.utc)

    # Julian Date (UTC ~ TDB here; fine for low-precision pointing)
    y = utc_datetime.year
    m = utc_datetime.month
    d = utc_datetime.day + (utc_datetime.hour + (utc_datetime.minute + utc_datetime.second/60)/60)/24
    if m <= 2:
        y -= 1
        m += 12
    A = math.floor(y/100)
    B = 2 - A + math.floor(A/4)
    JD = math.floor(365.25*(y+4716)) + math.floor(30.6001*(m+1)) + d + B - 1524.5
    T = (JD - 2451545.0) / 36525.0  # Julian centuries since J2000

    # Mean anomaly of the Sun (deg)
    M = (357.52911 + 35999.05029*T - 0.0001537*T*T) % 360.0
    # Mean longitude (deg)
    L0 = (280.46646 + 36000.76983*T + 0.0003032*T*T) % 360.0
    # Eccentricity of Earth's orbit
    e = 0.016708634 - 0.000042037*T - 0.0000001267*T*T
    Mr = math.radians(M)

    # Equation of center (deg)
    C = (1.914602 - 0.004817*T - 0.000014*T*T)*math.sin(Mr) \
        + (0.019993 - 0.000101*T)*math.sin(2*Mr) \
        + 0.000289*math.sin(3*Mr)

    # True longitude (deg)
    true_long = L0 + C
    # Apparent longitude (deg)
    omega = 125.04 - 1934.136*T
    lam = true_long - 0.00569 - 0.00478*math.sin(math.radians(omega))

    # Mean obliquity (deg) + nutation correction to get true obliquity
    eps0 = 23.439291 - 0.0130042*T - 1.64e-7*T*T + 5.04e-7*T*T*T
    eps = eps0 + 0.00256*math.cos(math.radians(omega))

    # Convert ecliptic lon to ICRF equatorial unit vector
    lam_r = math.radians(lam)
    eps_r = math.radians(eps)
    # Ecliptic unit vector (distance ignored)
    xe = math.cos(lam_r)
    ye = math.sin(lam_r)
    ze = 0.0
    # Rotate by obliquity to equatorial (ICRF)
    x = xe
    y = ye*math.cos(eps_r)
    z = ye*math.sin(eps_r)
    return _normalize(np.array([x, y, z], dtype=float))

# ---------- Build attitude with a chosen face pointing at Sun ----------
def _face_axes():
    return {
        "+X": np.array([1.0, 0.0, 0.0]),
        "-X": np.array([-1.0, 0.0, 0.0]),
        "+Y": np.array([0.0, 1.0, 0.0]),
        "-Y": np.array([0.0, -1.0, 0.0]),
        "+Z": np.array([0.0, 0.0, 1.0]),
        "-Z": np.array([0.0, 0.0, -1.0]),
    }

def _rotation_body_to_icrf_with_face(S_icrf, face_key="+Z", roll_rad=0.0):
    """
    Construct a proper rotation matrix R such that R * e_body = v_icrf,
    with the chosen face normal e_body aligned to S_icrf, and a roll about S_icrf.
    Columns of R are the ICRF directions of body +X,+Y,+Z.
    """
    S = _normalize(S_icrf)
    faces = _face_axes()
    if face_key not in faces:
        raise ValueError("face_key must be one of +X,-X,+Y,-Y,+Z,-Z")
    f = faces[face_key]  # body-face vector

    # Reference axis to define roll=0 (avoid near parallel with S)
    ref = np.array([0.0, 0.0, 1.0]) if abs(np.dot(S, [0,0,1])) < 0.95 else np.array([1.0, 0.0, 0.0])
    u0 = _normalize(ref - np.dot(ref, S)*S)      # in plane ⟂ S
    u = _rodrigues(u0, S, roll_rad)              # roll by ψ
    w = _normalize(np.cross(S, u))               # completes triad

    R = np.zeros((3,3))
    # Place columns so the chosen body axis equals ±S
    if face_key.endswith("X"):
        sign = np.sign(f[0])
        R[:,0] = sign * S
        R[:,1] = w
        R[:,2] = _normalize(np.cross(R[:,0], R[:,1]))
        R[:,1] = _normalize(np.cross(R[:,2], R[:,0]))
    elif face_key.endswith("Y"):
        sign = np.sign(f[1])
        R[:,1] = sign * S
        R[:,2] = u
        R[:,0] = _normalize(np.cross(R[:,1], R[:,2]))
        R[:,2] = _normalize(np.cross(R[:,0], R[:,1]))
    else:  # Z
        sign = np.sign(f[2])
        R[:,2] = sign * S
        R[:,0] = u
        R[:,1] = _normalize(np.cross(R[:,2], R[:,0]))
        R[:,0] = _normalize(np.cross(R[:,1], R[:,2]))

    # Ensure det +1
    if np.linalg.det(R) < 0:
        R[:,1] *= -1.0
    return R

# ---------- Main function ----------
def sun_pointing_quaternion_at_time(
    utc_datetime,
    face="+Z",
    roll_deg=0.0,
    convention="body_to_icrf",
    order="wxyz"
):
    """
    Return a unit quaternion for the attitude that points the chosen body face at the Sun at 'utc_datetime'.

    Parameters
    ----------
    utc_datetime : datetime (UTC; naive treated as UTC)
    face         : one of {+X,-X,+Y,-Y,+Z,-Z}, default +Z
    roll_deg     : roll about the Sun vector (deg), right-hand rule (ICRF)
    convention   : "body_to_icrf"  (quaternion rotates body->ICRF)
                   "icrf_to_body"  (inverse; for STK 'From File')
    order        : "wxyz" (scalar first) or "xyzw" (scalar last)

    Returns
    -------
    np.ndarray shape (4,) in requested component order and convention.
    """
    S = sun_vector_icrf(utc_datetime)
    R_b2i = _rotation_body_to_icrf_with_face(S, face_key=face, roll_rad=math.radians(roll_deg))
    q_wxyz = _rotmat_to_quat_wxyz(R_b2i)

    if convention == "icrf_to_body":
        # Inverse of unit quaternion is its conjugate
        w, x, y, z = q_wxyz
        q_wxyz = np.array([w, -x, -y, -z], dtype=float)

    if order.lower() == "wxyz":
        return q_wxyz
    elif order.lower() == "xyzw":
        w, x, y, z = q_wxyz
        return np.array([x, y, z, w], dtype=float)
    else:
        raise ValueError("order must be 'wxyz' or 'xyzw'.")

# -------------- Example --------------
if __name__ == "__main__":
    # Example: point body +Z at the Sun at a specific UTC time, with 30° roll about Sun
    t = datetime(2025, 1, 1, 12, 0, 0, tzinfo=timezone.utc)

    # For Python use (body->ICRF, scalar-first):
    q_b2i_wxyz = sun_pointing_quaternion_at_time(t, face="+Z", roll_deg=30.0,
                                                 convention="body_to_icrf", order="wxyz")
    print("Body->ICRF (w,x,y,z):", q_b2i_wxyz)

    # For STK 'From File' (ICRF->Body, scalar-last):
    q_icrf2body_xyzw = sun_pointing_quaternion_at_time(t, face="+Z", roll_deg=30.0,
                                                       convention="icrf_to_body", order="xyzw")
    print("ICRF->Body (x,y,z,w) for STK:", q_icrf2body_xyzw)
    
    


# S = sun_vector_icrf(t)
#         rolls = np.linspace(0.0, 2*np.pi, n_roll, endpoint=False)
#         quats = []
#         quats_by_face = {}
#         for psi in rolls:
#             R = complete_attitude_columns_for_face(S, "+Z", psi)
#             q = rotmat_to_quat(R)
#             quats.append(q)
#         quats_by_face["+Z"] = np.vstack(quats)

#         # Example: export +Z face
#         qZ = quats_by_face["+Z"]
#         df = pd.DataFrame(qZ, columns=["q0 (scalar)", "q1", "q2", "q3"])
#         df.to_csv("SunPointingQuaternions_Zface.csv", index=False)
#         print("Saved to SunPointingQuaternions_Zface.csv")

#         # Example: save +Z face sequence
#         save_stk_attitude(quats_by_face["+Z"], n_theta, n_roll, dt=10.0, filename="Zface_quaternions.a")
#         print("Saved attitude file: Zface_quaternions.a")