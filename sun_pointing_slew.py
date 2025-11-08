import numpy as np
import pandas as pd
# ------------------------------
# Quaternion utilities (scalar-first: w, x, y, z)
# ------------------------------
def rotmat_to_quat(R):
    """Convert a proper rotation matrix R (3x3) to quaternion (w,x,y,z)."""
    m = R
    tr = m[0,0] + m[1,1] + m[2,2]
    if tr > 0:
        S = np.sqrt(tr + 1.0) * 2.0
        w = 0.25 * S
        x = (m[2,1] - m[1,2]) / S
        y = (m[0,2] - m[2,0]) / S
        z = (m[1,0] - m[0,1]) / S
    else:
        if (m[0,0] > m[1,1]) and (m[0,0] > m[2,2]):
            S = np.sqrt(1.0 + m[0,0] - m[1,1] - m[2,2]) * 2.0
            w = (m[2,1] - m[1,2]) / S
            x = 0.25 * S
            y = (m[0,1] + m[1,0]) / S
            z = (m[0,2] + m[2,0]) / S
        elif (m[1,1] > m[2,2]):
            S = np.sqrt(1.0 + m[1,1] - m[0,0] - m[2,2]) * 2.0
            w = (m[0,2] - m[2,0]) / S
            x = (m[0,1] + m[1,0]) / S
            y = 0.25 * S
            z = (m[1,2] + m[2,1]) / S
        else:
            S = np.sqrt(1.0 + m[2,2] - m[0,0] - m[1,1]) * 2.0
            w = (m[1,0] - m[0,1]) / S
            x = (m[0,2] + m[2,0]) / S
            y = (m[1,2] + m[2,1]) / S
            z = 0.25 * S
    q = np.array([w, x, y, z], dtype=float)
    # normalize to be safe
    return q / np.linalg.norm(q)

def normalize(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("Zero-length vector.")
    return v / n

def rodrigues(v, k, ang):
    """Rotate vector v about axis k (unit) by angle ang."""
    k = normalize(k)
    v = np.asarray(v, dtype=float)
    c = np.cos(ang); s = np.sin(ang)
    return v*c + np.cross(k, v)*s + k*np.dot(k, v)*(1-c)

# ------------------------------
# Orbit geometry in ICRF (J2000)
# ------------------------------
def orbit_normal_icrf(incl_deg, raan_deg):
    """Orbital angular momentum unit vector ĥ in ICRF given inclination and RAAN (deg)."""
    i = np.deg2rad(incl_deg)
    O = np.deg2rad(raan_deg)
    # Standard expression: ĥ = [sin i * sin Ω, -sin i * cos Ω, cos i]
    return normalize(np.array([np.sin(i)*np.sin(O),
                               -np.sin(i)*np.cos(O),
                               np.cos(i)]))

def node_direction_icrf(raan_deg):
    """Ascending node direction n̂ in ICRF as a unit vector (deg)."""
    O = np.deg2rad(raan_deg)
    return np.array([np.cos(O), np.sin(O), 0.0], dtype=float)

# ------------------------------
# Sun vector ring for a given beta (β)
# ------------------------------
def sun_vectors_for_beta(beta_deg, h_hat, n_hat, n_theta=1):
    """
    Generate unit Sun vectors S(θ) that are at angle β from the orbital plane.
    β is the angle from the plane (so S·h_hat = sin β). θ runs around ĥ.
    """
    beta = np.deg2rad(beta_deg)
    # Orthonormal basis spanning the orbital plane:
    u = normalize(n_hat)                         # along line of nodes
    v = normalize(np.cross(h_hat, u))            # completes RHS basis with ĥ
    # Parametric circle around ĥ:
    # S(θ) = cosβ * (cosθ u + sinθ v) + sinβ * ĥ
    # thetas = np.linspace(298*np.pi/180, 2*np.pi, n_theta, endpoint=False)
    thetas = [298*np.pi/180]
    S_list = []
    for th in thetas:
        S = (np.cos(beta)*(np.cos(th)*u + np.sin(th)*v) +
             np.sin(beta)*h_hat)
        S_list.append(normalize(S))
    return np.vstack(S_list), thetas

# ------------------------------
# Attitudes with one face pointing at the Sun
# ------------------------------
def face_axes():
    """Return dict of face normals in body frame."""
    return {
        "+X": np.array([1.0, 0.0, 0.0]),
        "-X": np.array([-1.0, 0.0, 0.0]),
        "+Y": np.array([0.0, 1.0, 0.0]),
        "-Y": np.array([0.0, -1.0, 0.0]),
        "+Z": np.array([0.0, 0.0, 1.0]),
        "-Z": np.array([0.0, 0.0, -1.0]),
    }

def complete_attitude_columns_for_face(S, face_key, roll_rad):
    """
    Build a proper rotation matrix R (body->ICRF) such that the chosen face normal equals S,
    and the body is rolled by 'roll_rad' about S.
    Columns of R are the ICRF directions of body +X, +Y, +Z.
    """
    faces = face_axes()
    f = faces[face_key]  # face normal in body frame (±ex, ±ey, ±ez)

    # Choose a fixed reference vector not parallel to S to seed an orthonormal pair
    # (this defines roll=0 direction before applying 'roll_rad').
    ref = np.array([0.0, 0.0, 1.0]) if abs(np.dot(S, [0,0,1])) < 0.95 else np.array([1.0, 0.0, 0.0])
    u0 = normalize(ref - np.dot(ref, S)*S)      # perpendicular to S
    u = rodrigues(u0, S, roll_rad)              # roll about S by ψ
    w = normalize(np.cross(S, u))               # completes RHS triad (u, w, S)

    # Fill rotation matrix columns according to which face is bound to S
    R = np.zeros((3,3))
    # We want R * f_body = S; we’ll assign columns so that column of the bound axis equals S
    if face_key.endswith("X"):
        sign = np.sign(f[0])  # +1 or -1
        R[:,0] = sign * S
        # Assign remaining axes (ensure right-handedness)
        R[:,1] = w
        R[:,2] = np.cross(R[:,0], R[:,1])
        # Re-orthonormalize last in case of tiny drift
        R[:,2] = normalize(R[:,2])
        R[:,1] = normalize(np.cross(R[:,2], R[:,0]))
    elif face_key.endswith("Y"):
        sign = np.sign(f[1])
        R[:,1] = sign * S
        R[:,2] = u
        R[:,0] = np.cross(R[:,1], R[:,2])
        R[:,0] = normalize(R[:,0])
        R[:,2] = normalize(np.cross(R[:,0], R[:,1]))
    else:  # "Z"
        sign = np.sign(f[2])
        R[:,2] = sign * S
        R[:,0] = u
        R[:,1] = np.cross(R[:,2], R[:,0])
        R[:,1] = normalize(R[:,1])
        R[:,0] = normalize(np.cross(R[:,1], R[:,2]))

    # Ensure det +1 (proper rotation); fix if needed
    if np.linalg.det(R) < 0:
        R[:,1] = -R[:,1]

    return R

def sun_pointing_quaternions(beta_deg,
                             incl_deg=97.4959,
                             raan_deg=0.0,
                             n_theta=1,
                             n_roll=36,
                             faces=("+X","-X","+Y","-Y","+Z","-Z")):
    """
    Generate quaternions (w,x,y,z) in ICRF for all attitudes where *one* chosen body face points at the Sun.
    - beta_deg : Sun-plane angle β (deg), S·ĥ = sinβ (β>0 means Sun is on ĥ side of the plane).
    - incl_deg : orbit inclination (deg).
    - raan_deg : right ascension of ascending node (deg).
    - n_theta  : samples of Sun azimuth around ĥ (0..2π).
    - n_roll   : samples of free roll ψ around S (0..2π).
    - faces    : iterable of face keys from {±X, ±Y, ±Z}.
    Returns: dict(face -> array of shape [n_theta*n_roll, 4]) of quaternions (w,x,y,z),
             and the corresponding (theta, roll) arrays for reference.
    """
    h_hat = orbit_normal_icrf(incl_deg, raan_deg)
    n_hat = node_direction_icrf(raan_deg)
    S_list, thetas = sun_vectors_for_beta(beta_deg, h_hat, n_hat, n_theta=n_theta)
    rolls = np.linspace(0.0, 2*np.pi, n_roll, endpoint=False)

    out = {}
    for face in faces:
        quats = []
        for S in S_list:
            for psi in rolls:
                R = complete_attitude_columns_for_face(S, face, psi)
                q = rotmat_to_quat(R)
                quats.append(q)
        out[face] = np.vstack(quats)
    return out, thetas, rolls

def save_stk_attitude(quats, n_theta, n_roll, dt=1.0, filename="stk_attitude.a"):
    N = quats.shape[0]
    times = np.arange(N) * dt
    with open(filename, "w") as f:
        f.write("""stk.v.12.0

BEGIN Attitude

    NumberOfAttitudePoints		 {}

    BlockingFactor		 20

    InterpolationOrder		 2

    CentralBody		 Earth
    
    ScenarioEpoch		 3 Aug 2026 00:00:00.000000

    CoordinateAxes		 ICRF
    
    AttitudeTimeQuaternions\n""".format(n_theta*n_roll))
        for t, q in zip(times, quats):
            w, x, y, z = q
            f.write(f"{t:.3f} {x:.8f} {y:.8f} {z:.8f} {w:.8f}\n")
            # f.write(f"{t:.1f} {w:.8f} {x:.8f} {y:.8f} {z:.8f}\n")
        f.write("END Attitude\n")

# ------------------------------
# Example usage
# ------------------------------
if __name__ == "__main__":
    beta_deg = 50.0          # <-- example: angle between Sun and orbital plane
    incl_deg = 97.4959       # given
    raan_deg = 227.056            # set if known; 0° means ascending node along +X
    n_theta = 1             # Sun azimuth samples around orbit normal
    n_roll = 36              # roll samples about the Sun vector
    faces = ("+X","-X","+Y","-Y","+Z","-Z")

    quats_by_face, thetas, rolls = sun_pointing_quaternions(
        beta_deg, incl_deg, raan_deg, n_theta, n_roll, faces
    )
    
    
    # Example: export +Z face
    qZ = quats_by_face["+Z"]
    df = pd.DataFrame(qZ, columns=["q0 (scalar)", "q1", "q2", "q3"])
    df.to_csv("SunPointingQuaternions_Zface.csv", index=False)
    print("Saved to SunPointingQuaternions_Zface.csv")

    # Example: save +Z face sequence
    save_stk_attitude(quats_by_face["+Z"], n_theta, n_roll, dt=10.0, filename="Zface_quaternions.a")
    print("Saved attitude file: Zface_quaternions.a")