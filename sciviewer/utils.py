import numpy as np
import numpy.linalg as la

def angle_between(v1, v2):
    # Angle calculation A
#     amt = np.dot(v1, v2) / (la.norm(v1) * la.norm(v2))
#     if amt <= -1:
#         return -py5.PI
#     elif amt >= 1:
#         return 0
#     return np.arccos(amt)
    # Angle calculation B
    cosang = np.dot(v1, v2)
    sinang = la.norm(np.cross(v1, v2))
    return np.arctan2(sinang, cosang)