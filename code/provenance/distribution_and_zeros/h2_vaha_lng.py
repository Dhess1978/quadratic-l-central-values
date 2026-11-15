import math, numpy as np
def lngamma_c(z):
    s = np.zeros_like(z)
    for k in range(12): s = s - np.log(z + k)
    w = z + 12; B = [1/12, -1/360, 1/1260, -1/1680, 1/1188, -691/360360, 1/156, -3617/122400]
    t = (w - 0.5) * np.log(w) - w + 0.5 * math.log(2 * math.pi); iw = 1 / w; p = iw
    for b in B: t = t + b * p; p = p * iw * iw
    return s + t
