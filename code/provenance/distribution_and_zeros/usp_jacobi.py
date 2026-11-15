"""34.78: USp(2N) vážená Λ(1)^k přesně jako Jacobiho soubor (β = 2) v x = cos θ: váha (1 − x)^(k+½)(1 + x)^½.
Jednobodová hustota K_N(x, x), pravděpodobnost mezery E(t) = P(θ_min > t) = det(I − K) na (cos t, 1] (Bornemann, Gauss–Legendre),
E[x₁] = (N/π) ∫₀^π E(t) dt (konvence laboratoře x₁ = θ_min·N/π)."""
import numpy as np
from scipy.special import eval_jacobi, gammaln, roots_legendre
def _norm(n, a, b):   # ‖P_n^(a,b)‖² na [−1, 1] s vahou (1−x)^a (1+x)^b
    return np.exp((a + b + 1) * np.log(2) + gammaln(n + a + 1) + gammaln(n + b + 1) - np.log(2 * n + a + b + 1) - gammaln(n + a + b + 1) - gammaln(n + 1))
def phi(N, k, x):     # ortonormální funkce √w(x)·p_n(x), n < N; tvar (N, len(x))
    a, b = k + 0.5, 0.5; w = (1 - x) ** a * (1 + x) ** b
    return np.array([eval_jacobi(n, a, b, x) / np.sqrt(_norm(n, a, b)) for n in range(N)]) * np.sqrt(w)
def hustota_theta(N, k, th):  # hustota vlastních úhlů v θ ∈ [0, π] (integrál = N)
    x = np.cos(th); F = phi(N, k, x); return np.sum(F * F, axis=0) * np.sin(th)
def mezera(N, k, t, m=60):    # P(žádný úhel v [0, t]) = det(I − K) na x ∈ (cos t, 1]
    u, wq = roots_legendre(m); lo = np.cos(t); x = lo + (1 - lo) * (u + 1) / 2; wq = wq * (1 - lo) / 2
    F = phi(N, k, x) * np.sqrt(wq); K = F.T @ F
    return float(np.linalg.det(np.eye(m) - K))
def ex1(N, k, n=400):         # E[x₁] = (N/π)∫ E(t) dt
    t = np.linspace(0, np.pi, n + 1); E = np.array([mezera(N, k, tt) for tt in t]); 
    return float(N / np.pi * np.trapezoid(E, t)), t, E
if __name__ == '__main__':
    import json
    ref = json.load(open('/home/claude/br76/data/lh_ref.json'))['usp']
    for N in range(1, 7):
        m0, t, E = ex1(N, 0); th = np.linspace(1e-6, np.pi - 1e-6, 20001)
        print(N, 'E[x1] přesně %.5f, DPP %.5f' % (m0, ref[str(N)]['mx1']), ' ∫hustota k=0: %.6f, k=1: %.6f, k=2: %.6f' % tuple(np.trapezoid(hustota_theta(N, k, th), th) for k in (0, 1, 2)), ' E[x1] k=1 %.5f k=2 %.5f' % (ex1(N, 1)[0], ex1(N, 2)[0]))
