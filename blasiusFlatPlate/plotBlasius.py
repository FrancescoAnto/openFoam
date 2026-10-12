#!/usr/bin/env python3
"""Confronto OpenFOAM vs soluzione di Blasius (lastra piana laminare).

Legge direttamente constant/polyMesh e il campo U dell'ultimo tempo
(formato ascii), quindi funziona con qualsiasi versione di OpenFOAM senza
function object. Richiede numpy; matplotlib opzionale per il grafico.

Uso:  python3 plotBlasius.py [tempo]
"""
import os
import re
import sys
import numpy as np

U_INF, NU = 1.0, 1e-5
STATIONS = [0.25, 0.50, 0.75]


# ----------------------------------------------------------- lettura OpenFOAM
def _list_body(path):
    """Ritorna (n, testo) della prima lista OpenFOAM 'n ( ... )' nel file."""
    txt = open(path).read()
    txt = re.sub(r"/\*.*?\*/", "", txt, flags=re.S)
    txt = re.sub(r"//.*", "", txt)
    m = re.search(r"\n\s*(\d+)\s*\n\s*\(", txt)
    if m is None:
        raise ValueError(f"lista non trovata in {path} (serve writeFormat ascii)")
    return int(m.group(1)), txt[m.end():]


def read_vectors(path):
    n, body = _list_body(path)
    nums = re.findall(r"\(([^()]*)\)", body)[:n]
    return np.array([[float(v) for v in t.split()] for t in nums])


def read_labels(path):
    n, body = _list_body(path)
    return np.array(body[:body.index(")")].split()[:n], dtype=int)


def read_faces(path):
    n, body = _list_body(path)
    return [list(map(int, f.split())) for f in re.findall(r"\d+\(([^()]*)\)", body)[:n]]


def cell_centres(mesh="constant/polyMesh"):
    pts = read_vectors(os.path.join(mesh, "points"))
    faces = read_faces(os.path.join(mesh, "faces"))
    own = read_labels(os.path.join(mesh, "owner"))
    nei = read_labels(os.path.join(mesh, "neighbour"))
    ncell = own.max() + 1
    verts = [set() for _ in range(ncell)]
    for fi, c in enumerate(own):
        verts[c].update(faces[fi])
    for fi, c in enumerate(nei):
        verts[c].update(faces[fi])
    return np.array([pts[list(v)].mean(axis=0) for v in verts])


def read_U(time):
    txt = open(os.path.join(time, "U")).read()
    m = re.search(r"internalField\s+uniform\s+\(([^)]*)\)", txt)
    if m:
        raise ValueError("U uniforme: il caso non e' stato ancora risolto")
    n, body = _list_body(os.path.join(time, "U"))
    return read_vectors(os.path.join(time, "U"))


def latest_time():
    t = [d for d in os.listdir(".") if os.path.isdir(d) and re.fullmatch(r"[0-9.eE+-]+", d)
         and os.path.exists(os.path.join(d, "U"))]
    t = [d for d in t if float(d) > 0]
    if not t:
        sys.exit("Nessun tempo risolto trovato: lancia prima ./Allrun")
    return max(t, key=float)


# ------------------------------------------------------------------- Blasius
def blasius(eta_max=10.0, n=2001):
    """f''' + 0.5 f f'' = 0, f(0)=f'(0)=0, f'(inf)=1  (shooting + RK4)."""
    eta = np.linspace(0.0, eta_max, n)
    h = eta[1] - eta[0]
    rhs = lambda y: np.array([y[1], y[2], -0.5 * y[0] * y[2]])

    def integrate(fpp0):
        y = np.array([0.0, 0.0, fpp0]); out = [y]
        for _ in range(n - 1):
            k1 = rhs(y); k2 = rhs(y + h / 2 * k1); k3 = rhs(y + h / 2 * k2); k4 = rhs(y + h * k3)
            y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4); out.append(y)
        return np.array(out)

    a, b = 0.1, 1.0
    for _ in range(60):
        m = 0.5 * (a + b)
        a, b = (a, m) if integrate(m)[-1, 1] > 1.0 else (m, b)
    return eta, integrate(0.5 * (a + b))[:, 1], 0.5 * (a + b)


# ---------------------------------------------------------------------- main
def main():
    time = sys.argv[1] if len(sys.argv) > 1 else latest_time()
    print(f"Tempo analizzato: {time}")
    C = cell_centres()
    U = read_U(time)
    eta_b, fp_b, fpp0 = blasius()
    print(f"Blasius: f''(0) = {fpp0:.5f} (teorico 0.33206)\n")

    xs_unique = np.unique(np.round(C[:, 0], 9))
    results = []
    print(f"{'x [m]':>6} {'Re_x':>8} {'max|du/U|':>10} {'Cf CFD':>9} {'Cf Blas.':>9} "
          f"{'err Cf':>7} {'d99 CFD':>9} {'d99 Blas.':>9}")
    for x in STATIONS:
        xc = xs_unique[np.argmin(abs(xs_unique - x))]          # colonna di celle piu' vicina
        k = np.where(abs(C[:, 0] - xc) < 1e-9)[0]
        k = k[np.argsort(C[k, 1])]
        y, u = C[k, 1], U[k, 0]
        Rex = U_INF * xc / NU
        eta = y * np.sqrt(U_INF / (NU * xc))
        m = eta < 8
        err = np.max(np.abs(u[m] / U_INF - np.interp(eta[m], eta_b, fp_b)))
        cf = 2 * NU * (u[0] / y[0]) / U_INF**2                  # tau_w ~ mu*u1/y1
        cf_b = 0.664 / np.sqrt(Rex)
        d99 = np.interp(0.99, u[: np.argmax(u >= 0.99) + 1], y[: np.argmax(u >= 0.99) + 1])
        d99_b = 4.91 * xc / np.sqrt(Rex)
        print(f"{xc:6.3f} {Rex:8.0f} {err:10.4f} {cf:9.5f} {cf_b:9.5f} "
              f"{100 * (cf / cf_b - 1):+6.1f}% {1e3 * d99:7.2f}mm {1e3 * d99_b:7.2f}mm")
        results.append((xc, eta, u))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib non disponibile: niente grafico.")
        return
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(fp_b, eta_b, "k-", lw=2, label="Blasius")
    for mk, (xc, eta, u) in zip("osd", results):
        ax.plot(u / U_INF, eta, mk, ms=4, mfc="none", label=f"OpenFOAM x = {xc:.2f} m")
    ax.set(xlabel="u / U", ylabel=r"$\eta = y\,\sqrt{U/(\nu x)}$", ylim=(0, 8),
           title="Lastra piana: OpenFOAM vs Blasius")
    ax.grid(alpha=0.3); ax.legend(); fig.tight_layout()
    fig.savefig("blasius_comparison.png", dpi=150)
    print("\nGrafico salvato in blasius_comparison.png")


if __name__ == "__main__":
    main()
