"""
Write a PLUTO cooltable.dat (T [K], Lambda [erg cm^3 s^-1], n_e n_H norm.)
for metallicity Z (solar units) using cooling.py.

    python make_cooltable.py 0.3 cooltable.dat
"""
import sys
import cooling as C

if __name__ == "__main__":
    Z = float(sys.argv[1])
    out = sys.argv[2] if len(sys.argv) > 2 else "cooltable.dat"
    C.write_pluto_table(out, *C.make_table(Z))
    print("wrote", out)
