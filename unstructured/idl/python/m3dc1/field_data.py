from __future__ import annotations

from pathlib import Path

from .dimensions import dimensions


def field_data(name: str, itor: int = 1, filename: str | Path = "C1.h5") -> tuple[str, object]:
    """
    Python port of field_data.pro.
    Returns (symbol, units_dimension_vector).
    """
    n = name.lower()
    units = dimensions()

    if n in {"psi", "psi_i"}:
        return "psi", dimensions(b0=1, l0=1 + itor)
    if n in {"i", "i_i"}:
        return "I", dimensions(b0=1, l0=itor)
    if n in {"phi", "phi_i"}:
        return "phi", dimensions(v0=1, l0=1 - itor)
    if n in {"v", "v_i"}:
        return "V", dimensions(v0=1, l0=-itor)
    if n in {"chi", "chi_i"}:
        return "chi", dimensions(v0=1, l0=1 + 2 * itor)
    if n in {"eta", "eta_i"}:
        return "eta", dimensions(eta=1)
    if n in {"den", "den_i"}:
        return "$n_i$", dimensions(n0=1)
    if n in {"ne", "ne_i"}:
        return "$n_e$", dimensions(n0=1)
    if n in {"n_re", "n_re_i"}:
        return "$n_{RE}$", dimensions(n0=1)
    if n in {"p", "p_i"}:
        return "p", dimensions(p0=1)
    if n in {"pe", "pe_i"}:
        return "$p_e$", dimensions(p0=1)
    if n in {"te", "te_i"}:
        return "$T_e$", dimensions(temperature=1)
    if n in {"ti", "ti_i"}:
        return "$T_i$", dimensions(temperature=1)
    if n in {"sigma", "sigma_i"}:
        return "sigma", dimensions(n0=1, t0=-1)
    if n in {"force_phi", "force_phi_i", "pforce", "pforce_i"}:
        return "$F_p$", dimensions(p0=1, l0=-1)
    if n in {"heat_source", "heat_source_i"}:
        return "Q", dimensions(p0=1, t0=-1)
    if n == "kappa":
        return "kappa", dimensions(n0=1, l0=2, t0=-1)
    if n == "kappar":
        return "$kappa_{par}$", dimensions(n0=1, l0=2, t0=-1)
    if n == "denm":
        return "$D_{n}$", dimensions(l0=2, t0=-1)
    if n in {"visc", "visc_c", "visc_i", "visc_c_i"}:
        return "mu", dimensions(p0=1, t0=1)
    if n in {"jphi", "jphi_i"}:
        return "$j_{phi}$", dimensions(b0=1, l0=itor - 1)
    if n in {"vor", "vor_i"}:
        return "vor", dimensions(v0=1, l0=itor - 1)
    if n in {"com", "com_i"}:
        return "com", dimensions(v0=1, l0=-1)
    if n in {"torque_em", "torque_em_i"}:
        return "$tau_{EM}$", dimensions(p0=1)
    if n in {"e_r", "e_r_i"}:
        return "$E_R$", dimensions(potential=1, l0=-1)
    if n in {"e_phi", "e_phi_i"}:
        return "$E_{phi}$", dimensions(potential=1, l0=-1)
    if n in {"e_z", "e_z_i"}:
        return "$E_Z$", dimensions(potential=1, l0=-1)
    if n in {"e_par", "e_par_i"}:
        return "$E_{parallel}$", dimensions(potential=1, l0=-1)
    if n in {"potential", "potential_i"}:
        return "phi", dimensions(potential=1)
    if n == "frequency":
        return "omega", dimensions(t0=-1)
    if n == "kprad_sigma_i":
        return "$sigma_{KPRAD,i}$", dimensions(n0=1, t0=-1)
    if n == "kprad_sigma_e":
        return "$sigma_{KPRAD,e}$", dimensions(n0=1, t0=-1)
    if n == "kprad_totden":
        return "$n_{Z,KPRAD}$", dimensions(n0=1)
    if n == "kprad_rad":
        from .read_parameter import read_parameter

        version = float(read_parameter("version", filename=filename))
        if version < 22:
            return "$P_{KPRAD,total}$", dimensions(p0=1, t0=-1)
        return "$P_{KPRAD,line}$", dimensions(p0=1, t0=-1)
    if n == "kprad_brem":
        return "$P_{KPRAD,brem}$", dimensions(p0=1, t0=-1)
    if n == "kprad_ion":
        return "$P_{KPRAD,ion}$", dimensions(p0=1, t0=-1)
    if n == "kprad_reck":
        return "$P_{KPRAD,reck}$", dimensions(p0=1, t0=-1)
    if n == "kprad_recp":
        return "$P_{KPRAD,recp}$", dimensions(p0=1, t0=-1)
    if n == "rad_source":
        return "$Q_{rad}$", dimensions(p0=1, t0=-1)
    if n == "zeff":
        return "$Z_{eff}$", dimensions()

    if n.startswith("kprad_n"):
        from .read_parameter import read_parameter

        z = int(read_parameter("kprad_z", filename=filename))
        try:
            nz = int(name[8:10])
        except ValueError:
            nz = 0
        zstr = [
            "0",
            "H",
            "He",
            "Li",
            "Be",
            "B",
            "C",
            "N",
            "O",
            "F",
            "Ne",
            "Na",
            "Mg",
            "Al",
            "Si",
            "P",
            "S",
            "Cl",
            "Ar",
        ]
        nzstr = [
            "0",
            "I",
            "II",
            "III",
            "IV",
            "V",
            "VI",
            "VII",
            "VIII",
            "IX",
            "X",
            "XI",
            "XII",
            "XIII",
            "XIV",
            "XV",
            "XVI",
            "XVII",
            "XVIII",
            "XIX",
            "XX",
        ]
        z = max(0, min(z, len(zstr) - 1))
        nz = max(0, min(nz, len(nzstr) - 1))
        return f"{zstr[z]} {nzstr[nz]}", dimensions(n0=1)

    if n.startswith("kprad_particle_source"):
        from .read_parameter import read_parameter

        z = int(read_parameter("kprad_z", filename=filename))
        try:
            nz = int(name[22:24])
        except ValueError:
            nz = 0
        zstr = [
            "0",
            "H",
            "He",
            "Li",
            "Be",
            "B",
            "C",
            "N",
            "O",
            "F",
            "Ne",
            "Na",
            "Mg",
            "Al",
            "Si",
            "P",
            "S",
            "Cl",
            "Ar",
        ]
        nzstr = [
            "0",
            "I",
            "II",
            "III",
            "IV",
            "V",
            "VI",
            "VII",
            "VIII",
            "IX",
            "X",
            "XI",
            "XII",
            "XIII",
            "XIV",
            "XV",
            "XVI",
            "XVII",
            "XVIII",
            "XIX",
            "XX",
        ]
        z = max(0, min(z, len(zstr) - 1))
        nz = max(0, min(nz, len(nzstr) - 1))
        return f"{zstr[z]} {nzstr[nz]} Source", dimensions(n0=1, t0=-1)

    return name, units
