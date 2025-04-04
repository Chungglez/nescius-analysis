class Gas:
    def __init__(self, name: str, mass_fraction: float, molar_mass: float):
        self.name = name
        self.mass_fraction = mass_fraction
        self.M = molar_mass

        # Variables to be calculated later
        self.mol_fraction = None
        self.visc = None
        self.cp = None
        self.cond = None
        self.Pr = None
        self.rho = None

    def calc_cp(self, R, a, b, c, d, e, T):
        """
        Calculates the specific heat of the gas from a polynomial regression fit
        c_p/R = a + bT + cT^2 + dT^3 + eT^4

        :param R: gas constant  kJ/kg-K
        :param a: a
        :param b: b x 10^3      1/K
        :param c: c x 10^6      1/K^2
        :param d: d x 10^10     1/K^3
        :param e: e x 10^13
        :param T: gas temp      K
        """
        self.cp = R * (a + b * 1e-3 * T + c * 1e-6 * T ** 2 + d * 1e-10 * T ** 3 + e * 1e-13 * T ** 4)

    def calc_visc(self, a: float, b: float, c: float, T: float):
        self.visc = (a + b * T + c * T ** 2)*1e-7

    def calc_cond(self):
        self.cond = 1

    def calc_Pr(self):
        if None in [self.cp, self.visc, self.cond]:
            raise ValueError(
                "cp, viscosity, and thermal conductivity must all be defined before calculating Prandtl number.")
        self.Pr = self.cp * self.visc / self.cond

    def calc_mol_fraction(self, M_avg: float):
        self.mol_fraction = self.mass_fraction * M_avg / self.M


def mix_func(gas_i: Gas, gas_j: Gas) -> float:
    phi_ij_num = pow(1 + pow(gas_i.visc / gas_j.visc, 1 / 2) * pow(gas_j.M / gas_i.M, 1 / 4), 2)
    phi_ij_den = pow(8 * (1 + gas_i.M / gas_j.M), 1 / 2)
    return phi_ij_num / phi_ij_den
