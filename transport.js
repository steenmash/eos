import { getComponent } from "./thermo-data.js";

const UNIVERSAL_R = 8.314462618;

function mixtureAverages(components, composition) {
  let MW = 0;
  let Tc = 0;
  let Pc = 0;
  let Vc = 0;
  let omega = 0;
  let dipole = 0;
  let kappa = 0;
  for (let i = 0; i < components.length; i += 1) {
    const comp = getComponent(components[i]);
    const y = composition[i];
    MW += y * comp.molarMass;
    Tc += y * comp.criticalTemperature;
    Pc += y * comp.criticalPressure;
    Vc += y * comp.criticalVolume;
    omega += y * comp.acentricFactor;
    dipole += y * comp.dipoleMoment;
    kappa += y * comp.associationFactor;
  }
  return { MW, Tc, Pc, Vc, omega, dipole, kappa };
}

function zeroDensityViscosity(T, mixture) {
  const { MW, Tc, Vc, omega, dipole, kappa } = mixture;
  const epsilon = Tc / 1.2593;
  const muReduced =
    dipole === 0
      ? 0
      : (131.3 * dipole) /
        Math.sqrt(epsilon * Math.pow(Vc, 2 / 3));
  const Fc = 1 - 0.2756 * omega + 0.059035 * Math.pow(muReduced, 4) + kappa;
  const Tstar = 1.2593 * (T / Tc);
  const Omega =
    1.16145 / Math.pow(Tstar, 0.14874) +
    0.52487 * Math.exp(-0.7732 * Tstar) +
    2.16178 * Math.exp(-2.43787 * Tstar);
  const numerator = 40.785 * Fc * Math.sqrt(MW * T);
  const denom = Math.pow(Vc, 2 / 3) * Omega;
  const muMicroPoise = numerator / denom;
  return muMicroPoise * 1e-7; // Pa·s
}

function leeGonzalezEakinViscosity(T, density, mixture) {
  const { MW } = mixture;
  const temperatureRankine = T * (9 / 5);
  const gasDensity = density / 1000; // kg/m³ -> g/cm³
  const X = 3.448 + 986.4 / (temperatureRankine + 0.01) + 0.01009 * MW;
  const Y = 2.447 - 0.2224 * X;
  const viscosityCp = 1e-4 * Math.exp(X * Math.pow(gasDensity, Y));
  return viscosityCp * 1e-3; // Pa·s
}

export function transportProperties({
  temperature,
  pressure,
  density,
  composition,
  components,
  compressibility,
  heatCapacityMass,
}) {
  const mixture = mixtureAverages(components, composition);
  const mu0 = zeroDensityViscosity(temperature, mixture);
  const VcMix = mixture.Vc * 1e-6; // m³/mol
  const MWkg = mixture.MW / 1000; // kg/mol
  const rhoStar = density > 0 ? (density * VcMix) / MWkg : 0;
  let dynamicViscosity;
  if (compressibility > 0.85) {
    dynamicViscosity = leeGonzalezEakinViscosity(temperature, density, mixture);
  } else {
    const correction = 1 + 0.25 * Math.pow(Math.max(rhoStar, 0), 1.2);
    dynamicViscosity = mu0 * correction;
  }
  const kinematicViscosity = density > 0 ? dynamicViscosity / density : 0;
  const Rspec = UNIVERSAL_R / MWkg;
  const thermalConductivity = dynamicViscosity * (heatCapacityMass + 1.25 * Rspec);
  return {
    dynamicViscosity,
    kinematicViscosity,
    thermalConductivity,
  };
}
