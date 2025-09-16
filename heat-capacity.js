(function (global) {
const { getComponent } = global.ThermoData;

const REFERENCE_T = 298.15; // K

function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

function ensureRange(component, T) {
  const { idealCp } = component;
  const [min, max] = idealCp.range;
  return clamp(T, min, max);
}

function shomateCp(component, T) {
  const { idealCp } = component;
  const [A, B, C, D, E] = idealCp.coefficients;
  const t = T / 1000;
  return (
    A +
    B * t +
    C * t * t +
    D * t * t * t +
    (E !== 0 ? E / (t * t) : 0)
  );
}

function shomateIntegral(component, T1, T2) {
  const { idealCp } = component;
  const [A, B, C, D, E] = idealCp.coefficients;
  const toScaled = (T) => T / 1000;
  const integral = (T) => {
    const t = toScaled(T);
    const termA = A * T;
    const termB = (B / 2000) * T * T;
    const termC = (C / 3000000) * T * T * T;
    const termD = (D / 4000000000) * T * T * T * T;
    const termE = E === 0 ? 0 : -E * 1000000 / T;
    return termA + termB + termC + termD + termE;
  };
  return integral(T2) - integral(T1);
}

function idealGasHeatCapacity(componentId, temperature) {
  const component = getComponent(componentId);
  if (!component) throw new Error(`Unknown component: ${componentId}`);
  const T = ensureRange(component, temperature);
  return shomateCp(component, T);
}

function idealGasHeatCapacityMixture(components, composition, temperature) {
  let cp = 0;
  for (let i = 0; i < components.length; i += 1) {
    const comp = components[i];
    const frac = composition[i];
    cp += frac * idealGasHeatCapacity(comp.id, temperature);
  }
  return cp;
}

function idealGasEnthalpy(componentId, temperature) {
  const component = getComponent(componentId);
  if (!component) throw new Error(`Unknown component: ${componentId}`);
  const T = ensureRange(component, temperature);
  const reference = clamp(REFERENCE_T, component.idealCp.range[0], component.idealCp.range[1]);
  return shomateIntegral(component, reference, T);
}

function mixtureIdealEnthalpy(components, composition, temperature) {
  let h = 0;
  for (let i = 0; i < components.length; i += 1) {
    const comp = components[i];
    const frac = composition[i];
    h += frac * idealGasEnthalpy(comp.id, temperature);
  }
  return h;
}

const api = {
  idealGasHeatCapacity,
  idealGasHeatCapacityMixture,
  idealGasEnthalpy,
  mixtureIdealEnthalpy,
};

if (typeof module !== "undefined" && module.exports) {
  module.exports = api;
}

global.HeatCapacity = api;
})(typeof globalThis !== "undefined" ? globalThis : window);
