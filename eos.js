(function (global) {
const { GAS_CONSTANT, getComponent, binaryInteraction } = global.ThermoData;
const {
  idealGasHeatCapacityMixture,
  mixtureIdealEnthalpy,
} = global.HeatCapacity;
const { transportProperties } = global.Transport;

const SQRT2 = Math.SQRT2;
const R_SI = GAS_CONSTANT * 100; // J·L⁻¹ -> convert bar·L to J·mol⁻¹·K⁻¹ when multiplied by T

function wilsonK(component, T, P) {
  const Tc = component.criticalTemperature;
  const Pc = component.criticalPressure;
  const omega = component.acentricFactor;
  return (
    (Pc / P) *
    Math.exp(5.373 * (1 + omega) * (1 - Tc / T))
  );
}

function computePureParameters(component, T) {
  const Tc = component.criticalTemperature;
  const Pc = component.criticalPressure;
  const omega = component.acentricFactor;
  const Tr = Math.sqrt(T / Tc);
  const m = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
  const alpha = Math.pow(1 + m * (1 - Tr), 2);
  const a =
    0.45724 *
    ((GAS_CONSTANT * GAS_CONSTANT * Tc * Tc) / Pc) *
    alpha;
  const b = 0.0778 * ((GAS_CONSTANT * Tc) / Pc);
  const dalpha_dT =
    -m * (1 + m * (1 - Tr)) / (Math.sqrt(T * Tc));
  const da_dT =
    0.45724 *
    ((GAS_CONSTANT * GAS_CONSTANT * Tc * Tc) / Pc) *
    dalpha_dT;
  return { a, b, da_dT };
}

function cubicRoots(A, B) {
  const c2 = -1 + B;
  const c1 = A - 3 * B * B - 2 * B;
  const c0 = -A * B + B * B + B * B * B;
  const a = 1;
  const b = c2;
  const c = c1;
  const d = c0;
  const p = (3 * a * c - b * b) / (3 * a * a);
  const q =
    (27 * a * a * d - 9 * a * b * c + 2 * b * b * b) /
    (27 * a * a * a);
  const discriminant = (q * q) / 4 + (p * p * p) / 27;
  const roots = [];
  if (discriminant > 0) {
    const sqrtDisc = Math.sqrt(discriminant);
    const u = Math.cbrt(-q / 2 + sqrtDisc);
    const v = Math.cbrt(-q / 2 - sqrtDisc);
    roots.push(u + v - b / (3 * a));
  } else if (Math.abs(discriminant) < 1e-12) {
    if (Math.abs(q) < 1e-12) {
      roots.push(-b / (3 * a));
    } else {
      const u = Math.cbrt(-q / 2);
      roots.push(2 * u - b / (3 * a));
      roots.push(-u - b / (3 * a));
    }
  } else {
    const r = Math.sqrt(-p * p * p / 27);
    const phi = Math.acos(-q / (2 * Math.sqrt(-p * p * p / 27)));
    const m = 2 * Math.sqrt(-p / 3);
    roots.push(m * Math.cos(phi / 3) - b / (3 * a));
    roots.push(
      m * Math.cos((phi + 2 * Math.PI) / 3) - b / (3 * a)
    );
    roots.push(
      m * Math.cos((phi + 4 * Math.PI) / 3) - b / (3 * a)
    );
  }
  return roots.filter((root) => Number.isFinite(root));
}

function mixtureMolarMass(components, composition) {
  let mw = 0;
  for (let i = 0; i < components.length; i += 1) {
    mw += composition[i] * components[i].molarMass;
  }
  return mw;
}

function buildMixtureState(componentIds, composition, T, P) {
  const n = componentIds.length;
  const components = componentIds.map((id) => getComponent(id));
  const pure = components.map((comp) => computePureParameters(comp, T));
  const aij = Array.from({ length: n }, () => new Array(n).fill(0));
  const daij = Array.from({ length: n }, () => new Array(n).fill(0));
  const sumA = new Array(n).fill(0);
  let aMix = 0;
  let daMix_dT = 0;
  let bMix = 0;
  for (let i = 0; i < n; i += 1) {
    bMix += composition[i] * pure[i].b;
  }
  for (let i = 0; i < n; i += 1) {
    for (let j = 0; j < n; j += 1) {
      const kij = binaryInteraction(componentIds[i], componentIds[j]);
      const aijVal =
        (1 - kij) * Math.sqrt(pure[i].a * pure[j].a);
      aij[i][j] = aijVal;
      const term =
        0.5 * aijVal *
        ((pure[i].da_dT / pure[i].a) + (pure[j].da_dT / pure[j].a));
      daij[i][j] = term;
      aMix += composition[i] * composition[j] * aijVal;
      daMix_dT += composition[i] * composition[j] * term;
      sumA[i] += composition[j] * aijVal;
    }
  }
  const A = (aMix * P) / (Math.pow(GAS_CONSTANT * T, 2));
  const B = (bMix * P) / (GAS_CONSTANT * T);
  const roots = cubicRoots(A, B);
  const Zliquid = Math.min(...roots);
  const Zvapor = Math.max(...roots);
  const calcPhi = (Z) => {
    const logTerm = Math.log(
      (Z + (1 + SQRT2) * B) / (Z + (1 - SQRT2) * B)
    );
    const phi = new Array(n).fill(0);
    for (let i = 0; i < n; i += 1) {
      const term1 = (pure[i].b / bMix) * (Z - 1);
      const term2 = -Math.log(Math.max(Z - B, 1e-12));
      const term3 =
        (A / (2 * SQRT2 * B)) *
        (2 * sumA[i] / aMix - pure[i].b / bMix) *
        logTerm;
      phi[i] = Math.exp(term1 + term2 - term3);
    }
    return { phi, logTerm };
  };
  const liquidPhi = calcPhi(Zliquid);
  const vaporPhi = calcPhi(Zvapor);
  return {
    componentData: components,
    pure,
    composition,
    T,
    P,
    aMix,
    bMix,
    daMix_dT,
    A,
    B,
    Z: { liquid: Zliquid, vapor: Zvapor },
    phi: {
      liquid: liquidPhi.phi,
      vapor: vaporPhi.phi,
    },
    logTerm: {
      liquid: liquidPhi.logTerm,
      vapor: vaporPhi.logTerm,
    },
  };
}

function residualEnthalpy(state, phase, T) {
  const Z = state.Z[phase];
  const logTerm = state.logTerm[phase];
  const base = GAS_CONSTANT * T * (Z - 1);
  const part =
    ((T * state.daMix_dT - state.aMix) /
      (2 * SQRT2 * state.bMix)) *
    logTerm;
  const enthalpyBarL = base + part;
  return enthalpyBarL * 100; // convert (bar·L)/mol to J/mol
}

function totalEnthalpy(state, phase, composition, T) {
  const ideal = mixtureIdealEnthalpy(state.componentData, composition, T);
  const residual = residualEnthalpy(state, phase, T);
  return ideal + residual;
}

function heatCapacity(state, componentIds, composition, T, P, phase) {
  const delta = 0.5; // K
  const plusState = buildMixtureState(componentIds, composition, T + delta, P);
  const minusState = buildMixtureState(componentIds, composition, T - delta, P);
  const hPlus = totalEnthalpy(plusState, phase, composition, T + delta);
  const hMinus = totalEnthalpy(minusState, phase, composition, T - delta);
  return (hPlus - hMinus) / (2 * delta);
}

function solveRachfordRice(K, z) {
  const evalF = (V) => {
    let sum = 0;
    for (let i = 0; i < K.length; i += 1) {
      const denom = 1 + V * (K[i] - 1);
      if (denom <= 0) return Number.NaN;
      sum += (z[i] * (K[i] - 1)) / denom;
    }
    return sum;
  };
  const f0 = evalF(0);
  const f1 = evalF(1);
  if (Number.isNaN(f0) || Number.isNaN(f1)) {
    return { status: "invalid" };
  }
  if (f0 < 0 && f1 < 0) {
    return { status: "liquid" };
  }
  if (f0 > 0 && f1 > 0) {
    return { status: "vapor" };
  }
  let lower = 0;
  let upper = 1;
  let mid = 0.5;
  for (let iter = 0; iter < 100; iter += 1) {
    mid = 0.5 * (lower + upper);
    const fmid = evalF(mid);
    if (Math.abs(fmid) < 1e-10) break;
    if (fmid > 0) {
      lower = mid;
    } else {
      upper = mid;
    }
  }
  return { status: "two-phase", vaporFraction: mid };
}

function normaliseComposition(vector) {
  const sum = vector.reduce((acc, val) => acc + val, 0);
  return vector.map((val) => (sum === 0 ? 0 : val / sum));
}

function preparePhaseData(componentIds, composition, T, P, phase) {
  const state = buildMixtureState(componentIds, composition, T, P);
  const Z = state.Z[phase];
  const totalH = totalEnthalpy(state, phase, composition, T);
  const cpMolar = heatCapacity(state, componentIds, composition, T, P, phase);
  const cpMass = cpMolar / (mixtureMolarMass(state.componentData, composition) / 1000);
  const density =
    (P / (Z * GAS_CONSTANT * T)) *
    mixtureMolarMass(state.componentData, composition);
  const transport = transportProperties({
    temperature: T,
    pressure: P,
    density,
    composition,
    components: componentIds,
    compressibility: Z,
    heatCapacityMass: cpMass,
  });
  return {
    state,
    composition,
    componentIds,
    temperature: T,
    pressure: P,
    compressibility: Z,
    enthalpy: totalH,
    heatCapacityMolar: cpMolar,
    heatCapacityMass: cpMass,
    density,
    transport,
  };
}

function splitMixture(componentIds, z, K) {
  const solution = solveRachfordRice(K, z);
  if (solution.status === "liquid") {
    return { type: "liquid" };
  }
  if (solution.status === "vapor") {
    return { type: "vapor" };
  }
  if (solution.status !== "two-phase") {
    return { type: "invalid" };
  }
  const V = solution.vaporFraction;
  const x = new Array(z.length);
  const y = new Array(z.length);
  for (let i = 0; i < z.length; i += 1) {
    const denom = 1 + V * (K[i] - 1);
    x[i] = z[i] / denom;
    y[i] = K[i] * x[i];
  }
  const xNorm = normaliseComposition(x);
  const yNorm = normaliseComposition(y);
  return { type: "two-phase", vaporFraction: V, x: xNorm, y: yNorm };
}

function flashCalculation({ componentIds, composition, temperature, pressure }) {
  const T = temperature;
  const P = pressure;
  const components = componentIds.map((id) => getComponent(id));
  let K = components.map((comp) => Math.max(1e-6, wilsonK(comp, T, P)));
  let split = splitMixture(componentIds, composition, K);
  if (split.type !== "two-phase") {
    const phase = split.type === "liquid" ? "liquid" : "vapor";
    const phaseData = preparePhaseData(componentIds, composition, T, P, phase);
    return {
      vaporFraction: phase === "vapor" ? 1 : 0,
      liquidFraction: phase === "liquid" ? 1 : 0,
      phases: {
        [phase]: phaseData,
      },
    };
  }
  for (let iter = 0; iter < 100; iter += 1) {
    split = splitMixture(componentIds, composition, K);
    if (split.type !== "two-phase") break;
    const liquidState = buildMixtureState(componentIds, split.x, T, P);
    const vaporState = buildMixtureState(componentIds, split.y, T, P);
    const phiL = liquidState.phi.liquid;
    const phiV = vaporState.phi.vapor;
    const newK = new Array(K.length);
    let maxChange = 0;
    for (let i = 0; i < K.length; i += 1) {
      newK[i] = Math.max(1e-6, K[i] * (phiL[i] / phiV[i]));
      maxChange = Math.max(maxChange, Math.abs(newK[i] - K[i]));
    }
    K = newK;
    if (maxChange < 1e-7) {
      break;
    }
  }
  split = splitMixture(componentIds, composition, K);
  if (split.type !== "two-phase") {
    const phase = split.type === "liquid" ? "liquid" : "vapor";
    const phaseData = preparePhaseData(componentIds, composition, T, P, phase);
    return {
      vaporFraction: phase === "vapor" ? 1 : 0,
      liquidFraction: phase === "liquid" ? 1 : 0,
      phases: {
        [phase]: phaseData,
      },
    };
  }
  const liquidData = preparePhaseData(componentIds, split.x, T, P, "liquid");
  const vaporData = preparePhaseData(componentIds, split.y, T, P, "vapor");
  return {
    vaporFraction: split.vaporFraction,
    liquidFraction: 1 - split.vaporFraction,
    phases: {
      liquid: liquidData,
      vapor: vaporData,
    },
  };
}

const api = {
  flashCalculation,
};

if (typeof module !== "undefined" && module.exports) {
  module.exports = api;
}

global.EOS = api;
})(typeof globalThis !== "undefined" ? globalThis : window);
