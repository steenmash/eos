import { COMPONENTS, getComponent } from "./thermo-data.js";
import { flashCalculation } from "./eos.js";

const componentRows = document.getElementById("component-rows");
const addComponentBtn = document.getElementById("add-component");
const calculateBtn = document.getElementById("calculate");
const resultsContainer = document.getElementById("results");
const temperatureInput = document.getElementById("temperature");
const pressureInput = document.getElementById("pressure");
const massFlowInput = document.getElementById("mass-flow");

function createComponentRow(initialId = "", value = "") {
  const row = document.createElement("div");
  row.className = "component-row";

  const select = document.createElement("select");
  const emptyOption = document.createElement("option");
  emptyOption.value = "";
  emptyOption.textContent = "Компонент";
  select.appendChild(emptyOption);
  COMPONENTS.forEach((comp) => {
    const option = document.createElement("option");
    option.value = comp.id;
    option.textContent = `${comp.id} — ${comp.name}`;
    if (initialId === comp.id) option.selected = true;
    select.appendChild(option);
  });

  const input = document.createElement("input");
  input.type = "number";
  input.min = "0";
  input.step = "0.0001";
  input.placeholder = "Мольная доля";
  input.value = value;

  const removeBtn = document.createElement("button");
  removeBtn.type = "button";
  removeBtn.className = "remove-row";
  removeBtn.textContent = "×";
  removeBtn.addEventListener("click", () => {
    if (componentRows.children.length > 1) {
      row.remove();
    }
  });

  row.appendChild(select);
  row.appendChild(input);
  row.appendChild(removeBtn);
  componentRows.appendChild(row);
}

function initRows() {
  if (componentRows.children.length === 0) {
    createComponentRow("C1", 0.8);
    createComponentRow("C2", 0.1);
    createComponentRow("C3", 0.1);
  }
}

function parseInputs() {
  const rows = Array.from(componentRows.querySelectorAll(".component-row"));
  const componentIds = [];
  const fractions = [];
  rows.forEach((row) => {
    const select = row.querySelector("select");
    const input = row.querySelector("input");
    const id = select.value;
    const value = parseFloat(input.value);
    if (!id || Number.isNaN(value) || value < 0) {
      return;
    }
    componentIds.push(id);
    fractions.push(value);
  });
  if (componentIds.length === 0) {
    throw new Error("Необходимо задать хотя бы один компонент и его концентрацию.");
  }
  const sum = fractions.reduce((acc, val) => acc + val, 0);
  if (sum <= 0) {
    throw new Error("Сумма концентраций должна быть больше нуля.");
  }
  const composition = fractions.map((val) => val / sum);
  const temperatureC = parseFloat(temperatureInput.value);
  const pressure = parseFloat(pressureInput.value);
  const massFlow = parseFloat(massFlowInput.value);
  if (Number.isNaN(temperatureC) || Number.isNaN(pressure) || Number.isNaN(massFlow)) {
    throw new Error("Проверьте корректность условий (T, P, расход).");
  }
  if (pressure <= 0) {
    throw new Error("Давление должно быть положительным.");
  }
  if (massFlow < 0) {
    throw new Error("Массовый расход не может быть отрицательным.");
  }
  const temperature = temperatureC + 273.15;
  return {
    componentIds,
    composition,
    temperature,
    pressure,
    massFlow,
  };
}

function mixtureMolarMass(componentIds, composition) {
  let mw = 0;
  for (let i = 0; i < componentIds.length; i += 1) {
    const comp = getComponent(componentIds[i]);
    mw += composition[i] * comp.molarMass;
  }
  return mw;
}

function renderCompositionList(componentIds, composition) {
  const list = document.createElement("table");
  const tbody = document.createElement("tbody");
  componentIds.forEach((id, index) => {
    const row = document.createElement("tr");
    const nameCell = document.createElement("th");
    nameCell.textContent = id;
    const valueCell = document.createElement("td");
    valueCell.textContent = `${(composition[index] * 100).toFixed(3)} %`;
    row.appendChild(nameCell);
    row.appendChild(valueCell);
    tbody.appendChild(row);
  });
  list.appendChild(tbody);
  return list;
}

function formatNumber(value, digits = 3) {
  if (!Number.isFinite(value)) return "—";
  return Number(value).toFixed(digits);
}

function phaseHeading(name) {
  switch (name) {
    case "liquid":
      return "Жидкая фаза";
    case "vapor":
      return "Паровая фаза";
    default:
      return name;
  }
}

function renderPhaseBlock(name, data, flows) {
  const block = document.createElement("div");
  block.className = "phase-block";
  const title = document.createElement("h3");
  title.textContent = `${phaseHeading(name)} (${formatNumber(flows.molarFraction * 100, 2)} % мол.)`;
  block.appendChild(title);

  const table = document.createElement("table");
  const tbody = document.createElement("tbody");
  const addRow = (label, value, unit = "") => {
    const row = document.createElement("tr");
    const key = document.createElement("th");
    key.textContent = label;
    const val = document.createElement("td");
    val.textContent = `${value}${unit ? ` ${unit}` : ""}`;
    row.appendChild(key);
    row.appendChild(val);
    tbody.appendChild(row);
  };

  addRow("Компрессибилити Z", formatNumber(data.compressibility, 5));
  addRow("Плотность", formatNumber(data.density, 3), "кг/м³");
  addRow(
    "Динамическая вязкость",
    formatNumber(data.transport.dynamicViscosity * 1000, 4),
    "мПа·с"
  );
  addRow(
    "Кинематическая вязкость",
    formatNumber(data.transport.kinematicViscosity * 1e6, 3),
    "мм²/с"
  );
  addRow(
    "Теплопроводность",
    formatNumber(data.transport.thermalConductivity, 4),
    "Вт/(м·К)"
  );
  addRow(
    "Cp (мольная)",
    formatNumber(data.heatCapacityMolar, 2),
    "Дж/(моль·К)"
  );
  addRow(
    "Cp (массовая)",
    formatNumber(data.heatCapacityMass / 1000, 4),
    "кДж/(кг·К)"
  );
  addRow(
    "Энтальпия",
    formatNumber(data.enthalpy / 1000, 3),
    "кДж/моль"
  );
  addRow(
    "Мольный расход",
    formatNumber(flows.molar / 1000, 4),
    "кмоль/ч"
  );
  addRow("Массовый расход", formatNumber(flows.mass, 3), "кг/ч");
  addRow("Объёмный расход", formatNumber(flows.volume, 4), "м³/ч");

  table.appendChild(tbody);
  block.appendChild(table);

  const compositionTitle = document.createElement("h4");
  compositionTitle.textContent = "Состав (мольные доли)";
  block.appendChild(compositionTitle);
  block.appendChild(renderCompositionList(data.componentIds, data.composition));

  return block;
}

function renderResults(input, result) {
  resultsContainer.innerHTML = "";
  const summary = document.createElement("div");
  summary.className = "warning";
  summary.innerHTML = `
    <strong>Баланс:</strong> паровая доля = ${(result.vaporFraction * 100).toFixed(2)} %;
    жидкая доля = ${(result.liquidFraction * 100).toFixed(2)} %. Все доли указаны в мольном выражении.
  `;
  resultsContainer.appendChild(summary);

  const overallMW = mixtureMolarMass(input.componentIds, input.composition);
  const totalMolarFlow = overallMW > 0 ? (input.massFlow / (overallMW / 1000)) : 0;

  Object.entries(result.phases).forEach(([phaseName, data]) => {
    const fraction =
      phaseName === "vapor" ? result.vaporFraction : result.liquidFraction;
    const phaseMW = mixtureMolarMass(data.componentIds, data.composition);
    const molarFlow = totalMolarFlow * fraction;
    const massFlow = (molarFlow * phaseMW) / 1000;
    const volumeFlow = data.density > 0 ? massFlow / data.density : 0;
    const phaseBlock = renderPhaseBlock(phaseName, data, {
      molarFraction: fraction,
      molar: molarFlow,
      mass: massFlow,
      volume: volumeFlow,
    });
    resultsContainer.appendChild(phaseBlock);
  });
}

function renderError(message) {
  resultsContainer.innerHTML = "";
  const block = document.createElement("div");
  block.className = "warning";
  block.textContent = message;
  resultsContainer.appendChild(block);
}

addComponentBtn.addEventListener("click", () => {
  createComponentRow();
});

calculateBtn.addEventListener("click", () => {
  try {
    const input = parseInputs();
    const result = flashCalculation({
      componentIds: input.componentIds,
      composition: input.composition,
      temperature: input.temperature,
      pressure: input.pressure,
    });
    renderResults(input, result);
  } catch (error) {
    renderError(error.message);
  }
});

initRows();
