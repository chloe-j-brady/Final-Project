// the bar chart: for each group, how many times more often than chance
// a drug's look-alikes share its group. Data comes from scores.json.
import * as d3 from "d3";

const ROW_HEIGHT = 38;      // space for one group: its name, then its bar
const BAR_HEIGHT = 12;
const MARGIN = { top: 22, right: 44, bottom: 44, left: 4 };
const BAR_COLOR = "#BEBAB1";   // bars of groups that aren't picked

let bars;   // kept so colorBars() can recolor them later

// draws the chart once. onClick(group) runs when a bar is clicked.
export function drawBarChart(scores, onClick) {
  const data = [...scores.groups].sort((a, b) => b.times_chance - a.times_chance);   // biggest first

  const container = document.getElementById("barchart");
  const width = container.clientWidth;
  const height = MARGIN.top + data.length * ROW_HEIGHT + MARGIN.bottom;

  const svg = d3.select(container).append("svg")
    .attr("width", width)
    .attr("height", height)
    .attr("role", "img")
    .attr("aria-label", "Bar chart: times more often than chance that look-alike drugs share a group, for each of 14 groups");

  // x: times more than chance -> pixels
  const x = d3.scaleLinear()
    .domain([0, Math.ceil(d3.max(data, (d) => d.times_chance))])
    .range([MARGIN.left, width - MARGIN.right]);

  // y: where each group's row starts
  const rowTop = (i) => MARGIN.top + i * ROW_HEIGHT;

  // light vertical grid lines, behind everything
  svg.append("g")
    .selectAll("line")
    .data(x.ticks(6))
    .join("line")
    .attr("x1", (t) => x(t)).attr("x2", (t) => x(t))
    .attr("y1", MARGIN.top - 6).attr("y2", height - MARGIN.bottom)
    .attr("stroke", "#E2DFD8");

  // dashed line at 1× = what pure chance would give
  // (drawn before the rows, so the labels sit on top of it)
  svg.append("line")
    .attr("x1", x(1)).attr("x2", x(1))
    .attr("y1", MARGIN.top - 6).attr("y2", height - MARGIN.bottom)
    .attr("stroke", "#141413")
    .attr("stroke-dasharray", "4 3");
  svg.append("text")
    .attr("class", "chance-label")
    .attr("x", x(1) + 4)
    .attr("y", MARGIN.top - 10)
    .text("chance");

  // one row per group
  const rows = svg.append("g")
    .selectAll("g")
    .data(data)
    .join("g")
    .attr("class", "bar-row")
    .attr("transform", (d, i) => `translate(0, ${rowTop(i)})`)
    .style("cursor", "pointer")
    .on("click", (event, d) => onClick(d.group))
    .on("mousemove", (event, d) => showTooltip(event, d))
    .on("mouseleave", hideTooltip);

  // an invisible box over the whole row, so the click target is bigger than the bar
  rows.append("rect")
    .attr("x", 0).attr("y", 0)
    .attr("width", width).attr("height", ROW_HEIGHT)
    .attr("fill", "transparent");

  // group name above the bar
  rows.append("text")
    .attr("class", "bar-label")
    .attr("x", MARGIN.left)
    .attr("y", 12)
    .text((d) => d.name);

  // the bar itself
  bars = rows.append("rect")
    .attr("x", x(0))
    .attr("y", 17)
    .attr("height", BAR_HEIGHT)
    .attr("width", (d) => x(d.times_chance) - x(0))
    .attr("rx", 3)
    .attr("fill", BAR_COLOR);

  // the value at the end of the bar
  rows.append("text")
    .attr("class", "bar-value")
    .attr("x", (d) => x(d.times_chance) + 6)
    .attr("y", 27)
    .text((d) => d.times_chance.toFixed(1) + "×");

  // x axis along the bottom, with its label
  svg.append("g")
    .attr("class", "axis")
    .attr("transform", `translate(0, ${height - MARGIN.bottom})`)
    .call(d3.axisBottom(x).ticks(6).tickFormat((t) => t + "×").tickSizeOuter(0));
  svg.append("text")
    .attr("class", "axis-title")
    .attr("x", (MARGIN.left + width - MARGIN.right) / 2)
    .attr("y", height - 8)
    .attr("text-anchor", "middle")
    .text("Times more often than chance");
}

// picked: the 3 picked groups (or null), colors: their 3 colors
export function colorBars(picked, colors) {
  bars.attr("fill", (d) => {
    const slot = picked.indexOf(d.group);
    return slot === -1 ? BAR_COLOR : colors[slot];
  });
}

// ---------- tooltip ----------
const tooltip = document.getElementById("bar-tooltip");

function showTooltip(event, d) {
  const pct = (v) => Math.round(v * 100) + "%";
  tooltip.innerHTML = `
    <div class="name">${d.name}</div>
    <div><b>${pct(d.real)}</b> of look-alikes share this group</div>
    <div class="muted">${pct(d.chance)} by chance · ${d.n_drugs} drugs</div>`;
  const box = tooltip.parentElement.getBoundingClientRect();
  tooltip.style.left = Math.min(event.clientX - box.left + 12, box.width - 220) + "px";
  tooltip.style.top = event.clientY - box.top + 14 + "px";
  tooltip.hidden = false;
}

function hideTooltip() {
  tooltip.hidden = true;
}