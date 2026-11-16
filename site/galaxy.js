// the 3D galaxy: every drug is drawn as its own tiny ball-and-stick molecule,
// placed where t-SNE or PCA puts it
import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { drawBarChart, colorBars } from "./barchart.js";
import { initPanel, showDrug } from "./compare.js";

const GROUP_NAMES = {
  A: "Alimentary tract & metabolism", B: "Blood", C: "Cardiovascular",
  D: "Dermatological", G: "Genito-urinary & sex hormones", H: "Hormones (non-sex)",
  J: "Anti-infectives", L: "Cancer & immune", M: "Musculo-skeletal",
  N: "Nervous system", P: "Antiparasitic", R: "Respiratory",
  S: "Sensory organs", V: "Various",
};

// only 3 groups can be compared at once, so only 3 colors are needed.
// these 3 stay easy to tell apart, including for colorblind viewers.
const SLOT_COLORS = ["#2a78d6", "#eb6834", "#1baf7a"];   // blue, orange, green
const START_GROUPS = ["J", "C", "N"];                     // anti-infectives, cardiovascular, nervous system
const OTHER_COLOR = "#BEBAB1";                            // molecules not in a picked group

// classic chemistry colors for "Color by atom", by atomic number
const ATOM_COLORS = { 6: "#4D4D4D", 8: "#E0453A", 7: "#3A6FE0", 16: "#E3B53A", 9: "#3FB14A", 17: "#3FB14A", 35: "#3FB14A", 53: "#3FB14A" };
const OTHER_ATOM = "#B07AB0";
const BOND_COLOR = "#9C998F";

// molecule size settings (molecule files are in angstroms; the galaxy is 20 units across)
const MOL_SCALE = 0.02;      // 1 angstrom = 0.02 galaxy units, so a drug is about 0.2-0.4 units wide
const ATOM_RADIUS = 0.4;     // angstroms
const BOND_RADIUS = 0.13;    // angstroms

// ---------- load the data made by the pipeline ----------
const [drugs, scores, molecules, comparisons] = await Promise.all([
  fetch("data/drugs.json").then((r) => r.json()),
  fetch("data/scores.json").then((r) => r.json()),
  fetch("data/molecules.json").then((r) => r.json()),
  fetch("data/comparisons.json").then((r) => r.json()),
]);

const indexById = {};   // drug id -> position in the list (for the look-alike lines)
drugs.forEach((d, i) => (indexById[d.id] = i));

// ---------- which groups are picked ----------
// picked[0..2] holds a group code or null. A group keeps its slot (and color) until removed.
const picked = [null, null, null];
let colorMode = "group";   // "group" or "atoms"

function slotOf(group) { return picked.indexOf(group); }
function drugSlot(d) {
  for (let s = 0; s < 3; s++) if (picked[s] && d.groups.includes(picked[s])) return s;
  return -1;
}
function drugColor(d) {
  const s = drugSlot(d);
  return s === -1 ? OTHER_COLOR : SLOT_COLORS[s];
}

// ---------- scene, camera, lights ----------
const container = document.getElementById("galaxy");
const scene = new THREE.Scene();
scene.background = new THREE.Color("#EFECE6");
scene.fog = new THREE.Fog("#EFECE6", 16, 36);   // far-away molecules fade a little, which helps the 3D read

const camera = new THREE.PerspectiveCamera(50, container.clientWidth / container.clientHeight, 0.01, 1000);
camera.position.set(0, 0, 22);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(window.devicePixelRatio);
renderer.setSize(container.clientWidth, container.clientHeight);
container.appendChild(renderer.domElement);

// soft light from everywhere plus one light from above, so the balls look round
scene.add(new THREE.HemisphereLight("#ffffff", "#c9c4b8", 1.6));
const sun = new THREE.DirectionalLight("#ffffff", 1.4);
sun.position.set(5, 10, 8);
scene.add(sun);

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;   // smooth slow-down after you let go
controls.autoRotate = true;      // drift slowly until someone touches it
controls.autoRotateSpeed = 0.4;
controls.minDistance = 0.3;      // allow zooming right up to one molecule
controls.addEventListener("start", () => (controls.autoRotate = false));

// ---------- galaxy positions ----------
// t-SNE and PCA use very different number ranges, so scale each one
// to fit inside a sphere of radius 10, centered on 0
function scaled(layout) {
  const pts = drugs.map((d) => d[layout]);
  const mean = [0, 1, 2].map((k) => pts.reduce((sum, p) => sum + p[k], 0) / pts.length);
  const centered = pts.map((p) => p.map((v, k) => v - mean[k]));
  const maxDist = Math.max(...centered.map((p) => Math.hypot(...p)));
  return centered.map((p) => p.map((v) => (v / maxDist) * 10));
}
const layouts = { tsne: scaled("tsne"), pca: scaled("pca") };

// centers[i*3..i*3+2] = where drug i's molecule sits right now
const centers = new Float32Array(drugs.length * 3);
function centerOf(i) {
  return new THREE.Vector3(centers[i * 3], centers[i * 3 + 1], centers[i * 3 + 2]);
}

// ---------- build every molecule ----------
// 55,000 atoms and 59,000 bonds would be far too slow as separate objects,
// so ALL atoms share one "InstancedMesh" (one sphere shape, drawn many times)
// and ALL bonds share another (one cylinder shape, drawn many times)
const atomCount = molecules.reduce((n, m) => n + m.e.length, 0);
const bondCount = molecules.reduce((n, m) => n + m.b.length / 2, 0);

const atoms = new THREE.InstancedMesh(
  new THREE.IcosahedronGeometry(1, 1),   // a simple sphere
  new THREE.MeshStandardMaterial({ roughness: 0.55 }),
  atomCount);
const bonds = new THREE.InstancedMesh(
  new THREE.CylinderGeometry(1, 1, 1, 6, 1, true),   // a thin tube, 1 unit long
  new THREE.MeshStandardMaterial({ roughness: 0.6 }),
  bondCount);
scene.add(atoms, bonds);

// for each atom/bond: which drug it belongs to, and where it sits inside that molecule
const atomDrug = new Int32Array(atomCount);
const atomLocal = new Float32Array(atomCount * 3);
const atomElement = new Int32Array(atomCount);
const bondDrug = new Int32Array(bondCount);
const bondLocal = new Float32Array(bondCount * 3);   // the middle of the bond

const m4 = new THREE.Matrix4();
const q = new THREE.Quaternion();
const up = new THREE.Vector3(0, 1, 0);   // the cylinder's built-in direction
let a = 0, b = 0;
molecules.forEach((mol, i) => {
  const first = a;   // number of this molecule's first atom
  for (let k = 0; k < mol.e.length; k++, a++) {
    atomDrug[a] = i;
    atomElement[a] = mol.e[k];
    atomLocal.set([mol.p[k * 3] * MOL_SCALE, mol.p[k * 3 + 1] * MOL_SCALE, mol.p[k * 3 + 2] * MOL_SCALE], a * 3);
    // size never changes, so set it once now; only the position will move later
    m4.makeScale(ATOM_RADIUS * MOL_SCALE, ATOM_RADIUS * MOL_SCALE, ATOM_RADIUS * MOL_SCALE);
    atoms.setMatrixAt(a, m4);
  }
  for (let k = 0; k < mol.b.length; k += 2, b++) {
    const p1 = new THREE.Vector3().fromArray(atomLocal, (first + mol.b[k]) * 3);
    const p2 = new THREE.Vector3().fromArray(atomLocal, (first + mol.b[k + 1]) * 3);
    const dir = p2.clone().sub(p1);
    bondDrug[b] = i;
    bondLocal.set(p1.clone().add(p2).multiplyScalar(0.5).toArray(), b * 3);
    // turn the tube to point from atom 1 to atom 2, and stretch it to the bond length
    q.setFromUnitVectors(up, dir.clone().normalize());
    m4.compose(new THREE.Vector3(), q, new THREE.Vector3(BOND_RADIUS * MOL_SCALE, dir.length(), BOND_RADIUS * MOL_SCALE));
    bonds.setMatrixAt(b, m4);
  }
});

// move every atom and bond to (its molecule's center + its spot inside the molecule).
// a matrix stores position in slots 12, 13, 14, so only those need changing.
function placeMolecules() {
  const am = atoms.instanceMatrix.array;
  for (let k = 0; k < atomCount; k++) {
    const c = atomDrug[k] * 3;
    am[k * 16 + 12] = centers[c] + atomLocal[k * 3];
    am[k * 16 + 13] = centers[c + 1] + atomLocal[k * 3 + 1];
    am[k * 16 + 14] = centers[c + 2] + atomLocal[k * 3 + 2];
  }
  const bm = bonds.instanceMatrix.array;
  for (let k = 0; k < bondCount; k++) {
    const c = bondDrug[k] * 3;
    bm[k * 16 + 12] = centers[c] + bondLocal[k * 3];
    bm[k * 16 + 13] = centers[c + 1] + bondLocal[k * 3 + 1];
    bm[k * 16 + 14] = centers[c + 2] + bondLocal[k * 3 + 2];
  }
  atoms.instanceMatrix.needsUpdate = true;
  bonds.instanceMatrix.needsUpdate = true;
  hoverPoints.geometry.attributes.position.needsUpdate = true;
}

// ---------- colors ----------
const tmp = new THREE.Color();
function refreshColors() {
  for (let k = 0; k < atomCount; k++) {
    const hex = colorMode === "atoms"
      ? (ATOM_COLORS[atomElement[k]] || OTHER_ATOM)
      : drugColor(drugs[atomDrug[k]]);
    atoms.setColorAt(k, tmp.set(hex));
  }
  for (let k = 0; k < bondCount; k++) {
    const hex = colorMode === "atoms" ? BOND_COLOR : drugColor(drugs[bondDrug[k]]);
    bonds.setColorAt(k, tmp.set(hex));
  }
  if (atoms.instanceColor) atoms.instanceColor.needsUpdate = true;   // (empty if a file had no atoms/bonds)
  if (bonds.instanceColor) bonds.instanceColor.needsUpdate = true;
}

function setColorMode(mode) {
  colorMode = mode;
  document.getElementById("btn-group").classList.toggle("active", mode === "group");
  document.getElementById("btn-atoms").classList.toggle("active", mode === "atoms");
  document.getElementById("atom-key").hidden = mode !== "atoms";
  refreshColors();
}
document.getElementById("btn-group").addEventListener("click", () => setColorMode("group"));
document.getElementById("btn-atoms").addEventListener("click", () => setColorMode("atoms"));

// ---------- hovering ----------
// checking the mouse against 55,000 atoms would be slow, so hovering uses one
// invisible point at the center of each molecule instead
const hoverGeometry = new THREE.BufferGeometry();
hoverGeometry.setAttribute("position", new THREE.BufferAttribute(centers, 3));
const hoverPoints = new THREE.Points(hoverGeometry, new THREE.PointsMaterial({ visible: false }));
scene.add(hoverPoints);

// ---------- smooth animation between t-SNE and PCA ----------
let animation = null;   // { from, to, start } while the molecules are moving
const ANIMATION_MS = 1400;

function easeInOut(t) {
  // starts slow, speeds up, ends slow
  return t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2;
}

function setLayout(name, instant = false) {
  if (instant) {
    layouts[name].forEach((p, i) => centers.set(p, i * 3));
    placeMolecules();
    finishMove();
  } else {
    animation = { from: Float32Array.from(centers), to: layouts[name], start: performance.now() };
  }
  document.getElementById("btn-tsne").classList.toggle("active", name === "tsne");
  document.getElementById("btn-pca").classList.toggle("active", name === "pca");
  document.getElementById("axis-note").textContent = name === "pca"
    ? "Axes: " + scores.pca_axis_labels.join(" · ")
    : "t-SNE: axes have no units, only closeness matters";
}

function stepAnimation(now) {
  const t = Math.min((now - animation.start) / ANIMATION_MS, 1);   // 0 at the start, 1 at the end
  const e = easeInOut(t);
  for (let i = 0; i < drugs.length; i++) {
    for (let k = 0; k < 3; k++) {
      const from = animation.from[i * 3 + k];
      centers[i * 3 + k] = from + (animation.to[i][k] - from) * e;
    }
  }
  placeMolecules();
  if (t === 1) {
    animation = null;
    finishMove();
  }
}

// after moving, three.js needs to re-measure where everything is (for hovering and drawing)
function finishMove() {
  atoms.computeBoundingSphere();
  bonds.computeBoundingSphere();
  hoverGeometry.computeBoundingSphere();
}

document.getElementById("btn-tsne").addEventListener("click", () => setLayout("tsne"));
document.getElementById("btn-pca").addEventListener("click", () => setLayout("pca"));

// ---------- bar chart: click a bar to pick its group (up to 3) ----------
drawBarChart(scores, (g) => toggleGroup(g));

function toggleGroup(g) {
  const s = slotOf(g);
  if (s !== -1) {
    picked[s] = null;                     // already picked: remove it
  } else {
    const free = picked.indexOf(null);
    if (free === -1) return;              // already 3 picked: remove one first
    picked[free] = g;
  }
  if (colorMode === "atoms") setColorMode("group");   // picking a group means you want group colors
  else refreshColors();
  colorBars(picked, SLOT_COLORS);                      // the same colors in the bar chart
}

// ---------- look-alike lines ----------
// hover (or click) a molecule and thin lines connect it to its 5 closest look-alikes
const linePositions = new Float32Array(5 * 2 * 3);   // 5 lines, 2 ends, x y z
const lineGeometry = new THREE.BufferGeometry();
lineGeometry.setAttribute("position", new THREE.BufferAttribute(linePositions, 3));
const lines = new THREE.LineSegments(lineGeometry, new THREE.LineBasicMaterial({ color: "#141413" }));
lines.visible = false;
scene.add(lines);
let hoveredIndex = null;
let selectedIndex = null;   // the clicked molecule, shown in the panel on the right

function updateLines() {
  // show the hovered molecule's lines, or the clicked one's when not hovering anything
  const i = hoveredIndex ?? selectedIndex;
  if (i === null) { lines.visible = false; return; }
  const start = centerOf(i);
  drugs[i].lookalikes.forEach((la, k) => {
    const end = centerOf(indexById[la.id]);
    linePositions.set([start.x, start.y, start.z, end.x, end.y, end.z], k * 6);
  });
  lineGeometry.attributes.position.needsUpdate = true;
  lineGeometry.computeBoundingSphere();
  lines.visible = true;
}

// ---------- hover tooltip ----------
const raycaster = new THREE.Raycaster();
const mouse = new THREE.Vector2();
const tooltip = document.getElementById("tooltip");

function groupText(d) {
  return d.groups.length ? d.groups.map((g) => GROUP_NAMES[g]).join(", ") : "Unclassified";
}

// which molecule (if any) is under the mouse? returns its number, or null
function moleculeAt(event) {
  // convert the mouse position to the -1..1 range three.js uses
  const rect = container.getBoundingClientRect();
  mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
  mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;

  // how close the mouse must be grows with zoom distance, so it feels the same at any zoom
  raycaster.params.Points.threshold = Math.max(0.03, controls.getDistance() * 0.009);
  raycaster.setFromCamera(mouse, camera);
  const hits = raycaster.intersectObject(hoverPoints);
  return hits.length ? hits[0].index : null;
}

container.addEventListener("pointermove", (event) => {
  const rect = container.getBoundingClientRect();
  const found = moleculeAt(event);
  if (found === null) {
    tooltip.hidden = true;
    hoveredIndex = null;
    updateLines();
    return;
  }

  hoveredIndex = found;
  const d = drugs[hoveredIndex];
  const rows = d.lookalikes.slice(0, 4).map((la) => {
    const other = drugs[indexById[la.id]];
    return `<li><span class="dot" style="background:${drugColor(other)}"></span>${other.name}<span class="sim">${la.sim.toFixed(2)}</span></li>`;
  }).join("");

  tooltip.innerHTML = `
    <div class="name">${d.name}</div>
    <div class="groups">${groupText(d)}</div>
    <div class="lookalike-title">Closest look-alikes</div>
    <ul>${rows}</ul>`;

  // keep the tooltip inside the box
  const x = event.clientX - rect.left;
  const y = event.clientY - rect.top;
  tooltip.style.left = (x > rect.width - 270 ? x - 250 : x + 20) + "px";
  tooltip.style.top = Math.min(y + 12, rect.height - 170) + "px";
  tooltip.hidden = false;
  updateLines();
});
container.addEventListener("pointerleave", () => {
  tooltip.hidden = true;
  hoveredIndex = null;
  updateLines();
});

// ---------- click a molecule to open it in the panel ----------
// dragging to rotate also presses and releases the mouse, so only count it
// as a click if the mouse barely moved in between
let downX = 0, downY = 0;
container.addEventListener("pointerdown", (event) => { downX = event.clientX; downY = event.clientY; });
container.addEventListener("pointerup", (event) => {
  if (Math.hypot(event.clientX - downX, event.clientY - downY) > 5) return;   // it was a drag
  const found = moleculeAt(event);
  if (found === null) return;
  selectedIndex = found;
  updateLines();
  showDrug(found);
});

// the panel's "Back to overview" button clears the selection
initPanel(drugs, comparisons, GROUP_NAMES, () => {
  selectedIndex = null;
  updateLines();
});

// keep the picture the right size if the window changes
window.addEventListener("resize", () => {
  camera.aspect = container.clientWidth / container.clientHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(container.clientWidth, container.clientHeight);
});

// ---------- start ----------
setLayout("tsne", true);
START_GROUPS.forEach((g) => toggleGroup(g));

function animate(now) {
  requestAnimationFrame(animate);
  controls.update();
  if (animation) {
    stepAnimation(now);
    updateLines();
  }
  renderer.render(scene, camera);
}
requestAnimationFrame(animate);