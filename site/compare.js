// the panel on the right: shows the overall answer until you click a molecule,
// then that drug's story: do its look-alikes treat the same thing, and what do they share?
// molecules are drawn by RDKit.js (loaded in index.html), the same chemistry library as the Python side

const SHARED_COLOR = [0.62, 0.85, 0.79];   // light teal = the piece both molecules share
const DIFFERENT_COLOR = [0.97, 0.79, 0.68]; // light orange = the rest

const panel = document.getElementById("drug-panel");
const overviewHTML = panel.innerHTML;   // whatever index.html puts there first = the overview

let rdkit = null;
const rdkitReady = window.initRDKitModule().then((module) => (rdkit = module));

let drugs, comparisons, groupNames, onClose;
let current = null;   // the drug being shown, or null for the overview

export function initPanel(allDrugs, allComparisons, names, closeCallback) {
  drugs = allDrugs;
  comparisons = allComparisons;
  groupNames = names;
  onClose = closeCallback;
}

export function showOverview() {
  current = null;
  panel.innerHTML = overviewHTML;
  onClose();
}

// i = which drug, k = which of its look-alikes to compare against
export async function showDrug(i, k = null) {
  current = i;
  await rdkitReady;
  if (current !== i) return;   // another molecule was clicked while RDKit was loading

  const d = drugs[i];
  const rows = comparisons.pairs[i];

  // start with the closest look-alike that isn't identical (a 1.00 match is usually
  // a mirror-image form of the same molecule, which is all "shared" and not very interesting)
  if (k === null) {
    k = rows.findIndex((r) => r.shared < r.total);
    if (k === -1) k = 0;
  }
  const pair = rows[k];
  const identical = pair.shared === pair.total;
  const other = drugs[indexOf(pair.id)];

  panel.innerHTML = `
    <button class="back">← Back to overview</button>
    <div class="drug-name">${d.name}</div>
    <div class="drug-groups">${groupList(d)}</div>
    <p class="verdict">${verdict(d, rows)}</p>

    <div class="compare">
      <figure>${drawMolecule(i, pair.a, pair.ab)}<figcaption>${d.name}</figcaption></figure>
      <figure>${drawMolecule(indexOf(pair.id), pair.b, pair.bb)}<figcaption>${other.name}</figcaption></figure>
    </div>
    <div class="key"><span><i class="shared"></i>Shared part</span><span><i class="different"></i>Different</span></div>
    <p class="score">Similarity <b>${(pair.shared / pair.total).toFixed(2)}</b>
      = ${pair.shared} shared building blocks ÷ ${pair.total} in total</p>
    ${identical ? `<p class="note">A perfect 1.00 usually means a mirror-image form of the same molecule. The fingerprint can't tell left-handed from right-handed versions apart.</p>` : ""}

    <p class="list-title">Its 5 closest look-alikes <span>(click one to compare)</span></p>
    <ul class="la-list">
      ${rows.map((r, n) => lookalikeRow(d, r, n, n === k)).join("")}
    </ul>`;

  panel.querySelector(".back").addEventListener("click", showOverview);
  panel.querySelectorAll(".la-list button").forEach((button) => {
    button.addEventListener("click", () => showDrug(i, Number(button.dataset.k)));
  });
}

// ---------- pieces of the panel ----------

function indexOf(id) {
  return drugs.findIndex((d) => d.id === id);
}

function sharesGroup(a, b) {
  return a.groups.some((g) => b.groups.includes(g));
}

function groupList(d) {
  return d.groups.length ? d.groups.map((g) => groupNames[g]).join(", ") : "Unclassified";
}

// one sentence: how many of its look-alikes treat the same thing?
function verdict(d, rows) {
  const others = rows.map((r) => drugs[indexOf(r.id)]);

  if (d.groups.length === 0) {
    // no label, so it's left out of the test, but its look-alikes hint at what it might do
    const counts = {};
    others.forEach((o) => o.groups.forEach((g) => (counts[g] = (counts[g] || 0) + 1)));
    const top = Object.keys(counts).sort((a, b) => counts[b] - counts[a])[0];
    return top
      ? `Not labeled, so it's left out of the test. But <b>${counts[top]} of its 5</b> look-alikes are <b>${groupNames[top].toLowerCase()}</b> drugs, a hint at what it might do.`
      : "Not labeled, and neither are its look-alikes.";
  }

  const matches = others.filter((o) => sharesGroup(d, o)).length;
  const what = d.groups.length === 1 ? `also ${groupNames[d.groups[0]].toLowerCase()} drugs` : "share one of its groups";
  return `<b>${matches} of its 5</b> closest look-alikes are ${what}.`;
}

function lookalikeRow(d, r, n, isCurrent) {
  const o = drugs[indexOf(r.id)];
  const same = d.groups.length && sharesGroup(d, o);
  return `<li><button data-k="${n}" class="${isCurrent ? "current" : ""}">
    <span class="la-name">${o.name}<small>${groupList(o)}</small></span>
    <span class="la-match">${same ? "✓ same" : ""}</span>
    <span class="la-sim">${(r.shared / r.total).toFixed(2)}</span>
  </button></li>`;
}

// draw a molecule with the shared atoms/bonds in teal and everything else in orange
function drawMolecule(i, sharedAtoms, sharedBonds) {
  const mol = rdkit.get_mol(comparisons.smiles[i]);
  const atomCount = mol.get_num_atoms();
  const bondCount = JSON.parse(mol.get_json()).molecules[0].bonds.length;

  const atomColors = {}, bondColors = {};
  const sharedA = new Set(sharedAtoms), sharedB = new Set(sharedBonds);
  for (let a = 0; a < atomCount; a++) atomColors[a] = sharedA.has(a) ? SHARED_COLOR : DIFFERENT_COLOR;
  for (let b = 0; b < bondCount; b++) bondColors[b] = sharedB.has(b) ? SHARED_COLOR : DIFFERENT_COLOR;

  const svg = mol.get_svg_with_highlights(JSON.stringify({
    width: 220, height: 190,
    atoms: Object.keys(atomColors).map(Number),
    bonds: Object.keys(bondColors).map(Number),
    highlightAtomColors: atomColors,
    highlightBondColors: bondColors,
    highlightRadius: 0.4,                // bigger colored circles on the atoms
    highlightBondWidthMultiplier: 16,    // thicker colored bands along the bonds
    bondLineWidth: 1.5,
    backgroundColour: [1, 1, 1, 0],   // transparent
  }));
  mol.delete();   // RDKit.js needs molecules freed by hand
  return svg;
}