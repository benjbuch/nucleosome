// ── Mode toggle ───────────────────────────────────────────────────────────

const strictToggle = document.getElementById('strict-toggle');
const modeRow      = document.getElementById('mode-row');

function isStrict() { return strictToggle.checked; }

strictToggle.addEventListener('change', () => {
  lastKey = null;
  inputEl.dispatchEvent(new Event('input'));
});

// Show/hide the variant filter toggle based on whether the current parse
// contains at least one substitution.  Hidden when irrelevant.
function updateToggleVisibility(parsed) {
  if (!parsed) { modeRow.hidden = true; return; }
  const hasVariant = (function walk(node) {
    if (node.variant) return true;
    for (const s of node.slots ?? []) if (walk(s)) return true;
    for (const m of node.modifications ?? []) if (m.variant) return true;
    for (const n of node.nucleosomes ?? []) if (walk(n)) return true;
    return false;
  })(parsed);
  modeRow.hidden = !hasVariant;
}

// ── Interaction ───────────────────────────────────────────────────────────

const inputEl      = document.getElementById('notation-input');
const statusEl     = document.getElementById('parse-status');
const resultsEl    = document.getElementById('results');
const debugEl      = document.getElementById('debug-json');
const onboardingEl = document.getElementById('onboarding');

// Skeleton HTML keeps the two-column layout stable when nothing is rendered.
const pristineResults = resultsEl.innerHTML;

// Auto-complete paired brackets: ( → (), [ → [], { → {}
// With a selection, wraps the selected text: select "H3K27M", press [ → [H3K27M].
// Without a selection, inserts an empty pair and places the cursor inside.
const BRACKET_PAIRS = { '(': ')', '[': ']', '{': '}' };
inputEl.addEventListener('keydown', e => {
  const close = BRACKET_PAIRS[e.key];
  if (!close) return;
  e.preventDefault();
  const start = inputEl.selectionStart;
  const end   = inputEl.selectionEnd;
  const selected = inputEl.value.slice(start, end);
  inputEl.value = inputEl.value.slice(0, start) + e.key + selected + close + inputEl.value.slice(end);
  if (selected.length > 0) {
    // Keep the wrapped text selected (inside the brackets)
    inputEl.setSelectionRange(start + 1, start + 1 + selected.length);
  } else {
    inputEl.setSelectionRange(start + 1, start + 1);
  }
  inputEl.dispatchEvent(new Event('input'));
});

inputEl.addEventListener('input', async () => {
  const raw = inputEl.value;

  if (!raw.trim()) {
    inputEl.className = '';
    statusEl.innerHTML = '';
    statusEl.className = '';
    resultsEl.innerHTML = pristineResults;
    onboardingEl.hidden = false;
    debugEl.textContent = '';
    lastKey = null;
    updateToggleVisibility(null);
    return;
  }

  const parsed = tryParse(raw);
  const info   = describe(parsed);
  debugEl.textContent = parsed ? JSON.stringify(parsed, null, 2) : '';
  updateToggleVisibility(parsed);

  if (!info) {
    inputEl.className = 'invalid';
    statusEl.innerHTML = 'unrecognised notation';
    statusEl.className = 'err';
    return;
  }

  inputEl.className = 'valid';
  statusEl.innerHTML = `<code class="status-canonical">${info.canonical}</code><span class="status-detail">${info.detail}</span>`;
  statusEl.className = 'ok';

  const key = info.canonical;
  if (key === lastKey) return;
  lastKey = key;

  // Hide onboarding once we have a valid parse
  onboardingEl.hidden = true;

  if (!db) {
    resultsEl.innerHTML = '<p class="loading">Initialising database&hellip;</p>';
    return;
  }

  if (info.type === 'array') {
    await renderArray(parsed.nucleosomes, resultsEl, key);
  } else if (info.type === 'nucleosome') {
    await renderArray([parsed], resultsEl, key);
  } else {
    // Standalone histone or family: wrap in a virtual nucleosome so the
    // two-column grid structure is identical to a nucleosome query.
    // The card header label is hidden; the diagram column stays empty.
    await renderArray([{ slots: [parsed] }], resultsEl, key, { dimHeader: true });
  }
});

// ── Boot ──────────────────────────────────────────────────────────────────

// Hide toggle until a variant is typed
updateToggleVisibility(null);

initDuckDB().catch(err => {
  statusEl.innerHTML = `Database error: ${err.message}`;
  statusEl.className = 'err';
});
