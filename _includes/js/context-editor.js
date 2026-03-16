// ── Context editor ──────────────────────────────────────────────────────────

const CONTEXT_STORAGE_KEY = 'nucleosome-context';

// Mutable overrides — keyed by family or isoform token.
// Each entry can have: { uniprot_id, protein_name, species, taxon_id, modifications }
// Empty by default (no overrides — SQL resolves to human canonical).
let contextOverrides = {};

// Default YAML shown when no overrides are set — just a commented template.
const CONTEXT_TEMPLATE = `# Override materialization context per family/isoform.
# Unspecified fields fall back to human canonical (taxon 9606).
# Example:
#   H3:
#     taxon_id: 8355          # Xenopus laevis
#   H2A:
#     modifications:
#     - {position: 15, modification: ac}
#   H31:
#     modifications:
#     - {position:  96, variant: A}
#     - {position: 110, variant: A}
`;

// Serialize current overrides to YAML for the editor.
function overridesToYaml(ov) {
  if (!ov || Object.keys(ov).length === 0) return CONTEXT_TEMPLATE;
  // Clean empty modifications arrays before serializing
  const clean = {};
  for (const [key, entry] of Object.entries(ov)) {
    const c = {};
    for (const [f, v] of Object.entries(entry)) {
      if (f === 'modifications' && (!v || v.length === 0)) continue;
      c[f] = v;
    }
    if (Object.keys(c).length > 0) clean[key] = c;
  }
  if (Object.keys(clean).length === 0) return CONTEXT_TEMPLATE;
  return jsyaml.dump(clean, { lineWidth: -1, quotingType: "'", forceQuotes: false });
}

// Trigger a full re-render of current input.
function reRender() {
  lastKey = null;
  cachedSvg = null;
  clearContextCache();
  const el = document.getElementById('notation-input');
  if (el) el.dispatchEvent(new Event('input'));
}

// ── Apply / Reset ────────────────────────────────────────────────────────

function applyContext() {
  const textarea = document.getElementById('context-yaml');
  const errorEl  = document.getElementById('context-error');
  if (!textarea) return;

  try {
    const parsed = jsyaml.load(textarea.value);
    if (parsed == null) {
      // Empty or all comments — clear overrides
      contextOverrides = {};
      localStorage.removeItem(CONTEXT_STORAGE_KEY);
      errorEl.textContent = 'Context reset to defaults.';
      errorEl.className = 'context-applied';
  
      reRender();
      return;
    }
    if (typeof parsed !== 'object') {
      errorEl.textContent = 'YAML must be a mapping of family tokens to override entries.';
      errorEl.className = 'context-error';
      return;
    }
    contextOverrides = parsed;
    localStorage.setItem(CONTEXT_STORAGE_KEY, textarea.value);
    errorEl.textContent = 'Context applied.';
    errorEl.className = 'context-applied';

    reRender();
  } catch (e) {
    errorEl.textContent = `YAML error: ${e.message}`;
    errorEl.className = 'context-error';
  }
}

function resetContext() {
  contextOverrides = {};
  localStorage.removeItem(CONTEXT_STORAGE_KEY);
  const textarea = document.getElementById('context-yaml');
  const errorEl  = document.getElementById('context-error');
  if (textarea) textarea.value = CONTEXT_TEMPLATE;
  if (errorEl)  errorEl.textContent = '';
  reRender();
}

// ── Initialise ───────────────────────────────────────────────────────────

function initContextEditor() {
  const textarea = document.getElementById('context-yaml');
  if (!textarea) return;

  // Restore from localStorage if available
  const saved = localStorage.getItem(CONTEXT_STORAGE_KEY);
  if (saved) {
    try {
      const parsed = jsyaml.load(saved);
      if (parsed && typeof parsed === 'object' && Object.keys(parsed).length > 0) {
        contextOverrides = parsed;
        textarea.value = saved;
      } else {
        // All-comment YAML or empty object — treat as no overrides
        localStorage.removeItem(CONTEXT_STORAGE_KEY);
        textarea.value = CONTEXT_TEMPLATE;
      }
    } catch {
      // Corrupted storage — fall back to defaults
      localStorage.removeItem(CONTEXT_STORAGE_KEY);
      textarea.value = CONTEXT_TEMPLATE;
    }
  } else {
    textarea.value = CONTEXT_TEMPLATE;
  }

  document.getElementById('context-apply')?.addEventListener('click', applyContext);
  document.getElementById('context-reset')?.addEventListener('click', resetContext);
}

initContextEditor();
