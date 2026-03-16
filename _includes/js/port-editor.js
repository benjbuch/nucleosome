// ── Port-to context editor ──────────────────────────────────────────────────

const PORT_STORAGE_KEY = 'nucleosome-port';

// Mutable overrides for the port-to context — same shape as contextOverrides.
let portOverrides = {};

// Whether porting is active (port context has been applied with real entries).
function isPortActive() {
  return Object.keys(portOverrides).length > 0;
}

const PORT_TEMPLATE = `# Port-to context: same format as the main context.
# Define the target system for porting.
# Example — port to Xenopus:
#   H3:
#     taxon_id: 8355
#   H4:
#     taxon_id: 8355
`;

function portOverridesToYaml(ov) {
  if (!ov || Object.keys(ov).length === 0) return PORT_TEMPLATE;
  const clean = {};
  for (const [key, entry] of Object.entries(ov)) {
    const c = {};
    for (const [f, v] of Object.entries(entry)) {
      if (f === 'modifications' && (!v || v.length === 0)) continue;
      c[f] = v;
    }
    if (Object.keys(c).length > 0) clean[key] = c;
  }
  if (Object.keys(clean).length === 0) return PORT_TEMPLATE;
  return jsyaml.dump(clean, { lineWidth: -1, quotingType: "'", forceQuotes: false });
}

// Resolve a port-to context for a given family/isoform.
// Uses portOverrides instead of contextOverrides.
async function resolvePortContext(family, isoform) {
  const lookupKey = isoform ?? family;
  const userOv    = portOverrides[lookupKey] ?? portOverrides[family] ?? {};

  // If no port override applies to this family, fall back to the main context
  // (no diff to show for this slot).
  if (Object.keys(userOv).length === 0) return null;

  const taxon = userOv.taxon_id ?? 9606;

  const conn = await db.connect();
  try {
    const sql = `
      WITH pick AS (
        SELECT DISTINCT
          family, isoform, uniprot_id, protein_name,
          species, taxon_id, n_gene_loci, identity
        FROM read_parquet('${PROTEINS_URL}')
        WHERE family = $1
          AND ($2 IS NULL OR isoform = $2)
          AND ($3 IS NULL OR uniprot_id = $3)
          AND ($4 IS NULL OR protein_name = $4)
          AND ($5 IS NULL OR species = $5)
          AND taxon_id = $6
        ORDER BY n_gene_loci DESC, identity DESC
        LIMIT 1
      ),
      seq AS (
        SELECT
          p.uniprot_id,
          STRING_AGG(p.residue, '' ORDER BY p.protein_position) AS sequence
        FROM read_parquet('${PROTEINS_URL}') p
        JOIN pick ON p.uniprot_id = pick.uniprot_id
        WHERE p.protein_position >= 1
        GROUP BY p.uniprot_id
      ),
      mw AS (
        SELECT
          p.uniprot_id,
          SUM(rm.average)      + w.average      AS average_mw,
          SUM(rm.monoisotopic) + w.monoisotopic  AS monoisotopic_mw,
          COUNT(*)                               AS n_residues
        FROM read_parquet('${PROTEINS_URL}') p
        JOIN pick ON p.uniprot_id = pick.uniprot_id
        JOIN read_parquet('${MASSES_URL}') rm ON p.residue = rm.residue
        CROSS JOIN read_parquet('${WATER_URL}') w
        WHERE p.protein_position >= 1
        GROUP BY p.uniprot_id, w.average, w.monoisotopic
      )
      SELECT
        pick.*,
        seq.sequence,
        mw.average_mw,
        mw.monoisotopic_mw,
        mw.n_residues
      FROM pick
      JOIN seq ON pick.uniprot_id = seq.uniprot_id
      JOIN mw  ON pick.uniprot_id = mw.uniprot_id`;

    const stmt = await conn.prepare(sql);
    const result = await stmt.query(
      family,
      isoform ?? userOv.isoform ?? null,
      userOv.uniprot_id ?? null,
      userOv.protein_name ?? null,
      userOv.species ?? null,
      taxon
    );
    await stmt.close();

    const rows = result.toArray().map(r => r.toJSON());
    if (rows.length === 0) return null;

    const ctx = rows[0];
    ctx.modifications = userOv.modifications ?? [];
    return ctx;
  } finally {
    await conn.close();
  }
}

// ── Consensus-position alignment between source and port proteins ─────────
// Returns an object with:
//   srcToPort: Map<srcProteinPos, portProteinPos>  (mod position mapping)
//   portToSrc: Map<portProteinPos, srcProteinPos>  (reverse mapping)
//   driftByPort: Map<portProteinPos, { from, to }> (baseline residue diffs, keyed by port pos)
// Only aligned positions (both proteins have a consensus_position) are included.

const alignmentCache = new Map();

async function resolvePortAlignment(srcUniprotId, portUniprotId) {
  const key = `${srcUniprotId}→${portUniprotId}`;
  if (alignmentCache.has(key)) return alignmentCache.get(key);

  const conn = await db.connect();
  try {
    // Use family_position (aligned to the family reference) for cross-isoform
    // porting. This handles cases like H2A ↔ macroH2A where isoform-specific
    // consensus_position uses different coordinate systems, but family_position
    // is always relative to the same family reference (e.g. P0C0S8 for H2A).
    const sql = `
      WITH src AS (
        SELECT protein_position, family_position, residue
        FROM read_parquet('${PROTEINS_URL}')
        WHERE uniprot_id = $1 AND protein_position >= 1
          AND family_position IS NOT NULL
      ),
      port AS (
        SELECT protein_position, family_position, residue
        FROM read_parquet('${PROTEINS_URL}')
        WHERE uniprot_id = $2 AND protein_position >= 1
          AND family_position IS NOT NULL
      )
      SELECT
        src.protein_position  AS src_pos,
        src.residue           AS src_residue,
        port.protein_position AS port_pos,
        port.residue          AS port_residue
      FROM src
      JOIN port ON src.family_position = port.family_position
      ORDER BY src.protein_position`;
    const stmt = await conn.prepare(sql);
    const result = await stmt.query(srcUniprotId, portUniprotId);
    await stmt.close();

    const rows = result.toArray().map(r => r.toJSON());
    const srcToPort    = new Map();
    const portToSrc    = new Map();
    const driftByPort  = new Map();

    for (const r of rows) {
      srcToPort.set(r.src_pos, r.port_pos);
      portToSrc.set(r.port_pos, r.src_pos);
      if (r.src_residue !== r.port_residue) {
        driftByPort.set(r.port_pos, { from: r.src_residue, to: r.port_residue });
      }
    }

    const alignment = { srcToPort, portToSrc, driftByPort };
    alignmentCache.set(key, alignment);
    return alignment;
  } finally {
    await conn.close();
  }
}

function clearAlignmentCache() {
  alignmentCache.clear();
}

// ── Apply / Reset ────────────────────────────────────────────────────────

function applyPort() {
  const textarea = document.getElementById('port-yaml');
  const errorEl  = document.getElementById('port-error');
  if (!textarea) return;

  try {
    const parsed = jsyaml.load(textarea.value);
    if (parsed == null) {
      portOverrides = {};
      localStorage.removeItem(PORT_STORAGE_KEY);
      clearAlignmentCache();
      errorEl.textContent = 'Port context cleared.';
      errorEl.className = 'context-applied';
      reRender();
      return;
    }
    if (typeof parsed !== 'object') {
      errorEl.textContent = 'YAML must be a mapping of family tokens to override entries.';
      errorEl.className = 'context-error';
      return;
    }
    portOverrides = parsed;
    localStorage.setItem(PORT_STORAGE_KEY, textarea.value);
    clearAlignmentCache();
    errorEl.textContent = 'Port context applied.';
    errorEl.className = 'context-applied';
    reRender();
  } catch (e) {
    errorEl.textContent = `YAML error: ${e.message}`;
    errorEl.className = 'context-error';
  }
}

function resetPort() {
  portOverrides = {};
  localStorage.removeItem(PORT_STORAGE_KEY);
  clearAlignmentCache();
  const textarea = document.getElementById('port-yaml');
  const errorEl  = document.getElementById('port-error');
  if (textarea) textarea.value = PORT_TEMPLATE;
  if (errorEl)  errorEl.textContent = '';
  reRender();
}

// ── Initialise ───────────────────────────────────────────────────────────

function initPortEditor() {
  const textarea = document.getElementById('port-yaml');
  if (!textarea) return;

  const saved = localStorage.getItem(PORT_STORAGE_KEY);
  if (saved) {
    try {
      const parsed = jsyaml.load(saved);
      if (parsed && typeof parsed === 'object' && Object.keys(parsed).length > 0) {
        portOverrides = parsed;
        textarea.value = saved;
      } else {
        localStorage.removeItem(PORT_STORAGE_KEY);
        textarea.value = PORT_TEMPLATE;
      }
    } catch {
      localStorage.removeItem(PORT_STORAGE_KEY);
      textarea.value = PORT_TEMPLATE;
    }
  } else {
    textarea.value = PORT_TEMPLATE;
  }

  document.getElementById('port-apply')?.addEventListener('click', applyPort);
  document.getElementById('port-reset')?.addEventListener('click', resetPort);
}

initPortEditor();
