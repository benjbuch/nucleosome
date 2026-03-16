// ── DuckDB ────────────────────────────────────────────────────────────────

async function initDuckDB() {
  const bundles = duckdb.getJsDelivrBundles();
  const bundle  = await duckdb.selectBundle(bundles);
  const workerUrl = URL.createObjectURL(
    new Blob([`importScripts("${bundle.mainWorker}");`], { type: 'text/javascript' })
  );
  const worker = new Worker(workerUrl);
  db = new duckdb.AsyncDuckDB(new duckdb.ConsoleLogger('WARNING'), worker);
  await db.instantiate(bundle.mainModule, bundle.pthreadWorker);

  // Fetch parquet files and register as buffers.
  // HTTP range requests (registerFileURL) fail in Firefox on GitHub Pages,
  // so we fetch the full files up front and register them as local buffers.
  await Promise.all(
    [DATA_URL, PROTEINS_URL, MASSES_URL, PTM_URL, WATER_URL].map(async url => {
      const buf = await fetch(url).then(r => r.arrayBuffer());
      await db.registerFileBuffer(url, new Uint8Array(buf));
    })
  );
}

async function queryPosition(family, position, { variantAa = null, isoform = null } = {}) {
  const conn = await db.connect();
  try {
    // Build query with positional parameters to prevent injection.
    let sql = `
      SELECT
        uniprot_id,
        taxon_id,
        isoform,
        residue,
        variant_aa,
        ROUND(am_pathogenicity, 4) AS am_pathogenicity,
        am_class
      FROM read_parquet('${DATA_URL}')
      WHERE family = $1
        AND consensus_position = $2`;
    const params = [family, position];

    if (isoform) {
      params.push(isoform);
      sql += `\n            AND isoform = $${params.length}`;
    }
    if (variantAa) {
      params.push(variantAa);
      sql += `\n            AND variant_aa = $${params.length}`;
    }
    sql += `\n          ORDER BY variant_aa, am_pathogenicity DESC`;

    const stmt = await conn.prepare(sql);
    const result = await stmt.query(...params);
    await stmt.close();
    return result.toArray().map(r => r.toJSON());
  } finally {
    await conn.close();
  }
}

// ── Context resolution via proteins.parquet ──────────────────────────────

// Cache: contextKey → resolved context object
const contextCache = new Map();

function contextCacheKey(family, isoform, overrides) {
  return JSON.stringify([family, isoform, overrides]);
}

// Resolve a materialization context for a given family/isoform.
// `overrides` is an optional object with: { uniprot_id, protein_name, species, taxon_id }
// Returns: { family, isoform, uniprot_id, protein_name, species, taxon_id,
//            sequence, average_mw, monoisotopic_mw, n_residues, modifications }
// Returns null if no match.
async function resolveContext(family, isoform, overrides = {}) {
  const key = contextCacheKey(family, isoform, overrides);
  if (contextCache.has(key)) return contextCache.get(key);

  const conn = await db.connect();
  try {
    const lookupKey = isoform ?? family;
    const userOv    = contextOverrides[lookupKey] ?? contextOverrides[family] ?? {};
    const ov        = { ...userOv, ...overrides };

    // Default taxon to 9606 (human) if no constraint specified
    const taxon = ov.taxon_id ?? 9606;

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
      isoform ?? ov.isoform ?? null,
      ov.uniprot_id ?? null,
      ov.protein_name ?? null,
      ov.species ?? null,
      taxon
    );
    await stmt.close();

    const rows = result.toArray().map(r => r.toJSON());
    if (rows.length === 0) {
      contextCache.set(key, null);
      return null;
    }

    const ctx = rows[0];
    // Attach user-specified default modifications (from context editor).
    // Look up by: explicit overrides → input key → resolved isoform → family.
    // This ensures that an override on "H31" is picked up when querying "H3"
    // (which resolves to H3.1 / isoform H31).
    const mods = ov.modifications
      ?? contextOverrides[ctx.isoform]?.modifications
      ?? contextOverrides[ctx.family]?.modifications
      ?? [];
    ctx.modifications = mods;
    contextCache.set(key, ctx);
    return ctx;
  } finally {
    await conn.close();
  }
}

// Compute MW delta for a substitution: wt_residue → sub_residue
const subDeltaCache = new Map();

async function resolveSubstitutionDelta(wtResidue, subResidue) {
  const key = `${wtResidue}→${subResidue}`;
  if (subDeltaCache.has(key)) return subDeltaCache.get(key);

  const conn = await db.connect();
  try {
    const sql = `
      SELECT
        new.average       - old.average       AS average_delta,
        new.monoisotopic  - old.monoisotopic  AS monoisotopic_delta
      FROM read_parquet('${MASSES_URL}') old,
           read_parquet('${MASSES_URL}') new
      WHERE old.residue = $1
        AND new.residue = $2`;
    const stmt = await conn.prepare(sql);
    const result = await stmt.query(wtResidue, subResidue);
    await stmt.close();
    const rows = result.toArray().map(r => r.toJSON());
    const d = rows[0] ?? { average_delta: 0, monoisotopic_delta: 0 };
    subDeltaCache.set(key, d);
    return d;
  } finally {
    await conn.close();
  }
}

// Compute MW delta for a PTM
const ptmDeltaCache = new Map();

async function resolvePtmDelta(ptm) {
  if (ptmDeltaCache.has(ptm)) return ptmDeltaCache.get(ptm);

  const conn = await db.connect();
  try {
    const sql = `
      SELECT
        average       AS average_delta,
        monoisotopic  AS monoisotopic_delta
      FROM read_parquet('${PTM_URL}')
      WHERE ptm = $1`;
    const stmt = await conn.prepare(sql);
    const result = await stmt.query(ptm);
    await stmt.close();
    const rows = result.toArray().map(r => r.toJSON());
    const d = rows[0] ?? { average_delta: 0, monoisotopic_delta: 0 };
    ptmDeltaCache.set(ptm, d);
    return d;
  } finally {
    await conn.close();
  }
}

// Compute total MW delta for a list of modifications.
// Returns { average_delta, monoisotopic_delta }.
async function totalMwDelta(modifications) {
  let avg = 0, mono = 0;
  for (const m of modifications ?? []) {
    const wt  = m.residue;
    const sub = m.variant;
    const ptm = m.modification;
    if (sub && wt) {
      const d = await resolveSubstitutionDelta(wt, sub);
      avg  += d.average_delta;
      mono += d.monoisotopic_delta;
    }
    if (ptm) {
      const d = await resolvePtmDelta(ptm);
      avg  += d.average_delta;
      mono += d.monoisotopic_delta;
    }
  }
  return { average_delta: avg, monoisotopic_delta: mono };
}

// Flush context cache (called when user edits overrides).
function clearContextCache() {
  contextCache.clear();
}
