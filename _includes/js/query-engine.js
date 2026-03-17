// ── Universal query engine ──────────────────────────────────────────────────
// Layer registry: each layer defines a parquet URL, a SQL view, a render
// level, and a renderer.  Levels map to the parse-tree hierarchy:
//
//   array  →  rendered once above the nucleosome list
//   nucleosome  →  rendered once per nucleosome, above slots
//   slot  →  rendered once per histone slot, above mods
//   modification  →  rendered inside each mod-section (current AM behavior)
//
// Adding a new data layer requires only a new entry here + a render function.

const LAYER_LEVELS = ['array', 'nucleosome', 'slot', 'modification'];

const LAYERS = {
  alpha_missense: {
    level: 'modification',
    label: 'AlphaMissense',
    parquetUrl: DATA_URL,
    viewSql: (url) => `
      CREATE OR REPLACE VIEW layer_alpha_missense AS
      SELECT
        ROW_NUMBER() OVER () AS measurement_id,
        family, isoform, taxon_id, uniprot_id,
        consensus_position AS position,
        'variant' AS mod_type,
        variant_aa AS mod_value,
        'notation' AS source,
        NULL::SMALLINT AS nucleosome_position,
        NULL::SMALLINT AS histone_position,
        'alpha_missense' AS layer,
        family || ':' || residue || CAST(consensus_position AS TEXT) || variant_aa AS entry_key,
        ROUND(am_pathogenicity, 4) AS estimate,
        'pathogenicity' AS unit,
        ROUND(am_pathogenicity, 4) AS am_pathogenicity,
        am_class, residue, variant_aa
      FROM read_parquet('${url}')`,
    render: renderAmTable,
  }
};

// Return layer entries filtered to a given level.
function layersAt(level) {
  return Object.entries(LAYERS).filter(([, l]) => l.level === level);
}

async function initLayers() {
  const conn = await db.connect();
  try {
    for (const [name, layer] of Object.entries(LAYERS)) {
      await conn.query(layer.viewSql(layer.parquetUrl));
    }
  } finally {
    await conn.close();
  }
}

// Query layers at a given level for a family + position.
// Options: { modType, modValue, isoform }
// Returns: Map<layerName, rows[]>
async function queryLayers(level, family, position, { modType = null, modValue = null, isoform = null } = {}) {
  const layers = layersAt(level);
  if (layers.length === 0) return new Map();

  const conn = await db.connect();
  try {
    const results = new Map();
    for (const [name, layer] of layers) {
      let sql = `
        SELECT * FROM layer_${name}
        WHERE family = $1
          AND position = $2`;
      const params = [family, position];

      if (isoform) {
        params.push(isoform);
        sql += `\n          AND isoform = $${params.length}`;
      }
      if (modValue) {
        params.push(modValue);
        sql += `\n          AND mod_value = $${params.length}`;
      }
      sql += `\n        ORDER BY entry_key, estimate DESC`;

      const stmt = await conn.prepare(sql);
      const result = await stmt.query(...params);
      await stmt.close();
      results.set(name, result.toArray().map(r => r.toJSON()));
    }
    return results;
  } finally {
    await conn.close();
  }
}

// Render results from all layers into a container.
function renderLayerResults(resultsByLayer, info, container, opts = {}) {
  const entries = [...resultsByLayer.entries()];
  if (entries.length === 0) return;

  const multiLayer = entries.length > 1;

  for (const [name, rows] of entries) {
    const layer = LAYERS[name];
    if (!layer) continue;

    const section = document.createElement('div');
    section.className = 'layer-section';

    if (multiLayer) {
      const header = document.createElement('div');
      header.className = 'layer-header';
      header.textContent = layer.label;
      section.appendChild(header);
    }

    const body = document.createElement('div');
    section.appendChild(body);
    layer.render(rows, info, body, opts);

    container.appendChild(section);
  }
}
