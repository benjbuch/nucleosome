// ── Parser ────────────────────────────────────────────────────────────────

function tryParse(input) {
  try { return nucleosomeParser.parse(input.trim()); }
  catch { return null; }
}

// ── Schema vocabulary ─────────────────────────────────────────────────────
// PTM_NAMES is loaded from vocabulary.js (generated from vocabulary.yaml).
// family/isoform are normalized in the parse tree — no mapping needed.

// Canonical display string from parse tree fields.
//   { family:"H2A", isoform:"H2AZ" }  → "H2AZ"
//   { family:"H3",  isoform:"H33" }   → "H33"
//   { family:"H3" }                    → "H3"
function famStr(family, isoform) {
  if (isoform) return isoform;
  return family;
}

// Human-readable label for isoform context (null for plain family).
function isoformLabel(isoform) {
  if (isoform) return `${isoform} isoform`;
  return null;
}

// Describe a single modification entry → { canonical, type, pos, wt, ptm, sub }
function describeMod(m) {
  const pos    = m.position;
  const wt     = m.residue;
  const ptm    = m.modification;
  const sub    = m.variant;
  const prefix = wt ? `${wt}${pos}` : `${pos}`;

  if (ptm && sub) return { canonical: `${prefix}${sub}${ptm}`, type: 'mut+ptm',      pos, wt, ptm, sub };
  if (ptm)        return { canonical: `${prefix}${ptm}`,        type: 'modification', pos, wt, ptm };
  if (sub)        return { canonical: `${prefix}${sub}`,        type: 'substitution', pos, wt, sub };
  return           { canonical: prefix,                          type: 'position',     pos, wt };
}

// Canonical string for one nucleosome slot (used in nucleosome display).
// Certainty brackets are preserved: [exact], {native}, bare = unspecified.
function slotCanonical(slot) {
  const fs      = famStr(slot.family, slot.isoform);
  const mods    = (slot.modifications ?? []).map(m => describeMod(m).canonical);
  const content = mods.length ? `${fs}:${mods.join(':')}` : fs;
  if (slot.certainty === 'exact')  return `[${content}]`;
  if (slot.certainty === 'native') return `{${content}}`;
  return content;
}

// Returns a structured descriptor for every level of the schema.
// { type, canonical, detail, queryable, family, isoform, position, wt, variant, ptm }
function describe(parsed) {
  if (!parsed) return null;

  // ── Family ──
  if (parsed.level === 'family') {
    const fs    = famStr(parsed.family, parsed.isoform);
    const label = isoformLabel(parsed.isoform);
    return {
      type: 'family', canonical: fs,
      detail: label ?? 'all paralogs',
      family: parsed.family, isoform: parsed.isoform,
      queryable: false
    };
  }

  // ── Histone (standalone query) ──
  if (parsed.level === 'histone') {
    const mods = parsed.modifications ?? [];
    const fs    = famStr(parsed.family, parsed.isoform);
    const label = isoformLabel(parsed.isoform);

    if (mods.length === 0) {
      // Unmodified histone in certainty brackets, e.g. [H3]
      const cert = parsed.certainty === 'exact' ? 'exact' : parsed.certainty === 'native' ? 'native' : null;
      const canonical = cert === 'exact' ? `[${fs}]` : cert === 'native' ? `{${fs}}` : fs;
      const detail = label ?? 'all paralogs';
      return {
        type: 'family', canonical, detail,
        family: parsed.family, isoform: parsed.isoform,
        queryable: false
      };
    }

    const dm    = describeMod(mods[0]);
    const canonical = `${fs}:${mods.map(m => describeMod(m).canonical).join(':')}`;

    let detail;
    if      (dm.type === 'position')     detail = `position ${dm.pos} · all substitutions`;
    else if (dm.type === 'modification') detail = `position ${dm.pos} · ${PTM_NAMES[dm.ptm] ?? dm.ptm}`;
    else if (dm.type === 'substitution') detail = `position ${dm.pos} · ${dm.wt ?? '?'}\u2192${dm.sub}`;
    else if (dm.type === 'mut+ptm')      detail = `position ${dm.pos} · ${dm.sub} + ${PTM_NAMES[dm.ptm] ?? dm.ptm}`;
    if (label)           detail = `${label} · ${detail}`;
    if (mods.length > 1) detail += ` · ${mods.length} modifications`;

    return {
      type: dm.type, canonical, detail,
      family: parsed.family, isoform: parsed.isoform,
      position: dm.pos, wt: dm.wt, variant: dm.sub, ptm: dm.ptm,
      queryable: true
    };
  }

  // ── Nucleosome ──
  if (parsed.level === 'nucleosome') {
    const slots      = parsed.slots ?? [];
    const canonical  = `(${slots.map(slotCanonical).join(' \u0040 ')})`;
    const families   = [...new Set(slots.map(s => s.family))].join(' + ');
    const nModded    = slots.filter(s => (s.modifications ?? []).length > 0).length;
    const nSlots     = slots.length;
    let detail       = `nucleosome · ${nSlots} slot${nSlots !== 1 ? 's' : ''} · ${families}`;
    if (nModded) detail += ` · ${nModded} modified`;
    return { type: 'nucleosome', canonical, detail, queryable: false };
  }

  // ── Array ──
  if (parsed.level === 'array') {
    const n = parsed.nucleosomes?.length ?? 0;
    const canonical = (parsed.nucleosomes ?? [])
      .map(nuc => `(${(nuc.slots ?? []).map(slotCanonical).join(' \u0040 ')})`)
      .join(' ');
    return {
      type: 'array', canonical,
      detail: `${n}-nucleosome array`, queryable: false
    };
  }

  return null;
}
