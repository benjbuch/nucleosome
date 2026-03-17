// ── Render utilities ──────────────────────────────────────────────────────

function scoreColor(score) {
  if (score < 0.34)  return 'var(--accent-green)';
  if (score < 0.564) return 'var(--accent-gold)';
  return 'var(--accent-red)';
}

// Render an AM rows table into `container`. Shows meta line when showMeta=true.
function renderAmTable(rows, info, container, { showMeta = true } = {}) {
  if (rows.length === 0) {
    container.innerHTML = '<p class="loading">No AlphaMissense data.</p>';
    return;
  }

  const nProteins = new Set(rows.map(r => r.uniprot_id)).size;
  const nTaxa     = new Set(rows.map(r => r.taxon_id)).size;

  // Determine active measurement port for this slot (if any)
  const mPortKey = info.slotFamily
    ? `${info.slotFamily}:${info.slotIsoform ?? ''}`
    : null;
  const activeMPort = mPortKey ? measurementPort.get(mPortKey) : null;

  const rows_html = rows.map(r => {
    const score     = +r.am_pathogenicity;
    const pct       = Math.round(score * 100);
    const color     = scoreColor(score);
    const cls       = (r.am_class ?? '').toLowerCase();
    const isoform = r.isoform && r.isoform !== 'NA' ? r.isoform : null;
    const isoformDisplay = isoform ?? '&mdash;';
    const varLabel  = `${r.residue ?? '?'}${info.position}${r.variant_aa}`;
    const isActive  = activeMPort?.uniprot_id === r.uniprot_id;

    return `<tr data-uniprot="${r.uniprot_id}" data-taxon="${r.taxon_id}" data-isoform="${isoform ?? ''}"${isActive ? ' class="am-row-active"' : ''}>
      <td><code>${varLabel}</code></td>
      <td>${r.uniprot_id}</td>
      <td>${r.taxon_id}</td>
      <td>${isoformDisplay}</td>
      <td>
        <div class="score-wrap">
          <div class="score-bar-bg">
            <div class="score-bar-fill" style="width:${pct}%;background:${color}"></div>
          </div>
          <span class="score-num">${score.toFixed(3)}</span>
        </div>
      </td>
      <td><span class="badge ${cls}">${r.am_class}</span></td>
    </tr>`;
  }).join('');

  const meta = showMeta
    ? `<p class="result-meta">${rows.length} variant scores &middot; ${nProteins} protein${nProteins !== 1 ? 's' : ''} &middot; ${nTaxa} organism${nTaxa !== 1 ? 's' : ''}</p>`
    : '';

  container.innerHTML = `${meta}
    <div class="scrollable-table">
      <table style="width:100%">
        <thead>
          <tr>
            <th class="am-th">Variant</th>
            <th class="am-th">UniProt</th>
            <th class="am-th">Taxon</th>
            <th class="am-th">Isoform</th>
            <th class="am-th">Pathogenicity</th>
            <th class="am-th">Class</th>
          </tr>
        </thead>
        <tbody style="font-family:var(--font-mono);font-size:0.92rem">${rows_html}</tbody>
      </table>
    </div>`;

  // Click-to-port: clicking an AM row ports the slot's sequence display to
  // that measurement's protein.  Click again to dismiss.
  const tbody = container.querySelector('tbody');
  if (tbody && info.slotFamily) {
    tbody.addEventListener('click', (e) => {
      const tr = e.target.closest('tr[data-uniprot]');
      if (!tr) return;
      const uid   = tr.dataset.uniprot;
      const taxon = +(tr.dataset.taxon);
      const iso   = tr.dataset.isoform || null;
      const key   = `${info.slotFamily}:${info.slotIsoform ?? ''}`;
      if (measurementPort.get(key)?.uniprot_id === uid) {
        measurementPort.delete(key);  // toggle off
      } else {
        measurementPort.set(key, { uniprot_id: uid, taxon_id: taxon, isoform: iso });
      }
      lastKey = null;  // force re-render
      document.getElementById('notation-input').dispatchEvent(new Event('input'));
    });
  }
}

// Format full MW: average in kDa (1 dp) and monoisotopic in u with
// comma-thousands separator (1 dp). e.g. "15.3 kDa; 15,263.4 u"
function formatMW(average_mw, monoisotopic_mw) {
  const kDa  = (average_mw / 1000).toFixed(1) + '\u202fkDa';
  const mono = monoisotopic_mw.toLocaleString('en-US', {
    minimumFractionDigits: 1,
    maximumFractionDigits: 1
  }) + '\u202fu';
  return `${kDa}; ${mono}`;
}

// Format a MW delta: signed monoisotopic in u (3 dp) and average in Da (3 dp).
// Returns null when both deltas are zero (no modifications — nothing to show).
// e.g. "+42.0 Da; +42.011 u" for acetylation, "+2.946 u; +3.019 Da" for K→M
function formatMWDelta(avg_delta, mono_delta) {
  if (mono_delta === 0 && avg_delta === 0) return null;
  const sign  = n => n >= 0 ? '+' : '\u2212';
  const avg   = `${sign(avg_delta)}${Math.abs(avg_delta).toFixed(2)}\u202fDa`;
  const mono  = `${sign(mono_delta)}${Math.abs(mono_delta).toFixed(3)}\u202fu`;
  return `${avg}; ${mono}`;
}

// Format a protein sequence with modification positions highlighted.
// Substitutions replace the residue (K27M → "M"), PTMs add a superscript
// (K27ac → "K^ac"), combined shows both (K27Mac → "M^ac").
//
// When `drift` is provided (a Map of position → { from, to }), positions
// where the source and port-to contexts differ at baseline are annotated:
//   - drift only (no notation mod): outlined border
//   - drift + notation mod at same position: fill + outline
//   - notation mod neutralized (already native in port-to): gray fill + outline
function formatSequence(sequence, modifications, drift) {
  if (!sequence) return '';
  const modByPos = {};
  for (const m of modifications ?? []) modByPos[m.position] = m;
  return Array.from(sequence).map((aa, i) => {
    const pos = i + 1;
    const m = modByPos[pos];
    const d = drift?.get(pos);

    // Pure drift — no notation mod, just baseline difference
    if (d && !m) {
      const title = `${d.from}${pos} → ${d.to}${pos} (baseline drift)`;
      return `<mark class="seq-drift" title="${title}">${aa}</mark>`;
    }

    if (!m) return `<span title="${aa}${pos}">${aa}</span>`;

    const sub = m.variant;
    const ptm = m.modification;
    const display = sub ?? aa;
    const ptmHtml = ptm ? `<sup class="seq-ptm">${ptm}</sup>` : '';

    // Neutralized: notation mod is already the native state in port-to
    if (m._neutralized) {
      const title = `${m._sourceWt ?? '?'}${pos}→${display} (neutralized — native in target)`;
      return `<mark class="seq-neutralized" title="${title}">${display}${ptmHtml}</mark>`;
    }

    // Notation mod at a drift position: fill + outline
    const driftClass = d ? ' seq-drift-mod' : '';
    const title = `${aa}${pos}` + (sub ? `→${sub}` : '') + (ptm ? ` ${PTM_NAMES[ptm] ?? ptm}` : '')
      + (d ? ` (was ${d.from} in source context)` : '');
    return `<mark class="seq-highlight${sub ? ' seq-sub' : ''}${driftClass}" title="${title}">${display}${ptmHtml}</mark>`;
  }).join('');
}
