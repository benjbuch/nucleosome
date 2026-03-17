// ── Slot / nucleosome / array rendering ───────────────────────────────────

// Histone family colors — same palette as the SVG diagram trajectories.
const SLOT_COLOR = { H2A: '#b3b300', H2B: '#e31a1c', H3: '#1f78b4', H4: '#33a02c' };

// Persist open/closed state for sequence and mod-section <details>.
// Keys: "seq:<family>:<isoform>" or "mod:<family>:<pos>".
const detailsState = new Map();

// Measurement port state: clicking an AM row temporarily ports the slot's
// sequence display to that measurement's protein.
// Keys: "family:isoform" → { uniprot_id, taxon_id, isoform }
const measurementPort = new Map();

function detailsKey(prefix, family, isoform, pos) {
  return pos != null
    ? `${prefix}:${family}:${pos}`
    : `${prefix}:${family}:${isoform ?? ''}`;
}

// Render one nucleosome's slots as a <ul> into `container`.
// Each modification within a slot gets its own collapsible sub-section.
async function renderNucleosome(slots, container, key) {
  container.innerHTML = '';

  // ── Nucleosome-level layers (e.g. screen results) ──────────────────────
  if (db && layersAt('nucleosome').length > 0) {
    // TODO: nucleosome-level query needs the full slot set, not a single
    // family+position.  For now this is scaffolding — no layers use it yet.
  }

  const ul = document.createElement('ul');
  ul.className = 'nucleosome-list';
  container.appendChild(ul);

  for (const slot of slots) {
    if (key !== lastKey) return;

    let   mods = slot.modifications ?? [];
    const fs   = famStr(slot.family, slot.isoform);
    const head = mods.length ? slotCanonical(slot) : fs;
    const queryFamily = slot.family;
    const isoform     = slot.isoform;

    // Resolve context via SQL
    const ctx = await resolveContext(queryFamily, isoform);
    if (key !== lastKey) return;

    // Merge context modifications with notation modifications.
    //
    // A notation mod is "modifying" only if it specifies a variant or PTM.
    // Bare position queries (e.g. C96) are just queries — they don't
    // displace context mods at the same position.
    //
    // For modifying notation mods that overlap a context mod:
    //  - Redundant (same variant+PTM) → use the context mod (preserves
    //    the correct original-sequence wt residue for MW computation).
    //  - Different effect → override context, but pin wt residue to the
    //    original sequence residue so MW is always computed against raw_mw.
    //    Warn if the user's stated wt doesn't match the context state.
    if (ctx?.modifications?.length) {
      const seq = ctx.sequence ?? '';

      // Index context mods by position, filling in residue from sequence
      const ctxModByPos = new Map();
      for (const m of ctx.modifications) {
        ctxModByPos.set(m.position, {
          ...m,
          residue: m.residue ?? (seq[m.position - 1] || null),
        });
      }

      // Only modifying notation mods displace context mods
      const isModifying = m => m.variant != null || m.modification != null;
      const modifyingPositions = new Set(
        mods.filter(isModifying).map(m => m.position)
      );

      // Resolve each notation mod against context at the same position
      const resolvedNotation = mods.map(m => {
        const tagged = { ...m, source: 'notation' };
        const cm = ctxModByPos.get(m.position);
        if (!cm) return tagged;

        // Bare position query — coexists with context mod, doesn't override.
        // Still warn if stated wt doesn't match context state.
        if (!isModifying(m)) {
          const ctxResidue = cm.variant ?? cm.residue;
          if (m.residue && ctxResidue && m.residue !== ctxResidue) {
            tagged.contextConflict = ctxResidue;
          }
          return tagged;
        }

        // Redundant (same variant+PTM as context) → use context mod
        // with its correct original-sequence wt residue for MW.
        const sameVariant = (m.variant ?? null) === (cm.variant ?? null);
        const samePtm     = (m.modification ?? null) === (cm.modification ?? null);
        if (sameVariant && samePtm) {
          return { ...cm, source: 'context' };
        }

        // Different effect — override, but pin wt to the original sequence
        // residue so MW delta is always computed against raw_mw.
        const origResidue = seq[m.position - 1] || null;
        const ctxResidue  = cm.variant ?? cm.residue;
        if (m.residue && ctxResidue && m.residue !== ctxResidue) {
          tagged.contextConflict = ctxResidue;
        }
        return { ...tagged, residue: origResidue };
      });

      // Keep context mods at positions not displaced by modifying notation
      const contextMods = ctx.modifications
        .filter(m => !modifyingPositions.has(m.position))
        .map(m => ({
          ...m,
          residue: m.residue ?? (seq[m.position - 1] || null),
          source: 'context'
        }));

      mods = [...contextMods, ...resolvedNotation].sort((a, b) => a.position - b.position);
    }

    // ── Validate wt residues against the resolved sequence ──────────────
    // Runs for ALL notation mods, whether or not context modifications
    // exist.  The "effective" residue at a position is the context-modified
    // residue (if a context mod sits at that position) or the raw sequence
    // residue.  A mismatch means the user's stated wt is wrong.
    if (ctx?.sequence) {
      const seq = ctx.sequence;
      for (const m of mods) {
        if (m.source === 'context' || !m.residue || m.contextConflict) continue;
        const ctxMod = mods.find(cm =>
          cm.position === m.position && cm.source === 'context');
        const effective = ctxMod?.variant ?? seq[m.position - 1] ?? null;
        if (effective && m.residue !== effective) {
          m.contextConflict = effective;
        }
      }
    }

    // ── Port-to context resolution ────────────────────────────────────────
    // When porting is active, resolve the port-to context for this slot.
    // The display switches to the port-to sequence. Alignment is done via
    // consensus_position (not array index) so cross-isoform porting handles
    // Met cleavage differences and insertions/deletions correctly.
    //
    // We compute:
    //   - drift: positions where source and target baseline residues differ
    //     (keyed by port protein_position for display on the port sequence)
    //   - mapped mods: notation mod positions translated src→port via consensus
    //   - neutralized mods: notation mods already native in port-to
    let displayCtx = ctx;
    let drift = null; // Map<portPos, { from, to }>
    const mPortKey = `${queryFamily}:${isoform ?? ''}`;
    const mPort = measurementPort.get(mPortKey);

    if (mPort || isPortActive()) {
      let portCtx;
      if (mPort) {
        // Measurement port: resolve the specific measurement protein
        portCtx = await resolveContext(queryFamily, mPort.isoform, {
          uniprot_id: mPort.uniprot_id,
          taxon_id: mPort.taxon_id,
        });
      } else {
        portCtx = await resolvePortContext(queryFamily, isoform);
      }
      if (key !== lastKey) return;

      if (!portCtx && !mPort) {
        // Port context requested but failed to resolve for this slot.
        // Check if the user explicitly specified an override for this family.
        const lookupKey = isoform ?? queryFamily;
        const hasOverride = portOverrides[lookupKey] ?? portOverrides[queryFamily];
        if (hasOverride) {
          const errorEl = document.getElementById('port-error');
          if (errorEl) {
            const msg = `Could not resolve ${lookupKey} in port context.`;
            errorEl.textContent = msg;
            errorEl.className = 'context-error';
          }
        }
      }

      if (portCtx?.sequence && ctx?.uniprot_id) {
        displayCtx = portCtx;

        // Align via consensus_position JOIN
        const alignment = await resolvePortAlignment(ctx.uniprot_id, portCtx.uniprot_id);
        if (key !== lastKey) return;

        drift = alignment.driftByPort;
        const portSeq = portCtx.sequence;

        // Map each mod's position from source protein_position to port protein_position.
        // Mods that can't be mapped (position falls in an unaligned region) are flagged.
        for (const m of mods) {
          if (m.source === 'context') continue;
          const portPos = alignment.srcToPort.get(m.position);
          if (portPos == null) {
            // Position not alignable — flag it
            m._unmapped = true;
            continue;
          }

          // Store original position and remap
          m._srcPosition = m.position;
          m.position = portPos;

          // Check neutralization: variant already native in port-to?
          const portAa = portSeq[portPos - 1];
          if (portAa && m.variant && m.variant === portAa) {
            m._neutralized = true;
            m._sourceWt = ctx.sequence?.[m._srcPosition - 1] ?? m.residue;
          }
        }
      }
    }

    const li = document.createElement('li');
    li.className = 'nuc-slot';

    // Histone-family colored left border
    const color = SLOT_COLOR[queryFamily];
    if (color) li.style.setProperty('--slot-color', color);

    // Slot header: notation + context metadata on one line
    // When porting, show the port-to context metadata instead.
    const metaCtx = displayCtx ?? ctx;
    let headerHtml = `<div class="slot-header"><code>${head}</code>`;
    if (metaCtx) {
      if (mPort && displayCtx !== ctx) {
        headerHtml += `<code class="ctx-measurement-badge">source</code>`;
      } else if (isPortActive() && displayCtx !== ctx) {
        headerHtml += `<code class="ctx-port-badge">port</code>`;
      }
      headerHtml += `<span class="slot-context">`
        + `<span class="ctx-protein">${metaCtx.protein_name}</span>`
        + `<span class="ctx-id">${metaCtx.uniprot_id}</span>`
        + `<span class="ctx-taxon">${metaCtx.species} (${metaCtx.taxon_id})</span>`
        + `</span>`;
    }
    headerHtml += `</div>`;

    // Sequence display (open by default, state persisted across re-renders)
    if (displayCtx?.sequence) {
      const seqHtml = formatSequence(displayCtx.sequence, mods, drift);
      let mwStr = '';
      if (displayCtx.average_mw != null && displayCtx.monoisotopic_mw != null) {
        // Pin all residues to the raw sequence for MW computation.
        // MW is always relative to displayCtx.average_mw (raw sequence MW), so
        // deltas must be computed from the original residue, regardless of
        // what the user stated or what context modifications changed.
        const seq = displayCtx.sequence ?? '';
        const mwMods = mods.map(m => {
          const raw = seq[m.position - 1];
          return raw ? { ...m, residue: raw } : m;
        });
        const allDelta = await totalMwDelta(mwMods);
        if (key !== lastKey) return;
        const baseMW = formatMW(
          displayCtx.average_mw + allDelta.average_delta,
          displayCtx.monoisotopic_mw + allDelta.monoisotopic_delta
        );
        if (slot.certainty === 'native' || slot.certainty == null) {
          // Δ display: delta relative to context baseline, not raw sequence.
          // contextDelta accounts for mods the user chose as their baseline;
          // the visible Δ is only the additional effect of the notation.
          const ctxModsFull = (displayCtx.modifications ?? []).map(m => ({
            ...m, residue: m.residue ?? (seq[m.position - 1] || null)
          }));
          const ctxDelta = await totalMwDelta(ctxModsFull);
          if (key !== lastKey) return;
          const deltaStr = formatMWDelta(
            allDelta.average_delta    - ctxDelta.average_delta,
            allDelta.monoisotopic_delta - ctxDelta.monoisotopic_delta
          );
          mwStr = deltaStr ? ` (${baseMW}; ${deltaStr})` : ` (${baseMW})`;
        } else {
          mwStr = ` (${baseMW})`;
        }
      }

      const seqKey = detailsKey('seq', queryFamily, isoform);
      // Default to open; honour persisted state if user toggled it
      const seqOpen = detailsState.has(seqKey) ? detailsState.get(seqKey) : true;

      headerHtml += `<details class="slot-sequence" data-dk="${seqKey}"${seqOpen ? ' open' : ''}>`
        + `<summary>${displayCtx.sequence.length} aa${mwStr}</summary>`
        + `<code class="seq-display">${seqHtml}</code>`
        + `</details>`;
    }

    li.innerHTML = headerHtml;

    // Persist sequence toggle state
    const seqDetails = li.querySelector('.slot-sequence');
    if (seqDetails) {
      seqDetails.addEventListener('toggle', () => {
        detailsState.set(seqDetails.dataset.dk, seqDetails.open);
      });
    }

    ul.appendChild(li);

    // ── Slot-level layers (e.g. contacts, conservation) ────────────────
    if (db && layersAt('slot').length > 0) {
      // TODO: slot-level queries will pass all positions for this histone.
      // No layers use this level yet.
    }

    if (mods.length === 0) {
      const p = document.createElement('p');
      p.className = 'loading';
      p.textContent = 'No modification specified';
      li.appendChild(p);
      continue;
    }

    // One collapsible section per modification
    for (const mod of mods) {
      if (key !== lastKey) return;

      const dm   = describeMod(mod);
      const hint = dm.type === 'position'     ? `Position ${dm.pos} \u00b7 All substitutions`
                 : dm.type === 'modification' ? `Position ${dm.pos} \u00b7 ${PTM_NAMES[dm.ptm] ?? dm.ptm}`
                 : dm.type === 'substitution' ? `Position ${dm.pos} \u00b7 ${dm.wt ?? '?'}\u2192${dm.sub}`
                 :                              `Position ${dm.pos}`;

      const section = document.createElement('details');
      section.className = 'mod-section';

      // Persist mod-section open/closed state
      const modKey = detailsKey('mod', queryFamily, isoform, dm.pos);
      section.dataset.dk = modKey;
      // Default: collapsed.  Honour persisted state if user toggled it.
      if (detailsState.has(modKey) && detailsState.get(modKey)) {
        section.open = true;
      }
      section.addEventListener('toggle', () => {
        detailsState.set(section.dataset.dk, section.open);
      });

      let summaryHtml = hint;
      if (mod.contextConflict) {
        summaryHtml += `<span class="context-warning">Residue is ${mod.contextConflict} at this position</span>`;
      }
      const summary = document.createElement('summary');
      summary.innerHTML = summaryHtml;
      section.appendChild(summary);

      li.appendChild(section);

      // Skip AM table for context-source modifications — they represent
      // the baseline, not an active query.
      if (db && mod.source !== 'context') {
        const modBody = document.createElement('div');
        modBody.className = 'mod-body';
        section.appendChild(modBody);

        const modValue = (isStrict() && dm.type === 'substitution') ? dm.sub : null;
        const resultsByLayer = await queryLayers('modification', queryFamily, dm.pos, { modValue, isoform });
        if (key !== lastKey) return;
        renderLayerResults(resultsByLayer, { position: dm.pos, slotFamily: queryFamily, slotIsoform: isoform }, modBody, { showMeta: false });
      }
    }
  }
}

// Render an array of nucleosomes as a nested expandable list into `container`.
// Pass { dimHeader: true } when wrapping a standalone histone/family so the
// "Nucleosome N" label is hidden (the grid layout is still identical,
// preventing horizontal reflow on transitions).
async function renderArray(nucleosomes, container, key, { dimHeader = false } = {}) {
  // ── Array-level layers (e.g. dinucleosome data) ──────────────────────
  if (db && layersAt('array').length > 0) {
    // TODO: array-level query needs the full nucleosome set.
    // No layers use this level yet.
  }

  // Build the entire tree in a detached <ul>, then swap it into the
  // container in one shot.  This keeps old content visible while async
  // work (context resolution, MW, layer queries) runs — no flash of
  // empty content on measurement-port re-renders or fresh queries.
  const ul = document.createElement('ul');
  ul.className = 'array-list';

  for (let i = 0; i < nucleosomes.length; i++) {
    if (key !== lastKey) return;

    const nuc   = nucleosomes[i];
    const slots = nuc.slots ?? [];
    const label = slots.map(slotCanonical).join(' \u0040 ');

    // .nuc-item is a two-column grid: card (col 1) | diagram (col 2).
    // Keeping them as siblings — not parent/child — means the diagram
    // never participates in the card's internal layout, eliminating the
    // reflow that caused wobbling when content streamed in or the
    // sequence <details> was expanded.
    const li = document.createElement('li');
    li.className = 'nuc-item';

    const details = document.createElement('details');
    details.open = true;
    details.className = dimHeader ? 'nuc-details is-context' : 'nuc-details';

    const summary = document.createElement('summary');
    summary.innerHTML = `<code>Nucleosome ${i + 1}</code><span class="status-detail">${label}</span>`;
    details.appendChild(summary);

    const body = document.createElement('div');
    details.appendChild(body);
    li.appendChild(details);

    // Diagram panel is always present (col 2) so column widths never
    // change between query types.  SVG is only rendered for nucleosomes.
    const diagramPanel = document.createElement('object');
    diagramPanel.className = 'nuc-diagram-panel';
    li.appendChild(diagramPanel);
    if (!dimHeader) renderNucleosomeDiagram(slots, diagramPanel, key); // fire-and-forget

    ul.appendChild(li);

    await renderNucleosome(slots, body, key);
    if (key !== lastKey) return;
  }

  // Atomic swap: old content → new content in a single reflow.
  container.replaceChildren(ul);
}
