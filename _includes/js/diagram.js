    // ── Nucleosome diagram ────────────────────────────────────────────────────

    const HISTONE_LENGTH   = { H2A: 129, H2B: 125, H3: 135, H4: 102 };
    const TRAJECTORY_COLOR = { H2A: '#b3b300', H2B: '#e31a1c', H3: '#1f78b4', H4: '#33a02c' };

    const DOT_R           = 0.5;  // filled dot radius (user units)
    const STEM_LEN        = 3.0;  // radial stub from dot to polyline bend (user units)
    const LEADER_LEN      = 6.0;  // additional reach from stem tip to label (user units)
    const LEADER_SW       = 0.3;  // polyline stroke-width
    const LABEL_FS        = 3.0;  // font-size (unitless SVG user units, not CSS px)
    const MIN_ANG         = 0.20; // minimum angular gap between labels (radians)
    const MAX_ANG_DISP    = 0.45; // max angular displacement from natural angle (radians ~26°)
    const SPLIT_THRESHOLD = 14;   // split mirrored group when dots are farther apart (user units)
    const LABEL_MARGIN    = 0.5;  // keep labels this far inside the viewBox edge (user units)

    async function fetchNucleosomeSvg() {
      if (cachedSvg) return cachedSvg;
      const res  = await fetch(SVG_URL);
      const text = await res.text();
      cachedSvg  = new DOMParser().parseFromString(text, 'image/svg+xml').documentElement;
      return cachedSvg;
    }

    // Parse the translate(tx, ty) from the first <g> child's transform attribute.
    // Returns { tx, ty } or { tx: 0, ty: 0 } if none found.
    function parseGroupTranslate(svgEl) {
      const g = svgEl.querySelector('g[transform]');
      if (!g) return { tx: 0, ty: 0 };
      const m = g.getAttribute('transform').match(
        /translate\(\s*(-?[\d.]+)[\s,]+(-?[\d.]+)\s*\)/
      );
      return m ? { tx: parseFloat(m[1]), ty: parseFloat(m[2]) } : { tx: 0, ty: 0 };
    }

    // Derive nucleosome centre and viewBox bounds in group coordinate space.
    // group transform: svgCoord = groupCoord + (tx, ty)
    //   → groupCoord = svgCoord − (tx, ty)
    function svgGeometry(svgEl) {
      const vb = svgEl.getAttribute('viewBox').trim().split(/\s+/).map(Number);
      const { tx, ty } = parseGroupTranslate(svgEl);
      const center = { x: vb[0] + vb[2] / 2 - tx, y: vb[1] + vb[3] / 2 - ty };
      const bMinX  = vb[0] - tx + LABEL_MARGIN;
      const bMaxX  = vb[0] + vb[2] - tx - LABEL_MARGIN;
      const bMinY  = vb[1] - ty + LABEL_MARGIN;
      const bMaxY  = vb[1] + vb[3] - ty - LABEL_MARGIN;
      return { center, bMinX, bMaxX, bMinY, bMaxY };
    }

    // Find a trajectory path element for the given family and copy index.
    // Handles both naming conventions: "H2A:1:trajectory".
    function findTrajectory(svgEl, family, copy) {
      return svgEl.querySelector(`[id="${family}:${copy}:trajectory"]`);
    }

    // Return the {x, y} of node n (1-indexed) in a polyline path (m x0,y0 dx…).
    // Parses the d attribute directly — does not require the element to be in the DOM.
    function getPathNode(pathEl, n) {
      const nums = [...pathEl.getAttribute('d').matchAll(/-?[\d.]+(?:e[+-]?\d+)?/gi)].map(Number);
      let x = nums[0], y = nums[1];
      for (let i = 1; i < n; i++) { x += nums[2 * i]; y += nums[2 * i + 1]; }
      return { x, y };
    }

    // Spread angles so no two are closer than minGap.
    // Forward greedy sweep, re-centred to minimise max displacement.
    function spreadAngles(angles, minGap) {
      const n = angles.length;
      if (n <= 1) return [...angles];
      const ord = [...Array(n).keys()].sort((a, b) => angles[a] - angles[b]);
      const nat = ord.map(i => angles[i]);
      const adj = [...nat];
      for (let i = 1; i < n; i++) {
        if (adj[i] - adj[i - 1] < minGap) adj[i] = adj[i - 1] + minGap;
      }
      const excess = adj[n - 1] - nat[n - 1];
      for (let i = 0; i < n; i++) adj[i] -= excess / 2;
      const result = new Array(n);
      ord.forEach((origIdx, k) => { result[origIdx] = adj[k]; });
      return result;
    }

    // Render dot + angled polyline leader + text for every label.
    // Mirrored pairs (same mirrorGroup) share one label text; each dot gets its own leader.
    // labelList: [{family, x, y, position, mirrorGroup}] in group coordinate space.
    function renderDiagramLabels(svgEl, labelList) {
      if (!labelList.length) return;

      const ns    = 'http://www.w3.org/2000/svg';
      const group = svgEl.querySelector('g[transform]') || svgEl;
      const { center, bMinX, bMaxX, bMinY, bMaxY } = svgGeometry(svgEl);

      // ── Group mirrored pairs ──────────────────────────────────────────────
      // Items sharing a mirrorGroup share one label — unless their dots are so
      // far apart (opposite sides of the nucleosome) that a shared leader would
      // be longer than SPLIT_THRESHOLD; in that case each dot gets its own label.
      const groups = [];
      const seen   = new Set();

      for (let i = 0; i < labelList.length; i++) {
        if (seen.has(i)) continue;
        const it = labelList[i];
        if (it.mirrorGroup) {
          const grpIdx = [];
          for (let j = 0; j < labelList.length; j++) {
            if (!seen.has(j) && labelList[j].mirrorGroup === it.mirrorGroup) grpIdx.push(j);
          }
          grpIdx.forEach(j => seen.add(j));
          const dots = grpIdx.map(j => labelList[j]);
          // Split if any pair of dots exceeds the distance threshold.
          const tooFar = dots.length === 2 &&
            Math.hypot(dots[0].x - dots[1].x, dots[0].y - dots[1].y) > SPLIT_THRESHOLD;
          if (tooFar) {
            for (const d of dots) groups.push({ dots: [d], position: d.position, family: d.family });
          } else {
            groups.push({ dots, position: it.position, family: it.family });
          }
        } else {
          seen.add(i);
          groups.push({ dots: [it], position: it.position, family: it.family });
        }
      }

      // ── Natural angle = centroid direction from centre ──────────────────
      for (const g of groups) {
        const cx = g.dots.reduce((s, d) => s + d.x, 0) / g.dots.length;
        const cy = g.dots.reduce((s, d) => s + d.y, 0) / g.dots.length;
        g.naturalAngle = Math.atan2(cy - center.y, cx - center.x);
        g.rDot = Math.max(...g.dots.map(d => Math.hypot(d.x - center.x, d.y - center.y)));
      }

      // ── Spread angles (called once, before any rendering) ─────────────────
      const adjAngles = spreadAngles(groups.map(g => g.naturalAngle), MIN_ANG);

      // Clamp each adjusted angle to within MAX_ANG_DISP of its natural angle.
      // This prevents long diagonal leaders when the cluster is spread too far.
      for (let gi = 0; gi < adjAngles.length; gi++) {
        let delta = adjAngles[gi] - groups[gi].naturalAngle;
        // Normalise to (−π, π].
        while (delta >  Math.PI) delta -= 2 * Math.PI;
        while (delta < -Math.PI) delta += 2 * Math.PI;
        if (Math.abs(delta) > MAX_ANG_DISP) {
          adjAngles[gi] = groups[gi].naturalAngle + Math.sign(delta) * MAX_ANG_DISP;
        }
      }

      // ── Compute label positions, clamped to viewBox ───────────────────────
      for (let gi = 0; gi < groups.length; gi++) {
        const g   = groups[gi];
        const ang = adjAngles[gi];
        const ca  = Math.cos(ang), sa = Math.sin(ang);

        // Start at desired radius; pull inward if it falls outside the viewBox.
        let r = g.rDot + STEM_LEN + LEADER_LEN;
        while (r > g.rDot + STEM_LEN) {
          const lx = center.x + ca * r, ly = center.y + sa * r;
          if (lx >= bMinX && lx <= bMaxX && ly >= bMinY && ly <= bMaxY) break;
          r -= 0.5;
        }

        g.ang = ang;
        g.lx  = center.x + ca * r;
        g.ly  = center.y + sa * r;
      }

      // ── Render ────────────────────────────────────────────────────────────
      const frag = document.createDocumentFragment();

      for (const g of groups) {
        const color  = TRAJECTORY_COLOR[g.family] ?? '#888';
        const { lx, ly, ang } = g;
        const ca     = Math.cos(ang);
        const anchor = ca > 0.15 ? 'start' : ca < -0.15 ? 'end' : 'middle';

        // One leader per dot, all converging on the shared label point.
        for (const dot of g.dots) {
          const { x, y } = dot;
          const dotAng = Math.atan2(y - center.y, x - center.x);
          const bx     = x + Math.cos(dotAng) * STEM_LEN;
          const by     = y + Math.sin(dotAng) * STEM_LEN;

          const leader = document.createElementNS(ns, 'polyline');
          leader.setAttribute('points', `${x},${y} ${bx},${by} ${lx},${ly}`);
          leader.setAttribute('fill', 'none');
          leader.setAttribute('stroke', color);
          leader.setAttribute('stroke-width', LEADER_SW);
          leader.setAttribute('stroke-linecap', 'round');
          leader.setAttribute('stroke-linejoin', 'round');
          frag.appendChild(leader);

          const dotEl = document.createElementNS(ns, 'circle');
          dotEl.setAttribute('cx', x);  dotEl.setAttribute('cy', y);
          dotEl.setAttribute('r', DOT_R);
          dotEl.setAttribute('fill', color);
          frag.appendChild(dotEl);
        }

        // Single text label for the group.
        const text = document.createElementNS(ns, 'text');
        text.setAttribute('x', lx);              text.setAttribute('y', ly);
        text.setAttribute('text-anchor', anchor);
        text.setAttribute('dominant-baseline', 'central');
        text.setAttribute('font-size', LABEL_FS);
        text.setAttribute('font-family', 'sans-serif');
        // text.setAttribute('font-weight', 'bold');
        text.setAttribute('fill', color);
        text.setAttribute('paint-order', 'stroke');
        text.setAttribute('stroke', 'rgba(255,255,255,0.9)');
        text.setAttribute('stroke-width', '0.4');
        text.textContent = String(g.position);
        frag.appendChild(text);
      }

      group.appendChild(frag);
    }

    async function renderNucleosomeDiagram(slots, container, key) {
      const svgSrc = await fetchNucleosomeSvg();
      if (key !== lastKey) return;

      const svgEl = svgSrc.cloneNode(true);
      svgEl.style.overflow = 'visible';

      // Pass 1: count copies of each family to know whether to mirror.
      const copyCount = {};
      for (const slot of slots) copyCount[slot.family] = (copyCount[slot.family] || 0) + 1;

      // Pass 2: collect label positions.
      const copyIdx   = {};
      const placed    = new Set();
      const labelList = [];

      for (const slot of slots) {
        const fam = slot.family;
        copyIdx[fam] = (copyIdx[fam] || 0) + 1;
        const ci     = copyIdx[fam];
        const copies = copyCount[fam] === 1 ? [1, 2] : [ci];

        for (const mod of slot.modifications ?? []) {
          const pos    = mod.position;
          const maxLen = HISTONE_LENGTH[fam];
          if (!pos || !maxLen || pos < 1 || pos > maxLen) continue;

          for (const c of copies) {
            const uid = `${fam}:${c}:${pos}`;
            if (placed.has(uid)) continue;
            placed.add(uid);

            const pathEl = findTrajectory(svgEl, fam, c);
            if (!pathEl) continue;

            const { x, y } = getPathNode(pathEl, pos);
            // mirrorGroup tags a pair that should share one label text.
            const mirrorGroup = copies.length > 1 ? `${fam}:${pos}` : null;
            labelList.push({ family: fam, x, y, position: pos, mirrorGroup });
          }
        }
      }

      renderDiagramLabels(svgEl, labelList);

      // Remove trajectory paths from the rendered SVG.
      svgEl.querySelectorAll('[id$=":trajectory"]').forEach(el => el.remove());

      svgEl.setAttribute("xmlns", "http://www.w3.org/2000/svg");
      svgEl.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");

      // Expand viewBox to prevent label clipping (overflow:visible has no
      // effect inside <img>).  Pad by enough to cover leader + text extent.
      const pad = LABEL_FS;
      const vb  = svgEl.getAttribute('viewBox').trim().split(/\s+/).map(Number);
      svgEl.setAttribute('viewBox',
        `${vb[0] - pad} ${vb[1] - pad} ${vb[2] + 2 * pad} ${vb[3] + 2 * pad}`);

      // Serve as <img> so the browser offers "Save Image As…" on right-click.
      const blob = new Blob([new XMLSerializer().serializeToString(svgEl)],
                            { type: 'image/svg+xml' });
      const img  = document.createElement('img');
      img.src    = URL.createObjectURL(blob);
      img.style.width = '100%';
      container.appendChild(img);
    }
