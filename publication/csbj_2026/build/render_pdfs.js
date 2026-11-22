// Render supplementary HTML files to A4 PDF with the footer "iGEA · Supplementary Table N [· LLM interpretation] · page x of y".
// Usage: node render_pdfs.js <folder> [file.html ...]   (default: every "Supp Table*.html" and "Supp Note*.html" in the folder)
// Requires Node.js with Playwright (npm i playwright; npx playwright install chromium).
const { chromium } = require('playwright'); const fs = require('fs'); const path = require('path');
const LANDSCAPE = new Set(['05', '06', '08']);           // wide tables
const FIT = new Map([['01', 1373], ['02', 1300]]);       // one page sized to the content (width in px; no footer), as submitted
(async () => {
  const dir = path.resolve(process.argv[2] || '.');
  const files = process.argv.length > 3 ? process.argv.slice(3) :
    fs.readdirSync(dir).filter(f => /^Supp (Table|Note) .*\.html$/.test(f));
  const b = await chromium.launch(); const p = await b.newPage();
  for (const f of files) {
    const m = f.match(/^Supp (Table|Note) (\w+)/); const num = m[2].replace(/^0/, '');
    const label = (m[1] === 'Table' ? 'Supplementary Table ' : 'Supplementary Note ') + num + (f.includes('LLM answers') ? ' · LLM interpretation' : '');
    const land = m[1] === 'Table' && LANDSCAPE.has(m[2]);
    await p.goto('file://' + path.join(dir, f).split('/').map(encodeURIComponent).join('/').replace(/%3A/g, ':'));
    if (m[1] === 'Table' && FIT.has(m[2])) {
      const w = FIT.get(m[2]); await p.setViewportSize({ width: w, height: 900 }); await p.emulateMedia({ media: 'screen' });
      const h = await p.evaluate(() => document.documentElement.scrollHeight);
      await p.pdf({ path: path.join(dir, f.replace(/\.html$/, '.pdf')), width: w + 'px', height: (h + 40) + 'px', printBackground: true, pageRanges: '1' });
      await p.setViewportSize({ width: 1280, height: 720 }); await p.emulateMedia({ media: 'print' }); console.log('pdf (fit)', f); continue;
    }
    await p.pdf({ path: path.join(dir, f.replace(/\.html$/, '.pdf')), format: 'A4', landscape: land, printBackground: true,
      displayHeaderFooter: true, headerTemplate: '<span></span>',
      footerTemplate: `<div style="font:7pt Helvetica,Arial,sans-serif;color:#888;width:100%;text-align:center;">iGEA · ${label} · page <span class="pageNumber"></span> of <span class="totalPages"></span></div>`,
      margin: land ? { top: '14mm', bottom: '16mm', left: '12mm', right: '12mm' } : { top: '16mm', bottom: '18mm', left: '15mm', right: '15mm' } });
    console.log('pdf', f);
  }
  await b.close();
})();
