const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch(); const p = await b.newPage();
  for (const n of process.argv.slice(2)) {
    const m = n.match(/^Supp (Table|Note) (\w+)/); const llm = n.includes('LLM answers');
    const label = llm ? 'LLM interpretation' : (m[1]==='Table' ? 'Supplementary Table ' + m[2].replace(/^0/,'') : 'Supplementary Note ' + m[2]);
    await p.goto('file://' + process.cwd() + '/' + n + '.html');
    await p.pdf({ path: n + '.pdf', format: 'A4', landscape: process.env.LANDSCAPE==='1', printBackground: true, displayHeaderFooter: true, headerTemplate: '<span></span>',
      footerTemplate: `<div style="font:7pt Helvetica,Arial,sans-serif;color:#888;width:100%;text-align:center;">iGEA · ${label} · page <span class="pageNumber"></span> of <span class="totalPages"></span></div>`,
      margin: { top: '16mm', bottom: '18mm', left: '15mm', right: '15mm' } });
    console.log('pdf', n);
  }
  await b.close();
})();
