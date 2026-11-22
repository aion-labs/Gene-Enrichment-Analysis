const { chromium } = require('playwright');
(async () => {
  const [,, inDir, outDir, ...names] = process.argv;
  const b = await chromium.launch(); const p = await b.newPage();
  for (const n of names) {
    await p.goto('file://' + inDir + '/' + n + '.html');
    await p.pdf({ path: outDir + '/' + n + '.pdf', format: 'A4', printBackground: true, displayHeaderFooter: true,
      headerTemplate: '<span></span>',
      footerTemplate: '<div style="font:7pt Helvetica,Arial,sans-serif;color:#888;width:100%;text-align:center;">iGEA · LLM interpretation · page <span class="pageNumber"></span> of <span class="totalPages"></span></div>',
      margin: { top: '16mm', bottom: '18mm', left: '15mm', right: '15mm' } });
    console.log('pdf', n);
  }
  await b.close();
})();
