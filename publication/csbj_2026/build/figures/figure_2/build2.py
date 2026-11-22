import re,subprocess,sys
src=open('src.dot').read()
LAB={
'GOBP_GOLGI_VESICLE_TRANSPORT':'Golgi vesicle\\ntransport',
'GOBP_NUCLEAR_EXPORT':'Nuclear export',
'GOBP_DNA_TEMPLATED_TRANSCRIPTION_ELONGATION':'DNA-templated\\ntranscription elongation',
'GOBP_GLYCOLIPID_BIOSYNTHETIC_PROCESS':'Glycolipid\\nbiosynthetic process',
'GOBP_REGULATION_OF_PROTEIN_CONTAINING_COMPLEX_DISASSEMBLY':'Regulation of protein-containing\\ncomplex disassembly',
'GOBP_POSITIVE_REGULATION_OF_MOLECULAR_FUNCTION':'Positive regulation of\\nmolecular function',
'GOCC_NUCLEAR_PORE_OUTER_RING':'Nuclear pore\\nouter ring',
'GOCC_CORE_MEDIATOR_COMPLEX':'Core Mediator\\ncomplex',
'GOCC_GOLGI_TRANSPORT_COMPLEX':'Golgi transport\\ncomplex',
'GOCC_ENDOPLASMIC_RETICULUM_PROTEIN_CONTAINING_COMPLEX':'ER protein-containing\\ncomplex',
'GOCC_ENDOPLASMIC_RETICULUM_GOLGI_INTERMEDIATE_COMPARTMENT_MEMBRANE':'ER-Golgi intermediate\\ncompartment membrane',
'GOCC_NUCLEAR_EXOSOME_RNASE_COMPLEX':'Nuclear exosome\\n(RNase complex)',
'GOCC_CULLIN_RING_UBIQUITIN_LIGASE_COMPLEX':'Cullin-RING ubiquitin\\nligase complex',
'GOCC_NUCLEAR_PORE':'Nuclear pore',
'GOCC_CLATHRIN_COATED_VESICLE_MEMBRANE':'Clathrin-coated\\nvesicle membrane',
'GOCC_SPINDLE':'Spindle',
'GOCC_AUTOPHAGOSOME':'Autophagosome',
'GOMF_STRUCTURAL_CONSTITUENT_OF_NUCLEAR_PORE':'Structural constituent\\nof nuclear pore',
'GOMF_CATALYTIC_ACTIVITY_ACTING_ON_RNA':'Catalytic activity,\\nacting on RNA',
'GOMF_GTPASE_ACTIVATOR_ACTIVITY':'GTPase activator\\nactivity',
'GOMF_TRANSCRIPTION_COREGULATOR_ACTIVITY':'Transcription\\ncoregulator activity',
'GOMF_ENZYME_SUBSTRATE_ADAPTOR_ACTIVITY':'Enzyme-substrate\\nadaptor activity',
'GOMF_ALPHA_TUBULIN_BINDING':'Alpha-tubulin\\nbinding',
'GOMF_CARGO_ADAPTOR_ACTIVITY':'Cargo adaptor\\nactivity',
'GOMF_HEXOSYLTRANSFERASE_ACTIVITY':'Hexosyltransferase\\nactivity',
'GOMF_G_PROTEIN_ACTIVITY':'G protein activity'}
edges=re.findall(r'"g_([^"]+)" -- "t_([^"]+)"',src)
terms=list(LAB); genes=sorted({g for g,_ in edges})
L=float(sys.argv[1]); FS=float(sys.argv[2]); seed=sys.argv[3] if len(sys.argv)>3 else '7'; mode=sys.argv[4] if len(sys.argv)>4 else 'KK'
COL={'GOBP':('#F6B3AC','#C9413B'),'GOCC':('#B7D9F2','#2F7EC1'),'GOMF':('#DECBEA','#8A58A8')}
o=['graph G {',f'graph [layout=neato, mode={mode}, overlap=prism, overlap_scaling=-6, sep="+3", splines=false, outputorder=edgesfirst, bgcolor=white, pad=0.04, start={seed}];',
   f'node [fontname="Liberation Serif Bold", fontsize={FS}];','edge [penwidth=0.6];']
for t in terms:
    f,_=COL[t[:4]]
    o.append(f'"t_{t}" [label="{LAB[t]}", shape=box, style="filled", fillcolor="{f}", color="#4A4A4A", penwidth=0.6, margin="0.075,0.02", width=0, height=0];')
for g in genes: o.append(f'"g_{g}" [label="", shape=circle, style=filled, fillcolor="#1A1A1A", color="#1A1A1A", width=0.04, height=0.04, fixedsize=true];')
for g,t in edges: o.append(f'"g_{g}" -- "t_{t}" [color="{COL[t[:4]][1]}", len={L}];')
o.append('}')
open('fig2.dot','w').write('\n'.join(o))
subprocess.run(['neato','-Tsvg','fig2.dot','-o','fig2.svg'],check=True)
s=open('fig2.svg').read(); w,h=[float(x) for x in re.search(r'viewBox="[-\d.]+ [-\d.]+ ([\d.]+) ([\d.]+)"',s).groups()]
print(f'L={L} seed={seed} natural {w/72:.2f} x {h/72:.2f} in; at 7.48in: font {FS*538.6/w:.1f} pt; height {h/w*7.48:.2f} in')
