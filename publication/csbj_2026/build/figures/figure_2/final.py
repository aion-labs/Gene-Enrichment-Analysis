import re,subprocess,sys
L,FS,seed,mode=sys.argv[1:5]
subprocess.run(['python3','build2.py',L,FS,seed,mode],check=True)
subprocess.run(['neato','-Tdot','fig2.dot','-o','laid.dot'],check=True)
s=open('laid.dot').read()
W,H=map(float,re.search(r'bb="0,0,([\d.]+),([\d.]+)"',s).groups())
kx,ky=float(sys.argv[5]) if len(sys.argv)>5 else W-50, float(sys.argv[6]) if len(sys.argv)>6 else H-26
rows=[('#F6B3AC','GO Biological Process'),('#B7D9F2','GO Cellular Component'),('#DECBEA','GO Molecular Function')]
tr=''.join(f'<TR><TD WIDTH="10" HEIGHT="8" BGCOLOR="{c}"></TD><TD ALIGN="LEFT"><FONT POINT-SIZE="{FS}">{t}</FONT></TD></TR>' for c,t in rows)
tr+=f'<TR><TD WIDTH="10" HEIGHT="8"><FONT POINT-SIZE="9">●</FONT></TD><TD ALIGN="LEFT"><FONT POINT-SIZE="{FS}">Gene</FONT></TD></TR>'
key=f'\tkey [shape=plaintext, margin=0, pos="{kx},{ky}!", label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="1" CELLPADDING="1">{tr}</TABLE>>];\n'
i=s.rstrip().rfind('}'); s=s[:i]+key+'}\n'
s=s.replace('graph [','graph [size="7.48,7.48",',1)
open('final.dot','w').write(s)
subprocess.run(['neato','-n2','-Tsvg','final.dot','-o','f.svg'],check=True)
v=open('f.svg').read().replace('font-family="Liberation Serif Bold"','font-family="Liberation Serif" font-weight="bold"')
def rect(m):
    P=[tuple(map(float,p.split(','))) for p in m.group(3).split()]
    xs=[p[0] for p in P]; ys=[p[1] for p in P]
    return f'<rect x="{min(xs):.2f}" y="{min(ys):.2f}" width="{max(xs)-min(xs):.2f}" height="{max(ys)-min(ys):.2f}" rx="2" ry="2" fill="{m.group(1)}" stroke="{m.group(2)}" stroke-width="0.6"'
v=re.sub(r'<polygon fill="(#f6b3ac|#b7d9f2|#decbea)" stroke="(#4a4a4a)" stroke-width="[\d.]+" points="([^"]+)"',rect,v,flags=re.I)
open('f3.svg','w').write(v)
subprocess.run(['rsvg-convert','-w','1400','f3.svg','-b','white','-o',f'chk_{seed}.png'])
print(W,H,W/72,'in -> font',float(FS)*538.6/W)
