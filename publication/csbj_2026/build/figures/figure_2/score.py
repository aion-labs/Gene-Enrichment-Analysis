import subprocess,sys,re
def seg_rect(p,q,r):
    (x1,y1),(x2,y2)=p,q; cx,cy,w,h=r
    xmin,xmax,ymin,ymax=cx-w/2,cx+w/2,cy-h/2,cy+h/2
    # Liang-Barsky
    t0,t1=0,1; dx,dy=x2-x1,y2-y1
    for pp,qq in ((-dx,x1-xmin),(dx,xmax-x1),(-dy,y1-ymin),(dy,ymax-y1)):
        if pp==0:
            if qq<0: return False
        else:
            t=qq/pp
            if pp<0: t0=max(t0,t)
            else: t1=min(t1,t)
    return t0<t1
def score(seed,L,mode):
    subprocess.run(['python3','build2.py',str(L),FS,str(seed),mode],capture_output=True)
    out=subprocess.run(['neato','-Tplain','fig2.dot'],capture_output=True,text=True).stdout.split('\n')
    W,H=map(float,out[0].split()[2:4]); N={};E=[]
    for l in out:
        f=l.split()
        if f and f[0]=='node': N[f[1].strip('"')]=(float(f[2]),float(f[3]),float(f[4]),float(f[5]))
        if f and f[0]=='edge': E.append((f[1].strip('"'),f[2].strip('"')))
    bad=0
    for a,b in E:
        for n,r in N.items():
            if n.startswith('"t_') or n.startswith('t_'):
                if n in (a,b): continue
                if seg_rect(N[a][:2],N[b][:2],r): bad+=1
    # node-node overlap gene dot inside box
    dots=0
    for n,r in N.items():
        if n.startswith('t_'):
            for g,q in N.items():
                if g.startswith('g_') and abs(q[0]-r[0])<r[2]/2+0.03 and abs(q[1]-r[1])<r[3]/2+0.03: dots+=1
    return bad,dots,W,H
FS=sys.argv[1]
res=[]
for mode in ['KK','major']:
  for L in [0.22,0.25,0.28]:
    for s in range(1,41):
        b,d,W,H=score(s,L,mode); res.append((b+5*d+ (W*H)/2,b,d,round(W,2),round(H,2),s,L,mode))
res.sort()
for r in res[:10]: print(r)
