#!/usr/bin/env python3
"""Render ORIGINAL EFG numerical TM11/TM21 Hx/Hy samples to copyright-clean SVG."""
import argparse
import csv
import math
from pathlib import Path

def render(src,dest,title):
    with Path(src).open(encoding="utf-8",newline="") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=25*13:
        raise ValueError("Expected 25x13 exported field grid")
    points=[]
    for row in rows:
        x,y,hx,hy=[float(row[k]) for k in ("sample_x","sample_y","hx","hy")]
        if any(not math.isfinite(v) for v in (x,y,hx,hy)):
            raise ValueError("Nonfinite EFG H-vector data")
        points.append((x,y,hx,hy))
    selected=[p for p in points if round(p[0]*12)%2==0 and round(p[1]*12)%2==0]
    top=max(math.hypot(p[2],p[3]) for p in selected)
    if top<=0:raise ValueError("Zero EFG H field")
    # Scaling is arbitrary because eigenvectors have an undetermined norm
    # and global sign. Preserve RELATIVE arrow lengths to the maximum.
    W,H=960,540
    pad=60; sx=(W-2*pad)/2;sy=H-2*pad
    parts=[
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{W/2}" y="30" text-anchor="middle" font-size="20" '
        f'font-family="sans-serif">{title}</text>',
        f'<rect x="{pad}" y="{pad}" width="{W-2*pad}" height="{H-2*pad}" '
        'fill="none" stroke="#555" stroke-width="1.5"/>'
    ]
    for x,y,hx,hy in selected:
        X=pad+x*sx;Y=H-pad-y*sy
        # Arrow length is calibrated for the physical 2x1 rectangle.
        k=21/top
        dx=k*hx;dy=-k*hy
        ex=X+dx;ey=Y+dy
        length=math.hypot(dx,dy)
        if length<0.25:continue
        ux=dx/length;uy=dy/length
        lx=ex-6*ux-2.5*uy;ly=ey-6*uy+2.5*ux
        rx=ex-6*ux+2.5*uy;ry=ey-6*uy-2.5*ux
        parts.append(
            f'<path d="M {X:.3f},{Y:.3f} L {ex:.3f},{ey:.3f} '
            f'M {lx:.3f},{ly:.3f} L {ex:.3f},{ey:.3f} '
            f'L {rx:.3f},{ry:.3f}" fill="none" stroke="#2563eb" '
            'stroke-width="1.1" stroke-linecap="round"/>'
        )
    parts.extend([
        f'<text x="{W/2}" y="{H-17}" text-anchor="middle" font-size="14" '
        'font-family="sans-serif">x (0 a 2 m)</text>',
        f'<text x="18" y="{H/2}" transform="rotate(-90 18,270)" '
        'text-anchor="middle" font-size="14" font-family="sans-serif">y (0 a 1 m)</text>',
        '</svg>',
    ])
    dest=Path(dest)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text("\n".join(parts)+"\n",encoding="utf-8")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--title",required=True)
    args=p.parse_args()
    render(args.input,args.output,args.title)

if __name__=="__main__":
    main()
