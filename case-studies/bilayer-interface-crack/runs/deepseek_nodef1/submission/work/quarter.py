import numpy as np, fem

def quarter_point(mesh):
    coords=mesh['coords']; x=mesh['x']; y1=mesh['y1']; y2=mesh['y2']
    dx=x[1]  # first positive spacing from 0 (x index nxcorner=0? x[0]=-Lx, x[nxcorner]=0)
    # find index of x=0
    i0=np.argmin(np.abs(x))
    dxp=x[i0+1]-x[i0]; dxm=x[i0]-x[i0-1]
    dy1=y1[1]-y1[0]
    dy2=y2[1]-y2[0]
    # move nodes along interface (x=0) at y=+-dy/2
    for i,(xx,yy) in enumerate(coords):
        if abs(xx)<1e-14:
            if abs(abs(yy)-dy1/2)<1e-12:
                coords[i,1]=np.sign(yy)*dy1/4
            if abs(abs(yy)-dy2/2)<1e-12:
                coords[i,1]=np.sign(yy)*dy2/4
    # move nodes along y=0 at x=+-dx/2
    for i,(xx,yy) in enumerate(coords):
        if abs(yy)<1e-14:
            if abs(xx-dxp/2)<1e-12:
                coords[i,0]=dxp/4
            if abs(xx+dxm/2)<1e-12:
                coords[i,0]=-dxm/4
    return mesh
