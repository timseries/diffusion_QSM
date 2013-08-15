//landweber.cl: landweber iteration minimisation, multi-GPU ready

// Laplacian
#define Laplacian(LHS, RHS) \
if (r.x <= 1 && r.y <= 1 && r.z <= 1) \
{ \
  uint dist = r.x + r.y + r.z;  \
  if (dist==0)  \
    LHS -= RHS; \
  else if (dist==1) \
    LHS += C_10_96 * RHS; \
  else if (dist==2) \
    LHS += C_3_96 * RHS;  \
}

// Ax_b = A * x - b
// Calculate Ax and Laplacian Dx
__kernel void iterate1(__global Real* Ax_b,
                       __global Real* x,
                       __global Real* Dx,
                       __global int* FGindices,
                       __global Real* gx,
                       __global Real* gy,
                       __global Real* gz,
                       __global Real* ctr,
                       __global Real* sin2beta,
                       __global Real* mapx,
                       __global Real* mapy,
                       __global Real* mapz,
                       __global Real* skernel,
                       __global int* mask,
                       int dstart, int dend,
                       int fgstart, int fgend)
{
  const int p0 = get_global_id(0);
  const int P = p0 + dstart;
  if (P >= dend) return;
  int3 p;

  int size = dsize1 * dsize0;
  p.z = P / size;
  p.y = (P - p.z*size) / dsize0;
  p.x = P - p.z*size - p.y * dsize0;

  for (int O = fgstart; O < fgend; O++)
  {
    int mix_ = mask[FGindices[O]];
    int3 o;
    o.z = FGindices[O] / zoffset;
    o.y = (FGindices[O] - o.z * zoffset) / yoffset;
    o.x = FGindices[O] - o.y * yoffset - o.z * zoffset;

    uint3 r = abs(p-o);

    if (r.x > khalfsize || r.y > khalfsize || r.z > khalfsize) continue;

    // Linear system
    if (mix_ == -1) // spherical kernel
    {
      Ax_b[p0] += skernel[r.x+khalfsize + (r.y+khalfsize)*kyoffset + (r.z+khalfsize)*kzoffset] * x[O];
    }
    else if (r.x == 0 && r.y == 0 && r.z == 0)
    {
      Ax_b[p0] += ctr[mix_] * x[O];
    }
    else 
    {
      int3 r = p-o;
      Real rc = r.x*mapx[mix_] + r.y*mapy[mix_] + r.z*mapz[mix_];      
      if (fabs(rc) <= CYL2alpha)
      {
        Real r2 = 1.0f/(r.x*r.x + r.y*r.y + r.z*r.z - rc*rc); 
        Real rg = r.x*gx[mix_] + r.y*gy[mix_] + r.z*gz[mix_];
        
        Real h = CYL2alpha-fabs(rc);
        Real q = h*h * (CYL3alpha-h) * CYL4alpha3;
        
        Real val = CYLa * ( 2.0f*rg*rg*r2 - sin2beta[mix_]) * r2 * q; 
        //Another unnecessary branch?
        //if (fabs(val) >= thresholdB0) 
          Ax_b[p0] += val * x[O];
      }
    }

    // Laplacian
    //Rewritten laplacian to avoid calculating rx/ry/rz, slightly faster
    //Re-ordered laplacian to check z,y,x rather than x,y,z, much faster!
    //Replaced with distance function, even faster
    Laplacian(Dx[p0], x[O]);
  }
}

//Deltab
// Have Ax, Calculate Ax-b
__kernel void delta_b(__global Real* Ax_b,
                     __global Real* deltab,
                     int dstart, int dend)
{
  const int p0 = get_global_id(0);
  const int P = p0 + dstart;
  if (P >= dend) return;

  Ax_b[p0] -= deltab[p0];
}

// AtAx_b = A' * Ax_b
// Have Ax_b, Calculate A'Ax_b and Laplacian DtDx = D' * Dx
__kernel void iterate2(__global Real* AtAx_b,
                       __global Real* Ax_b,
                       __global Real* DtDx,
                       __global Real* Dx,
                       __global int* FGindices,
                       __global Real* gx,
                       __global Real* gy,
                       __global Real* gz,
                       __global Real* ctr,
                       __global Real* sin2beta,
                       __global Real* mapx,
                       __global Real* mapy,
                       __global Real* mapz,
                       __global Real* skernel,
                       __global int* mask,
                       int dstart, int dend, int doffset)
{
  const int O = get_global_id(0);

  if (O >= nFG) return;

  int mix_ = mask[FGindices[O]];

    Real3 g = (Real3)(gx[mix_], gy[mix_], gz[mix_]);
    Real sin2beta_ = sin2beta[mix_];
    Real ctr_ = ctr[mix_];
    Real3 map = (Real3)(mapx[mix_], mapy[mix_], mapz[mix_]);

  int oz = FGindices[O] / zoffset;
  int oy = (FGindices[O] - oz * zoffset) / yoffset;
  int ox = FGindices[O] - oy * yoffset - oz * zoffset;
   
  int pz = dstart / zoffset;
  int py = (dstart - pz * zoffset) / yoffset;
  int px = dstart - py * yoffset - pz * zoffset - 1;      

  uint3 r = (uint3)(abs(px-ox), abs(py-oy), abs(pz-oz));

  for (int P = dstart; P < dend; P++)
  {
    px++;
    if (px == dsize0)
    {
      px = 0;
      py++;
      if (py == dsize1)
      {
        py = 0;
        pz++;
        r.z = abs(pz-oz);
      }
      r.y = abs(py-oy);
    }
    r.x = abs(px-ox);

    if (r.x > khalfsize || r.y > khalfsize || r.z > khalfsize) continue;

    // Linear system
    if (mix_ == -1) // spherical kernel
    {
      AtAx_b[O] += skernel[r.x+khalfsize + (r.y+khalfsize)*kyoffset + (r.z+khalfsize)*kzoffset] * Ax_b[P-doffset];
    }
    else if (r.x == 0 && r.y == 0 && r.z == 0)
    {
      AtAx_b[O] += ctr_ * Ax_b[P-doffset];
    }
    else 
    {
      int3 r = (int3)(px-ox, py-oy, pz-oz);
      Real rc = r.x*map.x + r.y*map.y + r.z*map.z;  
      if (fabs(rc) <= CYL2alpha)
      {
        Real r2 = 1.0f/(r.x*r.x + r.y*r.y + r.z*r.z - rc*rc); 
        Real rg = r.x*g.x + r.y*g.y + r.z*g.z;
        
        Real h = CYL2alpha-fabs(rc);
        Real q = h*h * (CYL3alpha-h) * CYL4alpha3;
        
        Real val = CYLa * ( 2.0f*rg*rg*r2 - sin2beta_) * r2 * q; 
        //if (fabs(val) >= thresholdB0) 
          AtAx_b[O] += val * Ax_b[P-doffset];
      }
    }

    // Laplacian
    Laplacian(DtDx[O], Dx[P-doffset]);
  }
}

__kernel void zero(__global Real* Ax_b,
                   __global Real* Dx,
                   __global Real* AtAx_b,
                   __global Real* DtDx,
                   int dN)
{
//Multi-gpu version
  for (int i=get_global_id(0); i<N; i += get_global_size(0))
  {
    if (i < dN)
    {
      Ax_b[i] = 0.0f;
      Dx[i] = 0.0f;
    }
    if (i < nFG)
    {
      AtAx_b[i] = 0.0f;
      DtDx[i] = 0.0f;
    }
  }
}

//Calc new x and rms
__kernel void rms_new_x(__global Real* AtAx_b,
                        __global Real* DtDx,
                        __global Real* x,
                        __global struct OutputCL* out
                       )
{
  const int O = get_global_id(0);
  const int thread = get_local_id(0);
  const int group = get_group_id(0);

  __local Real rms_x[nthreads];
  __local Real rms_diff_x[nthreads];

  // Calculate new x and rms values
  if (O < nFG)
  {
    Real new_x = x[O] - alpha * tau * AtAx_b[O] - beta * tau * DtDx[O];
    rms_x[thread] = new_x * new_x;
    rms_diff_x[thread] = (new_x - x[O]) * (new_x - x[O]);
    x[O] = new_x;  //Copy result to x
  }
  else
  {
    rms_x[thread] = 0.0f;
    rms_diff_x[thread] = 0.0f;
  }

  //Reduction
  barrier(CLK_LOCAL_MEM_FENCE);
  for(unsigned int s=1; s < nthreads; s *= 2)
  {
    if ((thread % (2*s)) == 0 && thread+s < nthreads)
    {
      rms_x[thread] += rms_x[thread + s];
      rms_diff_x[thread] += rms_diff_x[thread + s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  // write result for this group to global mem
  if (thread == 0) 
  {
    out[group].rms_x = rms_x[0];
    out[group].rms_diff_x = rms_diff_x[0];
  }
}

__kernel void collect(__global struct OutputCL* out)
{
  if (get_global_id(0) == 0)
  {
    for (int i=1; i<get_num_groups(0); i++)
    {
      out[0].rms_x += out[i].rms_x;
      out[0].rms_diff_x += out[i].rms_diff_x;
    }
  }
}   
