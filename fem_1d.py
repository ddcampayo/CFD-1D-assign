# Implementation of "particle" class
# and "mesh" class

# Periodic boundary conditions

# Full projection assignment
# FLIP assignment

import numpy as np

#basic "particle" container

class parts(object):
  def __init__( self , N , do_quad = False ):
    self.do_quad = do_quad
    self.positions(N)
    self.volumes()

  def repos( self , r ):
    self.r = r
    self.N = self.r.size
    self.volumes()

  def displace( self , dr ):
    self.r += dr
    N = self.N

    for i in xrange(0,N):
      if (self.r[i] > 1) :
        self.r[i] -= 1
      elif (self.r[i] < 0) :
        self.r[i] += 1

    self.sort()
    #self.volumes()

  def randomize( self  ):
    #    self.r=np.zeros( N )

    N = self.N

    #include zero:
    self.r[0] = 0

    self.r[1:N] = np.sort(np.random.random( N-1 ) )

    self.N = self.r.size # just in case

    self.volumes()

    
  def sort( self  ):
    N = self.N

    for i in xrange(0,N):
      if (self.r[i] > 1) :
        self.r[i] -= 1
      elif (self.r[i] < 0) :
        self.r[i] += 1

    # sort positions, keeping data at each node
    idx   = np.argsort( self.r )
    self.r = self.r[idx];
    self.f = self.f[idx];
    self.volumes()

  def perturb( self , amount ):
    self.r

    N = self.N

    for i in xrange(0,N):

      hi=self.h[i]

      hm1=self.h[  (i - 1 + N) % N ]

      h = np.min(hi,hm1)
      
      ran = 2 * np.random.random() - 1 # in (-1 , 1)
      
      self.r[i] += ran * h * amount

    self.sort()

  def positions(self, N):

    #rp=np.zeros( Np )
    #include zero:
    #r[0] = 0

    # random.-
    #r[1:N] = np.sort(np.random.random( N-1 ) )

    # equispaced, cyclic

    self.r = np.linspace(0, 1 - 1.0/N , N )

    self.N = self.r.size # just in case


  def volumes(self):

    N = self.N

    self.h=np.zeros( N )
    self.h[0:N-1] = self.r[1:N] - self.r[0:N-1]

    self.h[N-1] = 1 + self.r[0] - self.r[N-1]

    self.h_mean = 1.0 / N

    self.vol=np.zeros( N )

#   particle volume:
    self.vol[ 1 : N ] = ( self.h[ 1 : N ] + self.h[0 : N-1] )/2.0
    self.vol[ 0 ] = ( self.h[ 0 ] + self.h[ N-1] )/2.0

    if self.do_quad:
      self.quad_compute()

  # def randomize_positions(self, alpha)
  #   N = self.N
  #   self.f = np.zeros( N )

  #   for ii in xrange( 0 , N ) :
  #     self.f[ii] = part.valueat( self.r[ii] )


  def quad(self):
    return self.do_quad

  def size(self):
    return self.N

  # def N(self):
  #   return self.N

  def assign_f( self , part ) :
    #    N = self.N
    #    self.f = np.zeros( N )

    #    for ii in xrange( 0 , N ) :
    #      self.f[ii] = part.valueat( self.r[ii] )
    self.f = part.valueat( self.r )

# Same thing .-

  def onto_delta( self , part ):
    self.assign_f( part )

  def create_f( self , f) :
    self.f = f

  def create_tophat_f( self ) :
    self.f = self.tophat_f()

  def tophat_f( self ) :
    N = self.N

    f = np.zeros( N )
    i = 0

    for x in self.r:
      if ((x>0.25) and (x<0.75)):
        f[i] = 1
      else:
        f[i] = 0
      i+=1

    return f

  def fidelity( self ) :
    f0 = self.tophat_f()
    return np.linalg.norm( self.f - f0) / np.linalg.norm(f0)

  
  def create_linear_f( self ) :
    self.f = self.r

  def create_quadratic_f( self ) :
    self.f = self.r**2

  def create_cos_f( self ) :
    self.f = np.cos( 2 * np.pi * self.r )
    
  def check_cos_f( self ) :
    dd = self.f - np.cos( 2 * np.pi * self.r )
    return np.linalg.norm(dd)
    
  def valueat( self , x ):
    r = self.r
    N = self.N

    #    ip1= ( np.searchsorted(r, x + 1e-16) ) % N
    ip1= ( np.searchsorted(r, x ) ) % N
    ii= (ip1 + N -1) % N
    xx=x-r[ii]
    if xx < 0:
      xx += 1
  
    h1 = 1 - xx/self.h[ii]
    h2=      xx/self.h[ii]

    f = self.f

    #    print h1,h2,f[ii] * h1 + f[ip1] * h2
    
    if self.do_quad:
      h12= h1*h2

      f12=0

      A = self.A

      for n in xrange( 0 , 4) :
        jj = (ii + n - 1 + N ) % N
        fj= f[ jj ]
        f12 += A[ ii , n ] * fj

      return f[ii] * h1 + f[ip1] * h2 + f12 * h12

    else:
      return f[ii] * h1 + f[ip1] * h2

  def int_f(self) :
    return np.dot( self.vol , self.f )

  def int_f2(self) :
    return np.dot( self.vol , np.power(self.f,2) )

  def quad_compute(self):
    N = self.N
    h = self.h

    self.A = np.zeros( shape = (N,4) )

    for ii in xrange( 0 , N ) :

      hi=h[ii]

      hm1=h[  (ii - 1 + N) % N ]

      hp1=h[  (ii+1) % N ]

      # minimum norm calculation

      M=np.matrix(
        [
          [           1 ,   1   , 1     , 1            ],
          [ -hm1-hi/2   , -hi/2 ,  hi/2 , hp1+hi/2     ],
          [ hm1*(hm1+hi),  0    ,  0    , hp1*(hp1+hi) ]
        ]
      )

      b=np.matrix( [ 0 , 0 , -hi**2 ] )

      # factor.-
      M[1,:] /= self.h_mean

      M[2,:] /= self.h_mean**2

      b /= self.h_mean**2

      aa = self.svdsolve( M , b.T )

      self.A[ii,:] = aa.T

  def svdsolve( self, a  , b):
#    import numpy as np
    u,s,v = np.linalg.svd(a,  full_matrices= False)
    c = np.dot( u.T , b )
    #    w = np.linalg.solve(np.diag(s),c)
    w = np.dot( np.diag(1.0/s) , c )
    x = np.dot( v.T , w )
    return x

  def onto_lumped( self , part ):
    self.onto_common( part )
    self.f /= self.vol

  def onto_common( self , part ):
    N = self.N
    ff = np.zeros( N )

    r = self.r
    h = self.h

    for ii in xrange( 0 , N ) :
      ip1 = (ii+1) % N
      r1 = r[ii]
      r2 = r[ip1]
      r12 = r1 + h[ii]/2
      if r12 > 1 :
        r12 -= 1
        
      f1 = part.valueat( r1 )
      f2 = part.valueat( r2 )
      f12 =part.valueat( r12)

      a = f1
      b = f2

      # Important integration choice:
      #if part.quad() :
      if True:
      #if False:
        c = 4*f12 - 2 * (a+b)
        ff[ii]  += ( a/3 + b/6 + c/12 ) * h[ii]
        ff[ip1] += ( a/6 + b/3 + c/12 ) * h[ii]
      else:
        ff[ii]  += ( a/3 + b/6 ) * h[ii]
        ff[ip1] += ( a/6 + b/3 ) * h[ii]

    self.f = ff 


class parts_at_random(parts):
#override positions function:
  def positions(self, N):

    self.r=np.zeros( N )
    #include zero:
    self.r[0] = 0

    self.r[1:N] = np.sort(np.random.random( N-1 ) )

    self.N = self.r.size # just in case



class mesh(parts):
#  "Derived 1D fem object"

#override parent constructor:
  def __init__( self , N , do_quad = False ):
    parts.__init__( self , N , do_quad )
    self.build_matrices()

  def repos( self , r ):
    parts.repos(self, r)
    self.build_matrices()

  def displace( self , dr ):
    parts.displace(self, dr)
    self.build_matrices()

  def randomize( self  ):    
    parts.randomize(self)
    self.build_matrices()

  def perturb( self , amount ):
    parts.perturb(self, amount)
    self.build_matrices()
   
    
  def build_matrices(self) :
    
    N = self.N
 
    self.s=np.zeros( shape=(N,N) )
    self.m=np.zeros( shape=(N,N) )
    self.l=np.zeros( shape=(N,N) )


    for ii in xrange( 0 , N ) :

      hi = self.h[ii]
      a = 1/hi

      jj = (ii +1 ) % N

      self.s[ii,ii] +=  a
      self.s[jj,ii]  = -a
      self.s[ii,jj]  = -a
      self.s[jj,jj] +=  a

      self.l[ii,ii] -=  0.5 
      self.l[ii,jj]  =  0.5 
      self.l[jj,ii]  = -0.5 
      self.l[jj,jj] +=  0.5 

      #lumped:
      #  m0(ii,ii)+=hi/2
      #  m0(jj,jj)+=hi/2

      #unlumped:
      self.m[ii,ii]+=hi/3
      self.m[jj,ii] =hi/6
      self.m[ii,jj] =hi/6
      self.m[jj,jj]+=hi/3

    if self.do_quad:

      for ii in xrange( 0 , N ) :

        hi=self.h[ii]

        for n in xrange(0,4) :
 
          jj= (ii - 1 + n + N ) % N

          AA=self.A[ii,n]

          aa=AA*hi/12

          self.m[ii  ,jj] += aa
          self.m[jj  ,ii] += aa

          self.m[ (ii +1)  % N , jj ] += aa
          self.m[ jj ,  (ii +1) % N ] += aa

          for n2 in xrange(0,4) :

            kk = (ii - 1 +n2 + N ) % N

            A2 = self.A[ii,n2]

            self.m[ jj , kk] += AA*A2*hi/30
            #      m(kk,jj)+=AA*A2*hi/30  ## it don't need this

            self.s[ jj , kk] += AA*A2/(3*hi)

            #      s(kk,jj)+=AA*A2/(3*hi)

            # s2[ jj , kk] += 4*AA*A2/(hi^3)

  def onto_full( self , part ):
    self.onto_common( part )

    N = self.N
    ff = np.zeros( N )

    ff = self.f

    self.f = np.linalg.solve( self.m , ff )

    
  def onto_full_old( self , part ):
    N = self.N
    ff = np.zeros( N )

    r = self.r
    h = self.h

    for ii in xrange( 0 , N ) :
      ip1 = (ii+1) % N
      r1 = r[ii]
      r2 = r[ip1]
      r12 = r1 + h[ii]/2
      if r12 > 1 :
        r12 -= 1
        
      f1 = part.valueat( r1 )
      f2 = part.valueat( r2 )
      f12 =part.valueat( r12)

      a = f1
      b = f2

      # Important integration choice:
      #if part.quad() :
      if True:
      #if False:
        c = 4*f12 - 2 * (a+b)
        ff[ii]  += ( a/3 + b/6 + c/12 ) * h[ii]
        ff[ip1] += ( a/6 + b/3 + c/12 ) * h[ii]
      else:
        ff[ii]  += ( a/3 + b/6 ) * h[ii]
        ff[ip1] += ( a/6 + b/3 ) * h[ii]
        
    self.f = np.linalg.solve( self.m , ff )


    

# nodal functions aka filters
# todo: early termination if outside suppor of node i

  def phi(self , i , r) :
    N = self.N

    ri = self.r[i]

    # distance :
    dd = r - ri

    if (dd >= 0.5 ) :
      dd =  dd - 1

    if (dd < -0.5 ) :
      dd =  dd + 1

    if (dd >= 0 ) :
      h = self.h[ i ]
      if (dd > h ) :
        return 0
      else :
        return 1 - dd/h

    if (dd < 0 ) :
      h = self.h[ (i - 1 ) % N ]
      if (dd < -h ) :
        return 0
      else :
        return 1 + dd/h



#   FLIP volume:
# TODO: extended functions

  def flip_volumes(self , part) :
    N = self.N
    r = self.r

    self.fvol=np.zeros( N )

    for p in xrange( 0 , part.size() ) :
    
      x = part.r[ p ]

      ip1= ( np.searchsorted( r , x ) ) % N
      ii= (ip1 + N -1) % N

      # self.fvol[i] += phi_i * part.f[ p ] * part.vol[ p ]

      self.fvol[ip1] += self.phi(ip1 , x ) * part.vol[ p ]
      self.fvol[ii ] += self.phi(ii  , x ) * part.vol[ p ]

      
    # in case no volume is gathered.-

    for i in xrange( 0 , N ) :
      if (self.fvol[i] < 1e-16) :
        self.fvol[i] = self.vol[i]

#   FLIP assign:
# TODO: extended functions

  def flip_assign(self , part) :

#    self.flip_volumes(part)

    N = self.N
    r = self.r

    self.f = np.zeros( N )

    for p in xrange( 0 , part.size() ) :
    
      x = part.r[ p ]

      ip1= ( np.searchsorted( r , x ) ) % N
      ii= (ip1 + N -1) % N

      # self.fvol[i] += phi_i * part.f[ p ] * part.vol[ p ]

      self.f[ip1] += self.phi(ip1 , x ) * part.f[ p ] * part.vol[ p ]
      self.f[ii ] += self.phi(ii  , x ) * part.f[ p ] * part.vol[ p ]
        
    self.f /= self.fvol


  def assign_f( self , part ) :
    N = self.N
    self.f = np.zeros( N )

    for ii in xrange( 0 , N ) :
      self.f[ii] = part.valueat( self.r[ii] )


class mesh_at_random(mesh):
#override positions function:
  def positions(self, N):

    self.r=np.zeros( N )
    #include zero:
    self.r[0] = 0

    self.r[1:N] = np.sort(np.random.random( N-1 ) )

    self.N = self.r.size # just in case

class mesh_with_pos(mesh):
#override __init__ function:

  def __init__( self , r , do_quad = False ):
    self.do_quad = do_quad
    self.r = r
    self.N = self.r.size
    self.volumes()


